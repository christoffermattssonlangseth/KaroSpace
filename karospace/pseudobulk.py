"""Pseudobulk differential expression utilities."""

from __future__ import annotations

import os
import re
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pandas.api.types import CategoricalDtype

from .console import log_detail, log_step, log_warning


_NUMPY_SLOGDET_WARNING_MESSAGE = r".*encountered in slogdet"
_NUMPY_SLOGDET_PYTHONWARNINGS_FILTERS = (
    "ignore:divide by zero encountered in slogdet:RuntimeWarning",
    "ignore:overflow encountered in slogdet:RuntimeWarning",
    "ignore:invalid value encountered in slogdet:RuntimeWarning",
)


def _install_numpy_slogdet_warning_filter() -> None:
    """Install a targeted slogdet warning filter for this process and child workers."""
    warnings.filterwarnings(
        "ignore",
        message=_NUMPY_SLOGDET_WARNING_MESSAGE,
        category=RuntimeWarning,
    )
    existing = os.environ.get("PYTHONWARNINGS", "")
    filters = [part.strip() for part in existing.split(",") if part.strip()]
    updated = False
    for filter_spec in _NUMPY_SLOGDET_PYTHONWARNINGS_FILTERS:
        if filter_spec not in filters:
            filters.append(filter_spec)
            updated = True
    if updated:
        os.environ["PYTHONWARNINGS"] = ",".join(filters)


_install_numpy_slogdet_warning_filter()


def _as_count_matrix(adata, counts_layer: Optional[str]) -> Tuple[Any, str, Optional[str]]:
    if counts_layer:
        if counts_layer in adata.layers:
            return adata.layers[counts_layer], str(counts_layer), None
        return adata.X, "X", f"counts layer '{counts_layer}' not found; using adata.X"
    return adata.X, "X", None


def _to_dense_counts(matrix) -> np.ndarray:
    dense = matrix.toarray() if sp.issparse(matrix) else np.asarray(matrix)
    dense = np.asarray(dense, dtype=np.float64)
    dense[~np.isfinite(dense)] = 0
    dense[dense < 0] = 0
    return np.rint(dense).astype(np.int64, copy=False)


def _cell_count_mask(matrix, min_cell_counts: int) -> np.ndarray:
    """Return cells meeting a raw total-count threshold without densifying cells x genes."""
    threshold = max(0, int(min_cell_counts))
    if threshold == 0:
        return np.ones(int(matrix.shape[0]), dtype=bool)
    totals = np.asarray(matrix.sum(axis=1)).ravel() if sp.issparse(matrix) else np.asarray(matrix).sum(axis=1)
    totals = np.asarray(totals, dtype=float)
    return np.isfinite(totals) & (totals >= threshold)


def _filter_pseudobulk_genes(
    counts: np.ndarray,
    metadata: pd.DataFrame,
    min_gene_counts: int,
) -> Tuple[np.ndarray, pd.DataFrame]:
    """Drop genes below a raw aggregate-count threshold for one DESeq2 fit."""
    threshold = max(0, int(min_gene_counts))
    if threshold == 0:
        return counts, metadata
    totals = np.asarray(counts, dtype=float).sum(axis=0)
    keep = np.isfinite(totals) & (totals >= threshold)
    filtered_counts = np.asarray(counts)[:, keep]
    filtered_meta = metadata.copy()
    gene_names = [str(g) for g in metadata.attrs.get("gene_names", [])]
    filtered_meta.attrs["gene_names"] = [
        gene for gene, include in zip(gene_names, keep) if include
    ]
    return filtered_counts, filtered_meta


def _positive_fraction(matrix, mask: np.ndarray, gene_indices: Sequence[int]) -> List[Optional[float]]:
    if not mask.any():
        return [None for _ in gene_indices]
    subset = matrix[mask]
    if sp.issparse(subset):
        subset = subset[:, list(gene_indices)]
        counts = np.asarray((subset > 0).sum(axis=0)).ravel()
    else:
        subset = np.asarray(subset)[:, list(gene_indices)]
        counts = np.count_nonzero(subset > 0, axis=0)
    denom = int(mask.sum())
    return [float(v) / denom for v in counts]


def _expression_prefilter_gene_names(
    expression_matrix,
    source_mask: np.ndarray,
    reference_mask: np.ndarray,
    var_names: Sequence[str],
    candidate_gene_names: Sequence[str],
    min_pct_expressed: float,
) -> List[str]:
    """Return candidate genes expressed in at least one side of a contrast."""
    min_pct = _normalize_pct_threshold(min_pct_expressed)
    candidate_names = [str(gene) for gene in candidate_gene_names]
    if min_pct <= 0 or not candidate_names:
        return candidate_names

    gene_to_idx = {str(gene): idx for idx, gene in enumerate(var_names)}
    valid_names: List[str] = []
    valid_indices: List[int] = []
    for gene in candidate_names:
        idx = gene_to_idx.get(gene)
        if idx is None:
            continue
        valid_names.append(gene)
        valid_indices.append(idx)
    if not valid_names:
        return []

    pct_source = _positive_fraction(expression_matrix, source_mask, valid_indices)
    pct_reference = _positive_fraction(expression_matrix, reference_mask, valid_indices)
    keep: List[str] = []
    for gene, source_value, reference_value in zip(valid_names, pct_source, pct_reference):
        source_pct = float(source_value) if source_value is not None and np.isfinite(source_value) else 0.0
        reference_pct = (
            float(reference_value)
            if reference_value is not None and np.isfinite(reference_value)
            else 0.0
        )
        if source_pct >= min_pct or reference_pct >= min_pct:
            keep.append(gene)
    return keep


def _subset_deseq2_dataset(dds: Any, gene_names: Optional[Sequence[str]]) -> Optional[Any]:
    """Subset a fitted PyDESeq2 dataset to contrast-testable genes."""
    if gene_names is None:
        return dds
    requested = [str(gene) for gene in gene_names]
    if not requested:
        return None
    try:
        dds_genes = [str(gene) for gene in dds.var_names]
    except Exception:
        return dds
    available = set(dds_genes)
    retained = [gene for gene in requested if gene in available]
    if not retained:
        return None
    try:
        subset = dds[:, retained].copy()
    except Exception:
        keep = np.asarray([gene in set(retained) for gene in dds_genes], dtype=bool)
        subset = dds[:, keep].copy()

    # AnnData slicing may return a plain AnnData object instead of preserving the
    # DeseqDataSet subclass. DeseqStats needs the fitted dataset methods and
    # runtime attributes, so restore them on the gene subset.
    try:
        if subset.__class__ is not dds.__class__:
            subset.__class__ = dds.__class__
    except Exception:
        pass

    for attr in (
        "refit_cooks",
        "low_memory",
        "formulaic_contrasts",
        "inference",
        "quiet",
        "fit_type",
        "min_replicates",
        "min_disp",
        "max_disp",
        "beta_tol",
        "min_mu",
        "max_iter",
        "n_cpus",
        "design_factors",
        "ref_level",
        "control_genes",
    ):
        if hasattr(dds, attr):
            try:
                setattr(subset, attr, getattr(dds, attr))
            except Exception:
                continue

    try:
        if "non_zero" in subset.var:
            non_zero = subset.var["non_zero"].to_numpy(dtype=bool)
        else:
            non_zero = np.ones(int(subset.n_vars), dtype=bool)
        subset.non_zero_idx = np.arange(int(subset.n_vars))[non_zero]
        subset.non_zero_genes = subset.var_names[non_zero]
    except Exception:
        pass

    try:
        original_new_zeroes = getattr(dds, "new_all_zeroes_genes", pd.Index([]))
        subset_var_names = {str(gene) for gene in subset.var_names}
        subset.new_all_zeroes_genes = pd.Index(
            [str(gene) for gene in original_new_zeroes if str(gene) in subset_var_names]
        )
    except Exception:
        subset.new_all_zeroes_genes = pd.Index([])
    return subset


def _adjust_pvalues(pvalues: np.ndarray, method: str = "fdr_bh") -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    adjusted = np.full(p.shape, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return adjusted
    idx = np.flatnonzero(finite)
    vals = np.clip(p[idx], 0, 1)
    method_norm = str(method or "fdr_bh").strip().lower().replace("-", "_")
    if method_norm in {"none", "raw", "pvalue", "pvalues"}:
        adjusted[idx] = vals
        adjusted[(adjusted == 0) & np.isfinite(adjusted)] = np.nextafter(0.0, 1.0)
        return adjusted
    if method_norm in {"bonferroni", "bonf"}:
        adjusted[idx] = np.clip(vals * len(vals), 0, 1)
        adjusted[(adjusted == 0) & np.isfinite(adjusted)] = np.nextafter(0.0, 1.0)
        return adjusted

    order = np.argsort(vals)
    ranked = vals[order]
    n = len(ranked)
    if method_norm in {"holm", "holm_bonferroni"}:
        adj = (n - np.arange(n, dtype=float)) * ranked
        adj = np.maximum.accumulate(adj)
    else:
        adj = ranked * n / (np.arange(n, dtype=float) + 1.0)
        adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    adjusted[idx[order]] = adj
    adjusted[(adjusted == 0) & np.isfinite(adjusted)] = np.nextafter(0.0, 1.0)
    return adjusted


def _normalize_pct_threshold(value: float) -> float:
    threshold = float(value or 0.0)
    if threshold > 1.0:
        threshold = threshold / 100.0
    return min(max(threshold, 0.0), 1.0)


def _json_float(value: float, digits: int = 6) -> Optional[float]:
    try:
        val = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(val):
        return None
    return float(round(val, digits))


@contextmanager
def _suppress_numpy_slogdet_warnings():
    """Mute expected singular-design warnings emitted inside PyDESeq2 stats."""
    _install_numpy_slogdet_warning_filter()
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=_NUMPY_SLOGDET_WARNING_MESSAGE,
            category=RuntimeWarning,
        )
        yield


def _format_duration(seconds: float) -> str:
    """Format a short, human-readable elapsed or estimated duration."""
    value = max(0.0, float(seconds))
    if value < 90.0:
        return f"{max(1, int(round(value)))} seconds"
    minutes, remainder = divmod(int(round(value)), 60)
    if minutes < 60:
        return f"{minutes} min" if remainder < 30 else f"{minutes} min {remainder} sec"
    hours, minutes = divmod(minutes, 60)
    return f"{hours} h {minutes} min"


def _estimate_shared_deseq2_fit_seconds(
    n_samples: int,
    n_genes: int,
    design_columns: int,
    n_cpus: int = 1,
) -> float:
    """Return a pre-fit runtime estimate.

    This is a workload heuristic, not a benchmark: DESeq2 convergence, BLAS,
    available memory, and CPU speed can substantially change the duration.
    """
    sample_gene_millions = (max(1, int(n_samples)) * max(1, int(n_genes))) / 1_000_000.0
    design_factor = 1.0 + 0.08 * max(0, int(design_columns) - 2)
    workers = max(1, int(n_cpus))
    parallel_speedup = 1.0 + 0.75 * (workers - 1)
    central_seconds = max(60.0, 6.0 * sample_gene_millions * design_factor / parallel_speedup)
    return central_seconds


def _compute_pseudobulk_sample_diagnostics(
    pair_counts: np.ndarray,
    pair_meta: pd.DataFrame,
    *,
    max_genes: int = 1000,
) -> Dict[str, Any]:
    counts = np.array(pair_counts, dtype=float, copy=True)
    counts[~np.isfinite(counts)] = 0
    counts[counts < 0] = 0
    if counts.ndim != 2 or counts.shape[0] < 2 or counts.shape[1] < 1:
        return {}

    library_sizes = counts.sum(axis=1)
    safe_libs = np.where(library_sizes > 0, library_sizes, 1.0)
    log_cpm = np.log1p((counts / safe_libs[:, None]) * 1_000_000.0)

    variances = np.var(log_cpm, axis=0)
    finite_var = np.isfinite(variances) & (variances > 0)
    if not finite_var.any():
        feature_idx = np.arange(log_cpm.shape[1])
    else:
        feature_idx = np.flatnonzero(finite_var)
    if feature_idx.size > int(max_genes):
        order = np.argsort(variances[feature_idx])[::-1][: int(max_genes)]
        feature_idx = feature_idx[order]

    matrix = log_cpm[:, feature_idx]
    matrix = matrix - np.nanmean(matrix, axis=0, keepdims=True)
    matrix[~np.isfinite(matrix)] = 0

    pca_points: List[List[Optional[float]]] = []
    pca_variance: List[Optional[float]] = [None, None]
    try:
        u, s, _vt = np.linalg.svd(matrix, full_matrices=False)
        scores = u[:, :2] * s[:2]
        if scores.shape[1] == 1:
            scores = np.column_stack([scores[:, 0], np.zeros(scores.shape[0])])
        eig = s ** 2
        total = float(eig.sum())
        if total > 0:
            pca_variance = [_json_float(eig[i] / total, 6) if i < eig.size else None for i in range(2)]
        pca_points = [[_json_float(x, 6), _json_float(y, 6)] for x, y in scores[:, :2]]
    except Exception:
        pca_points = [[None, None] for _ in range(counts.shape[0])]

    try:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import pdist, squareform

        condensed = pdist(matrix, metric="euclidean")
        dist = squareform(condensed)
        if condensed.size:
            linkage_matrix = linkage(condensed, method="average")
            order = leaves_list(linkage_matrix).astype(int).tolist()
        else:
            order = list(range(counts.shape[0]))
    except Exception:
        diff = matrix[:, None, :] - matrix[None, :, :]
        dist = np.sqrt(np.sum(diff * diff, axis=2))
        order = list(range(counts.shape[0]))

    ordered_dist = dist[np.ix_(order, order)]
    replicate_values = pair_meta["_pb_replicate"].astype(str).tolist()
    group_values = pair_meta["_pb_group"].astype(str).tolist()
    labels = [f"{replicate_values[i]} | {group_values[i]}" for i in range(len(replicate_values))]
    n_cells = pair_meta["n_cells"].to_numpy(dtype=float) if "n_cells" in pair_meta.columns else np.zeros(len(labels))

    return {
        "labels": labels,
        "replicates": replicate_values,
        "groups": group_values,
        "n_cells": [int(v) if np.isfinite(v) else None for v in n_cells],
        "library_size": [int(v) if np.isfinite(v) else None for v in library_sizes],
        "pca": pca_points,
        "pca_variance": pca_variance,
        "distance_order": [int(i) for i in order],
        "distance_labels": [labels[i] for i in order],
        "distance_groups": [group_values[i] for i in order],
        "distance_matrix": [
            [_json_float(v, 5) for v in row]
            for row in ordered_dist
        ],
        "distance_metric": "euclidean_log1p_cpm",
        "pca_features": int(feature_idx.size),
    }


def _compute_category_gene_means_from_aggregate(
    aggregate: np.ndarray,
    pb_meta: pd.DataFrame,
    categories: Sequence[str],
    gene_names: Sequence[str],
) -> Dict[str, Any]:
    """Summarize category-level means from replicate-level pseudobulk counts."""
    genes = [str(g) for g in gene_names]
    category_means: Dict[str, List[Optional[float]]] = {}
    category_cells: Dict[str, int] = {}
    aggregate = np.asarray(aggregate, dtype=float)
    n_cells = pb_meta["n_cells"].to_numpy(dtype=float) if "n_cells" in pb_meta else np.zeros(aggregate.shape[0])
    n_cells[~np.isfinite(n_cells)] = 0
    n_cells[n_cells < 0] = 0
    valid_rows = n_cells > 0
    per_sample_means = np.zeros_like(aggregate, dtype=float)
    if aggregate.size:
        np.divide(
            aggregate,
            n_cells[:, None],
            out=per_sample_means,
            where=valid_rows[:, None],
        )
    group_values = pb_meta["_pb_group"].astype(str).to_numpy()

    for category in [str(c) for c in categories]:
        mask = (group_values == category) & valid_rows
        cells = float(np.sum(n_cells[mask])) if mask.any() else 0.0
        category_cells[category] = int(cells)
        if mask.any():
            means = per_sample_means[mask].mean(axis=0)
            category_means[category] = [_json_float(v, 6) for v in means]
        else:
            category_means[category] = [0.0 for _ in genes]

    if "_pb_replicate" in pb_meta:
        replicate_values = pb_meta["_pb_replicate"].astype(str).to_numpy()
        replicate_means = []
        for replicate in sorted(set(replicate_values)):
            mask = (replicate_values == replicate) & valid_rows
            cells = float(np.sum(n_cells[mask])) if mask.any() else 0.0
            if cells > 0:
                replicate_means.append(np.asarray(aggregate[mask].sum(axis=0), dtype=float) / cells)
        if replicate_means:
            background = [_json_float(v, 6) for v in np.vstack(replicate_means).mean(axis=0)]
        else:
            background = [0.0 for _ in genes]
    elif valid_rows.any():
        background = [_json_float(v, 6) for v in per_sample_means[valid_rows].mean(axis=0)]
    else:
        background = [0.0 for _ in genes]

    return {
        "genes": genes,
        "categories": [str(c) for c in categories],
        "means": category_means,
        "background": background,
        "n_cells": category_cells,
        "source": "pseudobulk_aggregate",
    }


def _log2fc_from_pseudobulk_means(
    counts: np.ndarray,
    metadata: pd.DataFrame,
    source: str,
    reference: str,
) -> np.ndarray:
    """Compute source/reference log2FC from replicate-level aggregate per-cell means."""
    dense = np.asarray(counts, dtype=float)
    dense[~np.isfinite(dense)] = 0
    dense[dense < 0] = 0
    groups = metadata["_pb_group"].astype(str).to_numpy()
    if "n_cells" in metadata.columns:
        n_cells = metadata["n_cells"].to_numpy(dtype=float)
        n_cells[~np.isfinite(n_cells)] = 0
        n_cells[n_cells < 0] = 0
    else:
        n_cells = np.ones(dense.shape[0], dtype=float)

    source_mask = groups == str(source)
    reference_mask = groups == str(reference)
    n_vars = int(dense.shape[1])
    valid_rows = n_cells > 0
    per_sample_means = np.zeros_like(dense, dtype=float)
    if dense.size:
        np.divide(
            dense,
            n_cells[:, None],
            out=per_sample_means,
            where=valid_rows[:, None],
        )
    source_rows = source_mask & valid_rows
    reference_rows = reference_mask & valid_rows
    source_mean = (
        per_sample_means[source_rows].mean(axis=0)
        if source_rows.any()
        else np.zeros(n_vars, dtype=float)
    )
    reference_mean = (
        per_sample_means[reference_rows].mean(axis=0)
        if reference_rows.any()
        else np.zeros(n_vars, dtype=float)
    )
    tiny = np.nextafter(0.0, 1.0)
    log2fc = np.log2(np.maximum(source_mean, 0.0) + tiny) - np.log2(
        np.maximum(reference_mean, 0.0) + tiny
    )
    log2fc[~np.isfinite(log2fc)] = 0.0
    return np.asarray(log2fc, dtype=float)


def _fit_deseq2_pair(
    counts: np.ndarray,
    metadata: pd.DataFrame,
    source: str,
    reference: str,
    fit_type: str = "parametric",
    min_gene_counts: int = 0,
    n_cpus: int = 1,
) -> pd.DataFrame:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except Exception as exc:  # pragma: no cover - depends on optional runtime import
        raise RuntimeError(
            "PyDESeq2 is required for pseudobulk cluster DE. Install KaroSpace with "
            "pydeseq2 support or run `pip install pydeseq2`."
        ) from exc

    counts, metadata = _filter_pseudobulk_genes(counts, metadata, min_gene_counts)
    gene_names = [str(g) for g in metadata.attrs["gene_names"]]
    if not gene_names:
        return pd.DataFrame(columns=["log2FoldChange", "pvalue", "padj", "stat", "baseMean"])
    counts_df = pd.DataFrame(counts, index=metadata.index, columns=gene_names)
    meta = metadata[["_pb_replicate", "_pb_group"]].copy()
    meta["_pb_replicate"] = pd.Categorical(meta["_pb_replicate"])
    meta["_pb_group"] = pd.Categorical(
        meta["_pb_group"],
        categories=[str(reference), str(source)],
    )

    try:
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=meta,
            design="~ _pb_replicate + _pb_group",
            fit_type=fit_type,
            n_cpus=max(1, int(n_cpus)),
            quiet=True,
        )
    except TypeError:
        dds = DeseqDataSet(
            counts=counts_df,
            clinical=meta,
            design_factors=["_pb_replicate", "_pb_group"],
            fit_type=fit_type,
            refit_cooks=True,
            n_cpus=max(1, int(n_cpus)),
        )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        dds.deseq2()

    with _suppress_numpy_slogdet_warnings():
        try:
            stat_res = DeseqStats(
                dds,
                contrast=["_pb_group", str(source), str(reference)],
                quiet=True,
                n_cpus=max(1, int(n_cpus)),
            )
        except TypeError:
            stat_res = DeseqStats(
                dds,
                contrast=["_pb_group", str(source), str(reference)],
                n_cpus=max(1, int(n_cpus)),
            )
        stat_res.summary()
    result = stat_res.results_df.copy()
    gene_names = [str(g) for g in metadata.attrs["gene_names"]]
    mean_log2fc = _log2fc_from_pseudobulk_means(counts, metadata, source, reference)
    if len(mean_log2fc) == len(gene_names):
        result["log2FoldChange"] = pd.Series(mean_log2fc, index=gene_names).reindex(result.index)
    return result


def _fit_deseq2_shared_categories(
    counts: np.ndarray,
    metadata: pd.DataFrame,
    categories: Sequence[str],
    *,
    fit_type: str = "parametric",
    min_gene_counts: int = 0,
    n_cpus: int = 1,
) -> Tuple[Any, np.ndarray, pd.DataFrame]:
    """Fit one replicate-adjusted DESeq2 model for every retained category.

    The returned fitted dataset is deliberately reused for all pairwise and
    balanced-rest contrasts.  This keeps normalization and dispersion estimates
    common to the annotation column instead of re-estimating them per contrast.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
    except Exception as exc:  # pragma: no cover - depends on optional runtime import
        raise RuntimeError(
            "PyDESeq2 is required for pseudobulk category DE. Install KaroSpace with "
            "pydeseq2 support or run `pip install pydeseq2`."
        ) from exc

    fit_counts, fit_meta = _filter_pseudobulk_genes(counts, metadata, min_gene_counts)
    gene_names = [str(gene) for gene in fit_meta.attrs.get("gene_names", [])]
    if not gene_names:
        raise ValueError("no genes meet the raw pseudobulk count threshold")

    counts_df = pd.DataFrame(fit_counts, index=fit_meta.index, columns=gene_names)
    design_meta = fit_meta[["_pb_replicate", "_pb_group"]].copy()
    design_meta["_pb_replicate"] = pd.Categorical(design_meta["_pb_replicate"].astype(str))
    design_meta["_pb_group"] = pd.Categorical(
        design_meta["_pb_group"].astype(str),
        categories=[str(category) for category in categories],
    )
    try:
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=design_meta,
            design="~ _pb_replicate + _pb_group",
            fit_type=fit_type,
            n_cpus=max(1, int(n_cpus)),
            quiet=True,
        )
    except TypeError:
        dds = DeseqDataSet(
            counts=counts_df,
            clinical=design_meta,
            design_factors=["_pb_replicate", "_pb_group"],
            fit_type=fit_type,
            refit_cooks=True,
            n_cpus=max(1, int(n_cpus)),
        )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        dds.deseq2()
    return dds, fit_counts, fit_meta


def _deseq2_shared_contrast(
    dds: Any,
    contrast: np.ndarray,
    n_cpus: int = 1,
    gene_names: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Evaluate one numeric contrast from an already fitted DESeq2 dataset."""
    try:
        from pydeseq2.ds import DeseqStats
    except Exception as exc:  # pragma: no cover - depends on optional runtime import
        raise RuntimeError(
            "PyDESeq2 is required for pseudobulk category DE. Install KaroSpace with "
            "pydeseq2 support or run `pip install pydeseq2`."
        ) from exc
    vector = np.asarray(contrast, dtype=float)
    stats_dds = _subset_deseq2_dataset(dds, gene_names)
    if stats_dds is None:
        return pd.DataFrame(columns=["log2FoldChange", "pvalue", "padj", "stat", "baseMean"])
    with _suppress_numpy_slogdet_warnings():
        try:
            stats = DeseqStats(stats_dds, contrast=vector, quiet=True, n_cpus=max(1, int(n_cpus)))
        except TypeError:
            stats = DeseqStats(stats_dds, contrast=vector, n_cpus=max(1, int(n_cpus)))
        stats.summary()
    return stats.results_df.copy()


def _reverse_deseq2_result(result: pd.DataFrame) -> pd.DataFrame:
    """Return the opposite direction of a DESeq2 contrast result."""
    reversed_result = result.copy()
    for column in ("log2FoldChange", "stat"):
        if column in reversed_result.columns:
            reversed_result[column] = -pd.to_numeric(reversed_result[column], errors="coerce")
    return reversed_result


def _shared_group_contrast(dds: Any, source: str, reference: str) -> np.ndarray:
    """Return the model-matrix vector for source minus reference."""
    return np.asarray(
        dds.contrast(
            column="_pb_group",
            baseline=str(reference),
            group_to_compare=str(source),
        ),
        dtype=float,
    )


def _shared_category_design_rank(metadata: pd.DataFrame) -> Tuple[int, int]:
    """Return rank and column count for ``~ replicate + annotation``."""
    replicate = pd.get_dummies(metadata["_pb_replicate"].astype(str), drop_first=True, dtype=float)
    group = pd.get_dummies(metadata["_pb_group"].astype(str), drop_first=True, dtype=float)
    design = np.column_stack((np.ones(len(metadata), dtype=float), replicate.to_numpy(), group.to_numpy()))
    return int(np.linalg.matrix_rank(design)), int(design.shape[1])


# Dormant Complex design helpers. They have no supported caller while Complex
# design is in development and are retained only for later reactivation.
def _fit_deseq2_design(
    counts: np.ndarray,
    metadata: pd.DataFrame,
    formula: str,
    contrast: np.ndarray,
    *,
    fit_type: str = "parametric",
    dds: Any = None,
) -> Tuple[pd.DataFrame, Any]:
    """Fit one categorical DESeq2 design and evaluate a coefficient contrast."""
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except Exception as exc:  # pragma: no cover - depends on optional runtime import
        raise RuntimeError(
            "PyDESeq2 is required for complex pseudobulk DE. Install KaroSpace with "
            "pydeseq2 support or run `pip install pydeseq2`."
        ) from exc

    if dds is None:
        gene_names = [str(g) for g in metadata.attrs["gene_names"]]
        counts_df = pd.DataFrame(counts, index=metadata.index, columns=gene_names)
        design_metadata = metadata.drop(columns=["n_cells", "_pb_replicate", "_pb_group"], errors="ignore").copy()
        for column in design_metadata.columns:
            design_metadata[column] = pd.Categorical(design_metadata[column].astype(str))
        try:
            dds = DeseqDataSet(
                counts=counts_df,
                metadata=design_metadata,
                design=str(formula),
                fit_type=fit_type,
                n_cpus=1,
                quiet=True,
            )
        except TypeError:
            raise RuntimeError(
                "This PyDESeq2 version does not support formula-based Complex design fitting. "
                "Upgrade PyDESeq2 to a version with the `design` argument."
            )
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            dds.deseq2()
    with _suppress_numpy_slogdet_warnings():
        try:
            stat_res = DeseqStats(
                dds,
                contrast=np.asarray(contrast, dtype=float),
                quiet=True,
                n_cpus=1,
            )
        except TypeError:
            stat_res = DeseqStats(dds, contrast=np.asarray(contrast, dtype=float), n_cpus=1)
        stat_res.summary()
    return stat_res.results_df.copy(), dds


def _complex_coefficient_profile(
    coefficient: str,
    variables: Sequence[str],
    levels: Dict[str, Sequence[str]],
) -> Tuple[Dict[str, str], List[str]]:
    """Return the model-profile represented by a formulaic categorical coefficient."""
    profile = {
        str(variable): str((levels.get(str(variable)) or [""])[0])
        for variable in variables
    }
    matches = re.findall(r"([A-Za-z_][A-Za-z0-9_]*)\[T\.([^\]]+)\]", str(coefficient))
    changed: List[str] = []
    for variable, level in matches:
        if variable in profile:
            profile[variable] = str(level)
            changed.append(variable)
    return profile, list(dict.fromkeys(changed))


def compute_pseudobulk_complex_design_de(
    adata,
    design_audit: Optional[Dict[str, Any]],
    *,
    replicate: str,
    counts_layer: Optional[str] = "counts",
    min_cell_counts: int = 0,
    min_gene_counts: int = 0,
    min_cells: int = 20,
    min_replicates: int = 2,
    min_pct_expressed: float = 0.0,
    p_adjust_method: str = "fdr_bh",
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 0.5,
    fit_type: str = "parametric",
    max_coefficients: int = 32,
) -> Optional[Dict[str, Any]]:
    """Fit a formula-based categorical pseudobulk design and retain coefficient views.

    Each stored result is one estimable DESeq2 coefficient versus the reference
    level(s). This keeps a self-contained HTML export practical while supporting
    main effects and arbitrary-order categorical interaction terms.
    """
    if not isinstance(design_audit, dict) or not design_audit.get("valid"):
        return None
    formula = str(design_audit.get("formula") or "").strip()
    variables = [str(v) for v in (design_audit.get("variables") or [])]
    levels = {
        str(key): [str(value) for value in values]
        for key, values in (design_audit.get("levels") or {}).items()
    }
    sample_variables = [str(v) for v in (design_audit.get("sample_variables") or variables)]
    if not formula or not variables or not sample_variables:
        return None
    if replicate not in adata.obs.columns:
        print(f"  Warning: complex pseudobulk replicate '{replicate}' not found in obs.", flush=True)
        return None

    expression_matrix, counts_layer_used, warning = _as_count_matrix(adata, counts_layer)
    model_obs = adata.obs[sample_variables].copy()
    valid = model_obs.notna().all(axis=1).to_numpy(dtype=bool)
    valid &= _cell_count_mask(expression_matrix, min_cell_counts)
    if not valid.any():
        return None
    for variable in sample_variables:
        model_obs[variable] = model_obs[variable].astype(str)
    valid_indices = np.flatnonzero(valid)
    grouped_rows = model_obs.iloc[valid_indices].groupby(sample_variables, observed=True, sort=True)
    sample_keys = list(grouped_rows.groups.keys())
    if not sample_keys:
        return None
    key_to_row = {key if isinstance(key, tuple) else (key,): idx for idx, key in enumerate(sample_keys)}
    rows = np.empty(valid_indices.size, dtype=np.int64)
    for i, row in enumerate(model_obs.iloc[valid_indices][sample_variables].itertuples(index=False, name=None)):
        rows[i] = key_to_row[tuple(str(value) for value in row)]

    incidence = sp.csr_matrix(
        (np.ones(valid_indices.size, dtype=np.float64), (rows, valid_indices)),
        shape=(len(sample_keys), adata.n_obs),
    )
    aggregate = _to_dense_counts(incidence @ expression_matrix)
    cell_counts = np.bincount(rows, minlength=len(sample_keys)).astype(int)
    keep = cell_counts >= int(min_cells)
    if not keep.any():
        print("  Warning: complex pseudobulk DE was skipped: no pseudobulk sample meets the cell threshold.", flush=True)
        return None
    aggregate = aggregate[keep]
    retained_keys = [sample_keys[index] for index in np.flatnonzero(keep)]
    cell_counts = cell_counts[keep]
    pb_meta = pd.DataFrame(
        [dict(zip(sample_variables, key if isinstance(key, tuple) else (key,))) for key in retained_keys],
        index=[f"pb_complex_{index}" for index in range(len(retained_keys))],
    )
    for variable in variables:
        pb_meta[variable] = pd.Categorical(pb_meta[variable].astype(str), categories=levels.get(variable))
    pb_meta["n_cells"] = cell_counts
    pb_meta["_pb_replicate"] = pb_meta[replicate].astype(str)
    primary_variable = variables[-1]
    pb_meta["_pb_group"] = pb_meta[primary_variable].astype(str)
    pb_meta.attrs["gene_names"] = [str(g) for g in adata.var_names]
    aggregate, pb_meta = _filter_pseudobulk_genes(
        aggregate,
        pb_meta,
        min_gene_counts,
    )
    if aggregate.shape[1] == 0:
        print(
            "  Warning: complex pseudobulk DE was skipped: no genes meet the raw pseudobulk count threshold.",
            flush=True,
        )
        return None

    required_replicates = max(2, int(min_replicates))
    n_replicates = int(pb_meta["_pb_replicate"].nunique())
    if n_replicates < required_replicates:
        print(
            "  Warning: complex pseudobulk DE was skipped: "
            f"{n_replicates} biological replicate(s) are available; need >= {required_replicates}.",
            flush=True,
        )
        return None

    try:
        from formulaic import model_matrix

        formula_matrix = model_matrix(formula, pb_meta[variables], output="pandas")
        rank = int(np.linalg.matrix_rank(np.asarray(formula_matrix, dtype=float)))
        if rank < int(formula_matrix.shape[1]) or len(pb_meta) <= rank:
            print(
                "  Warning: complex pseudobulk DE was skipped: the filtered pseudobulk design "
                "is rank-deficient or has no residual degrees of freedom.",
                flush=True,
            )
            return None
        coefficient_names = [name for name in formula_matrix.columns if str(name).lower() != "intercept"]
    except Exception as exc:
        print(f"  Warning: complex pseudobulk DE was skipped: could not build design matrix ({exc}).", flush=True)
        return None

    if not coefficient_names:
        return None
    sample_diagnostics = _compute_pseudobulk_sample_diagnostics(aggregate, pb_meta)
    original_values = adata.obs[variables].copy()
    for variable in variables:
        original_values[variable] = original_values[variable].astype(str)
    valid_model_cells = valid.copy()
    baseline_profile = {variable: str((levels.get(variable) or [""])[0]) for variable in variables}
    results: List[Dict[str, Any]] = []
    fitted_dds = None
    skipped = max(0, len(coefficient_names) - int(max_coefficients))
    print(
        f"    - complex pseudobulk DE: fitting {len(pb_meta)} samples, {aggregate.shape[1]} genes, "
        f"{min(len(coefficient_names), int(max_coefficients))} coefficient contrast(s)",
        flush=True,
    )
    for coefficient in coefficient_names[: int(max_coefficients)]:
        contrast = np.zeros(len(formula_matrix.columns), dtype=float)
        contrast[list(formula_matrix.columns).index(coefficient)] = 1.0
        source_profile, changed_variables = _complex_coefficient_profile(coefficient, variables, levels)
        source_mask = valid_model_cells.copy()
        reference_mask = valid_model_cells.copy()
        for variable in variables:
            source_mask &= original_values[variable].to_numpy() == source_profile[variable]
            reference_mask &= original_values[variable].to_numpy() == baseline_profile[variable]
        n_source = int(np.count_nonzero(source_mask))
        n_reference = int(np.count_nonzero(reference_mask))
        try:
            fitted, fitted_dds = _fit_deseq2_design(
                aggregate,
                pb_meta,
                formula,
                contrast,
                fit_type=str(fit_type or "parametric"),
                dds=fitted_dds,
            )
            formatted = _format_result(
                fitted,
                source_mask=source_mask,
                reference_mask=reference_mask,
                expression_matrix=expression_matrix,
                var_names=adata.var_names,
                top_n=0,
                n_source=n_source,
                n_reference=n_reference,
                n_replicates=n_replicates,
                counts_layer_used=counts_layer_used,
                warning=warning,
                p_adjust_method=p_adjust_method,
                min_pct_expressed=min_pct_expressed,
                padj_cutoff=padj_cutoff,
                log2fc_cutoff=log2fc_cutoff,
                sample_diagnostics=sample_diagnostics,
                method="pseudobulk-deseq2-complex-design",
            )
            results.append({
                "id": f"coefficient-{len(results)}",
                "coefficient": str(coefficient),
                "variables": changed_variables,
                "primary_variable": changed_variables[-1] if changed_variables else primary_variable,
                "source_profile": source_profile,
                "reference_profile": baseline_profile,
                "result": formatted,
            })
        except Exception as exc:
            print(f"      - coefficient '{coefficient}': failed ({exc})", flush=True)
            results.append({
                "id": f"coefficient-{len(results)}",
                "coefficient": str(coefficient),
                "variables": changed_variables,
                "primary_variable": changed_variables[-1] if changed_variables else primary_variable,
                "source_profile": source_profile,
                "reference_profile": baseline_profile,
                "result": _empty_result(
                    "de_failed",
                    n_source=n_source,
                    n_reference=n_reference,
                    min_cells=int(min_cells),
                    min_replicates=required_replicates,
                    details=str(exc),
                ),
            })

    warnings: List[str] = []
    if skipped:
        warnings.append(
            f"Only the first {int(max_coefficients)} model coefficients were exported; {skipped} additional coefficient(s) were omitted."
        )
    return {
        "available": bool(results),
        "formula": formula,
        "baseline_profile": baseline_profile,
        "sample_count": int(len(pb_meta)),
        "design_rank": rank,
        "design_columns": int(formula_matrix.shape[1]),
        "coefficients": results,
        "warnings": warnings,
    }


def _empty_result(
    reason: str,
    *,
    n_source: int,
    n_reference: int,
    min_cells: int,
    min_replicates: int,
    details: Optional[str] = None,
) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        "available": False,
        "reason": reason,
        "method": "pseudobulk-deseq2",
        "genes": [],
        "logfoldchanges": [],
        "pvals": [],
        "pvals_adj": [],
        "scores": [],
        "pct_source": [],
        "pct_reference": [],
        "base_mean": [],
        "n_source": int(n_source),
        "n_reference": int(n_reference),
        "min_cells_required": int(min_cells),
        "min_replicates_required": int(min_replicates),
    }
    if details:
        result["details"] = details
    return result


def _empty_interaction_result(
    reason: str,
    *,
    n_contact: int,
    n_non_contact: int,
    min_cells: int,
    min_replicates: int,
    details: Optional[str] = None,
    pct_contact: float = 0.0,
    mean_target_neighbors_contact: float = 0.0,
    mean_target_neighbors_non_contact: float = 0.0,
    target_edge_count: float = 0.0,
    target_zscore: Optional[float] = None,
) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        "available": False,
        "reason": reason,
        "method": "pseudobulk-deseq2-contact",
        "genes": [],
        "logfoldchanges": [],
        "log2foldchanges": [],
        "pvals": [],
        "pvals_adj": [],
        "scores": [],
        "pct_source": [],
        "pct_reference": [],
        "base_mean": [],
        "n_contact": int(n_contact),
        "n_non_contact": int(n_non_contact),
        "min_cells_required": int(min_cells),
        "min_replicates_required": int(min_replicates),
        "pct_contact": float(pct_contact),
        "mean_target_neighbors_contact": float(mean_target_neighbors_contact),
        "mean_target_neighbors_non_contact": float(mean_target_neighbors_non_contact),
        "target_edge_count": float(target_edge_count),
        "target_zscore": target_zscore,
    }
    if details:
        result["details"] = details
    return result


def _format_result(
    df: pd.DataFrame,
    *,
    source_mask: np.ndarray,
    reference_mask: np.ndarray,
    expression_matrix,
    var_names: Sequence[str],
    top_n: int,
    n_source: int,
    n_reference: int,
    n_replicates: int,
    counts_layer_used: str,
    warning: Optional[str],
    p_adjust_method: str,
    min_pct_expressed: float,
    padj_cutoff: float,
    log2fc_cutoff: float,
    sample_diagnostics: Optional[Dict[str, Any]] = None,
    method: str = "pseudobulk-deseq2",
) -> Dict[str, Any]:
    p_adjust_method = str(p_adjust_method or "fdr_bh").strip().lower().replace("-", "_")
    min_pct_expressed = _normalize_pct_threshold(min_pct_expressed)
    padj_cutoff = min(max(float(padj_cutoff), 0.0), 1.0)
    log2fc_cutoff = max(float(log2fc_cutoff), 0.0)
    base_payload = {
        "method": str(method or "pseudobulk-deseq2"),
        "p_adjust_method": p_adjust_method,
        "min_pct_expressed": min_pct_expressed,
        "padj_cutoff": padj_cutoff,
        "log2fc_cutoff": log2fc_cutoff,
        "table_top_n": int(top_n),
        **({"pseudobulk_samples": sample_diagnostics} if sample_diagnostics else {}),
        "n_source": int(n_source),
        "n_reference": int(n_reference),
        "n_replicates": int(n_replicates),
        "counts_layer": counts_layer_used,
        **({"warning": warning} if warning else {}),
    }
    if df is None or df.empty:
        result = _empty_result("no_results", n_source=n_source, n_reference=n_reference, min_cells=0, min_replicates=0)
        result.update(base_payload)
        return result

    work = df.copy()
    if "log2FoldChange" not in work.columns:
        work["log2FoldChange"] = np.nan
    if "pvalue" not in work.columns:
        work["pvalue"] = np.nan
    work["padj"] = _adjust_pvalues(work["pvalue"].to_numpy(dtype=float), method=p_adjust_method)
    if "stat" not in work.columns:
        work["stat"] = np.nan
    if "baseMean" not in work.columns:
        work["baseMean"] = np.nan

    work["_gene"] = [str(idx) for idx in work.index]
    work = work[np.isfinite(work["log2FoldChange"].to_numpy(dtype=float))]
    if work.empty:
        return {
            "available": True,
            "genes": [],
            "log2foldchanges": [],
            "logfoldchanges": [],
            "pvals": [],
            "pvals_adj": [],
            "scores": [],
            "pct_source": [],
            "pct_reference": [],
            "base_mean": [],
            **base_payload,
        }

    gene_to_idx = {str(g): i for i, g in enumerate(var_names)}
    gene_indices = [gene_to_idx.get(g) for g in work["_gene"]]
    valid_positions = [i for i, idx in enumerate(gene_indices) if idx is not None]
    if len(valid_positions) != len(gene_indices):
        work = work.iloc[valid_positions]
        gene_indices = [gene_indices[i] for i in valid_positions]

    pct_source = _positive_fraction(expression_matrix, source_mask, gene_indices)
    pct_reference = _positive_fraction(expression_matrix, reference_mask, gene_indices)
    work["_pct_source"] = pct_source
    work["_pct_reference"] = pct_reference

    work["_padj_sort"] = work["padj"].fillna(np.inf)
    work["_pvalue_sort"] = work["pvalue"].fillna(np.inf)
    work["_abs_lfc"] = np.abs(work["log2FoldChange"].to_numpy(dtype=float))
    work = work.sort_values(
        ["_padj_sort", "_pvalue_sort", "_abs_lfc", "_gene"],
        ascending=[True, True, False, True],
    )
    pct_source = [work["_pct_source"].iloc[i] for i in range(len(work))]
    pct_reference = [work["_pct_reference"].iloc[i] for i in range(len(work))]

    def _finite_or_none(value):
        try:
            val = float(value)
        except (TypeError, ValueError):
            return None
        return val if np.isfinite(val) else None

    result = {
        "available": True,
        "genes": work["_gene"].astype(str).tolist(),
        "log2foldchanges": [_finite_or_none(v) for v in work["log2FoldChange"]],
        "logfoldchanges": [_finite_or_none(v) for v in work["log2FoldChange"]],
        "pvals": [_finite_or_none(v) for v in work["pvalue"]],
        "pvals_adj": [_finite_or_none(v) for v in work["padj"]],
        "scores": [_finite_or_none(v) for v in work["stat"]],
        "pct_source": pct_source,
        "pct_reference": pct_reference,
        "base_mean": [_finite_or_none(v) for v in work["baseMean"]],
        **base_payload,
    }
    return result


def _count_threshold_passing_genes(result: Dict[str, Any]) -> int:
    """Count DE genes satisfying the result's display/embed thresholds."""
    if not result or result.get("available") is False:
        return 0
    log2fc = result.get("log2foldchanges") or result.get("logfoldchanges") or []
    padj = result.get("pvals_adj") or []
    try:
        padj_cutoff = float(result.get("padj_cutoff", 0.05))
        log2fc_cutoff = float(result.get("log2fc_cutoff", 0.5))
    except (TypeError, ValueError):
        return 0
    count = 0
    for value, adjusted in zip(log2fc, padj):
        try:
            if float(adjusted) < padj_cutoff and abs(float(value)) >= log2fc_cutoff:
                count += 1
        except (TypeError, ValueError):
            continue
    return count


def _truncate_gene_result(result: Dict[str, Any], top_n: int) -> Dict[str, Any]:
    if int(top_n) < 1:
        return result
    n = min(int(top_n), len(result.get("genes") or []))
    fields = [
        "genes",
        "log2foldchanges",
        "logfoldchanges",
        "pvals",
        "pvals_adj",
        "scores",
        "pct_source",
        "pct_reference",
        "base_mean",
    ]
    for field in fields:
        if isinstance(result.get(field), list):
            result[field] = result[field][:n]
    result["table_top_n"] = int(top_n)
    return result


def _target_zscore_value(zscore: Optional[np.ndarray], source_idx: int, target_idx: int) -> Optional[float]:
    if not isinstance(zscore, np.ndarray):
        return None
    if source_idx >= zscore.shape[0] or target_idx >= zscore.shape[1]:
        return None
    value = float(zscore[source_idx, target_idx])
    return value if np.isfinite(value) else None


def _interaction_meta(
    *,
    target_neighbor_counts: np.ndarray,
    pos_mask: np.ndarray,
    neg_mask: np.ndarray,
    n_contact: int,
    n_non_contact: int,
    edge_count: float,
    target_zscore: Optional[float],
) -> Dict[str, Any]:
    return {
        "pct_contact": float((100.0 * n_contact) / max(1, n_contact + n_non_contact)),
        "mean_target_neighbors_contact": float(
            np.mean(target_neighbor_counts[pos_mask]) if n_contact > 0 else 0.0
        ),
        "mean_target_neighbors_non_contact": float(
            np.mean(target_neighbor_counts[neg_mask]) if n_non_contact > 0 else 0.0
        ),
        "target_edge_count": float(edge_count),
        "target_zscore": target_zscore,
    }


def compute_pseudobulk_interaction_markers(
    adata,
    annotation: str,
    *,
    replicate: str,
    graph,
    obs_idx: Sequence[int],
    labels: Sequence[int],
    categories: Sequence[str],
    neighbor_counts: np.ndarray,
    neighbor_zscore: Optional[np.ndarray] = None,
    neighbor_n_cells: Optional[Sequence[int]] = None,
    counts_layer: Optional[str] = "counts",
    min_cell_counts: int = 0,
    min_gene_counts: int = 0,
    top_targets: int = 8,
    top_genes: int = 20,
    min_cells: int = 30,
    min_neighbors: int = 1,
    min_replicates: int = 2,
    min_pct_expressed: float = 0.0,
    p_adjust_method: str = "fdr_bh",
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 0.5,
    fit_type: str = "parametric",
    n_cpus: int = 1,
) -> Optional[Dict[str, Dict[str, Dict[str, Any]]]]:
    """Compute replicate-aware contact-conditioned markers for one annotation.

    For each directional source -> target pair, source cells are split into
    contact+ and contact- within each replicate using the neighbor graph
    restricted to that replicate. Counts are then aggregated by
    replicate/contact status and tested with the same paired DESeq2 design used
    for category pseudobulk DE.
    """
    if replicate not in adata.obs.columns:
        log_warning(f"interaction pseudobulk replicate '{replicate}' not found in obs.", level=2)
        return None

    graph = graph.tocsr() if sp.issparse(graph) else sp.csr_matrix(graph)
    contact_graph = graph.copy()
    contact_graph.data = np.ones(contact_graph.nnz, dtype=np.float32)
    contact_graph.eliminate_zeros()
    obs_idx = np.asarray(obs_idx, dtype=np.int64)
    labels = np.asarray(labels, dtype=np.int32)
    categories = [str(cat) for cat in categories]
    if sp.issparse(neighbor_counts):
        counts = np.asarray(neighbor_counts.toarray(), dtype=float)
    else:
        counts = np.asarray(neighbor_counts, dtype=float)
    zscore = np.asarray(neighbor_zscore, dtype=float) if neighbor_zscore is not None else None
    if neighbor_n_cells is not None:
        n_cells = np.asarray(neighbor_n_cells, dtype=int)
    else:
        n_cells = np.bincount(labels[labels >= 0], minlength=len(categories)).astype(int)
    if graph.shape[0] != len(labels) or len(labels) != len(obs_idx):
        log_warning(f"interaction markers '{annotation}' graph/label size mismatch; skipping.", level=2)
        return None

    expression_matrix, counts_layer_used, warning = _as_count_matrix(adata, counts_layer)
    rep_values = adata.obs[replicate].astype(str).to_numpy()
    ctx_reps = rep_values[obs_idx]
    valid_reps = np.asarray(pd.notna(ctx_reps), dtype=bool)
    valid_reps &= _cell_count_mask(expression_matrix, min_cell_counts)[obs_idx]
    if not valid_reps.any():
        return None

    rep_to_positions: Dict[str, np.ndarray] = {
        str(rep): np.flatnonzero((ctx_reps == rep) & valid_reps)
        for rep in sorted(set(ctx_reps[valid_reps].astype(str)))
    }
    top_targets = int(top_targets)
    top_genes = int(top_genes)
    min_cells = int(min_cells)
    min_neighbors = int(min_neighbors)
    min_replicates = int(min_replicates)
    required_min_replicates = max(2, min_replicates)
    if len(rep_to_positions) < required_min_replicates:
        if len(rep_to_positions) == 1:
            reason = "only one biological replicate is available; pseudobulk analysis requires at least two"
        else:
            reason = (
                f"{len(rep_to_positions)} biological replicate(s) are available; "
                f"need >= {required_min_replicates}"
            )
        log_warning(
            f"pseudobulk interaction markers for '{annotation}' were skipped: {reason} "
            f"(replicate annotation: '{replicate}').",
            level=2,
        )
        return None
    rep_subgraph_cache: Dict[str, Any] = {
        rep: contact_graph[positions][:, positions]
        for rep, positions in rep_to_positions.items()
        if positions.size > 0
    }
    log_detail(
        f"{len(categories)} categories -> contact pseudobulk by {replicate}; "
        f"testing up to top {top_targets} target categories per source.",
        level=2,
    )

    results: Dict[str, Dict[str, Dict[str, Any]]] = {}
    target_neighbor_counts_cache: Dict[int, np.ndarray] = {}
    for source_idx, source_name in enumerate(categories):
        if source_idx >= len(n_cells) or int(n_cells[source_idx]) <= 0:
            continue
        source_mask = (labels == source_idx) & valid_reps
        if not source_mask.any():
            continue
        row = counts[source_idx] if source_idx < counts.shape[0] else None
        if row is None:
            continue

        candidate_targets = [
            t for t in range(len(categories))
            if t != source_idx and t < len(row) and float(row[t]) > 0
        ]
        if not candidate_targets:
            continue

        def _target_sort_key(tidx: int):
            zval = _target_zscore_value(zscore, source_idx, tidx)
            if zval is not None:
                return (0, -float(zval), -float(row[tidx]), categories[tidx])
            return (1, 0.0, -float(row[tidx]), categories[tidx])

        ranked_targets = sorted(candidate_targets, key=_target_sort_key)[:top_targets]
        source_result: Dict[str, Dict[str, Any]] = {}

        for target_idx in ranked_targets:
            target_name = categories[target_idx]
            target_neighbor_counts = target_neighbor_counts_cache.get(target_idx)
            if target_neighbor_counts is None:
                target_neighbor_counts = np.zeros(len(labels), dtype=float)
                target_mask = labels == target_idx

                # Classify contact status within each replicate using the induced
                # replicate subgraph. This prevents cross-section/mouse edges in a
                # global graph from defining contact status.
                for rep, positions in rep_to_positions.items():
                    target_vec = target_mask[positions].astype(np.float32, copy=False)
                    if not np.any(target_vec):
                        continue
                    subgraph = rep_subgraph_cache.get(rep)
                    if subgraph is None:
                        continue
                    target_neighbor_counts[positions] = np.asarray(subgraph.dot(target_vec)).ravel()
                target_neighbor_counts_cache[target_idx] = target_neighbor_counts
            target_mask = labels == target_idx

            pos_mask = source_mask & (target_neighbor_counts >= min_neighbors)
            neg_mask = source_mask & (target_neighbor_counts == 0)
            n_pos = int(np.count_nonzero(pos_mask))
            n_neg = int(np.count_nonzero(neg_mask))
            target_z = _target_zscore_value(zscore, source_idx, target_idx)
            meta = _interaction_meta(
                target_neighbor_counts=target_neighbor_counts,
                pos_mask=pos_mask,
                neg_mask=neg_mask,
                n_contact=n_pos,
                n_non_contact=n_neg,
                edge_count=float(row[target_idx]),
                target_zscore=target_z,
            )
            source_cell_mask = np.zeros(adata.n_obs, dtype=bool)
            reference_cell_mask = np.zeros(adata.n_obs, dtype=bool)
            source_cell_mask[obs_idx[pos_mask]] = True
            reference_cell_mask[obs_idx[neg_mask]] = True

            sample_keys: List[Tuple[str, str]] = []
            sample_index: Dict[Tuple[str, str], int] = {}
            selected_positions: List[int] = []
            selected_rows: List[int] = []
            for rep, positions in rep_to_positions.items():
                rep_pos = positions[pos_mask[positions]]
                rep_neg = positions[neg_mask[positions]]
                for status, status_positions in (("contact+", rep_pos), ("contact-", rep_neg)):
                    if status_positions.size < min_cells:
                        continue
                    key = (str(rep), status)
                    row_idx = sample_index.get(key)
                    if row_idx is None:
                        row_idx = len(sample_keys)
                        sample_index[key] = row_idx
                        sample_keys.append(key)
                    selected_positions.extend(status_positions.astype(int).tolist())
                    selected_rows.extend([row_idx] * int(status_positions.size))

            contact_reps = {rep for rep, status in sample_keys if status == "contact+"}
            non_contact_reps = {rep for rep, status in sample_keys if status == "contact-"}
            paired_reps = sorted(contact_reps & non_contact_reps)
            if len(paired_reps) < required_min_replicates:
                log_step(
                    f"{source_name} -> {target_name}: skipped, insufficient paired replicates "
                    f"({len(paired_reps)}; need >= {required_min_replicates})",
                    level=3,
                )
                source_result[target_name] = _empty_interaction_result(
                    "insufficient_replicates",
                    n_contact=n_pos,
                    n_non_contact=n_neg,
                    min_cells=min_cells,
                    min_replicates=required_min_replicates,
                    details=f"{len(paired_reps)} paired replicate(s) available",
                    **meta,
                )
                continue

            kept_sample_keys = [
                key
                for key in sample_keys
                if key[0] in paired_reps and key[1] in {"contact+", "contact-"}
            ]
            keep_sample = {key: idx for idx, key in enumerate(kept_sample_keys)}
            remapped_rows = []
            kept_obs_positions = []
            for position, row_idx in zip(selected_positions, selected_rows):
                key = sample_keys[row_idx]
                new_idx = keep_sample.get(key)
                if new_idx is None:
                    continue
                kept_obs_positions.append(int(obs_idx[position]))
                remapped_rows.append(int(new_idx))

            if not kept_obs_positions:
                source_result[target_name] = _empty_interaction_result(
                    "insufficient_cells",
                    n_contact=n_pos,
                    n_non_contact=n_neg,
                    min_cells=min_cells,
                    min_replicates=required_min_replicates,
                    **meta,
                )
                continue

            incidence = sp.csr_matrix(
                (
                    np.ones(len(kept_obs_positions), dtype=np.float64),
                    (np.asarray(remapped_rows, dtype=np.int64), np.asarray(kept_obs_positions, dtype=np.int64)),
                ),
                shape=(len(kept_sample_keys), adata.n_obs),
            )
            pair_counts = _to_dense_counts(incidence @ expression_matrix)
            cell_counts = np.bincount(np.asarray(remapped_rows, dtype=np.int64), minlength=len(keep_sample)).astype(int)
            ordered_keys = kept_sample_keys
            pair_meta = pd.DataFrame(
                {
                    "_pb_replicate": [key[0] for key in ordered_keys],
                    "_pb_group": [key[1] for key in ordered_keys],
                    "n_cells": cell_counts,
                },
                index=[f"pb_{i}" for i in range(len(ordered_keys))],
            )
            pair_meta.attrs["gene_names"] = [str(g) for g in adata.var_names]

            try:
                log_step(
                    f"{source_name} -> {target_name}: fitting contact+ vs contact- DESeq2 "
                    f"({len(paired_reps)} paired replicate"
                    f"{'s' if len(paired_reps) != 1 else ''}, "
                    f"{pair_counts.shape[0]} pseudobulk samples, "
                    f"{pair_counts.shape[1]} genes)",
                    level=3,
                )
                sample_diagnostics = _compute_pseudobulk_sample_diagnostics(pair_counts, pair_meta)
                pair_result = _fit_deseq2_pair(
                    pair_counts,
                    pair_meta,
                    "contact+",
                    "contact-",
                    fit_type=str(fit_type or "parametric"),
                    min_gene_counts=min_gene_counts,
                    n_cpus=n_cpus,
                )
                formatted = _format_result(
                    pair_result,
                    source_mask=source_cell_mask,
                    reference_mask=reference_cell_mask,
                    expression_matrix=expression_matrix,
                    var_names=adata.var_names,
                    top_n=top_genes,
                    n_source=n_pos,
                    n_reference=n_neg,
                    n_replicates=len(paired_reps),
                    counts_layer_used=counts_layer_used,
                    warning=warning,
                    p_adjust_method=p_adjust_method,
                    min_pct_expressed=min_pct_expressed,
                    padj_cutoff=padj_cutoff,
                    log2fc_cutoff=log2fc_cutoff,
                    sample_diagnostics=sample_diagnostics,
                )
                formatted["method"] = "pseudobulk-deseq2-contact"
                formatted["n_contact"] = int(n_pos)
                formatted["n_non_contact"] = int(n_neg)
                formatted["min_cells_required"] = int(min_cells)
                formatted["min_replicates_required"] = required_min_replicates
                formatted.update(meta)
                source_result[target_name] = _truncate_gene_result(formatted, top_genes)
                log_detail(
                    f"Stored top {min(top_genes, len(formatted.get('genes') or []))} marker genes "
                    "plus diagnostics for this source-target interaction.",
                    level=4,
                )
            except Exception as exc:
                log_step(f"{source_name} -> {target_name}: failed ({exc})", level=3)
                source_result[target_name] = _empty_interaction_result(
                    "de_failed",
                    n_contact=n_pos,
                    n_non_contact=n_neg,
                    min_cells=min_cells,
                    min_replicates=required_min_replicates,
                    details=str(exc),
                    **meta,
                )

        if source_result:
            results[source_name] = source_result

    if results:
        n_sources = len(results)
        n_targets = sum(len(targets) for targets in results.values() if isinstance(targets, dict))
        log_detail(
            f"Stored interaction-marker payload: {n_sources} source categories, "
            f"{n_targets} source-target tests.",
            level=2,
        )
    return results or None


def _compute_pseudobulk_group_de_shared(
    adata,
    groupby: str,
    *,
    replicate: str,
    pairwise_categories: Optional[Sequence[str]] = None,
    counts_layer: Optional[str] = "counts",
    min_cell_counts: int = 0,
    min_gene_counts: int = 0,
    min_cells: int = 20,
    min_replicates: int = 2,
    min_pct_expressed: float = 0.0,
    p_adjust_method: str = "fdr_bh",
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 0.5,
    fit_type: str = "parametric",
    n_cpus: int = 1,
) -> Optional[Dict[str, Dict[str, Dict[str, Any]]]]:
    """Fit one all-category model, then derive pairwise and balanced-rest DE.

    A category-versus-category test is the appropriate direct contrast from the
    shared ``~ replicate + annotation`` fit. A balanced-rest result is the
    source category minus the equally weighted mean of all
    other retained category coefficients. It is not a pooled-cell rest.
    """
    if groupby not in adata.obs.columns:
        log_warning(f"pseudobulk DE annotation '{groupby}' not found in obs.", level=2)
        return None
    if replicate not in adata.obs.columns:
        log_warning(f"pseudobulk DE replicate '{replicate}' not found in obs.", level=2)
        return None

    col = adata.obs[groupby]
    if pd.api.types.is_numeric_dtype(col):
        log_warning(f"pseudobulk DE annotation '{groupby}' is numeric; skipping.", level=2)
        return None
    if not isinstance(col.dtype, CategoricalDtype):
        col = col.astype("category")
    categories = [str(category) for category in col.cat.categories]
    if len(categories) < 2:
        return None
    if pairwise_categories is None:
        selected_pairwise_categories = list(categories)
    else:
        requested = [str(category) for category in pairwise_categories]
        selected_pairwise_categories = [category for category in categories if category in requested]
        missing = sorted(set(requested) - set(selected_pairwise_categories))
        if missing:
            log_warning(
                "requested Simple design pairwise categories not found in "
                f"'{groupby}': {', '.join(missing)}.",
                level=2,
            )

    expression_matrix, counts_layer_used, warning = _as_count_matrix(adata, counts_layer)
    rep_values = adata.obs[replicate].astype(str)
    group_values = col.astype(str)
    valid = rep_values.notna().to_numpy() & group_values.notna().to_numpy()
    valid &= np.asarray(col.cat.codes.to_numpy() >= 0, dtype=bool)
    valid &= _cell_count_mask(expression_matrix, min_cell_counts)
    if not valid.any():
        return None

    required_min_replicates = max(2, int(min_replicates))
    valid_indices = np.flatnonzero(valid)
    rep_valid = rep_values.iloc[valid_indices].to_numpy(dtype=str)
    group_valid = group_values.iloc[valid_indices].to_numpy(dtype=str)
    if len(set(rep_valid)) < required_min_replicates:
        available = len(set(rep_valid))
        reason = (
            "only one biological replicate is available; pseudobulk DE requires at least two"
            if available == 1
            else f"{available} biological replicate(s) are available; need >= {required_min_replicates}"
        )
        log_warning(
            f"pseudobulk DE for '{groupby}' was skipped: {reason} "
            f"(replicate annotation: '{replicate}').",
            level=2,
        )
        return None

    sample_keys: List[Tuple[str, str]] = []
    sample_index: Dict[Tuple[str, str], int] = {}
    rows = np.empty(valid_indices.size, dtype=np.int64)
    for index, (rep_value, group_value) in enumerate(zip(rep_valid, group_valid)):
        key = (str(rep_value), str(group_value))
        row = sample_index.get(key)
        if row is None:
            row = len(sample_keys)
            sample_index[key] = row
            sample_keys.append(key)
        rows[index] = row

    incidence = sp.csr_matrix(
        (np.ones(valid_indices.size, dtype=np.float64), (rows, valid_indices)),
        shape=(len(sample_keys), adata.n_obs),
    )
    aggregate = _to_dense_counts(incidence @ expression_matrix)
    cell_counts = np.bincount(rows, minlength=len(sample_keys)).astype(int)
    pb_meta = pd.DataFrame(
        {
            "_pb_replicate": [key[0] for key in sample_keys],
            "_pb_group": [key[1] for key in sample_keys],
            "n_cells": cell_counts,
        },
        index=[f"pb_{index}" for index in range(len(sample_keys))],
    )
    pb_meta.attrs["gene_names"] = [str(gene) for gene in adata.var_names]
    aggregate_summary = _compute_category_gene_means_from_aggregate(
        aggregate, pb_meta, categories, adata.var_names,
    )
    log_detail(
        f"Aggregated {len(pb_meta)} replicate x annotation pseudobulk samples from "
        f"{len(valid_indices)} retained cells; output is the count matrix used for DESeq2.",
        level=2,
    )

    sufficient_pb = pb_meta[pb_meta["n_cells"] >= int(min_cells)]
    retained_categories = [
        category
        for category in categories
        if sufficient_pb.loc[sufficient_pb["_pb_group"] == category, "_pb_replicate"].nunique()
        >= required_min_replicates
    ]
    omitted_categories = [category for category in categories if category not in retained_categories]
    if omitted_categories:
        log_detail(
            "Excluded categories with insufficient pseudobulk replicates from the shared fit: "
            + ", ".join(omitted_categories),
            level=2,
        )
    if len(retained_categories) < 2:
        log_warning(
            f"pseudobulk DE for '{groupby}' was skipped: fewer than two categories have "
            f">= {required_min_replicates} replicate pseudobulks with >= {int(min_cells)} cells.",
            level=2,
        )
        return None

    log_detail(
        f"Shared-fit cohort: {len(retained_categories)} retained categories; "
        f"requiring >= {required_min_replicates} replicate pseudobulks and >= {int(min_cells)} cells each",
        level=2,
    )

    model_mask = (
        pb_meta["_pb_group"].isin(retained_categories)
        & (pb_meta["n_cells"] >= int(min_cells))
    )
    model_positions = np.flatnonzero(model_mask.to_numpy())
    model_counts = aggregate[model_positions]
    model_meta = pb_meta.iloc[model_positions].copy()
    model_meta.attrs["gene_names"] = [str(gene) for gene in adata.var_names]
    design_rank, design_columns = _shared_category_design_rank(model_meta)
    residual_df = int(len(model_meta) - design_rank)
    if design_rank < design_columns or residual_df <= 0:
        log_warning(
            f"shared pseudobulk DE for '{groupby}' was skipped: the "
            f"~ {replicate} + {groupby} design is rank-deficient or has no residual "
            f"degrees of freedom (rank {design_rank}/{design_columns}; residual df {residual_df}).",
            level=2,
        )
        return None
    log_detail(
        f"Validated shared DESeq2 design ~ {replicate} + {groupby}: "
        f"rank {design_rank}/{design_columns}; residual df {residual_df}",
        level=2,
    )
    model_keys = set(
        zip(
            model_meta["_pb_replicate"].astype(str),
            model_meta["_pb_group"].astype(str),
        )
    )
    model_cell_mask = valid & np.fromiter(
        (
            (str(rep_value), str(group_value)) in model_keys
            for rep_value, group_value in zip(rep_values.to_numpy(), group_values.to_numpy())
        ),
        dtype=bool,
        count=adata.n_obs,
    )

    fit_counts, fit_meta = _filter_pseudobulk_genes(
        model_counts,
        model_meta,
        min_gene_counts,
    )
    if not fit_meta.attrs.get("gene_names"):
        log_warning(
            f"shared pseudobulk DE for '{groupby}' was skipped: no genes meet "
            f"min_gene_counts={max(0, int(min_gene_counts))}.",
            level=2,
        )
        return None
    log_detail(
        f"Shared-fit gene filter: min_gene_counts={max(0, int(min_gene_counts))}; "
        f"retained {fit_counts.shape[1]} of {model_counts.shape[1]} genes before DESeq2 fitting",
        level=2,
    )
    estimate = _estimate_shared_deseq2_fit_seconds(
        fit_counts.shape[0],
        fit_counts.shape[1],
        design_columns,
        n_cpus,
    )
    log_detail(
        "Fitting shared all-category DESeq2 model (normalization and dispersion estimation); "
        f"rough {max(1, int(n_cpus))}-CPU estimate: {_format_duration(estimate)}.",
        level=2,
    )
    fit_started = time.perf_counter()
    try:
        dds, fit_counts, fit_meta = _fit_deseq2_shared_categories(
            fit_counts,
            fit_meta,
            retained_categories,
            fit_type=str(fit_type or "parametric"),
            min_gene_counts=0,
            n_cpus=n_cpus,
        )
    except Exception as exc:
        log_warning(f"shared pseudobulk DE fit for '{groupby}' failed ({exc}).", level=2)
        return None

    fit_elapsed = time.perf_counter() - fit_started
    log_detail(
        f"Shared all-category DESeq2 model fit completed in {_format_duration(fit_elapsed)}.",
        level=2,
    )
    log_detail(
        f"Shared DESeq2 fit output: {len(retained_categories)} categories, "
        f"{fit_counts.shape[0]} pseudobulk samples, {fit_counts.shape[1]} genes "
        f"(~ {replicate} + {groupby}; design rank {design_rank}/{design_columns})",
        level=2,
    )
    model_info = {
        "model_formula": f"~ {replicate} + {groupby}",
        "model_categories": retained_categories,
        "contrast_reference": "balanced_rest",
    }
    results: Dict[str, Dict[str, Dict[str, Any]]] = {}

    def category_cell_mask(category: str) -> np.ndarray:
        return model_cell_mask & (group_values.to_numpy() == str(category))

    def format_shared_result(
        raw_result: pd.DataFrame,
        *,
        contrast_type: str,
        n_replicates: int,
        source_mask: np.ndarray,
        reference_mask: np.ndarray,
        fit_gene_count: int,
        test_gene_count: int,
    ) -> Dict[str, Any]:
        formatted = _format_result(
            raw_result,
            source_mask=source_mask,
            reference_mask=reference_mask,
            expression_matrix=expression_matrix,
            var_names=adata.var_names,
            top_n=0,
            n_source=int(np.count_nonzero(source_mask)),
            n_reference=int(np.count_nonzero(reference_mask)),
            n_replicates=n_replicates,
            counts_layer_used=counts_layer_used,
            warning=warning,
            p_adjust_method=p_adjust_method,
            min_pct_expressed=min_pct_expressed,
            padj_cutoff=padj_cutoff,
            log2fc_cutoff=log2fc_cutoff,
        )
        formatted.update(model_info)
        formatted["contrast_type"] = contrast_type
        formatted["min_pct_prefilter_gene_count"] = int(fit_gene_count)
        formatted["min_pct_retained_gene_count"] = int(test_gene_count)
        formatted["min_pct_removed_gene_count"] = int(max(0, fit_gene_count - test_gene_count))
        return formatted

    def evaluate_contrast(
        source: str,
        reference: str,
        contrast: np.ndarray,
        *,
        contrast_type: str,
        n_replicates: int,
        reference_mask: np.ndarray,
    ) -> Tuple[str, str, Dict[str, Any], str, Optional[str], str]:
        source_mask = category_cell_mask(source)
        n_source = int(np.count_nonzero(source_mask))
        n_reference = int(np.count_nonzero(reference_mask))
        label = f"{source} vs {reference}" if reference != "__rest__" else f"{source} vs balanced rest"
        try:
            fit_gene_names = [str(gene) for gene in fit_meta.attrs.get("gene_names", [])]
            test_gene_names = _expression_prefilter_gene_names(
                expression_matrix,
                source_mask,
                reference_mask,
                adata.var_names,
                fit_gene_names,
                min_pct_expressed,
            )
            raw_result = _deseq2_shared_contrast(
                dds,
                contrast,
                n_cpus=1,
                gene_names=test_gene_names,
            )
            formatted = format_shared_result(
                raw_result,
                contrast_type=contrast_type,
                n_replicates=n_replicates,
                source_mask=source_mask,
                reference_mask=reference_mask,
                fit_gene_count=len(fit_gene_names),
                test_gene_count=len(test_gene_names),
            )
            return source, reference, formatted, label, None, "shared DESeq2 contrast returned"
        except Exception as exc:
            failed = _empty_result(
                "de_failed",
                n_source=n_source,
                n_reference=n_reference,
                min_cells=int(min_cells),
                min_replicates=required_min_replicates,
                details=str(exc),
            )
            failed.update(model_info)
            failed["contrast_type"] = contrast_type
            return source, reference, failed, label, str(exc), "failed"

    def evaluate_pairwise_contrast(
        source: str,
        reference: str,
        contrast: np.ndarray,
        *,
        n_replicates: int,
        reference_mask: np.ndarray,
    ) -> List[Tuple[str, str, Dict[str, Any], str, Optional[str], str]]:
        source_mask = category_cell_mask(source)
        n_source = int(np.count_nonzero(source_mask))
        n_reference = int(np.count_nonzero(reference_mask))
        forward_label = f"{source} vs {reference}"
        reverse_label = f"{reference} vs {source}"
        try:
            fit_gene_names = [str(gene) for gene in fit_meta.attrs.get("gene_names", [])]
            test_gene_names = _expression_prefilter_gene_names(
                expression_matrix,
                source_mask,
                reference_mask,
                adata.var_names,
                fit_gene_names,
                min_pct_expressed,
            )
            raw_result = _deseq2_shared_contrast(
                dds,
                contrast,
                n_cpus=1,
                gene_names=test_gene_names,
            )
            forward = format_shared_result(
                raw_result,
                contrast_type="category_vs_category",
                n_replicates=n_replicates,
                source_mask=source_mask,
                reference_mask=reference_mask,
                fit_gene_count=len(fit_gene_names),
                test_gene_count=len(test_gene_names),
            )
            reverse = format_shared_result(
                _reverse_deseq2_result(raw_result),
                contrast_type="category_vs_category",
                n_replicates=n_replicates,
                source_mask=reference_mask,
                reference_mask=source_mask,
                fit_gene_count=len(fit_gene_names),
                test_gene_count=len(test_gene_names),
            )
            return [
                (source, reference, forward, forward_label, None, "shared DESeq2 contrast returned"),
                (reference, source, reverse, reverse_label, None, "derived reverse contrast returned"),
            ]
        except Exception as exc:
            forward_failed = _empty_result(
                "de_failed",
                n_source=n_source,
                n_reference=n_reference,
                min_cells=int(min_cells),
                min_replicates=required_min_replicates,
                details=str(exc),
            )
            forward_failed.update(model_info)
            forward_failed["contrast_type"] = "category_vs_category"
            reverse_failed = _empty_result(
                "de_failed",
                n_source=n_reference,
                n_reference=n_source,
                min_cells=int(min_cells),
                min_replicates=required_min_replicates,
                details=str(exc),
            )
            reverse_failed.update(model_info)
            reverse_failed["contrast_type"] = "category_vs_category"
            return [
                (source, reference, forward_failed, forward_label, str(exc), "failed"),
                (reference, source, reverse_failed, reverse_label, str(exc), "failed"),
            ]

    def log_returned_contrast(
        progress: str,
        label: str,
        formatted: Dict[str, Any],
        error: Optional[str],
        status: str,
    ) -> None:
        if error:
            log_step(f"[{progress}] {label}: failed ({error})", level=3)
            return
        log_step(f"[{progress}] {label}: {status}", level=3)
        if _normalize_pct_threshold(min_pct_expressed) > 0:
            retained_count = int(formatted.get("min_pct_retained_gene_count") or 0)
            prefilter_count = int(formatted.get("min_pct_prefilter_gene_count") or 0)
            log_detail(
                f"Minimum % cells retained {retained_count} "
                f"of {prefilter_count} fitted genes for DeseqStats",
                level=4,
            )
        log_detail(
            f"{_count_threshold_passing_genes(formatted)} "
            f"genes pass the DE thresholds (padj < {float(padj_cutoff):g} "
            f"and |log2FC| >= {float(log2fc_cutoff):g})",
            level=4,
        )

    def run_contrast_tasks(tasks: List[Dict[str, Any]], next_progress_fn) -> None:
        if not tasks:
            return
        workers = min(max(1, int(n_cpus)), len(tasks))
        log_detail(
            f"Running {len(tasks)} shared DESeq2 contrast"
            f"{'s' if len(tasks) != 1 else ''} with up to {workers} parallel worker"
            f"{'s' if workers != 1 else ''}.",
            level=2,
        )
        if workers == 1:
            for task in tasks:
                source, reference, formatted, label, error, status = evaluate_contrast(**task)
                results.setdefault(source, {})[reference] = formatted
                progress = next_progress_fn()
                log_returned_contrast(progress, label, formatted, error, status)
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = [executor.submit(evaluate_contrast, **task) for task in tasks]
                for future in as_completed(futures):
                    source, reference, formatted, label, error, status = future.result()
                    results.setdefault(source, {})[reference] = formatted
                    progress = next_progress_fn()
                    log_returned_contrast(progress, label, formatted, error, status)

    def run_pairwise_contrast_tasks(tasks: List[Dict[str, Any]], next_progress_fn) -> None:
        if not tasks:
            return
        workers = min(max(1, int(n_cpus)), len(tasks))
        log_detail(
            f"Running {len(tasks)} shared DESeq2 pairwise contrast task"
            f"{'s' if len(tasks) != 1 else ''} with up to {workers} parallel worker"
            f"{'s' if workers != 1 else ''}; each task returns both directions.",
            level=2,
        )
        if workers == 1:
            for task in tasks:
                for source, reference, formatted, label, error, status in evaluate_pairwise_contrast(**task):
                    results.setdefault(source, {})[reference] = formatted
                    progress = next_progress_fn()
                    log_returned_contrast(progress, label, formatted, error, status)
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = [executor.submit(evaluate_pairwise_contrast, **task) for task in tasks]
                for future in as_completed(futures):
                    for source, reference, formatted, label, error, status in future.result():
                        results.setdefault(source, {})[reference] = formatted
                        progress = next_progress_fn()
                        log_returned_contrast(progress, label, formatted, error, status)

    selected_retained = [category for category in selected_pairwise_categories if category in retained_categories]
    directed_pairs = len(selected_retained) * max(0, len(selected_retained) - 1)
    unique_pairs = [
        (source, reference)
        for source_index, source in enumerate(selected_retained)
        for reference in selected_retained[source_index + 1:]
    ]
    balanced_rest_total = len(retained_categories)
    pairwise_total = directed_pairs
    balanced_rest_completed = 0
    pairwise_completed = 0
    balanced_rest_tasks: List[Dict[str, Any]] = []
    pairwise_contrast_tasks: List[Dict[str, Any]] = []

    def next_balanced_rest_progress() -> str:
        nonlocal balanced_rest_completed
        balanced_rest_completed += 1
        return f"{balanced_rest_completed}/{balanced_rest_total}"

    def next_pairwise_progress() -> str:
        nonlocal pairwise_completed
        pairwise_completed += 1
        return f"{pairwise_completed}/{pairwise_total}"

    pair_diagnostics: Dict[str, Dict[str, Dict[str, Any]]] = {}

    rest_reference = "__rest__"
    log_detail(
        f"{len(retained_categories)} category-versus-balanced-rest contrasts from the shared fit; "
        "output feeds category-versus-rest marker summaries.",
        level=2,
    )
    for source in retained_categories:
        references = [category for category in retained_categories if category != source]
        source_reps = set(fit_meta.loc[fit_meta["_pb_group"] == source, "_pb_replicate"].astype(str))
        paired_reps = sorted(
            rep_value for rep_value in source_reps
            if any(
                rep_value in set(fit_meta.loc[fit_meta["_pb_group"] == reference, "_pb_replicate"].astype(str))
                for reference in references
            )
        )
        source_mask = category_cell_mask(source)
        reference_mask = model_cell_mask & (group_values.to_numpy() != source)
        if len(paired_reps) < required_min_replicates:
            progress = next_balanced_rest_progress()
            log_step(
                f"[{progress}] {source} vs balanced rest: "
                f"skipped, insufficient paired replicates "
                f"({len(paired_reps)}; need >= {required_min_replicates})",
                level=3,
            )
            empty = _empty_result(
                "insufficient_replicates",
                n_source=int(np.count_nonzero(source_mask)),
                n_reference=int(np.count_nonzero(reference_mask)),
                min_cells=int(min_cells),
                min_replicates=required_min_replicates,
                details=f"{len(paired_reps)} paired replicate(s) available",
            )
            empty.update(model_info)
            empty["contrast_type"] = "balanced_rest"
            results.setdefault(source, {})[rest_reference] = empty
            continue
        balanced_contrast = np.mean(
            [_shared_group_contrast(dds, source, reference) for reference in references],
            axis=0,
        )
        balanced_rest_tasks.append(
            {
                "source": source,
                "reference": rest_reference,
                "contrast": balanced_contrast,
                "contrast_type": "balanced_rest",
                "n_replicates": len(paired_reps),
                "reference_mask": reference_mask,
            }
        )

    run_contrast_tasks(balanced_rest_tasks, next_balanced_rest_progress)

    if unique_pairs:
        log_detail(
            f"Preparing {len(unique_pairs)} pairwise PCA/distance diagnostic payload"
            f"{'s' if len(unique_pairs) != 1 else ''} for the selected annotations; "
            "output feeds the Samples tab.",
            level=2,
        )
        for diagnostic_index, (source, reference) in enumerate(unique_pairs, start=1):
            keep_positions = np.flatnonzero(
                fit_meta["_pb_group"].isin([source, reference]).to_numpy()
            )
            log_step(
                f"Diagnostic [{diagnostic_index}/{len(unique_pairs)}] {source} vs {reference}",
                level=3,
            )
            pair_diagnostics.setdefault(source, {})[reference] = _compute_pseudobulk_sample_diagnostics(
                fit_counts[keep_positions],
                fit_meta.iloc[keep_positions].copy(),
            )
    log_detail(
        f"{len(selected_retained)} selected categories -> {directed_pairs} category-versus-category "
        f"contrasts from {len(unique_pairs)} fitted pairwise comparison"
        f"{'s' if len(unique_pairs) != 1 else ''}; output feeds Raw table, Genes, and Pathway Enrichment.",
        level=2,
    )
    for source, reference in unique_pairs:
        paired_reps = sorted(
            set(fit_meta.loc[fit_meta["_pb_group"] == source, "_pb_replicate"].astype(str))
            & set(fit_meta.loc[fit_meta["_pb_group"] == reference, "_pb_replicate"].astype(str))
        )
        source_mask = category_cell_mask(source)
        reference_mask = category_cell_mask(reference)
        if len(paired_reps) < required_min_replicates:
            for skipped_source, skipped_reference, skipped_source_mask, skipped_reference_mask in (
                (source, reference, source_mask, reference_mask),
                (reference, source, reference_mask, source_mask),
            ):
                progress = next_pairwise_progress()
                log_step(
                    f"[{progress}] {skipped_source} vs {skipped_reference} pair: "
                    f"skipped, insufficient paired replicates "
                    f"({len(paired_reps)}; need >= {required_min_replicates})",
                    level=3,
                )
                empty = _empty_result(
                    "insufficient_replicates",
                    n_source=int(np.count_nonzero(skipped_source_mask)),
                    n_reference=int(np.count_nonzero(skipped_reference_mask)),
                    min_cells=int(min_cells),
                    min_replicates=required_min_replicates,
                    details=f"{len(paired_reps)} paired replicate(s) available",
                )
                empty.update(model_info)
                empty["contrast_type"] = "category_vs_category"
                results.setdefault(skipped_source, {})[skipped_reference] = empty
            continue
        pairwise_contrast_tasks.append(
            {
                "source": source,
                "reference": reference,
                "contrast": _shared_group_contrast(dds, source, reference),
                "n_replicates": len(paired_reps),
                "reference_mask": reference_mask,
            }
        )

    run_pairwise_contrast_tasks(pairwise_contrast_tasks, next_pairwise_progress)

    if not results:
        return None
    stored_sources = [
        key for key, value in results.items()
        if not str(key).startswith("_") and isinstance(value, dict)
    ]
    stored_contrasts = sum(
        1
        for source in stored_sources
        for reference in (results.get(source) or {})
        if not str(reference).startswith("_")
    )
    log_detail(
        f"Stored pseudobulk result payload: {len(stored_sources)} source categories, "
        f"{stored_contrasts} contrasts, {'with' if pair_diagnostics else 'without'} pairwise diagnostics.",
        level=2,
    )
    results["_summary"] = {
        "category_gene_means": aggregate_summary,
        "replicate": str(replicate),
        "groupby": str(groupby),
        "counts_layer": counts_layer_used,
        "pairwise_categories": selected_pairwise_categories,
        "model_formula": f"~ {replicate} + {groupby}",
        "model_categories": retained_categories,
        "rest_definition": "balanced_equal_category_weight",
        "diagnostics": "pairwise",
        **({"pair_diagnostics": pair_diagnostics} if pair_diagnostics else {}),
    }
    return results


def compute_pseudobulk_group_de(
    adata,
    groupby: str,
    *,
    replicate: str,
    pairwise_categories: Optional[Sequence[str]] = None,
    counts_layer: Optional[str] = "counts",
    min_cell_counts: int = 0,
    min_gene_counts: int = 0,
    min_cells: int = 20,
    min_replicates: int = 2,
    min_pct_expressed: float = 0.0,
    p_adjust_method: str = "fdr_bh",
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 0.5,
    fit_type: str = "parametric",
    n_cpus: int = 1,
) -> Optional[Dict[str, Dict[str, Dict[str, Any]]]]:
    """Compute shared-fit category and balanced-rest pseudobulk DE.

    Kept as the public entry point used by the loader. The implementation lives
    above so the shared model path is explicit and independently testable.
    """
    return _compute_pseudobulk_group_de_shared(
        adata,
        groupby,
        replicate=replicate,
        pairwise_categories=pairwise_categories,
        counts_layer=counts_layer,
        min_cell_counts=min_cell_counts,
        min_gene_counts=min_gene_counts,
        min_cells=min_cells,
        min_replicates=min_replicates,
        min_pct_expressed=min_pct_expressed,
        p_adjust_method=p_adjust_method,
        padj_cutoff=padj_cutoff,
        log2fc_cutoff=log2fc_cutoff,
        fit_type=fit_type,
        n_cpus=n_cpus,
    )
