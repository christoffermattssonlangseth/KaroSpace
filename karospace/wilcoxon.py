"""Cell-level Wilcoxon marker utilities."""

from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from pandas.api.types import CategoricalDtype

from .console import log_detail, log_warning
from .pseudobulk import (
    _adjust_pvalues,
    _json_compact_float,
    _json_float,
    _normalize_pct_threshold,
    _positive_fraction,
)

_SCANPY_WILCOXON_FALLBACK_WARNED = False
_SCANPY_NORMALIZATION_FALLBACK_WARNED = False
_MISSING_WILCOXON_LAYER_WARNED: Set[str] = set()
_RAW_WILCOXON_NORMALIZATION_WARNED: Set[str] = set()
_WILCOXON_NORMALIZE_TARGET_SUM = 10000.0


def _copy_expression_matrix(matrix):
    return matrix.copy() if sp.issparse(matrix) else np.array(matrix, dtype=float, copy=True)


def _warn_raw_wilcoxon_normalized(source_label: str) -> None:
    if source_label in _RAW_WILCOXON_NORMALIZATION_WARNED:
        return
    log_warning(
        f"No 'normalized' layer found for Wilcoxon statistics; raw expression values from "
        f"{source_label} were library-size normalized to target_sum=10000 and log1p "
        "transformed before rank_genes_groups.",
        level=2,
    )
    _RAW_WILCOXON_NORMALIZATION_WARNED.add(source_label)


def _manual_log_normalize_raw_expression_matrix(matrix):
    if sp.issparse(matrix):
        normalized = matrix.astype(np.float64, copy=True).tocsr()
        normalized.data[~np.isfinite(normalized.data)] = 0.0
        row_sums = np.asarray(normalized.sum(axis=1)).ravel()
        scale = np.zeros_like(row_sums, dtype=np.float64)
        valid = np.isfinite(row_sums) & (row_sums > 0)
        scale[valid] = _WILCOXON_NORMALIZE_TARGET_SUM / row_sums[valid]
        normalized = sp.diags(scale).dot(normalized).tocsr()
        with np.errstate(invalid="ignore"):
            normalized.data = np.log1p(normalized.data)
        normalized.data[~np.isfinite(normalized.data)] = 0.0
        normalized.eliminate_zeros()
        return normalized

    normalized = np.array(matrix, dtype=np.float64, copy=True)
    normalized[~np.isfinite(normalized)] = 0.0
    row_sums = normalized.sum(axis=1)
    scale = np.zeros_like(row_sums, dtype=np.float64)
    valid = np.isfinite(row_sums) & (row_sums > 0)
    scale[valid] = _WILCOXON_NORMALIZE_TARGET_SUM / row_sums[valid]
    normalized *= scale[:, None]
    with np.errstate(invalid="ignore"):
        normalized = np.log1p(normalized)
    normalized[~np.isfinite(normalized)] = 0.0
    return normalized


def _scanpy_log_normalize_raw_expression_matrix(matrix):
    import scanpy as sc

    matrix_copy = _copy_expression_matrix(matrix)
    n_obs, n_vars = matrix_copy.shape
    obs = pd.DataFrame(index=[f"cell_{idx}" for idx in range(int(n_obs))])
    var = pd.DataFrame(index=[f"feature_{idx}" for idx in range(int(n_vars))])
    tmp = AnnData(X=matrix_copy, obs=obs, var=var)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Some cells have zero counts",
            category=UserWarning,
        )
        sc.pp.normalize_total(tmp, target_sum=_WILCOXON_NORMALIZE_TARGET_SUM)
    sc.pp.log1p(tmp)
    return tmp.X


def _log_normalize_raw_expression_matrix(matrix, *, source_name: str, source_label: str) -> Tuple[Any, str]:
    try:
        normalized = _scanpy_log_normalize_raw_expression_matrix(matrix)
    except Exception as exc:
        global _SCANPY_NORMALIZATION_FALLBACK_WARNED
        if not _SCANPY_NORMALIZATION_FALLBACK_WARNED:
            log_warning(
                f"Scanpy raw-expression normalization failed ({exc}); using an equivalent "
                "local normalize_total/log1p fallback for Wilcoxon statistics.",
                level=2,
            )
            _SCANPY_NORMALIZATION_FALLBACK_WARNED = True
        normalized = _manual_log_normalize_raw_expression_matrix(matrix)
    _warn_raw_wilcoxon_normalized(source_label)
    return normalized, f"{source_name}_log1p_normalized"


def _auto_wilcoxon_expression_matrix(adata) -> Tuple[Any, str]:
    layers = getattr(adata, "layers", None) or {}
    if "normalized" in layers:
        return layers["normalized"], "normalized"
    if "counts" in layers:
        return _log_normalize_raw_expression_matrix(
            layers["counts"],
            source_name="counts",
            source_label="adata.layers['counts']",
        )
    return _log_normalize_raw_expression_matrix(
        adata.X,
        source_name="X",
        source_label="adata.X",
    )


def resolve_wilcoxon_expression_matrix(adata, expression_layer: Optional[str] = None) -> Tuple[Any, str]:
    """Resolve the matrix used for cell-level Wilcoxon statistics."""
    layers = getattr(adata, "layers", None) or {}
    if expression_layer:
        if expression_layer in layers:
            if str(expression_layer) == "normalized":
                return layers[expression_layer], str(expression_layer)
            return _log_normalize_raw_expression_matrix(
                layers[expression_layer],
                source_name=str(expression_layer),
                source_label=f"adata.layers['{expression_layer}']",
            )
        matrix, layer_name = _auto_wilcoxon_expression_matrix(adata)
        warning_key = str(expression_layer)
        if warning_key not in _MISSING_WILCOXON_LAYER_WARNED:
            available = ", ".join(str(layer) for layer in layers.keys()) or "none"
            log_warning(
                f"wilcoxon_layer '{expression_layer}' is not present in adata.layers "
                f"(available: {available}); using {layer_name}.",
                level=2,
            )
            _MISSING_WILCOXON_LAYER_WARNED.add(warning_key)
        return matrix, layer_name
    return _auto_wilcoxon_expression_matrix(adata)


def _materialize_rows(matrix, row_mask: np.ndarray):
    subset = matrix[row_mask]
    return subset.copy() if sp.issparse(subset) else np.asarray(subset, dtype=float)


def _column_means(matrix, mask: np.ndarray) -> np.ndarray:
    if not mask.any():
        return np.zeros(int(matrix.shape[1]), dtype=float)
    subset = matrix[mask]
    means = np.asarray(subset.mean(axis=0), dtype=float).ravel() if sp.issparse(subset) else np.asarray(subset, dtype=float).mean(axis=0)
    means = np.asarray(means, dtype=float)
    means[~np.isfinite(means)] = 0.0
    return means


def _scanpy_wilcoxon_table(
    matrix,
    labels: np.ndarray,
    *,
    source: str,
    reference: str,
    feature_names: Sequence[str],
) -> pd.DataFrame:
    try:
        import scanpy as sc

        obs = pd.DataFrame(
            {"__group__": pd.Categorical(labels.astype(str))},
            index=[f"cell_{idx}" for idx in range(int(labels.shape[0]))],
        )
        var = pd.DataFrame(index=[str(feature) for feature in feature_names])
        tmp = AnnData(X=matrix, obs=obs, var=var)
        sc.tl.rank_genes_groups(
            tmp,
            "__group__",
            groups=[str(source)],
            reference=str(reference),
            method="wilcoxon",
            use_raw=False,
        )
        result = sc.get.rank_genes_groups_df(tmp, group=str(source))
        if result is None or result.empty:
            return pd.DataFrame(columns=["feature", "score", "pvalue"])
        return pd.DataFrame(
            {
                "feature": result["names"].astype(str),
                "score": pd.to_numeric(result.get("scores"), errors="coerce"),
                "pvalue": pd.to_numeric(result.get("pvals"), errors="coerce"),
            }
        )
    except Exception as exc:
        global _SCANPY_WILCOXON_FALLBACK_WARNED
        if not _SCANPY_WILCOXON_FALLBACK_WARNED:
            log_warning(f"Scanpy Wilcoxon failed ({exc}); using SciPy rank-sum fallback.", level=2)
            _SCANPY_WILCOXON_FALLBACK_WARNED = True
        return _scipy_ranksum_table(
            matrix,
            labels,
            source=source,
            reference=reference,
            feature_names=feature_names,
        )


def _matrix_to_dense(matrix) -> np.ndarray:
    dense = matrix.toarray() if sp.issparse(matrix) else np.asarray(matrix)
    dense = np.asarray(dense, dtype=float)
    dense[~np.isfinite(dense)] = 0.0
    return dense


def _scipy_ranksum_table(
    matrix,
    labels: np.ndarray,
    *,
    source: str,
    reference: str,
    feature_names: Sequence[str],
) -> pd.DataFrame:
    from scipy import stats

    labels = labels.astype(str)
    source_mask = labels == str(source)
    if str(reference) == "rest":
        reference_mask = ~source_mask
    else:
        reference_mask = labels == str(reference)
    dense = _matrix_to_dense(matrix)
    rows: List[Dict[str, Any]] = []
    for idx, feature in enumerate(feature_names):
        source_values = dense[source_mask, idx]
        reference_values = dense[reference_mask, idx]
        if source_values.size == 0 or reference_values.size == 0:
            score = 0.0
            pvalue = 1.0
        else:
            try:
                stat = stats.ranksums(source_values, reference_values)
                score = float(stat.statistic)
                pvalue = float(stat.pvalue)
            except Exception:
                score = 0.0
                pvalue = 1.0
        if not np.isfinite(pvalue):
            pvalue = 1.0
        rows.append({"feature": str(feature), "score": score, "pvalue": min(max(pvalue, 0.0), 1.0)})
    return pd.DataFrame(rows)


def _empty_result(
    reason: str,
    *,
    method: str,
    n_source: int,
    n_reference: int,
    min_cells: int,
    details: Optional[str] = None,
) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        "available": False,
        "reason": reason,
        "method": method,
        "features": [],
        "log2foldchanges": [],
        "pvals": [],
        "pvals_adj": [],
        "scores": [],
        "pct_source": [],
        "pct_reference": [],
        "base_mean": [],
        "n_source": int(n_source),
        "n_reference": int(n_reference),
        "min_cells_required": int(min_cells),
    }
    if details:
        result["details"] = details
    return result


def _format_wilcoxon_result(
    rank_table: pd.DataFrame,
    *,
    method: str,
    source_mask: np.ndarray,
    reference_mask: np.ndarray,
    expression_matrix,
    feature_names: Sequence[str],
    expression_layer_used: str,
    p_adjust_method: str,
    min_pct_expressed: float,
    padj_cutoff: float,
    log2fc_cutoff: float,
    top_n: int,
    log2fc_direction: str = "two_sided",
) -> Dict[str, Any]:
    n_source = int(source_mask.sum())
    n_reference = int(reference_mask.sum())
    min_pct = _normalize_pct_threshold(min_pct_expressed)
    p_adjust_method = str(p_adjust_method or "fdr_bh").strip().lower().replace("-", "_")
    padj_cutoff = min(max(float(padj_cutoff), 0.0), 1.0)
    log2fc_cutoff = max(float(log2fc_cutoff), 0.0)
    base = {
        "available": True,
        "method": method,
        "p_adjust_method": p_adjust_method,
        "min_pct_expressed": min_pct,
        "padj_cutoff": padj_cutoff,
        "log2fc_cutoff": log2fc_cutoff,
        "table_top_n": int(top_n),
        "n_source": n_source,
        "n_reference": n_reference,
        "n_replicates": 0,
        "counts_layer": expression_layer_used,
    }
    if rank_table is None or rank_table.empty:
        return {
            **base,
            "min_pct_feature_count": 0,
            "features": [],
            "log2foldchanges": [],
            "pvals": [],
            "pvals_adj": [],
            "scores": [],
            "pct_source": [],
            "pct_reference": [],
            "base_mean": [],
        }

    feature_to_idx = {str(feature): idx for idx, feature in enumerate(feature_names)}
    source_means = _column_means(expression_matrix, source_mask)
    reference_means = _column_means(expression_matrix, reference_mask)
    base_means = _column_means(expression_matrix, source_mask | reference_mask)
    eps = 1e-9

    work = rank_table.copy()
    work["_feature"] = work["feature"].astype(str)
    work["_feature_idx"] = [feature_to_idx.get(feature) for feature in work["_feature"]]
    work = work[work["_feature_idx"].notna()].copy()
    if work.empty:
        return {**base, "min_pct_feature_count": 0, "features": [], "log2foldchanges": [], "pvals": [], "pvals_adj": [], "scores": [], "pct_source": [], "pct_reference": [], "base_mean": []}
    indices = [int(idx) for idx in work["_feature_idx"]]
    pct_source = _positive_fraction(expression_matrix, source_mask, indices)
    pct_reference = _positive_fraction(expression_matrix, reference_mask, indices)
    log2fc = np.log2(source_means[indices] + eps) - np.log2(reference_means[indices] + eps)
    pvals = pd.to_numeric(work["pvalue"], errors="coerce").to_numpy(dtype=float).copy()
    pvals[~np.isfinite(pvals)] = 1.0
    pvals = np.clip(pvals, 0.0, 1.0)
    padj = np.asarray(_adjust_pvalues(pvals, p_adjust_method), dtype=float)
    padj[~np.isfinite(padj)] = 1.0

    work["_pct_source"] = [float(v) if v is not None and np.isfinite(v) else 0.0 for v in pct_source]
    work["_pct_reference"] = [float(v) if v is not None and np.isfinite(v) else 0.0 for v in pct_reference]
    work["_log2fc"] = log2fc
    work["_pvalue"] = pvals
    work["_padj"] = padj
    work["_base_mean"] = base_means[indices]
    work["_score"] = pd.to_numeric(work["score"], errors="coerce").fillna(0.0).to_numpy(dtype=float).copy()

    if min_pct > 0:
        work = work[
            (work["_pct_source"].to_numpy(dtype=float) >= min_pct)
            | (work["_pct_reference"].to_numpy(dtype=float) >= min_pct)
        ].copy()
    work = work[np.isfinite(work["_log2fc"].to_numpy(dtype=float))].copy()
    min_pct_feature_count = int(work.shape[0])
    if str(log2fc_direction or "two_sided") == "positive":
        work = work[work["_log2fc"].to_numpy(dtype=float) >= log2fc_cutoff].copy()
    if work.empty:
        return {**base, "min_pct_feature_count": min_pct_feature_count, "features": [], "log2foldchanges": [], "pvals": [], "pvals_adj": [], "scores": [], "pct_source": [], "pct_reference": [], "base_mean": []}

    work["_abs_lfc"] = np.abs(work["_log2fc"].to_numpy(dtype=float))
    work = work.sort_values(
        ["_padj", "_pvalue", "_abs_lfc", "_feature"],
        ascending=[True, True, False, True],
    )
    top_n = max(1, int(top_n))
    work = work.head(top_n)

    return {
        **base,
        "min_pct_feature_count": min_pct_feature_count,
        "features": work["_feature"].astype(str).tolist(),
        "log2foldchanges": [_json_compact_float(v, 6) for v in work["_log2fc"]],
        "pvals": [_json_compact_float(v, 6) for v in work["_pvalue"]],
        "pvals_adj": [_json_compact_float(v, 6) for v in work["_padj"]],
        "scores": [_json_compact_float(v, 6) for v in work["_score"]],
        "pct_source": [_json_compact_float(v, 5) for v in work["_pct_source"]],
        "pct_reference": [_json_compact_float(v, 5) for v in work["_pct_reference"]],
        "base_mean": [_json_compact_float(v, 6) for v in work["_base_mean"]],
    }


def _negate_json_number(value: Any) -> Optional[float]:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(numeric):
        return None
    return _json_compact_float(-numeric, 6)


def _invert_wilcoxon_pairwise_result(result: Dict[str, Any]) -> Dict[str, Any]:
    """Build the reverse pairwise contrast from a computed Wilcoxon result."""
    inverted = dict(result)
    inverted["n_source"] = int(result.get("n_reference", 0))
    inverted["n_reference"] = int(result.get("n_source", 0))
    inverted["log2foldchanges"] = [
        _negate_json_number(value) for value in (result.get("log2foldchanges") or [])
    ]
    inverted["scores"] = [
        _negate_json_number(value) for value in (result.get("scores") or [])
    ]
    inverted["pct_source"] = list(result.get("pct_reference") or [])
    inverted["pct_reference"] = list(result.get("pct_source") or [])
    return inverted


def _threshold_passing_feature_indices(result: Dict[str, Any], feature_to_idx: Dict[str, int]) -> List[int]:
    indices: List[int] = []
    padj_cutoff = float(result.get("padj_cutoff", 0.05))
    log2fc_cutoff = float(result.get("log2fc_cutoff", 0.0))
    for feature, padj, log2fc in zip(
        result.get("features") or [],
        result.get("pvals_adj") or [],
        result.get("log2foldchanges") or [],
    ):
        try:
            if float(padj) < padj_cutoff and abs(float(log2fc)) >= log2fc_cutoff:
                idx = feature_to_idx.get(str(feature))
                if idx is not None:
                    indices.append(idx)
        except (TypeError, ValueError):
            continue
    return indices


def _category_feature_mean_summary(
    expression_matrix,
    labels: np.ndarray,
    categories: Sequence[str],
    feature_names: Sequence[str],
    feature_indices: Sequence[int],
    *,
    source: str,
) -> Dict[str, Any]:
    ordered_indices = sorted(set(int(idx) for idx in feature_indices if int(idx) >= 0))
    all_mask = np.ones(int(labels.shape[0]), dtype=bool)
    background = _column_means(expression_matrix, all_mask)
    means: Dict[str, List[Optional[float]]] = {}
    n_cells: Dict[str, int] = {}
    for category in categories:
        mask = labels == str(category)
        n_cells[str(category)] = int(mask.sum())
        category_means = _column_means(expression_matrix, mask)
        means[str(category)] = [_json_float(category_means[i], 6) for i in ordered_indices]
    return {
        "features": [str(feature_names[i]) for i in ordered_indices],
        "categories": [str(category) for category in categories],
        "means": means,
        "background": [_json_float(background[i], 6) for i in ordered_indices],
        "n_cells": n_cells,
        "source": source,
    }


def compute_wilcoxon_group_de(
    adata,
    annotation_key: str,
    *,
    pairwise_categories: Optional[Sequence[str]] = None,
    expression_layer: Optional[str] = None,
    expression_matrix: Optional[Any] = None,
    expression_layer_used: Optional[str] = None,
    min_cells: int = 20,
    min_pct_expressed: float = 0.0,
    p_adjust_method: str = "fdr_bh",
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1,
    top_n_per_comparison: int = 300,
) -> Optional[Dict[str, Any]]:
    """Compute cell-level Wilcoxon category markers and pairwise comparisons."""
    if annotation_key not in adata.obs.columns:
        return None
    col = adata.obs[annotation_key]
    if pd.api.types.is_numeric_dtype(col):
        return None
    if not isinstance(col.dtype, CategoricalDtype):
        col = col.astype("category")
    labels = col.astype(str).to_numpy()
    categories = [str(category) for category in col.cat.categories]
    min_cells_eff = max(1, int(min_cells))
    retained_categories = [
        category for category in categories if int(np.count_nonzero(labels == category)) >= min_cells_eff
    ]
    if len(retained_categories) < 2:
        return None

    if expression_matrix is None:
        expression_matrix, expression_layer_used = resolve_wilcoxon_expression_matrix(adata, expression_layer)
    else:
        expression_layer_used = str(expression_layer_used or expression_layer or "provided")
    if expression_matrix is None or int(expression_matrix.shape[0]) != int(labels.shape[0]):
        return None
    feature_names = [str(feature) for feature in adata.var_names]
    retained_mask = np.isin(labels, retained_categories)
    retained_matrix = _materialize_rows(expression_matrix, retained_mask)
    retained_labels = labels[retained_mask]

    payload: Dict[str, Any] = {}
    feature_to_idx = {feature: idx for idx, feature in enumerate(feature_names)}
    summary_feature_indices: List[int] = []
    for category in retained_categories:
        source_mask = labels == category
        reference_mask = np.isin(labels, [c for c in retained_categories if c != category])
        if not source_mask.any() or int(reference_mask.sum()) < min_cells_eff:
            continue
        try:
            table = _scanpy_wilcoxon_table(
                retained_matrix,
                retained_labels,
                source=category,
                reference="rest",
                feature_names=feature_names,
            )
        except Exception as exc:
            log_warning(f"Wilcoxon markers for '{annotation_key}' category '{category}' failed ({exc}).", level=2)
            continue
        result = _format_wilcoxon_result(
            table,
            method="cell-wilcoxon-rest",
            source_mask=source_mask,
            reference_mask=reference_mask,
            expression_matrix=expression_matrix,
            feature_names=feature_names,
            expression_layer_used=expression_layer_used,
            p_adjust_method=p_adjust_method,
            min_pct_expressed=min_pct_expressed,
            padj_cutoff=padj_cutoff,
            log2fc_cutoff=log2fc_cutoff,
            top_n=top_n_per_comparison,
            log2fc_direction="positive",
        )
        payload.setdefault(category, {})["__rest__"] = result
        summary_feature_indices.extend(_threshold_passing_feature_indices(result, feature_to_idx))

    if pairwise_categories is None:
        pairwise = retained_categories
    else:
        requested_set = {str(category) for category in pairwise_categories}
        pairwise = [category for category in retained_categories if category in requested_set]
    for source_idx, source in enumerate(pairwise):
        for reference in pairwise[source_idx + 1 :]:
            source_mask = labels == source
            reference_mask = labels == reference
            n_source = int(source_mask.sum())
            n_reference = int(reference_mask.sum())
            if n_source < min_cells_eff or n_reference < min_cells_eff:
                result = _empty_result(
                    "too_few_cells",
                    method="cell-wilcoxon-pairwise",
                    n_source=n_source,
                    n_reference=n_reference,
                    min_cells=min_cells_eff,
                )
                payload.setdefault(source, {})[reference] = result
                payload.setdefault(reference, {})[source] = _invert_wilcoxon_pairwise_result(result)
                continue
            pair_mask = source_mask | reference_mask
            pair_labels = labels[pair_mask]
            pair_matrix = _materialize_rows(expression_matrix, pair_mask)
            try:
                table = _scanpy_wilcoxon_table(
                    pair_matrix,
                    pair_labels,
                    source=source,
                    reference=reference,
                    feature_names=feature_names,
                )
            except Exception as exc:
                result = _empty_result(
                    "wilcoxon_failed",
                    method="cell-wilcoxon-pairwise",
                    n_source=n_source,
                    n_reference=n_reference,
                    min_cells=min_cells_eff,
                    details=str(exc),
                )
                payload.setdefault(source, {})[reference] = result
                payload.setdefault(reference, {})[source] = _invert_wilcoxon_pairwise_result(result)
                continue
            result = _format_wilcoxon_result(
                table,
                method="cell-wilcoxon-pairwise",
                source_mask=source_mask,
                reference_mask=reference_mask,
                expression_matrix=expression_matrix,
                feature_names=feature_names,
                expression_layer_used=expression_layer_used,
                p_adjust_method=p_adjust_method,
                min_pct_expressed=min_pct_expressed,
                padj_cutoff=padj_cutoff,
                log2fc_cutoff=log2fc_cutoff,
                top_n=top_n_per_comparison,
            )
            reverse_result = _invert_wilcoxon_pairwise_result(result)
            payload.setdefault(source, {})[reference] = result
            payload.setdefault(reference, {})[source] = reverse_result
            summary_feature_indices.extend(_threshold_passing_feature_indices(result, feature_to_idx))
            summary_feature_indices.extend(_threshold_passing_feature_indices(reverse_result, feature_to_idx))

    if not any(not str(key).startswith("_") for key in payload):
        return None
    if not summary_feature_indices:
        for by_reference in payload.values():
            if not isinstance(by_reference, dict):
                continue
            for result in by_reference.values():
                if isinstance(result, dict):
                    summary_feature_indices.extend(
                        feature_to_idx[str(feature)]
                        for feature in (result.get("features") or [])[: min(10, len(result.get("features") or []))]
                        if str(feature) in feature_to_idx
                    )
    payload["_summary"] = {
        "category_feature_means": _category_feature_mean_summary(
            expression_matrix,
            labels,
            retained_categories,
            feature_names,
            summary_feature_indices,
            source="cell_wilcoxon",
        ),
        "source": "cell_wilcoxon",
    }
    log_detail(
        f"Stored Wilcoxon marker payload: {len([k for k in payload if not str(k).startswith('_')])} "
        f"categories, layer={expression_layer_used}.",
        level=2,
    )
    return payload


def _empty_interaction_result(
    reason: str,
    *,
    n_contact: int,
    n_non_contact: int,
    min_cells: int,
    pct_contact: float = 0.0,
    mean_target_neighbors_contact: float = 0.0,
    mean_target_neighbors_non_contact: float = 0.0,
    target_edge_count: float = 0.0,
    target_zscore: Optional[float] = None,
    details: Optional[str] = None,
) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        "available": False,
        "reason": reason,
        "method": "cell-wilcoxon-contact",
        "features": [],
        "log2foldchanges": [],
        "pvals": [],
        "pvals_adj": [],
        "scores": [],
        "pct_source": [],
        "pct_reference": [],
        "base_mean": [],
        "n_contact": int(n_contact),
        "n_non_contact": int(n_non_contact),
        "n_source": int(n_contact),
        "n_reference": int(n_non_contact),
        "min_cells_required": int(min_cells),
        "pct_contact": _json_compact_float(pct_contact, 6),
        "mean_target_neighbors_contact": _json_compact_float(mean_target_neighbors_contact, 6),
        "mean_target_neighbors_non_contact": _json_compact_float(mean_target_neighbors_non_contact, 6),
        "target_edge_count": _json_compact_float(target_edge_count, 6),
        "target_zscore": _json_compact_float(target_zscore, 6),
    }
    if details:
        result["details"] = details
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
        "pct_contact": _json_compact_float((100.0 * n_contact) / max(1, n_contact + n_non_contact), 6),
        "mean_target_neighbors_contact": _json_compact_float(
            np.mean(target_neighbor_counts[pos_mask]) if n_contact > 0 else 0.0, 6
        ),
        "mean_target_neighbors_non_contact": _json_compact_float(
            np.mean(target_neighbor_counts[neg_mask]) if n_non_contact > 0 else 0.0, 6
        ),
        "target_edge_count": _json_compact_float(edge_count, 6),
        "target_zscore": _json_compact_float(target_zscore, 6),
    }


def compute_wilcoxon_interaction_markers(
    adata,
    annotation_key: str,
    *,
    graph,
    obs_idx: np.ndarray,
    labels: np.ndarray,
    categories: Sequence[str],
    neighbor_counts: np.ndarray,
    neighbor_zscore: Optional[np.ndarray] = None,
    neighbor_n_cells: Optional[np.ndarray] = None,
    expression_layer: Optional[str] = None,
    expression_matrix: Optional[Any] = None,
    expression_layer_used: Optional[str] = None,
    top_targets: int = 5,
    top_features: int = 20,
    min_cells: int = 30,
    min_neighbors: int = 1,
    min_pct_expressed: float = 0.0,
    p_adjust_method: str = "fdr_bh",
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1,
) -> Optional[Dict[str, Any]]:
    """Compute contact-conditioned source-cell markers using cell-level Wilcoxon."""
    if annotation_key not in adata.obs.columns:
        return None
    labels = np.asarray(labels, dtype=np.int32)
    obs_idx = np.asarray(obs_idx, dtype=np.int64)
    categories = [str(category) for category in categories]
    if len(categories) < 2 or labels.size != obs_idx.size:
        return None
    if graph is None or neighbor_counts is None:
        return None
    if expression_matrix is None:
        expression_matrix, expression_layer_used = resolve_wilcoxon_expression_matrix(adata, expression_layer)
    else:
        expression_layer_used = str(expression_layer_used or expression_layer or "provided")
    if expression_matrix is None or int(expression_matrix.shape[0]) != int(adata.n_obs):
        return None
    feature_names = [str(feature) for feature in adata.var_names]
    top_targets = max(1, int(top_targets))
    top_features = max(1, int(top_features))
    min_cells_eff = max(1, int(min_cells))
    min_neighbors_eff = max(1, int(min_neighbors))
    neighbor_counts = np.asarray(neighbor_counts, dtype=float)
    if neighbor_counts.shape[0] != len(categories) or neighbor_counts.shape[1] != len(categories):
        return None

    payload: Dict[str, Any] = {}
    for source_idx, source_name in enumerate(categories):
        source_local = np.flatnonzero(labels == source_idx)
        if source_local.size < min_cells_eff * 2:
            continue
        target_scores: List[Tuple[float, float, int, str]] = []
        for target_idx, target_name in enumerate(categories):
            if target_idx == source_idx:
                continue
            zvalue = _target_zscore_value(neighbor_zscore, source_idx, target_idx)
            edge_count = float(neighbor_counts[source_idx, target_idx])
            primary = zvalue if zvalue is not None else edge_count
            target_scores.append((float(primary), edge_count, target_idx, target_name))
        target_scores.sort(key=lambda item: (-item[0], -item[1], item[3]))
        for _score, edge_count, target_idx, target_name in target_scores[:top_targets]:
            target_neighbor_counts = graph[obs_idx[source_local]][:, obs_idx[labels == target_idx]].sum(axis=1)
            target_neighbor_counts = np.asarray(target_neighbor_counts, dtype=float).ravel()
            pos_local = source_local[target_neighbor_counts >= min_neighbors_eff]
            neg_local = source_local[target_neighbor_counts < min_neighbors_eff]
            n_contact = int(pos_local.size)
            n_non_contact = int(neg_local.size)
            target_zscore = _target_zscore_value(neighbor_zscore, source_idx, target_idx)
            meta = _interaction_meta(
                target_neighbor_counts=target_neighbor_counts,
                pos_mask=target_neighbor_counts >= min_neighbors_eff,
                neg_mask=target_neighbor_counts < min_neighbors_eff,
                n_contact=n_contact,
                n_non_contact=n_non_contact,
                edge_count=edge_count,
                target_zscore=target_zscore,
            )
            if n_contact < min_cells_eff or n_non_contact < min_cells_eff:
                payload.setdefault(source_name, {})[target_name] = _empty_interaction_result(
                    "too_few_cells",
                    n_contact=n_contact,
                    n_non_contact=n_non_contact,
                    min_cells=min_cells_eff,
                    **meta,
                )
                continue
            source_obs_idx = obs_idx[pos_local]
            reference_obs_idx = obs_idx[neg_local]
            pair_obs_idx = np.concatenate([source_obs_idx, reference_obs_idx])
            pair_labels = np.asarray(["contact"] * n_contact + ["non_contact"] * n_non_contact, dtype=object)
            pair_matrix = expression_matrix[pair_obs_idx]
            if sp.issparse(pair_matrix):
                pair_matrix = pair_matrix.copy()
            else:
                pair_matrix = np.asarray(pair_matrix, dtype=float)
            source_mask = np.zeros(int(adata.n_obs), dtype=bool)
            reference_mask = np.zeros(int(adata.n_obs), dtype=bool)
            source_mask[source_obs_idx] = True
            reference_mask[reference_obs_idx] = True
            try:
                table = _scanpy_wilcoxon_table(
                    pair_matrix,
                    pair_labels,
                    source="contact",
                    reference="non_contact",
                    feature_names=feature_names,
                )
                result = _format_wilcoxon_result(
                    table,
                    method="cell-wilcoxon-contact",
                    source_mask=source_mask,
                    reference_mask=reference_mask,
                    expression_matrix=expression_matrix,
                    feature_names=feature_names,
                    expression_layer_used=expression_layer_used,
                    p_adjust_method=p_adjust_method,
                    min_pct_expressed=min_pct_expressed,
                    padj_cutoff=padj_cutoff,
                    log2fc_cutoff=log2fc_cutoff,
                    top_n=top_features,
                )
            except Exception as exc:
                result = _empty_interaction_result(
                    "wilcoxon_failed",
                    n_contact=n_contact,
                    n_non_contact=n_non_contact,
                    min_cells=min_cells_eff,
                    details=str(exc),
                    **meta,
                )
            result.update(meta)
            result["n_contact"] = n_contact
            result["n_non_contact"] = n_non_contact
            payload.setdefault(source_name, {})[target_name] = result

    if not payload:
        return None
    payload["_summary"] = {
        "source": "cell_wilcoxon_contact",
        "method": "cell-wilcoxon-contact",
        "top_targets": int(top_targets),
        "top_features": int(top_features),
        "min_cells": int(min_cells_eff),
        "min_neighbors": int(min_neighbors_eff),
    }
    log_detail(
        f"Stored Wilcoxon interaction-marker payload: {len([k for k in payload if not str(k).startswith('_')])} "
        f"source categories, layer={expression_layer_used}.",
        level=2,
    )
    return payload
