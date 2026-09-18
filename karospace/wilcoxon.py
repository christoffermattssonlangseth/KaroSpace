"""Cell-level Wilcoxon marker utilities."""

from __future__ import annotations

import math
import warnings
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

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
_MISSING_DISTRIBUTION_COUNTS_LAYER_WARNED: Set[str] = set()
_MISSING_DISTRIBUTION_NORMALIZED_LAYER_WARNED: Set[str] = set()
_MISSING_STATISTICS_FILTER_COUNTS_LAYER_WARNED: Set[str] = set()
_WILCOXON_NORMALIZE_TARGET_SUM = 10000.0
_WILCOXON_MODES = {"auto", "off", "force"}
_WILCOXON_RUNTIME_SAFETY_MULTIPLIER = 1.25


def normalize_wilcoxon_mode(mode: Optional[str]) -> str:
    value = str(mode or "auto").strip().lower()
    if value in {"", "true", "yes", "on", "1"}:
        return "auto"
    if value in {"false", "no", "none", "0"}:
        return "off"
    if value not in _WILCOXON_MODES:
        raise ValueError("wilcoxon must be one of: auto, off, force")
    return value


def parse_wilcoxon_runtime_limit(value: Any) -> float:
    if value is None:
        return 30.0 * 60.0
    if isinstance(value, (int, float)):
        seconds = float(value)
        if seconds < 0:
            raise ValueError("wilcoxon_runtime_limit must be >= 0")
        return seconds
    text = str(value).strip()
    if not text:
        return 30.0 * 60.0
    parts = text.split(":")
    try:
        if len(parts) == 1:
            seconds = float(parts[0])
        elif len(parts) == 3:
            hours, minutes, seconds_part = parts
            seconds = int(hours) * 3600.0 + int(minutes) * 60.0 + float(seconds_part)
        else:
            raise ValueError
    except ValueError as exc:
        raise ValueError("wilcoxon_runtime_limit must be seconds or HH:MM:SS") from exc
    if seconds < 0:
        raise ValueError("wilcoxon_runtime_limit must be >= 0")
    return float(seconds)


def format_wilcoxon_runtime(seconds: float) -> str:
    total = max(0, int(round(float(seconds))))
    hours, rem = divmod(total, 3600)
    minutes, secs = divmod(rem, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def decide_wilcoxon_runtime(
    mode: Optional[str],
    runtime_limit_seconds: float,
    estimate_items: Sequence[Dict[str, Any]],
    *,
    n_cpus: int = 1,
) -> Dict[str, Any]:
    mode_norm = normalize_wilcoxon_mode(mode)
    limit = max(0.0, float(runtime_limit_seconds))
    estimated = float(sum(float(item.get("estimated_seconds") or 0.0) for item in estimate_items))
    by_kind: Dict[str, float] = {}
    for item in estimate_items:
        kind = str(item.get("kind") or "wilcoxon")
        by_kind[kind] = by_kind.get(kind, 0.0) + float(item.get("estimated_seconds") or 0.0)
    if mode_norm == "off":
        should_run = False
        reason = "Wilcoxon disabled by mode=off."
    elif mode_norm == "force":
        should_run = True
        reason = "Wilcoxon mode=force; ignoring runtime limit."
    elif estimated > limit:
        should_run = False
        reason = (
            "Predicted Wilcoxon runtime exceeds runtime limit "
            f"({format_wilcoxon_runtime(estimated)} > {format_wilcoxon_runtime(limit)})."
        )
    else:
        should_run = True
        reason = (
            "Predicted Wilcoxon runtime is within runtime limit "
            f"({format_wilcoxon_runtime(estimated)} <= {format_wilcoxon_runtime(limit)})."
        )
    return {
        "mode": mode_norm,
        "should_run": bool(should_run),
        "reason": reason,
        "estimated_seconds": estimated,
        "runtime_limit_seconds": limit,
        "estimated_by_kind": by_kind,
        "n_cpus": max(1, int(n_cpus)),
        "items": list(estimate_items),
    }


def wilcoxon_runtime_decision_log_lines(decision: Mapping[str, Any]) -> List[Tuple[str, int, str]]:
    estimated = float(decision.get("estimated_seconds") or 0.0)
    limit = float(decision.get("runtime_limit_seconds") or 0.0)
    action = "run" if bool(decision.get("should_run")) else "skip"
    lines: List[Tuple[str, int, str]] = [
        (
            "step",
            0,
            "runtime_decision="
            f"{action}; estimated={format_wilcoxon_runtime(estimated)}; "
            f"limit={format_wilcoxon_runtime(limit)}",
        )
    ]
    reason = str(decision.get("reason") or "").strip()
    if reason:
        lines.append(("step", 0, reason))

    items = list(decision.get("items") or [])
    fixed_seconds = float(sum(float(item.get("estimated_fixed_seconds") or 0.0) for item in items))
    variable_seconds = float(sum(float(item.get("estimated_variable_seconds") or 0.0) for item in items))
    pre_safety_seconds = fixed_seconds + variable_seconds
    safety_multiplier = estimated / pre_safety_seconds if pre_safety_seconds > 0.0 else 1.0
    serial_fixed_seconds = float(
        sum(
            float(item.get("estimated_serial_fixed_seconds", item.get("estimated_fixed_seconds") or 0.0))
            for item in items
        )
    )
    serial_variable_seconds = float(
        sum(
            float(
                item.get(
                    "estimated_serial_variable_seconds",
                    item.get("estimated_variable_seconds") or 0.0,
                )
            )
            for item in items
        )
    )
    serial_seconds = serial_fixed_seconds + serial_variable_seconds
    comparison_count = int(sum(int(item.get("comparison_count") or 0) for item in items))
    work_units = float(sum(float(item.get("work_units") or 0.0) for item in items))
    if items and (fixed_seconds > 0.0 or variable_seconds > 0.0):
        n_cpus = max(1, int(decision.get("n_cpus") or 1))
        if n_cpus > 1:
            calculation = (
                "runtime_calculation="
                f"serial=(fixed={format_wilcoxon_runtime(serial_fixed_seconds)} + "
                f"variable={format_wilcoxon_runtime(serial_variable_seconds)}) = "
                f"{format_wilcoxon_runtime(serial_seconds)}; "
                f"parallel_adjusted=(fixed={format_wilcoxon_runtime(fixed_seconds)} + "
                f"variable={format_wilcoxon_runtime(variable_seconds)}) with n_cpus={n_cpus}; "
                f"* safety_multiplier={safety_multiplier:.3g} = {format_wilcoxon_runtime(estimated)}; "
                f"comparisons={comparison_count:,}; work_units={work_units:,.0f}"
            )
        else:
            calculation = (
                "runtime_calculation="
                f"(fixed={format_wilcoxon_runtime(fixed_seconds)} + "
                f"variable={format_wilcoxon_runtime(variable_seconds)}) * "
                f"safety_multiplier={safety_multiplier:.3g} = {format_wilcoxon_runtime(estimated)}; "
                f"comparisons={comparison_count:,}; work_units={work_units:,.0f}"
            )
        lines.append(
            (
                "detail",
                1,
                calculation,
            )
        )

    details_by_kind: Dict[str, Dict[str, float]] = {}
    for item in items:
        kind = str(item.get("kind") or "wilcoxon")
        detail = details_by_kind.setdefault(
            kind,
            {
                "fixed": 0.0,
                "variable": 0.0,
                "serial_fixed": 0.0,
                "serial_variable": 0.0,
                "comparisons": 0.0,
                "work_units": 0.0,
                "workers": 0.0,
            },
        )
        detail["fixed"] += float(item.get("estimated_fixed_seconds") or 0.0)
        detail["variable"] += float(item.get("estimated_variable_seconds") or 0.0)
        detail["serial_fixed"] += float(
            item.get("estimated_serial_fixed_seconds", item.get("estimated_fixed_seconds") or 0.0)
        )
        detail["serial_variable"] += float(
            item.get("estimated_serial_variable_seconds", item.get("estimated_variable_seconds") or 0.0)
        )
        detail["comparisons"] += float(int(item.get("comparison_count") or 0))
        detail["work_units"] += float(item.get("work_units") or 0.0)
        detail["workers"] = max(detail["workers"], float(item.get("estimated_parallel_workers") or 1.0))
    for kind, seconds in (decision.get("estimated_by_kind") or {}).items():
        detail = details_by_kind.get(str(kind))
        if detail is None:
            lines.append(("detail", 1, f"estimated_by_kind={kind}; runtime={format_wilcoxon_runtime(seconds)}"))
            continue
        serial_kind_seconds = detail["serial_fixed"] + detail["serial_variable"]
        worker_note = ""
        if max(1, int(decision.get("n_cpus") or 1)) > 1:
            worker_note = (
                f"serial={format_wilcoxon_runtime(serial_kind_seconds)}; "
                f"effective_workers={int(detail['workers'])}; "
            )
        lines.append(
            (
                "detail",
                1,
                f"estimated_by_kind={kind}; runtime={format_wilcoxon_runtime(seconds)}; "
                f"{worker_note}"
                f"fixed={format_wilcoxon_runtime(detail['fixed'])}; "
                f"variable={format_wilcoxon_runtime(detail['variable'])}; "
                f"comparisons={int(detail['comparisons']):,}; "
                f"work_units={detail['work_units']:,.0f}",
            )
        )
    return lines


def _copy_expression_matrix(matrix):
    return matrix.copy() if sp.issparse(matrix) else np.array(matrix, dtype=float, copy=True)


def _run_wilcoxon_tasks(tasks: Sequence[Any], worker, n_cpus: int) -> List[Any]:
    if not tasks:
        return []
    workers = min(max(1, int(n_cpus)), len(tasks))
    if workers <= 1:
        return [worker(task) for task in tasks]
    with ThreadPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(worker, tasks))


def library_size_normalize_expression_matrix(matrix, target_sum: float = _WILCOXON_NORMALIZE_TARGET_SUM):
    """Library-size normalize an expression matrix without log transformation."""
    target = float(target_sum)
    if sp.issparse(matrix):
        normalized = matrix.astype(np.float64, copy=True).tocsr()
        normalized.data[~np.isfinite(normalized.data)] = 0.0
        row_sums = np.asarray(normalized.sum(axis=1)).ravel()
        scale = np.zeros_like(row_sums, dtype=np.float64)
        valid = np.isfinite(row_sums) & (row_sums > 0)
        scale[valid] = target / row_sums[valid]
        normalized = sp.diags(scale).dot(normalized).tocsr()
        normalized.data[~np.isfinite(normalized.data)] = 0.0
        normalized.eliminate_zeros()
        return normalized

    normalized = np.array(matrix, dtype=np.float64, copy=True)
    normalized[~np.isfinite(normalized)] = 0.0
    row_sums = normalized.sum(axis=1)
    scale = np.zeros_like(row_sums, dtype=np.float64)
    valid = np.isfinite(row_sums) & (row_sums > 0)
    scale[valid] = target / row_sums[valid]
    normalized *= scale[:, None]
    normalized[~np.isfinite(normalized)] = 0.0
    return normalized


def log_normalize_expression_matrix(matrix, target_sum: float = _WILCOXON_NORMALIZE_TARGET_SUM):
    """Library-size normalize an expression matrix and apply log1p."""
    normalized = library_size_normalize_expression_matrix(matrix, target_sum=target_sum)
    if sp.issparse(normalized):
        normalized = normalized.tocsr(copy=True)
        with np.errstate(invalid="ignore"):
            normalized.data = np.log1p(normalized.data)
        normalized.data[~np.isfinite(normalized.data)] = 0.0
        normalized.eliminate_zeros()
        return normalized
    with np.errstate(invalid="ignore"):
        normalized = np.log1p(normalized)
    normalized[~np.isfinite(normalized)] = 0.0
    return normalized


def normalize_distribution_normalization(method: Optional[str]) -> str:
    text = str(method or "RC").strip().lower()
    if text in {"rc", "relative_counts", "relative-counts", "relativecounts"}:
        return "RC"
    if text in {"lognormalize", "log_normalize", "log-normalize", "lognorm"}:
        return "LogNormalize"
    raise ValueError("statistics_normalization must be one of: 'RC', 'LogNormalize'")


def resolve_distribution_expression_matrix(
    adata,
    counts_layer: Optional[str] = "counts",
    normalization: str = "RC",
    scale_factor: float = _WILCOXON_NORMALIZE_TARGET_SUM,
    normalized_layer: Optional[str] = None,
) -> Tuple[Any, str]:
    """Resolve display-scale expression values for Distribution panels."""
    layers = getattr(adata, "layers", None) or {}
    normalized_layer_name = str(normalized_layer or "").strip()
    if normalized_layer_name:
        if normalized_layer_name in layers:
            return layers[normalized_layer_name], f"{normalized_layer_name}_layer"
        if normalized_layer_name not in _MISSING_DISTRIBUTION_NORMALIZED_LAYER_WARNED:
            log_warning(
                f"statistics normalized layer '{normalized_layer_name}' was not found; "
                "falling back to the configured statistics count matrix and normalization.",
                level=2,
            )
            _MISSING_DISTRIBUTION_NORMALIZED_LAYER_WARNED.add(normalized_layer_name)

    method = normalize_distribution_normalization(normalization)
    source_matrix = adata.X
    source_name = "X"
    if counts_layer and counts_layer in layers:
        source_matrix = layers[counts_layer]
        source_name = str(counts_layer)
    elif counts_layer:
        missing_key = str(counts_layer)
        if missing_key not in _MISSING_DISTRIBUTION_COUNTS_LAYER_WARNED:
            log_warning(
                f"statistics counts layer '{counts_layer}' was not found; using adata.X for Distribution values.",
                level=2,
            )
            _MISSING_DISTRIBUTION_COUNTS_LAYER_WARNED.add(missing_key)

    if method == "RC":
        return (
            library_size_normalize_expression_matrix(source_matrix, target_sum=float(scale_factor)),
            f"{source_name}_library_normalized",
        )
    return (
        log_normalize_expression_matrix(source_matrix, target_sum=float(scale_factor)),
        f"{source_name}_log_normalized",
    )


def resolve_statistics_filter_counts_matrix(
    adata,
    counts_layer: Optional[str] = "counts",
) -> Tuple[Any, str]:
    """Resolve the raw-count matrix used for Statistics count filters."""
    layers = getattr(adata, "layers", None) or {}
    if counts_layer and counts_layer in layers:
        return layers[counts_layer], str(counts_layer)
    if counts_layer:
        missing_key = str(counts_layer)
        if missing_key not in _MISSING_STATISTICS_FILTER_COUNTS_LAYER_WARNED:
            log_warning(
                f"statistics counts layer '{counts_layer}' was not found; using adata.X for Statistics count filters.",
                level=2,
            )
            _MISSING_STATISTICS_FILTER_COUNTS_LAYER_WARNED.add(missing_key)
    return adata.X, "X"


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


def _subset_expression_matrix(matrix, row_mask: np.ndarray, feature_mask: np.ndarray):
    subset = matrix[row_mask]
    subset = subset[:, feature_mask]
    return subset.copy() if sp.issparse(subset) else np.asarray(subset, dtype=float)


def _matrix_axis_sum(matrix, axis: int) -> np.ndarray:
    values = np.asarray(matrix.sum(axis=axis)).ravel() if sp.issparse(matrix) else np.asarray(matrix, dtype=float).sum(axis=axis)
    values = np.asarray(values, dtype=float).ravel()
    values[~np.isfinite(values)] = 0.0
    return values


def _column_means(matrix, mask: np.ndarray) -> np.ndarray:
    if not mask.any():
        return np.zeros(int(matrix.shape[1]), dtype=float)
    subset = matrix[mask]
    means = np.asarray(subset.mean(axis=0), dtype=float).ravel() if sp.issparse(subset) else np.asarray(subset, dtype=float).mean(axis=0)
    means = np.asarray(means, dtype=float)
    means[~np.isfinite(means)] = 0.0
    return means


def _filtered_labels_and_feature_count_for_wilcoxon(
    labels: np.ndarray,
    expression_matrix,
    filter_expression_matrix=None,
    *,
    min_cell_counts: int = 0,
    min_feature_counts: int = 0,
) -> Tuple[np.ndarray, int]:
    if expression_matrix is None or int(expression_matrix.shape[0]) != int(labels.shape[0]):
        return labels[:0], 0
    if (
        filter_expression_matrix is None
        or int(filter_expression_matrix.shape[0]) != int(expression_matrix.shape[0])
        or int(filter_expression_matrix.shape[1]) != int(expression_matrix.shape[1])
    ):
        filter_expression_matrix = expression_matrix
    min_cell_counts_eff = max(0, int(min_cell_counts))
    min_feature_counts_eff = max(0, int(min_feature_counts))
    count_cell_mask = np.ones(int(labels.shape[0]), dtype=bool)
    if min_cell_counts_eff > 0:
        totals = _matrix_axis_sum(filter_expression_matrix, axis=1)
        count_cell_mask = np.isfinite(totals) & (totals >= min_cell_counts_eff)
    if not bool(count_cell_mask.any()):
        return labels[:0], 0
    n_features = int(expression_matrix.shape[1])
    if min_feature_counts_eff > 0:
        feature_totals = _matrix_axis_sum(filter_expression_matrix[count_cell_mask], axis=0)
        n_features = int(np.count_nonzero(np.isfinite(feature_totals) & (feature_totals >= min_feature_counts_eff)))
    return labels[count_cell_mask], n_features


def _wilcoxon_rank_work_units(cell_count: int, feature_count: int) -> float:
    cells = max(0, int(cell_count))
    features = max(0, int(feature_count))
    if cells == 0 or features == 0:
        return 0.0
    return float(cells * np.log2(max(2, cells)) * features)


def estimate_wilcoxon_group_workload(
    adata,
    annotation_key: str,
    *,
    expression_matrix,
    filter_expression_matrix=None,
    pairwise_categories: Optional[Sequence[str]] = None,
    min_cell_counts: int = 0,
    min_feature_counts: int = 0,
    min_cells: int = 20,
) -> List[Dict[str, Any]]:
    if annotation_key not in adata.obs.columns:
        return []
    col = adata.obs[annotation_key]
    if pd.api.types.is_numeric_dtype(col):
        return []
    if not isinstance(col.dtype, CategoricalDtype):
        col = col.astype("category")
    labels = col.astype(str).to_numpy()
    categories = [str(category) for category in col.cat.categories]
    labels, n_features = _filtered_labels_and_feature_count_for_wilcoxon(
        labels,
        expression_matrix,
        filter_expression_matrix,
        min_cell_counts=min_cell_counts,
        min_feature_counts=min_feature_counts,
    )
    if labels.size == 0 or n_features <= 0:
        return []
    min_cells_eff = max(1, int(min_cells))
    retained = [category for category in categories if int(np.count_nonzero(labels == category)) >= min_cells_eff]
    if len(retained) < 2:
        return []

    retained_mask = np.isin(labels, retained)
    retained_cells = int(retained_mask.sum())
    items: List[Dict[str, Any]] = [
        {
            "kind": "category_vs_rest",
            "annotation": str(annotation_key),
            "comparison_count": int(len(retained)),
            "cell_count": retained_cells,
            "feature_count": int(n_features),
            "work_units": _wilcoxon_rank_work_units(retained_cells, int(n_features)) * len(retained),
        }
    ]
    if pairwise_categories is None:
        pairwise = retained
    else:
        requested = {str(category) for category in pairwise_categories}
        pairwise = [category for category in retained if category in requested]
    pairwise_units = 0.0
    pairwise_count = 0
    for source_idx, source in enumerate(pairwise):
        source_n = int(np.count_nonzero(labels == source))
        for reference in pairwise[source_idx + 1 :]:
            reference_n = int(np.count_nonzero(labels == reference))
            if source_n >= min_cells_eff and reference_n >= min_cells_eff:
                pairwise_count += 1
                pairwise_units += _wilcoxon_rank_work_units(source_n + reference_n, int(n_features))
    if pairwise_count > 0:
        items.append(
            {
                "kind": "category_vs_category",
                "annotation": str(annotation_key),
                "comparison_count": int(pairwise_count),
                "cell_count": int(retained_cells),
                "feature_count": int(n_features),
                "work_units": float(pairwise_units),
            }
        )
    return items


def estimate_wilcoxon_interaction_workload(
    adata,
    annotation_key: str,
    *,
    expression_matrix,
    top_targets: int = 5,
    min_cells: int = 30,
) -> List[Dict[str, Any]]:
    if annotation_key not in adata.obs.columns:
        return []
    col = adata.obs[annotation_key]
    if pd.api.types.is_numeric_dtype(col):
        return []
    if not isinstance(col.dtype, CategoricalDtype):
        col = col.astype("category")
    labels = col.astype(str).to_numpy()
    categories = [str(category) for category in col.cat.categories]
    n_features = int(expression_matrix.shape[1]) if expression_matrix is not None else 0
    if n_features <= 0:
        return []
    min_cells_eff = max(1, int(min_cells))
    retained = [category for category in categories if int(np.count_nonzero(labels == category)) >= min_cells_eff * 2]
    if len(retained) < 2:
        return []
    target_count = max(1, int(top_targets))
    comparison_count = 0
    work_units = 0.0
    for source in retained:
        source_n = int(np.count_nonzero(labels == source))
        n_targets = min(target_count, max(0, len(retained) - 1))
        comparison_count += n_targets
        work_units += _wilcoxon_rank_work_units(source_n, n_features) * n_targets
    return [
        {
            "kind": "contact_conditioned",
            "annotation": str(annotation_key),
            "comparison_count": int(comparison_count),
            "cell_count": int(labels.size),
            "feature_count": int(n_features),
            "work_units": float(work_units),
        }
    ]


def _matrix_subset_rows_cols(matrix, row_indices: np.ndarray, col_indices: np.ndarray):
    subset = matrix[row_indices]
    subset = subset[:, col_indices]
    return subset.copy() if sp.issparse(subset) else np.asarray(subset, dtype=float)


def _fit_wilcoxon_runtime_model(
    samples: Sequence[Mapping[str, Any]],
    *,
    safety_multiplier: float = _WILCOXON_RUNTIME_SAFETY_MULTIPLIER,
) -> Dict[str, Any]:
    usable = [
        (
            float(sample.get("work_units") or 0.0),
            float(sample.get("elapsed_seconds") or 0.0),
        )
        for sample in samples
        if np.isfinite(float(sample.get("work_units") or 0.0))
        and np.isfinite(float(sample.get("elapsed_seconds") or 0.0))
        and float(sample.get("work_units") or 0.0) > 0.0
        and float(sample.get("elapsed_seconds") or 0.0) >= 0.0
    ]
    if not usable:
        return {
            "intercept_seconds": 0.0,
            "seconds_per_work_unit": 0.0,
            "safety_multiplier": max(1.0, float(safety_multiplier)),
        }
    x = np.asarray([item[0] for item in usable], dtype=float)
    y = np.asarray([item[1] for item in usable], dtype=float)
    if x.size >= 2 and float(np.ptp(x)) > 0.0:
        x_mean = float(np.mean(x))
        y_mean = float(np.mean(y))
        slope = float(np.sum((x - x_mean) * (y - y_mean)) / np.sum((x - x_mean) ** 2))
        intercept = y_mean - slope * x_mean
        if slope < 0.0:
            slope = 0.0
            intercept = y_mean
        if intercept < 0.0:
            intercept = 0.0
            slope = float(np.sum(x * y) / max(float(np.sum(x * x)), 1e-12))
    else:
        intercept = 0.0
        slope = float(y[0] / max(float(x[0]), 1.0))
    return {
        "intercept_seconds": float(max(0.0, intercept)),
        "seconds_per_work_unit": float(max(0.0, slope)),
        "safety_multiplier": max(1.0, float(safety_multiplier)),
    }


def _wilcoxon_runtime_sample_plan(
    n_cells: int,
    n_features: int,
    *,
    max_cells: int,
    max_features: int,
) -> List[Tuple[int, int]]:
    max_sample_cells = min(max(2, int(max_cells)), int(n_cells))
    if max_sample_cells % 2:
        max_sample_cells -= 1
    max_sample_cells = max(2, max_sample_cells)
    max_sample_features = min(max(1, int(max_features)), int(n_features))
    plan: List[Tuple[int, int]] = []
    for fraction in (0.25, 0.5, 1.0):
        sample_cells = min(max_sample_cells, max(2, int(round(max_sample_cells * fraction))))
        if sample_cells % 2:
            sample_cells -= 1
        sample_cells = max(2, sample_cells)
        sample_features = min(max_sample_features, max(1, int(round(max_sample_features * fraction))))
        dims = (sample_cells, sample_features)
        if dims not in plan:
            plan.append(dims)
    return plan


def _measure_wilcoxon_runtime_sample(
    matrix,
    *,
    rng: np.random.Generator,
    sample_cells: int,
    sample_features: int,
) -> Dict[str, Any]:
    n_cells = int(matrix.shape[0])
    n_features = int(matrix.shape[1])
    row_idx = np.sort(rng.choice(np.arange(n_cells), size=int(sample_cells), replace=False))
    col_idx = np.sort(rng.choice(np.arange(n_features), size=int(sample_features), replace=False))
    sample = _matrix_subset_rows_cols(matrix, row_idx, col_idx)
    sample_nnz = int(sample.nnz) if sp.issparse(sample) else int(np.count_nonzero(sample))
    sample_density = float(sample_nnz / max(int(sample_cells) * int(sample_features), 1))
    half = int(sample_cells) // 2
    labels = np.asarray(["calibration_a"] * half + ["calibration_b"] * (int(sample_cells) - half), dtype=object)
    feature_names = [f"feature_{idx}" for idx in range(int(sample_features))]
    fallback_before = bool(_SCANPY_WILCOXON_FALLBACK_WARNED)
    started = time.perf_counter()
    table = _scanpy_wilcoxon_table(
        sample,
        labels,
        source="calibration_a",
        reference="calibration_b",
        feature_names=feature_names,
    )
    source_mask = labels == "calibration_a"
    reference_mask = labels == "calibration_b"
    _format_wilcoxon_result(
        table,
        method="cell-wilcoxon-calibration",
        source_mask=source_mask,
        reference_mask=reference_mask,
        expression_matrix=sample,
        feature_names=feature_names,
        expression_layer_used="calibration_sample",
        p_adjust_method="fdr_bh",
        min_pct_expressed=0.0,
        padj_cutoff=0.05,
        log2fc_cutoff=0.0,
        top_n=int(sample_features),
    )
    elapsed = max(time.perf_counter() - started, 1e-9)
    if not fallback_before and bool(_SCANPY_WILCOXON_FALLBACK_WARNED):
        started = time.perf_counter()
        table = _scanpy_wilcoxon_table(
            sample,
            labels,
            source="calibration_a",
            reference="calibration_b",
            feature_names=feature_names,
        )
        _format_wilcoxon_result(
            table,
            method="cell-wilcoxon-calibration",
            source_mask=source_mask,
            reference_mask=reference_mask,
            expression_matrix=sample,
            feature_names=feature_names,
            expression_layer_used="calibration_sample",
            p_adjust_method="fdr_bh",
            min_pct_expressed=0.0,
            padj_cutoff=0.05,
            log2fc_cutoff=0.0,
            top_n=int(sample_features),
        )
        elapsed = max(time.perf_counter() - started, 1e-9)
    work_units = _wilcoxon_rank_work_units(int(sample_cells), int(sample_features))
    return {
        "sample_cells": int(sample_cells),
        "sample_features": int(sample_features),
        "sample_nnz": int(sample_nnz),
        "sample_density": float(sample_density),
        "sample_zero_fraction": float(1.0 - sample_density),
        "elapsed_seconds": float(elapsed),
        "work_units": work_units,
        "calibration_includes_formatting": True,
    }


def calibrate_wilcoxon_runtime(
    matrix,
    *,
    random_state: int = 0,
    max_cells: int = 2048,
    max_features: int = 512,
) -> Dict[str, Any]:
    n_cells = int(matrix.shape[0]) if matrix is not None else 0
    n_features = int(matrix.shape[1]) if matrix is not None and len(matrix.shape) > 1 else 0
    if n_cells < 2 or n_features < 1:
        return {
            "sample_cells": 0,
            "sample_features": 0,
            "sample_nnz": 0,
            "sample_density": 0.0,
            "sample_zero_fraction": 0.0,
            "elapsed_seconds": 0.0,
            "calibration_total_elapsed_seconds": 0.0,
            "work_units": 0.0,
            "intercept_seconds": 0.0,
            "seconds_per_work_unit": 0.0,
            "safety_multiplier": _WILCOXON_RUNTIME_SAFETY_MULTIPLIER,
            "calibration_samples": [],
            "calibration_includes_formatting": False,
        }
    rng = np.random.default_rng(int(random_state))
    sample_plan = _wilcoxon_runtime_sample_plan(
        n_cells,
        n_features,
        max_cells=max_cells,
        max_features=max_features,
    )
    warmup_cells, warmup_features = sample_plan[0]
    _measure_wilcoxon_runtime_sample(
        matrix,
        rng=rng,
        sample_cells=warmup_cells,
        sample_features=warmup_features,
    )
    calibration_samples = [
        _measure_wilcoxon_runtime_sample(
            matrix,
            rng=rng,
            sample_cells=sample_cells,
            sample_features=sample_features,
        )
        for sample_cells, sample_features in sample_plan
    ]
    largest = max(calibration_samples, key=lambda sample: float(sample.get("work_units") or 0.0))
    model = _fit_wilcoxon_runtime_model(calibration_samples)
    return {
        "sample_cells": int(largest.get("sample_cells") or 0),
        "sample_features": int(largest.get("sample_features") or 0),
        "sample_nnz": int(largest.get("sample_nnz") or 0),
        "sample_density": float(largest.get("sample_density") or 0.0),
        "sample_zero_fraction": float(largest.get("sample_zero_fraction") or 0.0),
        "elapsed_seconds": float(largest.get("elapsed_seconds") or 0.0),
        "calibration_total_elapsed_seconds": float(
            sum(float(sample.get("elapsed_seconds") or 0.0) for sample in calibration_samples)
        ),
        "work_units": float(largest.get("work_units") or 0.0),
        "calibration_samples": calibration_samples,
        "calibration_includes_formatting": True,
        **model,
    }


def apply_wilcoxon_runtime_calibration(
    items: Sequence[Dict[str, Any]],
    calibration: Mapping[str, Any],
    *,
    modality: str,
    n_cpus: int = 1,
) -> List[Dict[str, Any]]:
    seconds_per_unit = float(calibration.get("seconds_per_work_unit") or 0.0)
    intercept_seconds = float(calibration.get("intercept_seconds") or 0.0)
    safety_multiplier = float(calibration.get("safety_multiplier") or 1.0)
    estimated: List[Dict[str, Any]] = []
    requested_workers = max(1, int(n_cpus))
    for item in items:
        entry = dict(item)
        comparison_count = max(1, int(entry.get("comparison_count") or 1))
        effective_workers = min(requested_workers, comparison_count)
        serial_fixed_seconds = comparison_count * intercept_seconds
        serial_variable_seconds = float(entry.get("work_units") or 0.0) * seconds_per_unit
        fixed_seconds = math.ceil(comparison_count / effective_workers) * intercept_seconds
        variable_seconds = serial_variable_seconds / effective_workers
        entry["modality"] = str(modality)
        entry["estimated_fixed_seconds"] = float(fixed_seconds)
        entry["estimated_variable_seconds"] = float(variable_seconds)
        entry["estimated_serial_fixed_seconds"] = float(serial_fixed_seconds)
        entry["estimated_serial_variable_seconds"] = float(serial_variable_seconds)
        entry["estimated_parallel_workers"] = int(effective_workers)
        entry["estimated_requested_workers"] = int(requested_workers)
        entry["runtime_intercept_seconds"] = float(intercept_seconds)
        entry["runtime_seconds_per_work_unit"] = float(seconds_per_unit)
        entry["runtime_safety_multiplier"] = float(max(1.0, safety_multiplier))
        entry["estimated_seconds"] = float((fixed_seconds + variable_seconds) * max(1.0, safety_multiplier))
        estimated.append(entry)
    return estimated


def _scanpy_wilcoxon_table(
    matrix,
    labels: np.ndarray,
    *,
    source: str,
    reference: str,
    feature_names: Sequence[str],
) -> pd.DataFrame:
    global _SCANPY_WILCOXON_FALLBACK_WARNED
    if _SCANPY_WILCOXON_FALLBACK_WARNED:
        return _scipy_ranksum_table(
            matrix,
            labels,
            source=source,
            reference=reference,
            feature_names=feature_names,
        )
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
    summary_expression_matrix: Optional[Any] = None,
    summary_expression_source: Optional[str] = None,
    filter_expression_matrix: Optional[Any] = None,
    filter_expression_source: Optional[str] = None,
    min_cell_counts: int = 0,
    min_feature_counts: int = 0,
    min_cells: int = 20,
    min_pct_expressed: float = 0.0,
    p_adjust_method: str = "fdr_bh",
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1,
    top_n_per_comparison: int = 300,
    n_cpus: int = 1,
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

    if expression_matrix is None:
        expression_matrix, expression_layer_used = resolve_wilcoxon_expression_matrix(adata, expression_layer)
    else:
        expression_layer_used = str(expression_layer_used or expression_layer or "provided")
    if expression_matrix is None or int(expression_matrix.shape[0]) != int(labels.shape[0]):
        return None
    if (
        summary_expression_matrix is not None
        and (
            int(summary_expression_matrix.shape[0]) != int(labels.shape[0])
            or int(summary_expression_matrix.shape[1]) != int(expression_matrix.shape[1])
        )
    ):
        summary_expression_matrix = None
    feature_names = [str(feature) for feature in adata.var_names]
    if int(expression_matrix.shape[1]) != len(feature_names):
        return None

    min_cell_counts_eff = max(0, int(min_cell_counts))
    min_feature_counts_eff = max(0, int(min_feature_counts))
    count_filter_source = str(filter_expression_source or "provided")
    if (
        filter_expression_matrix is None
        or int(filter_expression_matrix.shape[0]) != int(expression_matrix.shape[0])
        or int(filter_expression_matrix.shape[1]) != int(expression_matrix.shape[1])
    ):
        filter_expression_matrix = expression_matrix
        count_filter_source = str(expression_layer_used or "expression_matrix")

    count_cell_mask = np.ones(int(labels.shape[0]), dtype=bool)
    if min_cell_counts_eff > 0:
        count_cell_totals = _matrix_axis_sum(filter_expression_matrix, axis=1)
        count_cell_mask = np.isfinite(count_cell_totals) & (count_cell_totals >= min_cell_counts_eff)
    if not bool(count_cell_mask.any()):
        return None

    feature_mask = np.ones(len(feature_names), dtype=bool)
    if min_feature_counts_eff > 0:
        count_filter_after_cells = filter_expression_matrix[count_cell_mask]
        feature_totals = _matrix_axis_sum(count_filter_after_cells, axis=0)
        feature_mask = np.isfinite(feature_totals) & (feature_totals >= min_feature_counts_eff)
    if not bool(feature_mask.any()):
        return None

    if not bool(count_cell_mask.all()) or not bool(feature_mask.all()):
        labels = labels[count_cell_mask]
        expression_matrix = _subset_expression_matrix(expression_matrix, count_cell_mask, feature_mask)
        if summary_expression_matrix is not None:
            summary_expression_matrix = _subset_expression_matrix(
                summary_expression_matrix,
                count_cell_mask,
                feature_mask,
            )
        feature_names = [feature for feature, keep in zip(feature_names, feature_mask) if bool(keep)]

    retained_categories = [
        category for category in categories if int(np.count_nonzero(labels == category)) >= min_cells_eff
    ]
    if len(retained_categories) < 2:
        return None
    retained_mask = np.isin(labels, retained_categories)
    retained_matrix = _materialize_rows(expression_matrix, retained_mask)
    retained_labels = labels[retained_mask]
    n_cpus_eff = max(1, int(n_cpus))

    payload: Dict[str, Any] = {}

    def _compute_rest(category: str) -> Tuple[str, Optional[Dict[str, Any]]]:
        source_mask = labels == category
        reference_mask = np.isin(labels, [c for c in retained_categories if c != category])
        if not source_mask.any() or int(reference_mask.sum()) < min_cells_eff:
            return category, None
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
            return category, None
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
        return category, result

    for category, result in _run_wilcoxon_tasks(retained_categories, _compute_rest, n_cpus_eff):
        if result is None:
            continue
        payload.setdefault(category, {})["__rest__"] = result

    if pairwise_categories is None:
        pairwise = retained_categories
    else:
        requested_set = {str(category) for category in pairwise_categories}
        pairwise = [category for category in retained_categories if category in requested_set]
    pairwise_tasks = [
        (source, reference)
        for source_idx, source in enumerate(pairwise)
        for reference in pairwise[source_idx + 1 :]
    ]

    def _compute_pairwise(task: Tuple[str, str]) -> Tuple[str, str, Dict[str, Any], Dict[str, Any]]:
        source, reference = task
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
            return source, reference, result, _invert_wilcoxon_pairwise_result(result)
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
            return source, reference, result, _invert_wilcoxon_pairwise_result(result)
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
        return source, reference, result, _invert_wilcoxon_pairwise_result(result)

    for source, reference, result, reverse_result in _run_wilcoxon_tasks(pairwise_tasks, _compute_pairwise, n_cpus_eff):
        payload.setdefault(source, {})[reference] = result
        payload.setdefault(reference, {})[source] = reverse_result

    if not any(not str(key).startswith("_") for key in payload):
        return None
    payload["_summary"] = {
        "category_feature_means": _category_feature_mean_summary(
            summary_expression_matrix if summary_expression_matrix is not None else expression_matrix,
            labels,
            retained_categories,
            feature_names,
            range(len(feature_names)),
            source=str(summary_expression_source or "cell_wilcoxon"),
        ),
        "source": "cell_wilcoxon",
        "filter_expression_source": count_filter_source,
        "min_cell_counts": int(min_cell_counts_eff),
        "min_feature_counts": int(min_feature_counts_eff),
        "n_cells_after_count_filter": int(labels.shape[0]),
        "n_features_after_count_filter": int(len(feature_names)),
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
    n_cpus: int = 1,
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
    tasks: List[Tuple[int, str, int, str, float]] = []
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
            tasks.append((source_idx, source_name, target_idx, target_name, edge_count))

    def _compute_interaction(task: Tuple[int, str, int, str, float]) -> Tuple[str, str, Dict[str, Any]]:
        source_idx, source_name, target_idx, target_name, edge_count = task
        source_local = np.flatnonzero(labels == source_idx)
        target_neighbor_counts = graph[obs_idx[source_local]][:, obs_idx[labels == target_idx]].sum(axis=1)
        target_neighbor_counts = np.asarray(target_neighbor_counts, dtype=float).ravel()
        pos_mask = target_neighbor_counts >= min_neighbors_eff
        neg_mask = target_neighbor_counts < min_neighbors_eff
        pos_local = source_local[pos_mask]
        neg_local = source_local[neg_mask]
        n_contact = int(pos_local.size)
        n_non_contact = int(neg_local.size)
        target_zscore = _target_zscore_value(neighbor_zscore, source_idx, target_idx)
        meta = _interaction_meta(
            target_neighbor_counts=target_neighbor_counts,
            pos_mask=pos_mask,
            neg_mask=neg_mask,
            n_contact=n_contact,
            n_non_contact=n_non_contact,
            edge_count=edge_count,
            target_zscore=target_zscore,
        )
        if n_contact < min_cells_eff or n_non_contact < min_cells_eff:
            result = _empty_interaction_result(
                "too_few_cells",
                n_contact=n_contact,
                n_non_contact=n_non_contact,
                min_cells=min_cells_eff,
                **meta,
            )
            return source_name, target_name, result
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
        return source_name, target_name, result

    for source_name, target_name, result in _run_wilcoxon_tasks(tasks, _compute_interaction, n_cpus):
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
