"""
Data loading utilities for spatial transcriptomics data.

Handles loading h5ad files with scanpy and extracting spatial coordinates,
gene expression, and metadata for visualization.
"""

import json
import os
import re
import shutil
import tempfile
from copy import deepcopy
from dataclasses import dataclass, field
from itertools import combinations, product
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from pandas.api.types import CategoricalDtype
from scipy.spatial import cKDTree
from scipy.sparse import issparse

from .console import log_detail, log_step, log_warning


sc = ad  # Compatibility alias; this module only needs AnnData/read_h5ad, not Scanpy.
COMPANION_ANALYTICS_STORAGE = "json-string-v1"
COMPANION_ANALYTICS_JSON_FIELDS = {
    "cluster_de_json": "cluster_de",
    "pseudobulk_de_json": "pseudobulk_de",
    "neighbor_stats_json": "neighbor_stats",
    "interaction_markers_json": "interaction_markers",
    "gene_correlations_json": "gene_correlations",
    "spatial_variable_genes_json": "spatial_variable_genes",
    "cluster_gene_means_json": "cluster_gene_means",
}


def _clean_pseudobulk_category_list(value: Any, option_name: str) -> Optional[List[str]]:
    if value is None:
        return None
    if isinstance(value, str):
        categories = [item.strip() for item in value.split(",") if item.strip()]
        return categories or None
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        categories = [str(item).strip() for item in value if str(item).strip()]
        return categories or None
    raise ValueError(f"{option_name} values must be category strings or lists of category strings")


def normalize_pseudobulk_simple_constrast_categories(
    value: Any,
    annotation_columns: Sequence[str],
    *,
    option_name: str = "pseudobulk_simple_constrast_categories",
) -> Optional[Dict[str, Optional[List[str]]]]:
    """Normalize optional Simple design category filters by annotation column.

    A flat category list is only meaningful when a single annotation is analyzed.
    With multiple pseudobulk annotation columns, callers must provide either a
    mapping keyed by annotation name or a nested list whose order matches
    ``annotation_columns``.
    """
    columns = [str(column).strip() for column in annotation_columns if str(column).strip()]
    if not columns or value is None:
        return None
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return None
        if text[0] in "[{":
            try:
                value = json.loads(text)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{option_name} must be valid JSON ({exc})") from exc
        else:
            if len(columns) > 1:
                raise ValueError(
                    f"{option_name} is ambiguous with multiple pseudobulk annotations. "
                    "Use a JSON object keyed by annotation name or a nested JSON list in "
                    "the order: " + ", ".join(columns) + "."
                )
            return {columns[0]: _clean_pseudobulk_category_list(text, option_name)}
    if isinstance(value, Mapping):
        unknown = sorted(str(key) for key in value.keys() if str(key) not in set(columns))
        if unknown:
            raise ValueError(
                f"{option_name} contains annotation column(s) not requested for pseudobulk DE: "
                + ", ".join(unknown)
            )
        return {
            column: _clean_pseudobulk_category_list(value.get(column), option_name)
            for column in columns
        }
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        nested = any(
            isinstance(item, (Mapping, list, tuple, set, np.ndarray))
            and not isinstance(item, (str, bytes, bytearray))
            for item in value
        )
        if not nested:
            if len(columns) > 1:
                raise ValueError(
                    f"{option_name} is ambiguous with multiple pseudobulk annotations. "
                    "Use a JSON object keyed by annotation name or a nested JSON list in "
                    "the order: " + ", ".join(columns) + "."
                )
            return {columns[0]: _clean_pseudobulk_category_list(value, option_name)}
        if len(value) != len(columns):
            raise ValueError(
                f"{option_name} nested list must contain one category list per pseudobulk "
                f"annotation ({len(columns)} expected: {', '.join(columns)})."
            )
        normalized: Dict[str, Optional[List[str]]] = {}
        for column, item in zip(columns, value):
            if isinstance(item, Mapping):
                if "categories" in item:
                    item = item.get("categories")
                elif column in item:
                    item = item.get(column)
                else:
                    raise ValueError(
                        f"{option_name} nested mapping for '{column}' must contain "
                        "'categories' or the annotation name."
                    )
            normalized[column] = _clean_pseudobulk_category_list(item, option_name)
        return normalized
    raise ValueError(
        f"{option_name} must be a category list, a JSON object keyed by annotation, "
        "or a nested list matching the pseudobulk annotation order."
    )


# Common conventions for spatial coordinates in adata.obsm, tried in order when
# the requested key is absent. "spatial" is the squidpy/scanpy standard.
_SPATIAL_KEY_FALLBACKS = (
    "spatial",
    "X_spatial",
    "spatial_coords",
    "X_spatial_coords",
    "Spatial",
    "spatialcoords",
)

# Dormant Complex design support. It is intentionally not exposed through the
# API or CLI, and no export path invokes it while the feature is in development.
_PSEUDOBULK_MODEL_FACTOR = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _build_pseudobulk_model_design_audit(
    obs: pd.DataFrame,
    formula: Optional[str],
    pseudobulk_replicate_annotation: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """Parse a safe categorical DESeq2 design subset and describe its pseudobulk grid."""
    text = str(formula or "").strip()
    if not text:
        return None
    audit: Dict[str, Any] = {
        "formula": text,
        "valid": False,
        "terms": [],
        "variables": [],
        "errors": [],
        "warnings": [],
        "sample_count": 0,
        "design_rank": None,
        "design_columns": 0,
        "levels": {},
        "sample_preview": [],
        "candidate_contrasts": [],
        "sample_variables": [],
        "pseudobulk_replicate_annotation": None,
    }
    if not text.startswith("~"):
        audit["errors"].append("The model must start with '~'.")
        return audit
    expression = text[1:].strip()
    if not expression:
        audit["errors"].append("The model has no terms after '~'.")
        return audit
    if any(token in expression for token in ("(", ")", "|", "-", "/")):
        audit["errors"].append(
            "Only categorical main effects and ':' or '*' interactions are supported."
        )
        return audit

    terms: List[Tuple[str, ...]] = []
    for raw_term in (item.strip() for item in expression.split("+")):
        if not raw_term:
            audit["errors"].append("Empty model term.")
            continue
        separators = "*" if "*" in raw_term else ":"
        factors = tuple(item.strip() for item in raw_term.split(separators))
        if not factors or any(not _PSEUDOBULK_MODEL_FACTOR.match(item) for item in factors):
            audit["errors"].append(f"Unsupported model term '{raw_term}'.")
            continue
        expanded = (
            [combo for size in range(1, len(factors) + 1) for combo in combinations(factors, size)]
            if separators == "*"
            else [factors]
        )
        for term in expanded:
            if term not in terms:
                terms.append(term)
    if audit["errors"]:
        return audit

    variables = list(dict.fromkeys(factor for term in terms for factor in term))
    audit["terms"] = [":".join(term) for term in terms]
    audit["variables"] = variables
    missing = [factor for factor in variables if factor not in obs.columns]
    if missing:
        audit["errors"].append("Unknown obs annotation(s): " + ", ".join(missing) + ".")
        return audit
    numeric = [factor for factor in variables if pd.api.types.is_numeric_dtype(obs[factor])]
    if numeric:
        audit["errors"].append(
            "Numeric covariates are not supported in Complex design yet: " + ", ".join(numeric) + "."
        )
        return audit

    replicate = str(pseudobulk_replicate_annotation or "").strip()
    if replicate and replicate not in obs.columns:
        audit["errors"].append(
            f"Pseudobulk replicate annotation '{replicate}' is not an obs annotation."
        )
        return audit
    sample_variables = list(variables)
    if replicate and replicate not in sample_variables:
        sample_variables.insert(0, replicate)
    audit["sample_variables"] = sample_variables
    audit["pseudobulk_replicate_annotation"] = replicate or None

    model_obs = obs[sample_variables].dropna().copy()
    if model_obs.empty:
        audit["errors"].append("No cells have complete values for every model annotation.")
        return audit
    for factor in variables:
        model_obs[factor] = model_obs[factor].astype(str)
    levels = {factor: sorted(model_obs[factor].unique().tolist()) for factor in variables}
    audit["levels"] = levels
    single_level = [factor for factor, values in levels.items() if len(values) < 2]
    if single_level:
        audit["errors"].append("Annotation(s) have fewer than two levels: " + ", ".join(single_level) + ".")
        return audit

    grouped = model_obs.groupby(sample_variables, observed=True, sort=True).size().reset_index(name="n_cells")
    audit["sample_count"] = int(len(grouped))
    audit["sample_preview"] = [
        {"values": {factor: str(row[factor]) for factor in sample_variables}, "n_cells": int(row["n_cells"])}
        for _, row in grouped.head(200).iterrows()
    ]
    if len(grouped) > 200:
        audit["warnings"].append("Only the first 200 pseudobulk samples are shown in the balance preview.")

    basis: Dict[str, np.ndarray] = {}
    for factor in variables:
        one_hot = pd.get_dummies(grouped[factor], drop_first=True, dtype=float).to_numpy(dtype=float)
        basis[factor] = one_hot
    columns = [np.ones(len(grouped), dtype=float)]
    for term in terms:
        arrays = [basis[factor] for factor in term]
        if any(array.shape[1] == 0 for array in arrays):
            continue
        for chosen in product(*(range(array.shape[1]) for array in arrays)):
            value = np.ones(len(grouped), dtype=float)
            for array, index in zip(arrays, chosen):
                value *= array[:, index]
            columns.append(value)
    matrix = np.column_stack(columns)
    audit["design_columns"] = int(matrix.shape[1])
    audit["design_rank"] = int(np.linalg.matrix_rank(matrix))
    if audit["design_rank"] < audit["design_columns"]:
        audit["errors"].append(
            "The design is rank-deficient: one or more model effects are confounded."
        )
        return audit
    audit["candidate_contrasts"] = [":".join(term) for term in terms if len(term) > 1]
    audit["valid"] = True
    return audit


def _rgb_nums_to_hex(nums) -> Optional[str]:
    """Convert a 3/4-number RGB(A) sequence to a ``#RRGGBB`` hex string.

    Accepts matplotlib-style 0..1 floats or 0..255 ints. Returns ``None`` if the
    sequence isn't a usable RGB triple.
    """
    try:
        vals = [float(v) for v in nums]
    except (TypeError, ValueError):
        return None
    if len(vals) < 3:
        return None
    r, g, b = vals[0], vals[1], vals[2]
    # Heuristic: matplotlib RGBA are 0..1 floats; scale those to 0..255.
    if max(r, g, b) <= 1.0:
        r, g, b = r * 255.0, g * 255.0, b * 255.0
    clamp = lambda x: max(0, min(255, int(round(x))))
    return "#%02X%02X%02X" % (clamp(r), clamp(g), clamp(b))


def _to_css_color(c) -> Optional[str]:
    """Normalize one ``uns[*_colors]`` entry to a CSS-usable color string.

    Handles hex/``rgb()``/named strings (passed through), bytes, numpy RGB(A)
    float/int arrays, and stringified arrays like ``"[0.88 0.47 0.72]"`` (which
    is what ``str(np.array([...]))`` produces — the source of cells rendering
    with no color). Returns ``None`` if it can't be made into a color.
    """
    if isinstance(c, bytes):
        c = c.decode("utf-8", "replace")
    if isinstance(c, np.str_):
        c = str(c)
    if isinstance(c, str):
        s = c.strip()
        if s.startswith("[") and s.endswith("]"):
            parts = s[1:-1].replace(",", " ").split()
            hexed = _rgb_nums_to_hex(parts)
            return hexed if hexed is not None else None
        return s or None
    # numpy scalar / 0-d
    if isinstance(c, np.generic):
        return str(c)
    # array / sequence of channel values
    return _rgb_nums_to_hex(c)


def _normalize_uns_palette(uns_palette) -> Optional[List[str]]:
    """Convert a full ``uns[*_colors]`` palette to CSS strings.

    Returns ``None`` if any entry can't be converted, so the caller can fall back
    to the viewer's default palette rather than emit invalid colors.
    """
    out = []
    for c in uns_palette:
        css = _to_css_color(c)
        if css is None:
            return None
        out.append(css)
    return out


def _resolve_spatial_key(adata, spatial_key: str = "spatial") -> str:
    """Return a usable spatial-coordinate key in ``adata.obsm``.

    Uses ``spatial_key`` when present; otherwise falls back to common naming
    conventions (``X_spatial``, ``spatial_coords``, ...). Raises a clear error
    listing the available keys when nothing usable is found.
    """
    def _is_usable(key: str) -> bool:
        if key not in adata.obsm:
            return False
        arr = np.asarray(adata.obsm[key])
        return arr.ndim == 2 and arr.shape[1] >= 2

    if spatial_key in adata.obsm:
        if not _is_usable(spatial_key):
            raise ValueError(
                f"adata.obsm['{spatial_key}'] is not a 2D array with at least two columns."
            )
        return spatial_key

    for candidate in _SPATIAL_KEY_FALLBACKS:
        if candidate != spatial_key and _is_usable(candidate):
            print(
                f"  Note: '{spatial_key}' not in adata.obsm; using '{candidate}' for spatial coordinates."
            )
            return candidate

    available = ", ".join(sorted(adata.obsm.keys())) or "(none)"
    raise ValueError(
        f"Spatial coordinates not found in adata.obsm['{spatial_key}']. "
        f"Available obsm keys: {available}. "
        f"Pass spatial_key=... to select one."
    )


def _set_spatial_from_obs_columns(
    adata,
    spatial_columns: Tuple[str, str],
    spatial_key: str = "spatial",
) -> str:
    """Create ``adata.obsm[spatial_key]`` from two numeric ``adata.obs`` columns."""
    if len(spatial_columns) != 2:
        raise ValueError("spatial_columns must contain exactly two obs column names")

    x_col, y_col = [str(col).strip() for col in spatial_columns]
    if not x_col or not y_col:
        raise ValueError("spatial_columns entries must be non-empty obs column names")

    missing = [col for col in (x_col, y_col) if col not in adata.obs.columns]
    if missing:
        available = ", ".join(map(str, adata.obs.columns)) or "(none)"
        raise ValueError(
            f"Spatial coordinate column(s) not found in adata.obs: {', '.join(missing)}. "
            f"Available obs columns: {available}."
        )

    coords = []
    for col in (x_col, y_col):
        values = pd.to_numeric(adata.obs[col], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(
                f"Spatial coordinate column '{col}' must contain only finite numeric values."
            )
        coords.append(values)

    adata.obsm[spatial_key] = np.column_stack(coords)
    print(
        f"  Using obs columns '{x_col}' and '{y_col}' as adata.obsm['{spatial_key}'] spatial coordinates."
    )
    return spatial_key


def _numeric_category_perm(categories: List) -> Optional[List[int]]:
    """Return a permutation ordering ``categories`` numerically, or ``None``.

    Many cluster annotations are stored as string categories ("0", "1", "10", ...),
    which pandas orders lexicographically so "10" sorts before "2". When *every*
    label parses as a number, return ``order`` such that ``order[new_idx] = old_idx``
    gives an ascending-numeric ordering. Returns ``None`` for non-numeric labels
    (preserving any intentional categorical order) or when already sorted.
    """
    strs = [str(c) for c in categories]
    nums: List[float] = []
    for s in strs:
        try:
            nums.append(float(s))
        except (TypeError, ValueError):
            return None
    if len(nums) <= 1:
        return None
    order = sorted(range(len(strs)), key=lambda i: (nums[i], strs[i]))
    if order == list(range(len(strs))):
        return None
    return order


def _compute_positive_fraction(
    matrix,
    mask: np.ndarray,
    gene_positions: List[Optional[int]],
) -> List[Optional[float]]:
    out: List[Optional[float]] = [None] * len(gene_positions)
    if matrix is None:
        return out
    mask = np.asarray(mask, dtype=bool)
    if matrix.shape[0] == 0 or not mask.any():
        return [0.0 if pos is not None and pos >= 0 else None for pos in gene_positions]

    valid = [(idx, pos) for idx, pos in enumerate(gene_positions) if pos is not None and pos >= 0]
    if not valid:
        return out

    valid_indices = [int(pos) for _, pos in valid]
    subset = matrix[mask][:, valid_indices]
    if issparse(subset):
        pct = np.asarray((subset > 0).mean(axis=0)).ravel()
    else:
        pct = np.mean(np.asarray(subset) > 0, axis=0)
    for (out_idx, _), value in zip(valid, np.asarray(pct).ravel()):
        out[out_idx] = float(value)
    return out


def _format_h5ad_removed_paths(paths: List[str], limit: int = 8) -> str:
    if len(paths) <= limit:
        return ", ".join(paths)
    return ", ".join(paths[:limit]) + f", ... ({len(paths) - limit} more)"


def _h5ad_attr_text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    if isinstance(value, np.bytes_):
        return value.item().decode("utf-8", "replace")
    return str(value)


def _copy_h5ad_for_repair(src_path: str) -> str:
    fd, tmp_path = tempfile.mkstemp(suffix=".h5ad")
    os.close(fd)
    shutil.copy2(src_path, tmp_path)
    return tmp_path


def _strip_null_encoded_h5ad_entries(src_path: str) -> Tuple[str, List[str]]:
    """Copy an h5ad file and remove null-encoded optional metadata entries."""
    import h5py

    tmp_path = _copy_h5ad_for_repair(src_path)

    removed_paths: List[str] = []
    with h5py.File(tmp_path, "r+") as handle:
        if "uns" not in handle:
            return tmp_path, removed_paths

        def _walk(group, prefix: str = "") -> None:
            for key in list(group.keys()):
                obj = group[key]
                path = f"{prefix}/{key}" if prefix else f"/{key}"
                if _h5ad_attr_text(obj.attrs.get("encoding-type", "")).lower() == "null":
                    removed_paths.append(path)
                    del group[key]
                    continue
                if isinstance(obj, h5py.Group):
                    _walk(obj, path)

        _walk(handle["uns"], "/uns")

    return tmp_path, removed_paths


def _strip_optional_uns_h5ad_entries(src_path: str) -> Tuple[str, List[str]]:
    """Copy an h5ad file and clear optional ``uns`` metadata for legacy read retries."""
    import h5py

    tmp_path = _copy_h5ad_for_repair(src_path)

    removed_paths: List[str] = []
    with h5py.File(tmp_path, "r+") as handle:
        if "uns" not in handle:
            return tmp_path, removed_paths

        removed_paths = [f"/uns/{key}" for key in handle["uns"].keys()]
        del handle["uns"]
        uns = handle.create_group("uns")
        uns.attrs["encoding-type"] = "dict"
        uns.attrs["encoding-version"] = "0.1.0"

    return tmp_path, removed_paths


def _is_h5ad_encoding_error(exc: Exception) -> bool:
    messages = []
    current: Optional[BaseException] = exc
    while current is not None:
        messages.append(str(current))
        current = current.__cause__ or current.__context__
    text = "\n".join(messages)
    return any(
        needle in text
        for needle in (
            "No read method registered for IOSpec",
            "IOSpec(",
            "encoding_type=",
            "encoding-type",
        )
    )


def _read_h5ad_with_fallback(path: str) -> sc.AnnData:
    """Read h5ad, retrying with legacy optional metadata stripped if needed."""
    try:
        return sc.read_h5ad(path)
    except PermissionError as exc:
        message = f"Unable to read h5ad file: {path}."
        if str(path).startswith("/Volumes/"):
            message += (
                " macOS denied access to the mounted volume. Grant your terminal or IDE access "
                "to external/removable volumes or Full Disk Access, or copy the file to a local "
                "path such as ~/Downloads and retry."
            )
        else:
            message += (
                " The OS denied access. Check file permissions for the current process or copy "
                "the file to a readable local path and retry."
            )
        raise PermissionError(message) from exc
    except Exception as exc:
        if not _is_h5ad_encoding_error(exc):
            raise

        repair_attempts = []
        if "encoding_type='null'" in str(exc) or 'encoding_type="null"' in str(exc):
            repair_attempts.append(
                (
                    "unsupported null-encoded H5AD metadata",
                    _strip_null_encoded_h5ad_entries,
                    "null metadata field(s)",
                )
            )
        repair_attempts.append(
            (
                "legacy H5AD optional metadata encoding",
                _strip_optional_uns_h5ad_entries,
                "optional uns entry/entries",
            )
        )

        for description, repair, label in repair_attempts:
            temp_path = None
            try:
                temp_path, removed_paths = repair(path)
                if not removed_paths:
                    continue
                print(f"  Detected {description}; retrying with a sanitized temporary copy...")
                print(f"  Removed {len(removed_paths)} {label}: {_format_h5ad_removed_paths(removed_paths)}")
                return sc.read_h5ad(temp_path)
            except Exception:
                continue
            finally:
                if temp_path and os.path.exists(temp_path):
                    os.remove(temp_path)

        raise


def _looks_like_spatialdata(obj: Any) -> bool:
    """Return whether ``obj`` behaves like a SpatialData container."""
    return hasattr(obj, "tables") and not isinstance(obj, sc.AnnData)


def _read_spatialdata_zarr(path: str) -> Any:
    """Read a SpatialData Zarr store, keeping spatialdata as an optional dependency."""
    try:
        import spatialdata as sd
    except ImportError as exc:
        raise ImportError(
            "Reading SpatialData input requires the optional 'spatialdata' package. "
            "Install it with `pip install spatialdata` or `pip install karospace[spatialdata]`."
        ) from exc

    try:
        return sd.read_zarr(path)
    except AttributeError:
        from spatialdata import read_zarr

        return read_zarr(path)


def _coerce_spatialdata_table(
    sdata: Any,
    spatialdata_table: Optional[str] = None,
) -> Tuple[sc.AnnData, str]:
    tables = getattr(sdata, "tables", None)
    if tables is None or not hasattr(tables, "keys"):
        raise ValueError("SpatialData input does not expose a tables mapping.")

    table_keys = list(tables.keys())
    if not table_keys:
        raise ValueError("SpatialData input contains no AnnData tables.")

    chosen_key: Any
    if spatialdata_table:
        matches = [key for key in table_keys if str(key) == str(spatialdata_table)]
        if not matches:
            available = ", ".join(map(str, table_keys))
            raise ValueError(
                f"SpatialData table '{spatialdata_table}' not found. Available tables: {available}"
            )
        chosen_key = matches[0]
    elif len(table_keys) == 1:
        chosen_key = table_keys[0]
    elif "table" in table_keys:
        chosen_key = "table"
    else:
        available = ", ".join(map(str, table_keys))
        raise ValueError(
            "SpatialData input contains multiple AnnData tables. "
            f"Pass spatialdata_table=... or --spatialdata-table. Available tables: {available}"
        )

    table = tables[chosen_key]
    if not isinstance(table, sc.AnnData):
        if hasattr(table, "to_adata"):
            table = table.to_adata()
        else:
            raise TypeError(
                f"SpatialData table '{chosen_key}' is not an AnnData object and cannot be converted."
            )
    return table, str(chosen_key)


def _coerce_input_to_anndata(
    data: Any,
    spatialdata_table: Optional[str] = None,
) -> Tuple[sc.AnnData, str, Optional[str]]:
    if isinstance(data, sc.AnnData):
        return data, "AnnData object", None

    if _looks_like_spatialdata(data):
        adata, table_key = _coerce_spatialdata_table(data, spatialdata_table)
        return adata, f"SpatialData object table '{table_key}'", table_key

    if isinstance(data, (str, os.PathLike)):
        path = str(data)
        path_obj = os.fspath(data)
        lower = path.lower()
        if lower.endswith(".zarr") or os.path.isdir(path_obj):
            sdata = _read_spatialdata_zarr(path)
            adata, table_key = _coerce_spatialdata_table(sdata, spatialdata_table)
            return adata, f"SpatialData store {path} table '{table_key}'", table_key
        return _read_h5ad_with_fallback(path), path, None

    raise TypeError(
        "Input must be a path to an .h5ad file, a SpatialData .zarr store, "
        "an AnnData object, or a SpatialData object."
    )


def _normalize_spatialdata_region_list(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, np.str_)):
        text = value.decode("utf-8") if isinstance(value, bytes) else str(value)
        return [text]
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, (pd.Index, list, tuple, set)):
        return [str(item) for item in value if item is not None]
    return [str(value)]


def _get_spatialdata_attrs(adata: sc.AnnData) -> Dict[str, Any]:
    attrs = adata.uns.get("spatialdata_attrs", {})
    if hasattr(attrs, "get"):
        return attrs
    try:
        return dict(attrs)
    except Exception:
        return {}


def _resolve_groupby_for_spatialdata(
    adata: sc.AnnData,
    groupby: str,
    spatialdata_table_key: Optional[str],
) -> str:
    if spatialdata_table_key is None or groupby in adata.obs.columns:
        return groupby

    attrs = _get_spatialdata_attrs(adata)
    region_key_raw = attrs.get("region_key")
    region_key = str(region_key_raw).strip() if region_key_raw is not None else ""
    if groupby == "sample_id" and region_key and region_key in adata.obs.columns:
        print(
            f"  Note: groupby 'sample_id' not found; using SpatialData region_key "
            f"'{region_key}' for sections."
        )
        return region_key

    if groupby == "sample_id":
        regions = _normalize_spatialdata_region_list(attrs.get("region"))
        fallback_col = "_karospace_spatialdata_region"
        if len(regions) == 1:
            value = regions[0]
        else:
            value = str(spatialdata_table_key)
            if regions:
                print(
                    "  Warning: SpatialData table has multiple regions but no per-cell region_key; "
                    f"exporting it as one section named '{value}'."
                )
        adata.obs[fallback_col] = pd.Categorical([value] * adata.n_obs)
        print(
            f"  Note: groupby 'sample_id' not found; using '{fallback_col}' for one SpatialData section."
        )
        return fallback_col

    return groupby


def inspect_input_file(data: Any, spatialdata_table: Optional[str] = None) -> Dict[str, Any]:
    """Read an input object/file and summarize its available cell metadata.

    This intentionally performs no coordinate validation, section construction,
    downsampling, or analytical calculation. It is suitable for choosing CLI
    metadata, annotation, and pseudobulk design arguments before export.
    """
    max_examples = 10
    adata, source_label, table_key = _coerce_input_to_anndata(data, spatialdata_table)
    metadata: List[Dict[str, Any]] = []
    for column in adata.obs.columns:
        series = adata.obs[column]
        observed = series.dropna()
        examples = [str(value) for value in observed.drop_duplicates().head(max_examples).tolist()]
        if isinstance(series.dtype, CategoricalDtype):
            value_type = "categorical"
        elif pd.api.types.is_bool_dtype(series):
            value_type = "boolean"
        elif pd.api.types.is_numeric_dtype(series):
            value_type = "numeric"
        else:
            value_type = "text"
        metadata.append(
            {
                "name": str(column),
                "type": value_type,
                "dtype": str(series.dtype),
                "n_unique": int(observed.nunique()),
                "n_missing": int(series.isna().sum()),
                "examples": examples,
            }
        )
    return {
        "path": str(source_label),
        "spatialdata_table": table_key,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "metadata": metadata,
    }


def _normalize_uns_scalar(value: Any) -> Any:
    """Normalize AnnData uns scalars/arrays to plain Python values when possible."""
    while isinstance(value, np.ndarray) and value.ndim == 0:
        value = value.item()
    if isinstance(value, np.ndarray):
        if value.size == 1:
            return value.reshape(()).item()
        return value
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, memoryview):
        return value.tobytes()
    if isinstance(value, bytearray):
        return bytes(value)
    return value


def _normalize_uns_text(value: Any) -> Optional[str]:
    value = _normalize_uns_scalar(value)
    if value is None:
        return None
    if isinstance(value, np.ndarray):
        if value.size != 1:
            return None
        value = value.reshape(()).item()
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _normalize_uns_text_list(value: Any) -> List[str]:
    value = _normalize_uns_scalar(value)
    if value is None:
        return []
    if isinstance(value, np.ndarray):
        items = value.ravel().tolist()
    elif isinstance(value, (list, tuple)):
        items = list(value)
    else:
        items = [value]

    normalized: List[str] = []
    for item in items:
        text = _normalize_uns_text(item)
        if text is None:
            continue
        text = text.strip()
        if text:
            normalized.append(text)
    return normalized


def _load_companion_analytics(adata) -> Dict[str, Any]:
    """Load precomputed viewer analytics emitted by KaroSpaceCompanion from adata.uns."""
    uns = getattr(adata, "uns", None)
    if uns is None or "karospace_companion" not in uns:
        return {}

    companion = uns.get("karospace_companion")
    if companion is None:
        return {}
    if not hasattr(companion, "get"):
        try:
            companion = dict(companion)
        except Exception:
            return {}

    storage = _normalize_uns_text(companion.get("analytics_storage"))
    if storage != COMPANION_ANALYTICS_STORAGE:
        return {}

    analytics: Dict[str, Any] = {}
    if "analytics_columns" in companion:
        analytics["analytics_columns"] = _normalize_uns_text_list(companion.get("analytics_columns"))

    for uns_key, analytics_key in COMPANION_ANALYTICS_JSON_FIELDS.items():
        if uns_key not in companion:
            continue
        raw = _normalize_uns_text(companion.get(uns_key))
        if raw is None or not raw.strip():
            continue
        try:
            analytics[analytics_key] = json.loads(raw)
        except Exception as exc:
            print(
                f"  Warning: Could not parse KaroSpaceCompanion analytics field "
                f"'{uns_key}': {exc}"
            )

    return analytics


def _strip_category_pseudobulk_sample_diagnostics(payload: Any) -> Any:
    """Drop legacy all-category PCA/distance diagnostics from category DE payloads.

    Pairwise diagnostics live under ``_summary.pair_diagnostics`` and are kept.
    Contact-marker diagnostics are handled separately in the interaction marker
    payload.
    """
    if not isinstance(payload, dict):
        return payload
    cleaned = deepcopy(payload)

    def _walk(node: Any) -> None:
        if not isinstance(node, dict):
            return
        node.pop("pseudobulk_samples", None)
        for value in node.values():
            if isinstance(value, dict):
                _walk(value)
            elif isinstance(value, list):
                for item in value:
                    _walk(item)

    _walk(cleaned)
    return cleaned


@dataclass
class SectionData:
    """Data for a single tissue section."""
    section_id: str
    coordinates: np.ndarray  # (n_cells, 2) array of x, y coordinates
    metadata: Dict[str, str] = field(default_factory=dict)

    @property
    def n_cells(self) -> int:
        return self.coordinates.shape[0]

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """Return (xmin, xmax, ymin, ymax)."""
        return (
            float(self.coordinates[:, 0].min()),
            float(self.coordinates[:, 0].max()),
            float(self.coordinates[:, 1].min()),
            float(self.coordinates[:, 1].max()),
        )


@dataclass
class Modality:
    """One modality (RNA, protein, etc.) — n_cells × n_features matrix + var."""
    name: str
    matrix: Any
    var: pd.DataFrame
    layers: Dict[str, Any] = field(default_factory=dict)
    value_kind: str = "counts"  # 'counts' (RNA-like) or 'intensity' (protein-like)
    label: Optional[str] = None

    @property
    def feature_names(self) -> List[str]:
        return list(self.var.index)

    @property
    def n_features(self) -> int:
        return int(self.matrix.shape[1])

    def feature_position(self, name: str) -> int:
        return self.feature_names.index(name)

    def get_feature_vector(self, name: str) -> np.ndarray:
        idx = self.feature_position(name)
        source = self.layers.get("normalized")
        if source is None:
            source = self.matrix
        col = source[:, idx]
        if issparse(col):
            return np.asarray(col.toarray()).ravel()
        return np.asarray(col).ravel()


@dataclass
class SpatialDataset:
    """Container for spatial transcriptomics dataset."""
    adata: sc.AnnData
    sections: List[SectionData]
    groupby: str
    obs_columns: List[str]
    var_names: List[str]
    metadata_section: List[str] = field(default_factory=list)
    metadata_section_extra: List[str] = field(default_factory=list)
    metadata_value_order: Optional[Dict[str, List[str]]] = None
    modalities: Dict[str, Modality] = field(default_factory=dict)
    default_modality: str = "rna"
    spatial_key: str = "spatial"

    def __post_init__(self) -> None:
        def _clean(values: Optional[List[str]]) -> List[str]:
            cleaned = []
            for value in values or []:
                if value is None:
                    continue
                text = str(value).strip()
                if text:
                    cleaned.append(text)
            return list(dict.fromkeys(cleaned))

        if self.metadata_section:
            self.metadata_section = _clean(self.metadata_section)
        else:
            self.metadata_section = []
        self.metadata_section_extra = _clean(self.metadata_section_extra)

    @property
    def n_sections(self) -> int:
        return len(self.sections)

    @property
    def n_cells(self) -> int:
        return self.adata.n_obs

    @property
    def has_umap(self) -> bool:
        """Check if UMAP coordinates are available."""
        return "X_umap" in self.adata.obsm

    def get_companion_analytics(self) -> Dict[str, Any]:
        """Return normalized KaroSpaceCompanion analytics embedded in adata.uns."""
        return _load_companion_analytics(self.adata)

    def _resolve_modality(self, modality: Optional[str]) -> Optional[Modality]:
        if not self.modalities:
            return None
        name = modality or self.default_modality
        return self.modalities.get(name)

    def get_annotation_data(
        self,
        annotation: str,
        modality: Optional[str] = None,
    ) -> Tuple[np.ndarray, bool, Optional[List[str]]]:
        """
        Get annotation or gene values for all cells.

        Parameters
        ----------
        annotation : str
            Column in obs or feature name within the active modality
        modality : str, optional
            Modality name to look up the feature in. Defaults to ``default_modality``.

        Returns
        -------
        values : np.ndarray
            Numeric values for each cell
        is_continuous : bool
            Whether data is continuous
        categories : list or None
            Category names if categorical, else None
        """
        if annotation in self.adata.obs.columns:
            col = self.adata.obs[annotation]
            if isinstance(col.dtype, CategoricalDtype):
                categories = list(col.cat.categories)
                values = col.cat.codes.to_numpy().astype(float)
                # Handle NaN codes (-1)
                values[values < 0] = np.nan
                return values, False, categories
            elif pd.api.types.is_numeric_dtype(col):
                values = col.to_numpy(dtype=float)
                return values, True, None
            else:
                # Convert to categorical
                cat = col.astype("category")
                categories = list(cat.cat.categories)
                values = cat.cat.codes.to_numpy().astype(float)
                values[values < 0] = np.nan
                return values, False, categories

        mod = self._resolve_modality(modality)
        if mod is not None and annotation in mod.feature_names:
            return mod.get_feature_vector(annotation), True, None

        if annotation in self.adata.var_names:
            # Back-compat fallback when modality registry is unpopulated.
            gene_idx = self.adata.var_names.get_loc(annotation)
            expr_layer = None
            if "normalized" in self.adata.layers:
                expr_layer = self.adata.layers["normalized"]
            x = expr_layer[:, gene_idx] if expr_layer is not None else self.adata.X[:, gene_idx]
            if issparse(x):
                values = np.asarray(x.toarray()).ravel()
            else:
                values = np.asarray(x).ravel()
            return values, True, None

        raise KeyError(f"{annotation!r} not found in obs columns or modality {modality or self.default_modality!r}")

    def get_section_indices(self) -> Dict[str, np.ndarray]:
        """Get cell indices for each section."""
        indices = {}
        gvals = self.adata.obs[self.groupby].astype(str).to_numpy()
        for section in self.sections:
            indices[section.section_id] = np.flatnonzero(gvals == section.section_id)
        return indices

    def get_metadata_filters(self) -> Dict[str, List[str]]:
        """Get unique values for filterable metadata columns."""
        filters = {}
        for col in self.metadata_section:
            if col in self.adata.obs.columns:
                unique_vals = list(self.adata.obs[col].dropna().astype(str).unique())
                custom_order = None
                if self.metadata_value_order and col in self.metadata_value_order:
                    custom_order = [str(v) for v in self.metadata_value_order[col]]
                if custom_order:
                    custom_set = set(custom_order)
                    ordered = [v for v in custom_order if v in unique_vals]
                    remaining = [v for v in unique_vals if v not in custom_set]
                    if col == "last_day":
                        def _sort_key(v):
                            try:
                                return (0, float(v))
                            except ValueError:
                                return (1, v)
                        remaining = sorted(remaining, key=_sort_key)
                    else:
                        remaining = sorted(remaining)
                    filters[col] = ordered + remaining
                elif col == "last_day":
                    def _sort_key(v):
                        try:
                            return (0, float(v))
                        except ValueError:
                            return (1, v)
                    filters[col] = sorted(unique_vals, key=_sort_key)
                else:
                    filters[col] = sorted(unique_vals)
        return filters

    def _get_export_section_indices(self, downsample: Optional[int] = None) -> Dict[str, np.ndarray]:
        """Return per-section obs indices used for export, including deterministic downsampling."""
        section_indices = self.get_section_indices()
        export_indices: Dict[str, np.ndarray] = {}
        for section in self.sections:
            idx = section_indices[section.section_id]
            if downsample and len(idx) > downsample:
                rng = np.random.default_rng(42)
                idx = rng.choice(idx, size=downsample, replace=False)
                idx = np.sort(idx)
            export_indices[section.section_id] = np.asarray(idx)
        return export_indices

    def _collect_gene_data(
        self,
        genes: Optional[List[str]] = None,
        modality: Optional[str] = None,
    ) -> Dict[str, Dict[str, Union[np.ndarray, float]]]:
        """Collect feature vectors and ranges for export within a modality."""
        gene_data: Dict[str, Dict[str, Union[np.ndarray, float]]] = {}
        mod = self._resolve_modality(modality)
        feature_set = set(mod.feature_names) if mod is not None else set(self.adata.var_names)
        for gene in genes or []:
            if gene not in feature_set:
                continue
            try:
                vals, _, _ = self.get_annotation_data(gene, modality=modality)
                finite = np.isfinite(vals)
                gene_vmin = float(np.nanmin(vals[finite])) if finite.any() else 0.0
                gene_vmax = float(np.nanmax(vals[finite])) if finite.any() else 1.0
                gene_data[gene] = {
                    "values": vals,
                    "vmin": gene_vmin,
                    "vmax": gene_vmax,
                }
            except Exception as e:
                print(f"  Warning: Could not load feature '{gene}': {e}")
        return gene_data

    @staticmethod
    def _resolve_gene_encodings(
        gene_data: Dict[str, Dict[str, Union[np.ndarray, float]]],
        gene_encoding: str,
        gene_sparse_zero_threshold: float,
    ) -> Dict[str, str]:
        gene_encodings: Dict[str, str] = {}
        for gene, gdata in gene_data.items():
            if gene_encoding == "dense":
                gene_encodings[gene] = "dense"
            elif gene_encoding == "sparse":
                gene_encodings[gene] = "sparse"
            else:
                vals = np.asarray(gdata["values"])
                finite = np.isfinite(vals)
                nonzero = finite & (vals != 0)
                zero_frac = 1.0
                if vals.size:
                    zero_frac = 1.0 - (float(np.count_nonzero(nonzero)) / float(vals.size))
                gene_encodings[gene] = "sparse" if zero_frac >= float(gene_sparse_zero_threshold) else "dense"
        return gene_encodings

    @staticmethod
    def _validate_gene_export_options(
        gene_encoding: str,
        gene_sparse_zero_threshold: float,
        gene_sparse_pack_min_nnz: int,
    ) -> None:
        gene_encoding = str(gene_encoding or "auto").lower()
        if gene_encoding not in {"auto", "dense", "sparse"}:
            raise ValueError("gene_encoding must be one of: 'auto', 'dense', 'sparse'")
        if not (0.0 <= float(gene_sparse_zero_threshold) <= 1.0):
            raise ValueError("gene_sparse_zero_threshold must be between 0 and 1")
        if int(gene_sparse_pack_min_nnz) < 0:
            raise ValueError("gene_sparse_pack_min_nnz must be >= 0")

    @staticmethod
    def _validate_gene_value_encoding(gene_value_encoding: str) -> str:
        normalized = str(gene_value_encoding or "uint16").lower()
        if normalized not in {"uint16", "uint8"}:
            raise ValueError("gene_value_encoding must be one of: 'uint16', 'uint8'")
        return normalized

    @staticmethod
    def _quantize_to_uint(values: np.ndarray, vmin: float, vmax: float, max_code: int, dtype) -> np.ndarray:
        arr = np.asarray(values, dtype=np.float32)
        finite = np.isfinite(arr)
        out = np.zeros(arr.shape, dtype=dtype)
        if not finite.any():
            return out
        if not np.isfinite(vmin):
            vmin = float(np.nanmin(arr[finite])) if finite.any() else 0.0
        if not np.isfinite(vmax):
            vmax = float(np.nanmax(arr[finite])) if finite.any() else vmin
        if vmax <= vmin:
            return out
        max_code = int(max_code)
        scale = np.float32(max_code / (float(vmax) - float(vmin)))
        clipped = np.clip(np.rint((arr[finite] - np.float32(vmin)) * scale), 0, max_code)
        out[finite] = clipped.astype(dtype, copy=False)
        return out

    @classmethod
    def _quantize_values(cls, values: np.ndarray, vmin: float, vmax: float, value_encoding: str) -> np.ndarray:
        if value_encoding == "uint16":
            return cls._quantize_to_uint(values, vmin, vmax, 65535, np.uint16)
        if value_encoding == "uint8":
            return cls._quantize_to_uint(values, vmin, vmax, 255, np.uint8)
        raise ValueError(f"Unsupported quantized gene value encoding: {value_encoding}")

    @staticmethod
    def _estimate_packed_dense_bytes(n_values: int, nan_count: int = 0, value_bytes: int = 4) -> int:
        return max(0, int(n_values)) * max(1, int(value_bytes)) + max(0, int(nan_count)) * 4

    @staticmethod
    def _estimate_packed_sparse_bytes(nnz: int, nan_count: int = 0, value_bytes: int = 4) -> int:
        return max(0, int(nnz)) * (4 + max(1, int(value_bytes))) + max(0, int(nan_count)) * 4

    @staticmethod
    def _gene_value_encoding_bytes(gene_value_encoding: str) -> int:
        if gene_value_encoding == "uint16":
            return 2
        if gene_value_encoding == "uint8":
            return 1
        return 4

    @classmethod
    def _resolve_sidecar_gene_encoding_mode(
        cls,
        resolved_gene_encoding: str,
        n_values: int,
        nonzero_count: int,
        nan_count: int,
        gene_sparse_zero_threshold: float,
        gene_value_encoding: str,
    ) -> str:
        if resolved_gene_encoding == "dense":
            return "dense"
        if resolved_gene_encoding == "sparse":
            return "sparse"

        zero_frac = 1.0
        if n_values:
            zero_frac = 1.0 - (float(nonzero_count) / float(n_values))

        # Sidecar payloads compare packed dense quantized buffers against packed
        # sparse index/value buffers. Prefer sparse whenever it is explicitly
        # requested via zero-threshold or whenever it is estimated to be smaller.
        value_bytes = cls._gene_value_encoding_bytes(gene_value_encoding)
        dense_bytes = cls._estimate_packed_dense_bytes(n_values, nan_count, value_bytes=value_bytes)
        sparse_bytes = cls._estimate_packed_sparse_bytes(nonzero_count, nan_count, value_bytes=value_bytes)
        if zero_frac >= float(gene_sparse_zero_threshold) or sparse_bytes <= dense_bytes:
            return "sparse"
        return "dense"

    @staticmethod
    def _serialize_gene_section_values(
        section_vals: np.ndarray,
        mode: str,
        gene_sparse_pack: bool,
        gene_sparse_pack_min_nnz: int,
        b64_encoder,
    ) -> Dict:
        if mode == "sparse":
            finite = np.isfinite(section_vals)
            nonzero = finite & (section_vals != 0)
            nz_idx = np.flatnonzero(nonzero).astype(np.uint32)
            nz_vals = np.asarray(section_vals[nonzero], dtype=np.float32)
            if bool(gene_sparse_pack) and int(nz_idx.size) >= int(gene_sparse_pack_min_nnz):
                sparse_entry = {
                    "ib64": b64_encoder(np.asarray(nz_idx, dtype="<u4")),
                    "vb64": b64_encoder(np.asarray(nz_vals, dtype="<f4")),
                }
            else:
                sparse_entry = {
                    "i": nz_idx.astype(int).tolist(),
                    "v": nz_vals.astype(float).tolist(),
                }
            nan_idx = np.flatnonzero(np.isnan(section_vals)).astype(int)
            if nan_idx.size:
                sparse_entry["nan"] = nan_idx.tolist()
            return {"sparse": sparse_entry}

        return {
            "dense": [float(v) if np.isfinite(v) else None for v in section_vals]
        }

    def to_gene_sidecar_data(
        self,
        genes: Optional[List[str]] = None,
        downsample: Optional[int] = None,
        export_indices: Optional[Dict[str, np.ndarray]] = None,
        gene_encoding: str = "auto",
        gene_value_encoding: str = "uint16",
        gene_sparse_zero_threshold: float = 0.8,
        gene_sparse_pack: bool = True,
        gene_sparse_pack_min_nnz: int = 256,
        modality: Optional[str] = None,
    ) -> Dict:
        """Export a gene-major sidecar payload for downstream viewer loading."""
        self._validate_gene_export_options(gene_encoding, gene_sparse_zero_threshold, gene_sparse_pack_min_nnz)
        resolved_value_encoding = self._validate_gene_value_encoding(gene_value_encoding)
        export_indices = export_indices or self._get_export_section_indices(downsample=downsample)

        mod = self._resolve_modality(modality)
        if mod is not None:
            feature_pos = {n: i for i, n in enumerate(mod.feature_names)}
            source_matrix = mod.layers.get("normalized", mod.matrix)
        else:
            feature_pos = {n: i for i, n in enumerate(self.adata.var_names)}
            source_matrix = self.adata.layers["normalized"] if "normalized" in self.adata.layers else self.adata.X

        unique_genes = []
        seen = set()
        for gene in genes or []:
            if gene in seen or gene not in feature_pos:
                continue
            seen.add(gene)
            unique_genes.append(gene)

        def _b64(arr: np.ndarray) -> str:
            import base64

            carr = np.ascontiguousarray(arr)
            return base64.b64encode(carr.tobytes(order="C")).decode("ascii")

        genes_payload = {}
        genes_meta = {}
        gene_encodings: Dict[str, str] = {}
        gene_value_encodings: Dict[str, str] = {}
        if unique_genes:
            gene_indices = [feature_pos[gene] for gene in unique_genes]
            expr_layer = source_matrix
            batch = expr_layer[:, gene_indices]
            resolved_gene_encoding = str(gene_encoding or "auto").lower()

            if issparse(batch):
                batch_csr = batch.tocsr()
                batch_csc = batch.tocsc()
                section_batches = {
                    section.section_id: batch_csr[export_indices[section.section_id], :].tocsc()
                    for section in self.sections
                }
                n_obs = int(batch.shape[0])
                for gene_pos, gene in enumerate(unique_genes):
                    col = batch_csc.getcol(gene_pos)
                    raw_vals = np.asarray(col.data, dtype=float).ravel()
                    finite_vals = raw_vals[np.isfinite(raw_vals)]
                    if finite_vals.size:
                        raw_min = float(np.min(finite_vals))
                        raw_max = float(np.max(finite_vals))
                        has_implicit_zeros = int(np.count_nonzero(col.data != 0)) < n_obs
                        gene_vmin = min(0.0, raw_min) if has_implicit_zeros else raw_min
                        gene_vmax = max(0.0, raw_max) if has_implicit_zeros else raw_max
                    else:
                        gene_vmin = 0.0
                        gene_vmax = 0.0

                    if resolved_gene_encoding == "dense":
                        mode = "dense"
                    elif resolved_gene_encoding == "sparse":
                        mode = "sparse"
                    else:
                        finite_mask = np.isfinite(raw_vals)
                        nonzero_mask = finite_mask & (raw_vals != 0)
                        mode = self._resolve_sidecar_gene_encoding_mode(
                            resolved_gene_encoding=resolved_gene_encoding,
                            n_values=n_obs,
                            nonzero_count=int(np.count_nonzero(nonzero_mask)),
                            nan_count=int(raw_vals.size - np.count_nonzero(finite_mask)),
                            gene_sparse_zero_threshold=gene_sparse_zero_threshold,
                            gene_value_encoding=resolved_value_encoding,
                        )
                    gene_encodings[gene] = mode
                    gene_value_encodings[gene] = resolved_value_encoding

                    sections_payload = {}
                    for section in self.sections:
                        sec_col = section_batches[section.section_id].getcol(gene_pos)
                        if mode == "sparse":
                            section_vals = np.asarray(sec_col.data, dtype=float).ravel()
                            section_idx = np.asarray(sec_col.indices, dtype=np.int64).ravel()
                            finite = np.isfinite(section_vals)
                            nonzero = finite & (section_vals != 0)
                            nz_idx = section_idx[nonzero].astype(np.uint32, copy=False)
                            nz_vals = section_vals[nonzero].astype(np.float32, copy=False)
                            if resolved_value_encoding in {"uint16", "uint8"}:
                                quantized = self._quantize_values(nz_vals, gene_vmin, gene_vmax, resolved_value_encoding)
                                sparse_entry = {
                                    "ib64": _b64(np.asarray(nz_idx, dtype="<u4")),
                                    ("vq16b64" if resolved_value_encoding == "uint16" else "vq8b64"): _b64(
                                        np.asarray(
                                            quantized,
                                            dtype="<u2" if resolved_value_encoding == "uint16" else "u1",
                                        )
                                    ),
                                }
                            else:
                                if bool(gene_sparse_pack) and int(nz_idx.size) >= int(gene_sparse_pack_min_nnz):
                                    sparse_entry = {
                                        "ib64": _b64(np.asarray(nz_idx, dtype="<u4")),
                                        "vb64": _b64(np.asarray(nz_vals, dtype="<f4")),
                                    }
                                else:
                                    sparse_entry = {
                                        "i": nz_idx.astype(int).tolist(),
                                        "v": nz_vals.astype(float).tolist(),
                                    }
                            nan_idx = section_idx[np.isnan(section_vals)].astype(int, copy=False)
                            if nan_idx.size:
                                sparse_entry["nan"] = nan_idx.tolist()
                            sections_payload[section.section_id] = {"sparse": sparse_entry}
                        else:
                            dense_vals = np.asarray(sec_col.toarray(), dtype=np.float32).ravel()
                            finite_mask = np.isfinite(dense_vals)
                            if resolved_value_encoding in {"uint16", "uint8"}:
                                quantized = self._quantize_values(
                                    np.where(finite_mask, dense_vals, 0.0),
                                    gene_vmin,
                                    gene_vmax,
                                    resolved_value_encoding,
                                )
                                dense_entry = {
                                    ("dq16b64" if resolved_value_encoding == "uint16" else "dq8b64"): _b64(
                                        np.asarray(
                                            quantized,
                                            dtype="<u2" if resolved_value_encoding == "uint16" else "u1",
                                        )
                                    )
                                }
                            else:
                                dense_entry = {
                                    "db64": _b64(
                                        np.asarray(
                                            np.where(finite_mask, dense_vals, 0.0),
                                            dtype="<f4",
                                        )
                                    )
                                }
                            if not finite_mask.all():
                                dense_entry["nan"] = np.flatnonzero(~finite_mask).astype(int).tolist()
                            sections_payload[section.section_id] = dense_entry

                    genes_payload[gene] = {"sections": sections_payload}
                    genes_meta[gene] = {"vmin": gene_vmin, "vmax": gene_vmax}
            else:
                batch_dense = np.asarray(batch)
                if batch_dense.ndim == 1:
                    batch_dense = batch_dense.reshape(-1, 1)
                for gene_pos, gene in enumerate(unique_genes):
                    vals = np.asarray(batch_dense[:, gene_pos], dtype=float).ravel()
                    finite = np.isfinite(vals)
                    gene_vmin = float(np.nanmin(vals[finite])) if finite.any() else 0.0
                    gene_vmax = float(np.nanmax(vals[finite])) if finite.any() else 0.0
                    if resolved_gene_encoding == "dense":
                        mode = "dense"
                    elif resolved_gene_encoding == "sparse":
                        mode = "sparse"
                    else:
                        nonzero = finite & (vals != 0)
                        mode = self._resolve_sidecar_gene_encoding_mode(
                            resolved_gene_encoding=resolved_gene_encoding,
                            n_values=int(vals.size),
                            nonzero_count=int(np.count_nonzero(nonzero)),
                            nan_count=int(vals.size - np.count_nonzero(finite)),
                            gene_sparse_zero_threshold=gene_sparse_zero_threshold,
                            gene_value_encoding=resolved_value_encoding,
                        )
                    gene_encodings[gene] = mode
                    gene_value_encodings[gene] = resolved_value_encoding

                    sections_payload = {}
                    for section in self.sections:
                        idx = export_indices[section.section_id]
                        section_vals = np.asarray(batch_dense[idx, gene_pos], dtype=np.float32).ravel()
                        if mode == "dense":
                            finite_mask = np.isfinite(section_vals)
                            if resolved_value_encoding in {"uint16", "uint8"}:
                                quantized = self._quantize_values(
                                    np.where(finite_mask, section_vals, 0.0),
                                    gene_vmin,
                                    gene_vmax,
                                    resolved_value_encoding,
                                )
                                dense_entry = {
                                    ("dq16b64" if resolved_value_encoding == "uint16" else "dq8b64"): _b64(
                                        np.asarray(
                                            quantized,
                                            dtype="<u2" if resolved_value_encoding == "uint16" else "u1",
                                        )
                                    )
                                }
                            else:
                                dense_entry = {
                                    "db64": _b64(
                                        np.asarray(
                                            np.where(finite_mask, section_vals, 0.0),
                                            dtype="<f4",
                                        )
                                    )
                                }
                            if not finite_mask.all():
                                dense_entry["nan"] = np.flatnonzero(~finite_mask).astype(int).tolist()
                            sections_payload[section.section_id] = dense_entry
                        else:
                            finite_vals = np.isfinite(section_vals)
                            nonzero = finite_vals & (section_vals != 0)
                            nz_idx = np.flatnonzero(nonzero).astype(np.uint32)
                            nz_vals = np.asarray(section_vals[nonzero], dtype=np.float32)
                            if resolved_value_encoding in {"uint16", "uint8"}:
                                quantized = self._quantize_values(nz_vals, gene_vmin, gene_vmax, resolved_value_encoding)
                                sparse_entry = {
                                    "ib64": _b64(np.asarray(nz_idx, dtype="<u4")),
                                    ("vq16b64" if resolved_value_encoding == "uint16" else "vq8b64"): _b64(
                                        np.asarray(
                                            quantized,
                                            dtype="<u2" if resolved_value_encoding == "uint16" else "u1",
                                        )
                                    ),
                                }
                            else:
                                sections_payload[section.section_id] = self._serialize_gene_section_values(
                                    section_vals=section_vals,
                                    mode=mode,
                                    gene_sparse_pack=gene_sparse_pack,
                                    gene_sparse_pack_min_nnz=gene_sparse_pack_min_nnz,
                                    b64_encoder=_b64,
                                )
                                continue
                            nan_idx = np.flatnonzero(np.isnan(section_vals)).astype(int)
                            if nan_idx.size:
                                sparse_entry["nan"] = nan_idx.tolist()
                            sections_payload[section.section_id] = {"sparse": sparse_entry}
                    genes_payload[gene] = {"sections": sections_payload}
                    genes_meta[gene] = {"vmin": gene_vmin, "vmax": gene_vmax}

        return {
            "format": "karospace-gene-sidecar-v1",
            "modality": (mod.name if mod is not None else "rna"),
            "genes_meta": genes_meta,
            "gene_encodings": gene_encodings,
            "gene_value_encodings": gene_value_encodings,
            "genes": genes_payload,
        }

    def to_json_data(
        self,
        annotation: str,
        downsample: Optional[int] = None,
        cells_annotations: Optional[List[str]] = None,
        genes: Optional[List[str]] = None,
        gene_encoding: str = "auto",
        gene_sparse_zero_threshold: float = 0.8,
        gene_sparse_pack: bool = True,
        gene_sparse_pack_min_nnz: int = 256,
        pseudobulk_de_groupby: Optional[List[str]] = None,
        pseudobulk_replicate_annotation: Optional[str] = None,
        pseudobulk_simple_constrast_categories: Any = None,
        pseudobulk_counts_layer: Optional[str] = "counts",
        pseudobulk_min_cell_counts: int = 0,
        pseudobulk_min_gene_counts: int = 0,
        pseudobulk_min_replicates: int = 2,
        pseudobulk_min_pct_expressed: float = 0.0,
        pseudobulk_p_adjust_method: str = "fdr_bh",
        pseudobulk_padj_cutoff: float = 0.05,
        pseudobulk_log2fc_cutoff: float = 0.5,
        pseudobulk_deseq2_fit_type: str = "parametric",
        pseudobulk_n_cpus: int = 1,
        pseudobulk_embed_top_n_per_comparison: int = 20,
        interaction_markers_groupby: Optional[List[str]] = None,
        neighbor_stats_groupby: Optional[List[str]] = None,
        neighbor_stats_permutations: int = 0,
        neighbor_stats_seed: int = 0,
        interaction_markers_top_targets: int = 8,
        interaction_markers_top_genes: int = 20,
        interaction_markers_min_cells: int = 30,
        interaction_markers_min_neighbors: int = 1,
        section_rotations: Optional[Dict[str, float]] = None,
        deconvolutions: Optional[Dict[str, str]] = None,
    ) -> Dict:
        """
        Export dataset to JSON-serializable format for the HTML viewer.

        Parameters
        ----------
        annotation : str
            Initial cell annotation column or gene
        downsample : int, optional
            If set, randomly downsample to this many cells per section
        cells_annotations : list, optional
            Additional cell obs columns to include for annotation switching.
        genes : list, optional
            Gene names to include for expression visualization
        gene_encoding : str
            "dense", "sparse", or "auto" (default: "auto"). When "sparse" (or when
            "auto" decides to use sparse), per-section gene vectors are stored as
            (index, value) pairs for non-zero entries to reduce HTML size for
            zero-inflated expression matrices.
        gene_sparse_zero_threshold : float
            Only used when gene_encoding="auto". Use sparse encoding when the
            fraction of zeros is >= this threshold (default: 0.8).
        gene_sparse_pack : bool
            When using sparse gene encoding, store indices/values as base64 typed arrays
            (smaller + faster JSON parse for large datasets). Default: True.
        gene_sparse_pack_min_nnz : int
            Only pack sparse arrays when non-zero entries in a section are >= this value.
            Default: 256.
        pseudobulk_de_groupby : list, optional
            Internal list of obs columns to compute shared-fit pseudobulk DE for.
            Exporter supplies the initial annotation and pseudobulk additional annotations.
        pseudobulk_replicate_annotation : str, optional
            Obs annotation used as the biological replicate for pseudobulk analyses.
            Defaults to the dataset groupby annotation.
        pseudobulk_simple_constrast_categories : list or dict, optional
            Categories to include in Simple design category-versus-category
            contrasts. Use a flat list only when one annotation is analyzed.
            With multiple pseudobulk annotation columns, pass a dict keyed by
            annotation name or a nested list matching the annotation order. All
            retained categories remain in the shared fit and receive a
            balanced-rest contrast.
        pseudobulk_counts_layer : str, optional
            AnnData layer containing raw counts for pseudobulk aggregation.
            Defaults to "counts" when present, otherwise adata.X.
        pseudobulk_min_cell_counts : int
            Exclude cells below this total raw-count threshold before pseudobulk
            aggregation. Zero disables filtering.
        pseudobulk_min_gene_counts : int
            Exclude genes below this total raw pseudobulk-count threshold in the
            shared DESeq2 fit. Zero disables filtering.
        pseudobulk_min_replicates : int
            Minimum paired replicates required for a reported group-vs-group
            contrast.
            Pseudobulk DE always requires at least two replicates.
        pseudobulk_min_pct_expressed : float
            Minimum fraction of cells expressing a gene required in at least one
            compared group before DeseqStats is run. Values > 1 are interpreted
            as percentages.
        pseudobulk_p_adjust_method : str
            Multiple-testing correction method for pseudobulk p-values.
        pseudobulk_padj_cutoff : float
            Adjusted p-value cutoff used by the viewer volcano/table.
        pseudobulk_log2fc_cutoff : float
            Absolute log2 fold-change cutoff used by the viewer volcano/table.
        pseudobulk_deseq2_fit_type : str
            PyDESeq2 dispersion trend fit type, "parametric" or "mean".
        pseudobulk_n_cpus : int
            Number of CPU workers used by PyDESeq2 pseudobulk fitting and
            contrast evaluation. Must be at least one.
        pseudobulk_embed_top_n_per_comparison : int
            Maximum significant DE genes to auto-embed per category or contact
            comparison. Explicitly requested ``genes`` are always embedded.
        interaction_markers_groupby : list, optional
            Internal list of obs columns to compute contact-conditioned
            pseudobulk interaction markers for. Empty/None disables them.
        neighbor_stats_groupby : list, optional
            Obs columns to compute neighbor composition stats for (categorical only)
        neighbor_stats_permutations : int
            Number of permutations for neighbor enrichment z-scores (0 disables)
        neighbor_stats_seed : int
            Random seed used for neighbor permutations
        interaction_markers_top_targets : int
            Number of target cell types to evaluate per source (ranked by z-score or edge count).
        interaction_markers_top_genes : int
            Number of top genes to keep per source-target interaction.
        interaction_markers_min_cells : int
            Minimum cells required per replicate in both contact+ and contact-
            pseudobulk samples.
        interaction_markers_min_neighbors : int
            Minimum number of target neighbors for a source cell to be labeled contact+.
        Returns
        -------
        dict
            JSON-serializable data structure
        """
        coords = np.asarray(self.adata.obsm[self.spatial_key])[:, :2]
        export_section_indices = self._get_export_section_indices(downsample=downsample)

        # Get UMAP coordinates if available
        umap_coords = None
        umap_bounds = None
        if self.has_umap:
            umap_coords = np.asarray(self.adata.obsm["X_umap"])[:, :2]
            finite_umap_mask = np.isfinite(umap_coords).all(axis=1)
            if finite_umap_mask.any():
                finite_umap = umap_coords[finite_umap_mask]
                # Compute global UMAP bounds from finite rows only so a few bad
                # embedding values do not poison the entire exported viewer.
                umap_bounds = {
                    "xmin": float(finite_umap[:, 0].min()),
                    "xmax": float(finite_umap[:, 0].max()),
                    "ymin": float(finite_umap[:, 1].min()),
                    "ymax": float(finite_umap[:, 1].max()),
                }
            else:
                umap_coords = None

        # Get neighborhood graph if available
        neighbor_graph = None
        neighbor_graph_key = None
        for key in ("spatial_connectivities", "connectivities", "neighbors", "neighbor_graph"):
            if key in self.adata.obsp:
                neighbor_graph = self.adata.obsp[key]
                neighbor_graph_key = key
                break
        if neighbor_graph is not None:
            if not issparse(neighbor_graph):
                neighbor_graph = sp.csr_matrix(neighbor_graph)
            else:
                neighbor_graph = neighbor_graph.tocsr()
            # Some writers store CSR with mismatched index dtypes (e.g. squidpy
            # spatial_connectivities: int64 indptr + int32 indices). scipy fancy
            # indexing then raises "Output dtype not compatible with inputs", so
            # canonicalize indptr/indices to a single consistent dtype.
            if neighbor_graph.indptr.dtype != neighbor_graph.indices.dtype:
                maxval = max(
                    neighbor_graph.shape[0],
                    neighbor_graph.shape[1],
                    int(neighbor_graph.nnz),
                    1,
                )
                idx_dtype = np.int64 if maxval > np.iinfo(np.int32).max else np.int32
                neighbor_graph.indptr = neighbor_graph.indptr.astype(idx_dtype, copy=False)
                neighbor_graph.indices = neighbor_graph.indices.astype(idx_dtype, copy=False)

        # Get initial annotation data
        values, is_continuous, categories = self.get_annotation_data(annotation)

        # Compute global bounds for initial annotation
        if is_continuous:
            finite_mask = np.isfinite(values)
            if finite_mask.any():
                global_vmin = float(np.nanmin(values[finite_mask]))
                global_vmax = float(np.nanmax(values[finite_mask]))
            else:
                global_vmin, global_vmax = 0.0, 1.0
        else:
            global_vmin, global_vmax = None, None

        # Build list of all annotations to export
        all_colors = [annotation]
        if cells_annotations:
            all_colors.extend([c for c in cells_annotations if c != annotation and c in self.obs_columns])

        # Pre-compute all color data
        color_data = {}
        for col in all_colors:
            try:
                vals, is_cont, cats = self.get_annotation_data(col)
                cat_perm = None
                if not is_cont and cats:
                    cat_perm = _numeric_category_perm(cats)
                    if cat_perm is not None:
                        old_to_new = np.empty(len(cats), dtype=np.int64)
                        for new_idx, old_idx in enumerate(cat_perm):
                            old_to_new[old_idx] = new_idx
                        cats = [cats[old_idx] for old_idx in cat_perm]
                        finite = np.isfinite(vals)
                        remapped = vals.copy()
                        remapped[finite] = old_to_new[vals[finite].astype(np.int64)]
                        vals = remapped
                if is_cont:
                    finite = np.isfinite(vals)
                    col_vmin = float(np.nanmin(vals[finite])) if finite.any() else 0.0
                    col_vmax = float(np.nanmax(vals[finite])) if finite.any() else 1.0
                else:
                    col_vmin, col_vmax = None, None
                color_data[col] = {
                    "values": vals,
                    "is_continuous": is_cont,
                    "categories": cats,
                    "vmin": col_vmin,
                    "vmax": col_vmax,
                    "_cat_perm": cat_perm,
                }
            except Exception as e:
                    log_warning(f"Could not load annotation '{col}': {e}", level=1)

        # Pre-compute deconvolution proportion matrices (per-cell × K cell types)
        # deconvolutions: {display_name: obsm_key} where adata.obsm[obsm_key] is (N, K)
        decon_data: Dict[str, Dict[str, Any]] = {}
        if deconvolutions:
            for decon_name, obsm_key in deconvolutions.items():
                if not isinstance(obsm_key, str) or obsm_key not in self.adata.obsm:
                    log_warning(f"deconvolution obsm key '{obsm_key}' not found; skipping '{decon_name}'.", level=1)
                    continue
                raw = self.adata.obsm[obsm_key]
                if hasattr(raw, "columns"):
                    cat_labels = [str(c) for c in raw.columns]
                    matrix = np.asarray(raw.to_numpy(), dtype=np.float32)
                else:
                    matrix = np.asarray(raw, dtype=np.float32)
                    uns_cols_key = f"{obsm_key}_columns"
                    if uns_cols_key in self.adata.uns:
                        cat_labels = [str(c) for c in self.adata.uns[uns_cols_key]]
                    else:
                        cat_labels = [f"comp_{i}" for i in range(matrix.shape[1])]
                if matrix.ndim != 2 or matrix.shape[0] != self.adata.n_obs:
                    log_warning(f"deconvolution '{decon_name}' has bad shape {matrix.shape}; skipping.", level=1)
                    continue
                if len(cat_labels) != matrix.shape[1]:
                    cat_labels = [f"comp_{i}" for i in range(matrix.shape[1])]
                # Row-normalize so wedges sum to 1; tolerate already-normalized inputs.
                row_sums = matrix.sum(axis=1, keepdims=True)
                row_sums = np.where(np.isfinite(row_sums) & (row_sums > 0), row_sums, 1.0)
                matrix = matrix / row_sums
                # Compute dominant component per cell (categorical code, NaN where row is all-zero/NaN)
                dominant = np.argmax(matrix, axis=1).astype(np.float32)
                row_max = matrix.max(axis=1)
                dominant[~np.isfinite(row_max) | (row_max <= 0)] = np.nan
                decon_data[decon_name] = {
                    "matrix": matrix.astype(np.float32, copy=False),
                    "categories": cat_labels,
                    "k": int(matrix.shape[1]),
                    "dominant": dominant,
                }

        gene_encoding = str(gene_encoding or "auto").lower()
        self._validate_gene_export_options(gene_encoding, gene_sparse_zero_threshold, gene_sparse_pack_min_nnz)
        if int(pseudobulk_embed_top_n_per_comparison) < 0:
            raise ValueError("pseudobulk_embed_top_n_per_comparison must be >= 0")
        if int(interaction_markers_top_targets) < 1:
            raise ValueError("interaction_markers_top_targets must be >= 1")
        if int(interaction_markers_top_genes) < 1:
            raise ValueError("interaction_markers_top_genes must be >= 1")
        if int(interaction_markers_min_cells) < 1:
            raise ValueError("interaction_markers_min_cells must be >= 1")
        if int(interaction_markers_min_neighbors) < 1:
            raise ValueError("interaction_markers_min_neighbors must be >= 1")

        def _b64(arr: np.ndarray) -> str:
            import base64
            carr = np.ascontiguousarray(arr)
            return base64.b64encode(carr.tobytes(order="C")).decode("ascii")

        # Prepare float32 views to avoid per-section dtype conversions when packing arrays.
        coords_f4 = np.asarray(coords, dtype=np.float32, order="C")
        umap_f4 = None
        if umap_coords is not None:
            umap_f4 = np.asarray(umap_coords, dtype=np.float32, order="C")

        for col, cdata in color_data.items():
            if cdata.get("is_continuous"):
                cdata["_values_f4"] = np.asarray(cdata["values"], dtype=np.float32, order="C")
            else:
                # Categorical codes are already numeric; keep float32 for compatibility (NaN for missing).
                cdata["_values_f4"] = np.asarray(cdata["values"], dtype=np.float32, order="C")

        # Get metadata filters
        metadata_filters = self.get_metadata_filters()
        companion_analytics = self.get_companion_analytics()

        def _prepare_neighbor_stats_groupby(
            groupby_name: str,
            precomputed_entry: Optional[dict] = None,
        ) -> Tuple[Optional[dict], Optional[dict]]:
            if groupby_name not in self.adata.obs.columns:
                log_warning(f"neighbor stats annotation '{groupby_name}' not found in obs.", level=1)
                return None, None

            col = self.adata.obs[groupby_name]
            if pd.api.types.is_numeric_dtype(col):
                log_warning(f"neighbor stats annotation '{groupby_name}' is numeric; skipping.", level=1)
                return None, None
            if not isinstance(col.dtype, CategoricalDtype):
                col = col.astype("category")

            categories = [str(cat) for cat in col.cat.categories]
            codes = col.cat.codes.to_numpy()
            valid_mask = codes >= 0
            if not valid_mask.any():
                log_warning(f"neighbor stats annotation '{groupby_name}' has no valid categories.", level=1)
                return None, None

            if valid_mask.all():
                graph = neighbor_graph
                labels = codes
                obs_idx = np.arange(self.adata.n_obs, dtype=np.int64)
            else:
                obs_idx = np.flatnonzero(valid_mask).astype(np.int64)
                graph = neighbor_graph[obs_idx][:, obs_idx]
                labels = codes[valid_mask]

            n_cells = np.bincount(labels, minlength=len(categories)).astype(int)
            if graph is None or graph.shape[0] == 0:
                log_warning(f"neighbor stats annotation '{groupby_name}' has an empty graph.", level=1)
                return None, None

            def _build_context(entry_counts, entry_zscore, entry_n_cells, entry_mean_degree):
                return {
                    "categories": categories,
                    "labels": labels.astype(np.int32, copy=False),
                    "graph": graph.tocsr(),
                    "obs_idx": obs_idx,
                    "counts": entry_counts,
                    "zscore": entry_zscore,
                    "n_cells": entry_n_cells,
                    "mean_degree": entry_mean_degree,
                }

            if precomputed_entry is not None:
                companion_categories = [
                    str(cat) for cat in (precomputed_entry.get("categories") or [])
                ]
                if companion_categories and companion_categories != categories:
                    log_warning(
                        f"companion neighbor stats '{groupby_name}' category mismatch; recomputing.",
                        level=1,
                    )
                else:
                    try:
                        counts = np.asarray(precomputed_entry.get("counts"), dtype=float)
                        if counts.shape != (len(categories), len(categories)):
                            raise ValueError("counts shape mismatch")

                        n_cells_ctx = n_cells
                        if precomputed_entry.get("n_cells") is not None:
                            n_cells_ctx = np.asarray(precomputed_entry.get("n_cells"), dtype=int)
                            if n_cells_ctx.shape != (len(categories),):
                                raise ValueError("n_cells shape mismatch")

                        if precomputed_entry.get("mean_degree") is None:
                            mean_degree = np.zeros(len(categories), dtype=float)
                            valid_cells = n_cells_ctx > 0
                            mean_degree[valid_cells] = (
                                counts.sum(axis=1)[valid_cells] / n_cells_ctx[valid_cells]
                            )
                        else:
                            mean_degree = np.asarray(
                                precomputed_entry.get("mean_degree"),
                                dtype=float,
                            )
                            if mean_degree.shape != (len(categories),):
                                raise ValueError("mean_degree shape mismatch")

                        zscore = None
                        if precomputed_entry.get("zscore") is not None:
                            zscore = np.asarray(precomputed_entry.get("zscore"), dtype=float)
                            if zscore.shape != counts.shape:
                                raise ValueError("zscore shape mismatch")

                        entry = {
                            "categories": categories,
                            "counts": counts.tolist(),
                            "n_cells": n_cells_ctx.astype(int).tolist(),
                            "mean_degree": mean_degree.tolist(),
                        }
                        if precomputed_entry.get("perm_n") is not None:
                            entry["perm_n"] = int(precomputed_entry.get("perm_n"))
                        if zscore is not None:
                            entry["zscore"] = zscore.tolist()

                        return entry, _build_context(counts, zscore, n_cells_ctx, mean_degree)
                    except Exception as exc:
                        log_warning(
                            f"companion neighbor stats '{groupby_name}' malformed; recomputing ({exc}).",
                            level=1,
                        )

            onehot = sp.csr_matrix(
                (np.ones(len(labels), dtype=float), (np.arange(len(labels)), labels)),
                shape=(len(labels), len(categories)),
            )
            counts = onehot.T.dot(graph).dot(onehot)
            if issparse(counts):
                counts = counts.toarray()
            counts = np.asarray(counts, dtype=float)
            row_sums = counts.sum(axis=1)
            mean_degree = np.zeros(len(categories), dtype=float)
            valid_cells = n_cells > 0
            mean_degree[valid_cells] = row_sums[valid_cells] / n_cells[valid_cells]

            zscore = None
            entry = {
                "categories": categories,
                "counts": counts.tolist(),
                "n_cells": n_cells.tolist(),
                "mean_degree": mean_degree.tolist(),
            }
            if neighbor_stats_permutations and neighbor_stats_permutations > 0:
                rng = np.random.default_rng(int(neighbor_stats_seed))
                perm_mean = np.zeros_like(counts, dtype=float)
                perm_m2 = np.zeros_like(counts, dtype=float)
                for i in range(int(neighbor_stats_permutations)):
                    perm_labels = rng.permutation(labels)
                    perm_onehot = sp.csr_matrix(
                        (
                            np.ones(len(perm_labels), dtype=float),
                            (np.arange(len(perm_labels)), perm_labels),
                        ),
                        shape=(len(perm_labels), len(categories)),
                    )
                    perm_counts = perm_onehot.T.dot(graph).dot(perm_onehot)
                    if issparse(perm_counts):
                        perm_counts = perm_counts.toarray()
                    perm_counts = np.asarray(perm_counts, dtype=float)
                    delta = perm_counts - perm_mean
                    perm_mean += delta / (i + 1)
                    perm_m2 += delta * (perm_counts - perm_mean)
                if neighbor_stats_permutations > 1:
                    perm_var = perm_m2 / (neighbor_stats_permutations - 1)
                else:
                    perm_var = np.zeros_like(counts, dtype=float)
                perm_std = np.sqrt(perm_var)
                zscore = np.zeros_like(counts, dtype=float)
                valid_std = perm_std > 0
                zscore[valid_std] = (
                    counts[valid_std] - perm_mean[valid_std]
                ) / perm_std[valid_std]
                entry["perm_n"] = int(neighbor_stats_permutations)
                entry["zscore"] = zscore.tolist()

            return entry, _build_context(counts, zscore, n_cells, mean_degree)

        dispersion_baseline_cache: Dict[str, Optional[Dict[str, Any]]] = {}

        def _section_dispersion_baseline(section_id: str, xy_valid: np.ndarray) -> Optional[Dict[str, Any]]:
            cached = dispersion_baseline_cache.get(section_id)
            if cached is not None or section_id in dispersion_baseline_cache:
                return cached

            n_xy = int(xy_valid.shape[0])
            if n_xy < 2:
                dispersion_baseline_cache[section_id] = None
                return None

            try:
                tree = cKDTree(xy_valid)
                max_k = min(256, n_xy - 1)
                sample_size = min(5000, n_xy)
                if n_xy > sample_size:
                    seed = 42 + sum(bytearray(str(section_id).encode("utf-8")))
                    rng = np.random.default_rng(seed)
                    sample_idx = np.sort(rng.choice(n_xy, size=sample_size, replace=False))
                    query_coords = xy_valid[sample_idx]
                else:
                    query_coords = xy_valid

                distances, _indices = tree.query(query_coords, k=max_k + 1)
                distances = np.asarray(distances, dtype=float)
                if distances.ndim == 1:
                    distances = distances.reshape(-1, 1)
                neighbor_distances = distances[:, 1:]
                neighbor_distances[~np.isfinite(neighbor_distances)] = np.nan
                with np.errstate(invalid="ignore"):
                    mean_by_rank = np.nanmean(neighbor_distances, axis=0)
                mean_by_rank = mean_by_rank[np.isfinite(mean_by_rank)]
                if mean_by_rank.size == 0:
                    dispersion_baseline_cache[section_id] = None
                    return None
                baseline = {
                    "n": n_xy,
                    "mean_by_rank": mean_by_rank.astype(float, copy=False),
                }
                dispersion_baseline_cache[section_id] = baseline
                return baseline
            except Exception:
                dispersion_baseline_cache[section_id] = None
                return None

        def _expected_same_type_nn_from_baseline(
            baseline: Optional[Dict[str, Any]],
            category_fraction: float,
        ) -> Optional[float]:
            if baseline is None or category_fraction <= 0:
                return None
            mean_by_rank = np.asarray(baseline.get("mean_by_rank"), dtype=float)
            mean_by_rank = mean_by_rank[np.isfinite(mean_by_rank)]
            if mean_by_rank.size == 0:
                return None
            p = float(np.clip(category_fraction, 1e-12, 1.0))
            ranks = np.arange(mean_by_rank.size, dtype=float)
            weights = p * np.power(1.0 - p, ranks)
            expected = float(np.sum(weights * mean_by_rank))
            tail = max(0.0, 1.0 - float(np.sum(weights)))
            if tail > 0:
                expected += tail * float(mean_by_rank[-1])
            return expected if expected > 0 else None

        def _compute_full_dispersion_for_labels(
            *,
            categories: List[str],
            labels_full: np.ndarray,
            graph_full: Optional[Any],
        ) -> Optional[List[Dict[str, Any]]]:
            if not categories:
                return None
            labels_full = np.asarray(labels_full, dtype=float)
            if labels_full.shape[0] != self.adata.n_obs:
                return None
            valid = np.isfinite(labels_full) & (labels_full >= 0)
            if not np.any(valid):
                return None

            labels_int = np.full(labels_full.shape[0], -1, dtype=np.int32)
            labels_int[valid] = np.rint(labels_full[valid]).astype(np.int32)
            labels_int[(labels_int < 0) | (labels_int >= len(categories))] = -1

            section_indices_full = self.get_section_indices()
            agg = [
                {"catName": str(cat), "n": 0, "saiSum": 0.0, "saiW": 0, "nniSum": 0.0, "nniW": 0}
                for cat in categories
            ]
            rng = np.random.default_rng(42)

            for section in self.sections:
                idx = np.asarray(section_indices_full.get(section.section_id, []), dtype=np.int64)
                if idx.size == 0:
                    continue
                section_coords = coords_f4[idx]
                finite_xy = np.isfinite(section_coords).all(axis=1)
                if not np.any(finite_xy):
                    continue
                section_labels = labels_int[idx]
                if not np.any(section_labels >= 0):
                    continue

                section_graph = None
                degree = None
                if graph_full is not None:
                    try:
                        section_graph = graph_full[idx][:, idx].tocsr()
                        section_graph = section_graph.maximum(section_graph.T).tocsr()
                        section_graph.setdiag(0)
                        section_graph.eliminate_zeros()
                        degree = np.diff(section_graph.indptr).astype(float)
                    except Exception:
                        section_graph = None
                        degree = None

                xy_valid = section_coords[finite_xy]
                n_section_xy = int(xy_valid.shape[0])
                section_baseline = _section_dispersion_baseline(section.section_id, xy_valid)

                for cat_idx, cat_name in enumerate(categories):
                    local_mask = section_labels == cat_idx
                    n_cat = int(np.count_nonzero(local_mask))
                    if n_cat <= 0:
                        continue
                    agg[cat_idx]["n"] += n_cat

                    if section_graph is not None and degree is not None:
                        type_vec = local_mask.astype(np.float32, copy=False)
                        same_counts = np.asarray(section_graph.dot(type_vec)).ravel()
                        type_degree = degree[local_mask]
                        frac = np.zeros(n_cat, dtype=float)
                        valid_degree = type_degree > 0
                        if np.any(valid_degree):
                            frac[valid_degree] = same_counts[local_mask][valid_degree] / type_degree[valid_degree]
                        agg[cat_idx]["saiSum"] += float(np.sum(frac))
                        agg[cat_idx]["saiW"] += n_cat

                    xy_mask = local_mask & finite_xy
                    n_xy = int(np.count_nonzero(xy_mask))
                    if n_xy >= 3:
                        cat_coords = section_coords[xy_mask]
                        try:
                            tree = cKDTree(cat_coords)
                            if n_xy > 4000:
                                sample_size = min(3000, n_xy)
                                sample_idx = np.sort(rng.choice(n_xy, size=sample_size, replace=False))
                                query_coords = cat_coords[sample_idx]
                            else:
                                query_coords = cat_coords
                            distances, _indices = tree.query(query_coords, k=2)
                            nearest = np.asarray(distances[:, 1], dtype=float)
                            nearest = nearest[np.isfinite(nearest)]
                            if nearest.size:
                                observed_mean = float(np.mean(nearest))
                                category_fraction = n_xy / n_section_xy if n_section_xy > 0 else 0.0
                                expected_mean = _expected_same_type_nn_from_baseline(
                                    section_baseline,
                                    category_fraction,
                                )
                                if expected_mean is not None and expected_mean > 0:
                                    nni = observed_mean / expected_mean
                                    agg[cat_idx]["nniSum"] += float(nni) * n_xy
                                    agg[cat_idx]["nniW"] += n_xy
                        except Exception:
                            pass

            rows: List[Dict[str, Any]] = []
            for row in agg:
                n = int(row["n"])
                if n <= 0:
                    continue
                sai = (float(row["saiSum"]) / float(row["saiW"])) if int(row["saiW"]) > 0 else None
                nni = (float(row["nniSum"]) / float(row["nniW"])) if int(row["nniW"]) > 0 else None
                rows.append(
                    {
                        "catName": str(row["catName"]),
                        "n": n,
                        "sai": sai,
                        "nni": nni,
                    }
                )
            return rows or None

        def _compute_full_dispersion_stats(log_as_neighbor_child: bool = False) -> Dict[str, List[Dict[str, Any]]]:
            dispersion: Dict[str, List[Dict[str, Any]]] = {}
            candidates: Dict[str, Tuple[List[str], np.ndarray]] = {}
            for col, cdata in color_data.items():
                if cdata.get("is_continuous"):
                    continue
                cats = [str(cat) for cat in (cdata.get("categories") or [])]
                if cats:
                    candidates[col] = (cats, np.asarray(cdata.get("values"), dtype=float))

            if not candidates:
                return dispersion
            message = (
                f"Computing full-cell spatial dispersion for {len(candidates)} annotation"
                f"{'s' if len(candidates) != 1 else ''}; output feeds Neighbors > Dispersion."
            )
            log_detail(message, level=1)
            for color_col, (cats, labels) in candidates.items():
                rows = _compute_full_dispersion_for_labels(
                    categories=cats,
                    labels_full=labels,
                    graph_full=neighbor_graph,
                )
                if rows:
                    dispersion[color_col] = rows
                    log_detail(
                        f"{color_col}: stored dispersion rows for {len(rows)} categor"
                        f"{'ies' if len(rows) != 1 else 'y'}.",
                        level=2,
                    )
            return dispersion

        replicate_override = str(pseudobulk_replicate_annotation or "").strip()
        pseudobulk_replicate_name = replicate_override or str(self.groupby)
        if pseudobulk_replicate_name not in self.adata.obs.columns:
            raise ValueError(
                "pseudobulk_replicate_annotation "
                f"'{pseudobulk_replicate_name}' is not an obs column"
            )
        marker_genes = {}
        pseudobulk_de = {}
        requested_pseudobulk_de_groupby = list(pseudobulk_de_groupby or [])
        pseudobulk_simple_categories_by_groupby = normalize_pseudobulk_simple_constrast_categories(
            pseudobulk_simple_constrast_categories,
            requested_pseudobulk_de_groupby,
        )
        pending_pseudobulk_de_groupby = list(requested_pseudobulk_de_groupby)
        companion_pseudobulk_de = (
            companion_analytics.get("pseudobulk_de")
            or companion_analytics.get("cluster_de")
        )
        if (
            not replicate_override
            and pending_pseudobulk_de_groupby
            and isinstance(companion_pseudobulk_de, dict)
        ):
            reused_pseudobulk_de_groupby = [
                groupby_name
                for groupby_name in pending_pseudobulk_de_groupby
                if groupby_name in companion_pseudobulk_de
            ]
            if reused_pseudobulk_de_groupby:
                log_step(
                    f"Reusing KaroSpaceCompanion pseudobulk DE for "
                    f"{len(reused_pseudobulk_de_groupby)} annotation column"
                    f"{'s' if len(reused_pseudobulk_de_groupby) != 1 else ''}."
                )
                for groupby_name in reused_pseudobulk_de_groupby:
                    pseudobulk_de[groupby_name] = _strip_category_pseudobulk_sample_diagnostics(
                        companion_pseudobulk_de[groupby_name]
                    )
            pending_pseudobulk_de_groupby = [
                groupby_name
                for groupby_name in pending_pseudobulk_de_groupby
                if groupby_name not in pseudobulk_de
            ]
        if pending_pseudobulk_de_groupby:
            log_step(
                f"Computing pseudobulk differential expression for "
                f"{len(pending_pseudobulk_de_groupby)} annotation column"
                f"{'s' if len(pending_pseudobulk_de_groupby) != 1 else ''}; "
                "output feeds Insights > Compare > Per sample."
            )
            pseudobulk_min_cells_n = 20
            pseudobulk_min_cell_counts_n = int(pseudobulk_min_cell_counts)
            pseudobulk_min_gene_counts_n = int(pseudobulk_min_gene_counts)
            pseudobulk_min_rep_n = int(pseudobulk_min_replicates)
            pseudobulk_min_pct_n = float(pseudobulk_min_pct_expressed)
            pseudobulk_padj_cutoff_n = float(pseudobulk_padj_cutoff)
            pseudobulk_log2fc_cutoff_n = float(pseudobulk_log2fc_cutoff)
            pseudobulk_n_cpus_n = max(1, int(pseudobulk_n_cpus))

            from .pseudobulk import compute_pseudobulk_group_de

            for groupby_name in pending_pseudobulk_de_groupby:
                groupby_pairwise_categories = (
                    (pseudobulk_simple_categories_by_groupby or {}).get(groupby_name)
                    if pseudobulk_simple_categories_by_groupby
                    else None
                )
                log_step(
                    f"Pseudobulk DE: annotation column {groupby_name}; replicate={pseudobulk_replicate_name}",
                    level=1,
                )
                pairwise_categories_label = (
                    ", ".join(str(value) for value in groupby_pairwise_categories)
                    if groupby_pairwise_categories
                    else "all categories"
                )
                log_detail("Parameters:", level=2)
                log_detail(
                    f"model=shared_all_category(~ {pseudobulk_replicate_name} + {groupby_name})",
                    level=3,
                )
                log_detail("rest=balanced_equal_category_weight", level=3)
                log_detail(f"counts_layer={pseudobulk_counts_layer or 'X'}", level=3)
                log_detail(f"min_cell_counts={pseudobulk_min_cell_counts_n}", level=3)
                log_detail(f"min_gene_counts={pseudobulk_min_gene_counts_n}", level=3)
                log_detail(f"min_cells_per_pseudobulk={pseudobulk_min_cells_n}", level=3)
                log_detail(f"min_replicates={max(2, pseudobulk_min_rep_n)}", level=3)
                log_detail(f"min_pct_expressed={pseudobulk_min_pct_n:g}", level=3)
                log_detail(f"p_adjust={pseudobulk_p_adjust_method}", level=3)
                log_detail(f"padj_cutoff={pseudobulk_padj_cutoff_n:g}", level=3)
                log_detail(f"log2fc_cutoff={pseudobulk_log2fc_cutoff_n:g}", level=3)
                log_detail(f"fit_type={pseudobulk_deseq2_fit_type}", level=3)
                log_detail(f"n_cpus={pseudobulk_n_cpus_n}", level=3)
                log_detail("diagnostics=pairwise", level=3)
                log_detail(f"reported_pairwise_categories={pairwise_categories_label}", level=3)
                groupby_results = compute_pseudobulk_group_de(
                    self.adata,
                    groupby_name,
                    replicate=pseudobulk_replicate_name,
                    pairwise_categories=groupby_pairwise_categories,
                    counts_layer=pseudobulk_counts_layer,
                    min_cell_counts=pseudobulk_min_cell_counts_n,
                    min_gene_counts=pseudobulk_min_gene_counts_n,
                    min_cells=pseudobulk_min_cells_n,
                    min_replicates=pseudobulk_min_rep_n,
                    min_pct_expressed=pseudobulk_min_pct_n,
                    p_adjust_method=pseudobulk_p_adjust_method,
                    padj_cutoff=pseudobulk_padj_cutoff_n,
                    log2fc_cutoff=pseudobulk_log2fc_cutoff_n,
                    fit_type=pseudobulk_deseq2_fit_type,
                    n_cpus=max(1, int(pseudobulk_n_cpus)),
                )
                if groupby_results:
                    pseudobulk_de[groupby_name] = groupby_results

        # Compute neighbor composition stats
        neighbor_stats = {}
        neighbor_stats_context = {}
        requested_neighbor_stats_groupby = list(neighbor_stats_groupby or [])
        requested_neighbor_stats_groupby_set = set(requested_neighbor_stats_groupby)
        companion_neighbor_stats = companion_analytics.get("neighbor_stats")
        if requested_neighbor_stats_groupby and isinstance(companion_neighbor_stats, dict):
            reused_neighbor_groupby = [
                groupby_name
                for groupby_name in requested_neighbor_stats_groupby
                if groupby_name in companion_neighbor_stats
            ]
            if reused_neighbor_groupby:
                log_step(
                    f"Reusing KaroSpaceCompanion neighbor stats for {len(reused_neighbor_groupby)} "
                    f"annotation column{'s' if len(reused_neighbor_groupby) != 1 else ''}; "
                    "output feeds Neighbors > Enrichment/Interactions."
                )
                for groupby_name in reused_neighbor_groupby:
                    neighbor_stats[groupby_name] = companion_neighbor_stats[groupby_name]

        pending_neighbor_stats_groupby = [
            groupby_name
            for groupby_name in requested_neighbor_stats_groupby
            if groupby_name not in neighbor_stats
        ]
        printed_neighbor_stats_header = False
        if neighbor_graph is not None and pending_neighbor_stats_groupby:
            for groupby_name in pending_neighbor_stats_groupby:
                log_step(f"Neighbor stats: annotation column {groupby_name}")
                log_detail(
                    "Computing observed neighbor composition; output feeds Neighbors > Enrichment "
                    "and Neighbors > Interactions.",
                    level=1,
                )
                printed_neighbor_stats_header = True
                entry, context = _prepare_neighbor_stats_groupby(groupby_name)
                if entry is None or context is None:
                    continue
                neighbor_stats[groupby_name] = entry
                neighbor_stats_context[groupby_name] = context
                log_detail(
                    f"Stored neighbor matrix for {len(entry.get('categories') or [])} categories.",
                    level=1,
                )

        # Compute contact-conditioned interaction markers:
        # for source S and target T, compare source cells contacting T vs source cells not contacting T.
        interaction_markers = {}
        requested_interaction_markers_groupby = list(interaction_markers_groupby or [])
        pending_interaction_markers_groupby = list(requested_interaction_markers_groupby)
        companion_interaction_markers = companion_analytics.get("interaction_markers")
        if (
            not replicate_override
            and pending_interaction_markers_groupby
            and isinstance(companion_interaction_markers, dict)
        ):
            reused_interaction_groupby = [
                groupby_name
                for groupby_name in pending_interaction_markers_groupby
                if groupby_name in companion_interaction_markers
            ]
            if reused_interaction_groupby:
                log_step(
                    f"Reusing KaroSpaceCompanion interaction markers for "
                    f"{len(reused_interaction_groupby)} "
                    f"annotation column{'s' if len(reused_interaction_groupby) != 1 else ''}."
                )
                for groupby_name in reused_interaction_groupby:
                    interaction_markers[groupby_name] = companion_interaction_markers[groupby_name]
            pending_interaction_markers_groupby = [
                groupby_name
                for groupby_name in pending_interaction_markers_groupby
                if groupby_name not in interaction_markers
            ]
        if neighbor_graph is not None and pending_interaction_markers_groupby:
            log_step(
                f"Computing contact-conditioned pseudobulk interaction markers for "
                f"{len(pending_interaction_markers_groupby)} annotation column"
                f"{'s' if len(pending_interaction_markers_groupby) != 1 else ''}; "
                "output feeds Neighbors > Interactions."
            )
            top_targets = int(interaction_markers_top_targets)
            top_genes = int(interaction_markers_top_genes)
            min_cells = int(interaction_markers_min_cells)
            min_neighbors = int(interaction_markers_min_neighbors)
            interaction_replicate_name = pseudobulk_replicate_name
            interaction_min_rep_n = int(pseudobulk_min_replicates)
            from .pseudobulk import compute_pseudobulk_interaction_markers

            for groupby_name in pending_interaction_markers_groupby:
                log_step(
                    f"Interaction markers: annotation column {groupby_name}; replicate={interaction_replicate_name}",
                    level=1,
                )
                if groupby_name not in neighbor_stats_context:
                    companion_neighbor_entry = None
                    if isinstance(companion_neighbor_stats, dict):
                        companion_neighbor_entry = companion_neighbor_stats.get(groupby_name)
                    entry, context = _prepare_neighbor_stats_groupby(
                        groupby_name,
                        precomputed_entry=companion_neighbor_entry,
                    )
                    if context is not None:
                        neighbor_stats_context[groupby_name] = context
                    if (
                        entry is not None
                        and groupby_name in requested_neighbor_stats_groupby_set
                        and groupby_name not in neighbor_stats
                    ):
                        neighbor_stats[groupby_name] = entry

                ctx = neighbor_stats_context.get(groupby_name)
                if ctx is None:
                    log_warning(
                        f"interaction markers '{groupby_name}' unavailable "
                        "(missing neighbor stats for this annotation).",
                        level=2,
                    )
                    continue

                categories = [str(c) for c in ctx["categories"]]
                labels = np.asarray(ctx["labels"], dtype=np.int32)
                graph = ctx["graph"]
                obs_idx = np.asarray(ctx["obs_idx"], dtype=np.int64)
                counts = np.asarray(ctx["counts"], dtype=float)
                zscore = ctx.get("zscore")
                n_cells = np.asarray(ctx["n_cells"], dtype=int)

                group_interactions = compute_pseudobulk_interaction_markers(
                    self.adata,
                    groupby_name,
                    replicate=interaction_replicate_name,
                    graph=graph,
                    obs_idx=obs_idx,
                    labels=labels,
                    categories=categories,
                    neighbor_counts=counts,
                    neighbor_zscore=zscore,
                    neighbor_n_cells=n_cells,
                    counts_layer=pseudobulk_counts_layer,
                    min_cell_counts=int(pseudobulk_min_cell_counts),
                    min_gene_counts=int(pseudobulk_min_gene_counts),
                    top_targets=top_targets,
                    top_genes=top_genes,
                    min_cells=min_cells,
                    min_neighbors=min_neighbors,
                    min_replicates=interaction_min_rep_n,
                    min_pct_expressed=pseudobulk_min_pct_expressed,
                    p_adjust_method=pseudobulk_p_adjust_method,
                    padj_cutoff=pseudobulk_padj_cutoff,
                    log2fc_cutoff=pseudobulk_log2fc_cutoff,
                    fit_type=pseudobulk_deseq2_fit_type,
                    n_cpus=max(1, int(pseudobulk_n_cpus)),
                )
                if group_interactions:
                    interaction_markers[groupby_name] = group_interactions

        dispersion_stats = _compute_full_dispersion_stats(
            log_as_neighbor_child=printed_neighbor_stats_header
        )

        def _significant_de_genes(
            payload: Any,
            padj_threshold: float,
            log2fc_threshold: float,
            *,
            exclude_category_vs_rest: bool = False,
            limit_per_comparison: int = 20,
        ) -> List[str]:
            found: List[str] = []
            limit = int(limit_per_comparison)
            if limit == 0:
                return found

            def _walk(node: Any, key: Optional[str] = None) -> None:
                if not isinstance(node, dict):
                    return
                if exclude_category_vs_rest and key == "__rest__":
                    return
                genes_list = node.get("genes")
                padj_list = node.get("pvals_adj")
                log2fc_list = node.get("log2foldchanges")
                if not isinstance(log2fc_list, list):
                    log2fc_list = node.get("logfoldchanges")
                if (
                    isinstance(genes_list, list)
                    and isinstance(padj_list, list)
                    and isinstance(log2fc_list, list)
                ):
                    ranked = []
                    for idx, (gene, padj, log2fc) in enumerate(zip(genes_list, padj_list, log2fc_list)):
                        try:
                            padj_value = float(padj)
                            log2fc_value = float(log2fc)
                        except (TypeError, ValueError):
                            continue
                        if (
                            np.isfinite(padj_value)
                            and np.isfinite(log2fc_value)
                            and padj_value < padj_threshold
                            and abs(log2fc_value) >= log2fc_threshold
                        ):
                            ranked.append((padj_value, -abs(log2fc_value), idx, str(gene)))
                    ranked.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
                    if limit > 0:
                        ranked = ranked[:limit]
                    found.extend(gene for *_score, gene in ranked)
                for child_key, value in node.items():
                    if isinstance(value, dict):
                        _walk(value, str(child_key))

            _walk(payload)
            return found

        def _pseudobulk_de_marker_genes(
            payload: Any,
            padj_threshold: float,
            log2fc_threshold: float,
            limit_per_category: int = 50,
        ) -> Dict[str, Dict[str, List[str]]]:
            markers: Dict[str, Dict[str, List[str]]] = {}
            if not isinstance(payload, dict):
                return markers
            rest_reference_key = "__rest__"
            for color_name, by_source in payload.items():
                if not isinstance(by_source, dict):
                    continue
                color_markers: Dict[str, List[str]] = {}
                for source, by_reference in by_source.items():
                    if str(source).startswith("_") or not isinstance(by_reference, dict):
                        continue
                    ranked: Dict[str, Tuple[float, float]] = {}
                    references = (
                        [(rest_reference_key, by_reference[rest_reference_key])]
                        if rest_reference_key in by_reference
                        else list(by_reference.items())
                    )
                    for reference, result in references:
                        if (
                            str(reference).startswith("_")
                            and str(reference) != rest_reference_key
                        ) or not isinstance(result, dict):
                            continue
                        genes_list = result.get("genes")
                        padj_list = result.get("pvals_adj")
                        log2fc_list = result.get("log2foldchanges") or result.get("logfoldchanges")
                        if (
                            not isinstance(genes_list, list)
                            or not isinstance(padj_list, list)
                            or not isinstance(log2fc_list, list)
                        ):
                            continue
                        for idx, (gene, padj, log2fc) in enumerate(zip(genes_list, padj_list, log2fc_list)):
                            try:
                                padj_value = float(padj)
                                log2fc_value = float(log2fc)
                            except (TypeError, ValueError):
                                continue
                            if (
                                np.isfinite(padj_value)
                                and np.isfinite(log2fc_value)
                                and padj_value < padj_threshold
                                and log2fc_value >= log2fc_threshold
                            ):
                                key = str(gene)
                                prev = ranked.get(key)
                                score = (padj_value, -abs(log2fc_value))
                                if prev is None or score < prev:
                                    ranked[key] = score
                    if ranked:
                        color_markers[str(source)] = [
                            gene
                            for gene, _score in sorted(
                                ranked.items(),
                                key=lambda item: (item[1][0], item[1][1], item[0]),
                            )[: int(limit_per_category)]
                        ]
                if color_markers:
                    markers[str(color_name)] = color_markers
            return markers

        requested_genes = list(genes or [])
        de_embed_limit = int(pseudobulk_embed_top_n_per_comparison)
        de_gene_candidates = []
        de_gene_candidates.extend(
            _significant_de_genes(
                pseudobulk_de,
                float(pseudobulk_padj_cutoff),
                float(pseudobulk_log2fc_cutoff),
                limit_per_comparison=de_embed_limit,
            )
        )
        de_gene_candidates.extend(
            _significant_de_genes(
                interaction_markers,
                float(pseudobulk_padj_cutoff),
                float(pseudobulk_log2fc_cutoff),
                limit_per_comparison=de_embed_limit,
            )
        )
        export_genes = list(
            dict.fromkeys(
                gene
                for gene in [*requested_genes, *de_gene_candidates]
                if gene in self.adata.var_names
            )
        )
        if de_gene_candidates:
            log_step("Preparing expression viewer gene payload")
            log_detail(
                f"Embedding {len(export_genes)} requested/significant DE gene"
                f"{'s' if len(export_genes) != 1 else ''} in the HTML expression viewer "
                f"(padj < {float(pseudobulk_padj_cutoff):g}, "
                f"|log2FC| >= {float(pseudobulk_log2fc_cutoff):g}; "
                f"automatic cap={de_embed_limit} per comparison).",
                level=1,
            )
        marker_genes = _pseudobulk_de_marker_genes(
            pseudobulk_de,
            float(pseudobulk_padj_cutoff),
            float(pseudobulk_log2fc_cutoff),
        )

        gene_data = self._collect_gene_data(export_genes)
        gene_encodings = self._resolve_gene_encodings(gene_data, gene_encoding, gene_sparse_zero_threshold)

        # Build section data with all color layers
        sections_data = []
        for section in self.sections:
            idx = export_section_indices[section.section_id]
            rotation_deg = float((section_rotations or {}).get(section.section_id, 0.0))

            section_coords = coords_f4[idx]

            # Get UMAP coordinates for this section if available
            section_umap = None
            if umap_f4 is not None:
                section_umap = umap_f4[idx]

            # Build color values for this section
            section_colors = {}
            section_colors_b64 = {}
            for col, cdata in color_data.items():
                section_vals = cdata.get("_values_f4", cdata["values"])[idx]
                section_colors_b64[col] = _b64(section_vals.astype("<f4", copy=False))

            # Deconvolution proportions: store per-section as flat row-major float32,
            # plus the dominant-component code array (used for legend/filter pipeline).
            section_proportions_b64: Dict[str, Dict[str, Any]] = {}
            for decon_name, ddata in decon_data.items():
                sec_matrix = ddata["matrix"][idx]  # (n_section, K)
                sec_dominant = ddata["dominant"][idx]
                # Dominant code feeds the existing categorical annotation pipeline.
                section_colors_b64[decon_name] = _b64(sec_dominant.astype("<f4", copy=False))
                section_proportions_b64[decon_name] = {
                    "k": ddata["k"],
                    "b64": _b64(np.ascontiguousarray(sec_matrix, dtype="<f4")),
                }

            # Build gene expression values for this section
            section_genes_dense = {}
            section_genes_sparse = {}
            for gene, gdata in gene_data.items():
                section_vals = gdata["values"][idx]
                mode = gene_encodings.get(gene, "dense")
                payload = self._serialize_gene_section_values(
                    section_vals=np.asarray(section_vals),
                    mode=mode,
                    gene_sparse_pack=gene_sparse_pack,
                    gene_sparse_pack_min_nnz=gene_sparse_pack_min_nnz,
                    b64_encoder=_b64,
                )
                if "sparse" in payload:
                    section_genes_sparse[gene] = payload["sparse"]
                else:
                    section_genes_dense[gene] = payload["dense"]

            section_entry = {
                "id": section.section_id,
                "metadata": section.metadata,
                "n_cells": int(len(idx)),
                "rotation_deg": rotation_deg,
                "x": None,
                "y": None,
                "xb64": None,
                "yb64": None,
                "obs_idx": None,
                "obs_idxb64": None,
                "colors": section_colors,
                "colors_b64": section_colors_b64,
                "proportions_b64": section_proportions_b64,
                "genes": section_genes_dense,
                "genes_sparse": section_genes_sparse,
                "bounds": {
                    "xmin": float(section_coords[:, 0].min()) if len(idx) > 0 else 0,
                    "xmax": float(section_coords[:, 0].max()) if len(idx) > 0 else 0,
                    "ymin": float(section_coords[:, 1].min()) if len(idx) > 0 else 0,
                    "ymax": float(section_coords[:, 1].max()) if len(idx) > 0 else 0,
                }
            }

            # Add UMAP coordinates if available
            if section_umap is not None:
                section_entry["umap_x"] = None
                section_entry["umap_y"] = None
                section_entry["umap_xb64"] = None
                section_entry["umap_yb64"] = None
            else:
                section_entry["umap_x"] = None
                section_entry["umap_y"] = None
                section_entry["umap_xb64"] = None
                section_entry["umap_yb64"] = None

            section_entry["xb64"] = _b64(section_coords[:, 0].astype("<f4", copy=False))
            section_entry["yb64"] = _b64(section_coords[:, 1].astype("<f4", copy=False))
            section_entry["obs_idxb64"] = _b64(np.asarray(idx, dtype="<u4"))
            if section_umap is not None:
                section_entry["umap_xb64"] = _b64(section_umap[:, 0].astype("<f4", copy=False))
                section_entry["umap_yb64"] = _b64(section_umap[:, 1].astype("<f4", copy=False))

            if neighbor_graph is not None:
                subgraph = neighbor_graph[idx][:, idx]
                if issparse(subgraph) and subgraph.nnz > 0:
                    upper = sp.triu(subgraph, k=1).tocoo()
                    rows = np.asarray(upper.row, dtype=np.uint32)
                    cols = np.asarray(upper.col, dtype=np.uint32)
                    pairs = np.empty(rows.size * 2, dtype=np.uint32)
                    pairs[0::2] = rows
                    pairs[1::2] = cols
                    section_entry["edges"] = None
                    section_entry["edges_b64"] = _b64(pairs.astype("<u4", copy=False))
                else:
                    section_entry["edges"] = []
                    section_entry["edges_b64"] = None

            sections_data.append(section_entry)

        # Build color metadata
        colors_meta = {}
        uns = getattr(self.adata, "uns", {}) or {}
        for col, cdata in color_data.items():
            meta = {
                "is_continuous": cdata["is_continuous"],
                "categories": cdata["categories"],
                "vmin": cdata["vmin"],
                "vmax": cdata["vmax"],
            }
            if not cdata["is_continuous"]:
                uns_palette = uns.get(f"{col}_colors")
                cats = cdata.get("categories") or []
                if uns_palette is not None and len(uns_palette) == len(cats) and len(cats) > 0:
                    try:
                        # Normalize to CSS colors. uns palettes may be hex strings
                        # OR matplotlib RGB(A) float arrays — the latter str()'d to
                        # "[0.88 0.47 ...]", which the canvas can't parse (cells
                        # render with no color). Fall back to the default palette
                        # if anything fails to convert.
                        palette_list = _normalize_uns_palette(uns_palette)
                        if palette_list is not None:
                            # uns palettes are aligned to the original category order;
                            # follow the same numeric reorder applied to the categories.
                            perm = cdata.get("_cat_perm")
                            if perm is not None and len(perm) == len(palette_list):
                                palette_list = [palette_list[old_idx] for old_idx in perm]
                            meta["palette"] = palette_list
                    except Exception:
                        pass
            colors_meta[col] = meta
        # Deconvolution proportion entries: appear as additional categorical annotation modes.
        for decon_name, ddata in decon_data.items():
            colors_meta[decon_name] = {
                "is_continuous": False,
                "is_proportions": True,
                "categories": list(ddata["categories"]),
                "n_components": ddata["k"],
                "vmin": None,
                "vmax": None,
            }

        # Build gene metadata
        genes_meta = {}
        for gene, gdata in gene_data.items():
            genes_meta[gene] = {
                "vmin": gdata["vmin"],
                "vmax": gdata["vmax"],
            }

        return {
            "initial_color": annotation,
            "groupby": self.groupby,
            "pseudobulk_replicate_annotation": pseudobulk_replicate_name,
            "pseudobulk_settings": {
                "min_replicates": max(2, int(pseudobulk_min_replicates)),
                "min_cell_counts": int(pseudobulk_min_cell_counts),
                "min_gene_counts": int(pseudobulk_min_gene_counts),
                "n_cpus": max(1, int(pseudobulk_n_cpus)),
                "diagnostics": "pairwise",
                "p_adjust_method": str(pseudobulk_p_adjust_method or "fdr_bh"),
                "min_pct_expressed": float(pseudobulk_min_pct_expressed),
                "padj_cutoff": float(pseudobulk_padj_cutoff),
                "log2fc_cutoff": float(pseudobulk_log2fc_cutoff),
                "embed_top_n_per_comparison": int(pseudobulk_embed_top_n_per_comparison),
            },
            "colors_meta": colors_meta,
            "genes_meta": genes_meta,
            "gene_encodings": gene_encodings,
            "metadata_filters": metadata_filters,
            "metadata_section": list(self.metadata_section),
            "metadata_section_extra": list(self.metadata_section_extra),
            "n_sections": len(sections_data),
            "total_cells": sum(s["n_cells"] for s in sections_data),
            "loaded_genes": len(genes_meta),
            "sections": sections_data,
            "available_colors": list(color_data.keys()) + list(decon_data.keys()),
            "available_deconvolutions": list(decon_data.keys()),
            "available_genes": list(gene_data.keys()),
            "marker_genes": marker_genes,
            "pseudobulk_de": pseudobulk_de,
            "has_umap": umap_coords is not None,
            "umap_bounds": umap_bounds,
            "has_neighbors": neighbor_graph is not None,
            "neighbors_key": neighbor_graph_key,
            "neighbor_stats": neighbor_stats,
            "interaction_markers": interaction_markers,
            "dispersion_stats": dispersion_stats,
        }


_MODALITY_OBSM_SKIP_KEYS = {"spatial", "deconvolutions"}
_MODALITY_OBSM_SKIP_PREFIXES = ("X_",)


def _coerce_modality_var(uns_var: Any, expected_rows: int) -> Optional[pd.DataFrame]:
    if uns_var is None:
        return None
    if isinstance(uns_var, pd.DataFrame):
        df = uns_var.copy()
    else:
        try:
            df = pd.DataFrame(uns_var)
        except Exception:
            return None
    if df.shape[0] != expected_rows:
        return None
    if df.index.name is None:
        for cand in ("protein", "gene", "feature", "name"):
            if cand in df.columns:
                df = df.set_index(cand)
                break
    df.index = df.index.map(str)
    return df


def _detect_modalities(adata: sc.AnnData) -> Dict[str, Modality]:
    modalities: Dict[str, Modality] = {}
    rna_var = adata.var.copy()
    rna_layers: Dict[str, Any] = {}
    if "normalized" in adata.layers:
        rna_layers["normalized"] = adata.layers["normalized"]
    if "counts" in adata.layers:
        rna_layers["counts"] = adata.layers["counts"]
    modalities["rna"] = Modality(
        name="rna",
        matrix=adata.X,
        var=rna_var,
        layers=rna_layers,
        value_kind="counts",
        label="RNA",
    )

    for key in list(adata.obsm.keys()):
        if key in _MODALITY_OBSM_SKIP_KEYS:
            continue
        if any(key.startswith(p) for p in _MODALITY_OBSM_SKIP_PREFIXES):
            continue
        matrix = adata.obsm[key]
        shape = getattr(matrix, "shape", None)
        if shape is None or len(shape) != 2 or shape[1] <= 0:
            continue
        var_df = _coerce_modality_var(adata.uns.get(f"{key}_var"), int(shape[1]))
        if var_df is None:
            continue
        modalities[key] = Modality(
            name=key,
            matrix=matrix,
            var=var_df,
            layers={},
            value_kind="intensity",
            label=key.replace("_", " "),
        )
    return modalities


def load_spatial_data(
    path: Any,
    groupby: str = "sample_id",
    spatial_key: str = "spatial",
    spatial_columns: Optional[Tuple[str, str]] = None,
    group_order: Optional[List[str]] = None,
    metadata_section: Optional[List[str]] = None,
    metadata_section_extra: Optional[List[str]] = None,
    metadata_value_order: Optional[Dict[str, List[str]]] = None,
    metadata_max_columns: Optional[int] = None,
    spatialdata_table: Optional[str] = None,
) -> SpatialDataset:
    """
    Load spatial transcriptomics data from an AnnData/H5AD or SpatialData input.

    Parameters
    ----------
    path : str, AnnData, or SpatialData
        Path to an .h5ad file, path to a SpatialData .zarr store, an AnnData
        object, or a SpatialData object.
    groupby : str
        Column in obs to group sections by
    spatial_key : str
        Key in obsm containing spatial coordinates
    spatial_columns : tuple, optional
        Two obs columns (x, y) used to create ``adata.obsm[spatial_key]`` before
        loading. Use this when coordinates are stored as separate metadata
        columns rather than an obsm matrix.
    group_order : list, optional
        Custom order for sections
    metadata_section : list, optional
        Obs columns to use for section metadata and visual filter chips.
    metadata_section_extra : list, optional
        Additional obs columns to store as section metadata without visual filter chips.
    metadata_value_order : dict, optional
        Custom ordering for metadata values per column (e.g. {"course": ["A", "B"]})
        If group_order is not provided, the first key in this dict is used to order sections
        by that metadata column (unknowns last, then section_id sort).
    metadata_max_columns : int, optional
        Limit the number of metadata columns used (order preserved)
    spatialdata_table : str, optional
        AnnData table key to use when the input is a SpatialData object/store.
        Required when a SpatialData input contains multiple tables and no table
        named "table" exists.

    Returns
    -------
    SpatialDataset
        Loaded dataset ready for visualization
    """
    adata, source_label, spatialdata_table_key = _coerce_input_to_anndata(path, spatialdata_table)
    log_step(f"Loading input data from {source_label}")
    log_detail(f"Loaded AnnData table with {adata.n_obs:,} cells x {adata.n_vars:,} genes.")

    groupby = _resolve_groupby_for_spatialdata(adata, groupby, spatialdata_table_key)

    if spatial_columns is not None:
        spatial_key = _set_spatial_from_obs_columns(adata, spatial_columns, spatial_key)
    else:
        spatial_key = _resolve_spatial_key(adata, spatial_key)

    if groupby not in adata.obs.columns:
        raise ValueError(f"Groupby column '{groupby}' not found in adata.obs")

    # Determine section order
    gser = adata.obs[groupby]
    gser_str = gser.astype(str)
    if group_order is not None:
        section_ids = [str(g) for g in group_order if str(g) in gser_str.unique()]
    else:
        order_by_meta = None
        if metadata_value_order:
            order_by_meta = next(iter(metadata_value_order.keys()), None)
        if order_by_meta and order_by_meta in adata.obs.columns:
            desired_order = [str(v) for v in metadata_value_order.get(order_by_meta, [])]
            desired_index = {v: i for i, v in enumerate(desired_order)}
            section_ids = []
            for sid in gser_str.unique():
                mask = gser_str == str(sid)
                vals = adata.obs.loc[mask, order_by_meta].dropna().astype(str).unique()
                meta_value = vals[0] if len(vals) == 1 else "mixed"
                section_ids.append((str(sid), meta_value))
            def _order_key(item):
                sid, meta_value = item
                if meta_value in desired_index:
                    return (0, desired_index[meta_value], sid)
                return (1, meta_value, sid)
            section_ids = [sid for sid, _ in sorted(section_ids, key=_order_key)]
        elif isinstance(gser.dtype, CategoricalDtype) and gser.cat.ordered:
            section_ids = [str(c) for c in gser.cat.categories if str(c) in gser_str.unique()]
        else:
            section_ids = sorted(gser_str.unique())

    log_detail(
        f"Resolved {len(section_ids)} section{'s' if len(section_ids) != 1 else ''} "
        f"from obs column '{groupby}'."
    )

    def _clean_column_list(values: Optional[List[str]]) -> List[str]:
        cleaned = []
        for col in values or []:
            if col is None:
                continue
            text = str(col).strip()
            if text:
                cleaned.append(text)
        return list(dict.fromkeys(cleaned))

    # Determine section metadata columns.
    if metadata_section is None:
        metadata_section = ["course", "region", "condition", "timepoint", "last_score", "last_day"]
    else:
        metadata_section = _clean_column_list(metadata_section)
    metadata_section_extra = _clean_column_list(metadata_section_extra)
    if metadata_max_columns is not None:
        if metadata_max_columns < 0:
            raise ValueError("metadata_max_columns must be >= 0")
        metadata_section = metadata_section[:metadata_max_columns]
    section_metadata_columns = list(dict.fromkeys([*metadata_section, *metadata_section_extra]))
    if section_metadata_columns:
        log_detail(
            "Section metadata columns: "
            + ", ".join(section_metadata_columns)
            + (
                f" ({len(metadata_section)} shown in the visual params bar; "
                f"{len(metadata_section_extra)} stored as section-only metadata)."
            )
        )
    else:
        log_detail("No section metadata columns were requested.")

    # Build section data
    coords = np.asarray(adata.obsm[spatial_key])[:, :2]
    gvals = gser.astype(str).to_numpy()

    sections = []
    for sid in section_ids:
        mask = gvals == sid
        section_coords = coords[mask]

        # Extract metadata
        metadata = {}
        for meta_col in section_metadata_columns:
            if meta_col in adata.obs.columns:
                vals = adata.obs.loc[mask, meta_col].dropna().astype(str).unique()
                if len(vals) == 1:
                    metadata[meta_col] = vals[0]
                elif len(vals) > 1:
                    metadata[meta_col] = "mixed"

        sections.append(SectionData(
            section_id=sid,
            coordinates=section_coords,
            metadata=metadata,
        ))

    # Get available columns for coloring
    obs_columns = [
        col for col in adata.obs.columns
        if isinstance(adata.obs[col].dtype, CategoricalDtype)
        or pd.api.types.is_numeric_dtype(adata.obs[col])
    ]

    modalities = _detect_modalities(adata)
    default_modality = "rna" if "rna" in modalities else next(iter(modalities))
    extra_modalities = [m for m in modalities if m != default_modality]
    if extra_modalities:
        log_detail(f"Modalities: {default_modality} (default) + {', '.join(extra_modalities)}.")

    return SpatialDataset(
        adata=adata,
        sections=sections,
        groupby=groupby,
        obs_columns=obs_columns,
        var_names=list(adata.var_names),
        metadata_section=metadata_section,
        metadata_section_extra=metadata_section_extra,
        metadata_value_order=metadata_value_order,
        modalities=modalities,
        default_modality=default_modality,
        spatial_key=spatial_key,
    )
