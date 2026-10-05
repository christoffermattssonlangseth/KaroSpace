#!/usr/bin/env python
"""Compare notebook and KaroSpace shared DESeq2 fit paths.

Run from the KaroSpace repository root with the KaroSpace conda environment:

    conda run -n karospace-pseudobulk python compare_deseq2_shared_fit.py

or:

    /Users/bastien.herve/miniconda3/envs/karospace-pseudobulk/bin/python compare_deseq2_shared_fit.py
"""

from __future__ import annotations

import argparse
import inspect
import pickle
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_ARGS_PATH = Path(
    "/Users/bastien.herve/Library/CloudStorage/OneDrive-KarolinskaInstitutet/"
    "Documents/VS Code/KI/Tools/KaroSpace/tests/"
    "karospace-deseq2-shared-fit-args-leiden_rna-1789571685897474000.pkl"
)


def fit_deseq2_shared(counts, metadata, feature_names, retained_categories):
    from pydeseq2.dds import DeseqDataSet

    counts_df = pd.DataFrame(counts, index=metadata.index, columns=feature_names)
    design_meta = metadata[["_pb_replicate", "_pb_group"]].copy()
    design_meta["_pb_replicate"] = pd.Categorical(design_meta["_pb_replicate"].astype(str))
    design_meta["_pb_group"] = pd.Categorical(
        design_meta["_pb_group"].astype(str),
        categories=[str(c) for c in retained_categories],
    )
    try:
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=design_meta,
            design="~ _pb_replicate + _pb_group",
            fit_type=pseudobulk_fit_type,
            n_cpus=max(1, int(pseudobulk_n_cpus)),
            quiet=True,
        )
    except TypeError:
        dds = DeseqDataSet(
            counts=counts_df,
            clinical=design_meta,
            design_factors=["_pb_replicate", "_pb_group"],
            fit_type=pseudobulk_fit_type,
            refit_cooks=True,
            n_cpus=max(1, int(pseudobulk_n_cpus)),
        )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        dds.deseq2()
    return dds


def run_deseq2_pairwise(dds, source, reference):
    from pydeseq2.ds import DeseqStats

    contrast = np.asarray(
        dds.contrast(
            column="_pb_group",
            baseline=str(reference),
            group_to_compare=str(source),
        ),
        dtype=float,
    )
    try:
        kwargs = {
            "dds": dds,
            "contrast": contrast,
            "quiet": True,
            "n_cpus": 1,
            "independent_filter": False,
        }
        if "independent_filter" not in inspect.signature(DeseqStats.__init__).parameters:
            kwargs.pop("independent_filter")
        stat_res = DeseqStats(**kwargs)
    except TypeError:
        stat_res = DeseqStats(dds, contrast=contrast, n_cpus=1)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        stat_res.summary()
    return stat_res.results_df.copy()


def _array_max_abs_diff(left: Any, right: Any) -> float:
    left_arr = np.asarray(left)
    right_arr = np.asarray(right)
    if left_arr.shape != right_arr.shape:
        return float("nan")
    if not (np.issubdtype(left_arr.dtype, np.number) and np.issubdtype(right_arr.dtype, np.number)):
        return 0.0 if np.array_equal(left_arr, right_arr) else float("nan")
    if left_arr.size == 0:
        return 0.0
    return float(np.nanmax(np.abs(left_arr.astype(float) - right_arr.astype(float))))


def _compare_frame(left: pd.DataFrame, right: pd.DataFrame, name: str) -> dict[str, Any]:
    row: dict[str, Any] = {
        "object": name,
        "type": "DataFrame",
        "left_shape": str(left.shape),
        "right_shape": str(right.shape),
        "index_equal": bool(left.index.equals(right.index)),
        "columns_equal": bool(left.columns.equals(right.columns)),
        "max_abs_numeric_diff": 0.0,
        "n_non_numeric_diff_columns": 0,
        "status": "match",
    }
    if left.shape != right.shape or not left.index.equals(right.index) or not left.columns.equals(right.columns):
        row["status"] = "diff"
        return row

    max_abs = 0.0
    non_numeric_diff_columns = 0
    for column in left.columns:
        left_col = left[column]
        right_col = right[column]
        if pd.api.types.is_numeric_dtype(left_col) and pd.api.types.is_numeric_dtype(right_col):
            diff = _array_max_abs_diff(left_col.to_numpy(), right_col.to_numpy())
            if np.isfinite(diff):
                max_abs = max(max_abs, diff)
        elif not left_col.equals(right_col):
            non_numeric_diff_columns += 1

    row["max_abs_numeric_diff"] = max_abs
    row["n_non_numeric_diff_columns"] = non_numeric_diff_columns
    if max_abs > 0 or non_numeric_diff_columns:
        row["status"] = "diff"
    return row


def _compare_mapping(left: Any, right: Any, name: str) -> dict[str, Any]:
    left_keys = set(left.keys())
    right_keys = set(right.keys())
    return {
        "object": name,
        "type": "mapping",
        "left_shape": len(left_keys),
        "right_shape": len(right_keys),
        "index_equal": left_keys == right_keys,
        "columns_equal": None,
        "max_abs_numeric_diff": np.nan,
        "n_non_numeric_diff_columns": len(left_keys.symmetric_difference(right_keys)),
        "status": "match" if left_keys == right_keys else "diff",
    }


def compare_dds_objects(dds_from_notebook_fit, dds_from_karospace_fit) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    for attr_name in ["shape", "n_obs", "n_vars", "obs_names", "var_names", "X"]:
        left_value = getattr(dds_from_notebook_fit, attr_name)
        right_value = getattr(dds_from_karospace_fit, attr_name)
        if isinstance(left_value, pd.Index):
            status = "match" if left_value.equals(right_value) else "diff"
            max_abs = np.nan
            left_shape = len(left_value)
            right_shape = len(right_value)
        elif attr_name == "X":
            max_abs = _array_max_abs_diff(left_value, right_value)
            status = "match" if max_abs == 0 else "diff"
            left_shape = np.asarray(left_value).shape
            right_shape = np.asarray(right_value).shape
        else:
            status = "match" if left_value == right_value else "diff"
            max_abs = np.nan
            left_shape = left_value
            right_shape = right_value
        rows.append(
            {
                "object": f"dds.{attr_name}",
                "type": type(left_value).__name__,
                "left_shape": str(left_shape),
                "right_shape": str(right_shape),
                "index_equal": status == "match",
                "columns_equal": None,
                "max_abs_numeric_diff": max_abs,
                "n_non_numeric_diff_columns": 0 if status == "match" else 1,
                "status": status,
            }
        )

    for attr_name in ["obs", "var", "design_matrix"]:
        left_value = getattr(dds_from_notebook_fit, attr_name, None)
        right_value = getattr(dds_from_karospace_fit, attr_name, None)
        if isinstance(left_value, pd.DataFrame) and isinstance(right_value, pd.DataFrame):
            rows.append(_compare_frame(left_value, right_value, f"dds.{attr_name}"))

    for attr_name in ["layers", "obsm", "varm", "obsp", "varp", "uns"]:
        left_value = getattr(dds_from_notebook_fit, attr_name, None)
        right_value = getattr(dds_from_karospace_fit, attr_name, None)
        if hasattr(left_value, "keys") and hasattr(right_value, "keys"):
            rows.append(_compare_mapping(left_value, right_value, f"dds.{attr_name}"))
            for key in sorted(set(left_value.keys()) & set(right_value.keys()), key=str):
                left_item = left_value[key]
                right_item = right_value[key]
                if isinstance(left_item, pd.DataFrame) and isinstance(right_item, pd.DataFrame):
                    rows.append(_compare_frame(left_item, right_item, f"dds.{attr_name}[{key!r}]"))
                else:
                    max_abs = _array_max_abs_diff(left_item, right_item)
                    same = bool(np.array_equal(np.asarray(left_item), np.asarray(right_item))) if np.isfinite(max_abs) else False
                    rows.append(
                        {
                            "object": f"dds.{attr_name}[{key!r}]",
                            "type": type(left_item).__name__,
                            "left_shape": str(getattr(left_item, "shape", None)),
                            "right_shape": str(getattr(right_item, "shape", None)),
                            "index_equal": None,
                            "columns_equal": None,
                            "max_abs_numeric_diff": max_abs,
                            "n_non_numeric_diff_columns": 0 if same else 1,
                            "status": "match" if same else "diff",
                        }
                    )

    for attr_name in [
        "design_factors",
        "continuous_factors",
        "refit_cooks",
        "min_replicates",
        "min_disp",
        "max_disp",
        "ref_level",
        "fit_type",
        "quiet",
        "n_cpus",
        "size_factors_fit_type",
        "control_genes",
    ]:
        if hasattr(dds_from_notebook_fit, attr_name) or hasattr(dds_from_karospace_fit, attr_name):
            left_value = getattr(dds_from_notebook_fit, attr_name, None)
            right_value = getattr(dds_from_karospace_fit, attr_name, None)
            status = "match" if repr(left_value) == repr(right_value) else "diff"
            rows.append(
                {
                    "object": f"dds.{attr_name}",
                    "type": type(left_value).__name__,
                    "left_shape": repr(left_value),
                    "right_shape": repr(right_value),
                    "index_equal": status == "match",
                    "columns_equal": None,
                    "max_abs_numeric_diff": np.nan,
                    "n_non_numeric_diff_columns": 0 if status == "match" else 1,
                    "status": status,
                }
            )

    return pd.DataFrame(rows)


def compare_result_tables(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    common_columns = [
        column
        for column in ["baseMean", "log2FoldChange", "pvalue", "padj", "stat"]
        if column in left.columns and column in right.columns
    ]
    rows = []
    for column in common_columns:
        delta = pd.to_numeric(left[column], errors="coerce") - pd.to_numeric(right[column], errors="coerce")
        delta_values = delta.to_numpy(dtype=float)
        rows.append(
            {
                "column": column,
                "max_abs_diff": float(np.nanmax(np.abs(delta_values))),
                "n_nonzero_diff": int(np.count_nonzero(np.nan_to_num(delta_values, nan=0.0))),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare notebook fit_deseq2_shared against KaroSpace _fit_deseq2_shared_categories."
    )
    parser.add_argument("--args-pkl", type=Path, default=DEFAULT_ARGS_PATH)
    parser.add_argument("--source", default="1")
    parser.add_argument("--reference", default="2")
    parser.add_argument("--fit-type", default="mean", choices=["parametric", "mean"])
    parser.add_argument("--n-cpus", type=int, default=1)
    parser.add_argument("--show-features", default="CYBB,CD209,SIGLEC1")
    parser.add_argument("--out-prefix", type=Path, default=None)
    args = parser.parse_args()

    global pseudobulk_fit_type, pseudobulk_n_cpus
    pseudobulk_fit_type = str(args.fit_type)
    pseudobulk_n_cpus = max(1, int(args.n_cpus))

    repo_root = Path(__file__).resolve().parent
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from karospace.pseudobulk import _fit_deseq2_shared_categories

    with args.args_pkl.expanduser().open("rb") as handle:
        shared_fit_args = pickle.load(handle)

    fit_counts = shared_fit_args["fit_counts"]
    fit_meta = shared_fit_args["fit_meta"]
    retained_categories = [str(category) for category in shared_fit_args["retained_categories"]]
    feature_names = [str(feature) for feature in fit_meta.attrs["feature_names"]]

    print(f"Python: {sys.executable}")
    try:
        import pydeseq2

        print(f"PyDESeq2: {getattr(pydeseq2, '__version__', 'unknown')}")
    except Exception as exc:
        print(f"PyDESeq2 version unavailable: {exc}")
    print(f"Args pickle: {args.args_pkl}")
    print(f"fit_counts: shape={fit_counts.shape}, dtype={fit_counts.dtype}")
    print(f"fit_meta: shape={fit_meta.shape}, dtypes={fit_meta.dtypes.astype(str).to_dict()}")
    print(f"features: {len(feature_names):,}")
    print(f"retained_categories: {retained_categories}")
    print(f"fit_type={pseudobulk_fit_type}; n_cpus={pseudobulk_n_cpus}; contrast={args.source} vs {args.reference}")

    print("\nFitting notebook fit_deseq2_shared...")
    dds_notebook = fit_deseq2_shared(fit_counts, fit_meta, feature_names, retained_categories)

    print("Fitting KaroSpace _fit_deseq2_shared_categories...")
    dds_karospace, karospace_fit_counts, karospace_fit_meta = _fit_deseq2_shared_categories(
        fit_counts,
        fit_meta,
        retained_categories,
        fit_type=pseudobulk_fit_type,
        min_feature_counts=0,
        n_cpus=pseudobulk_n_cpus,
    )

    args_report = pd.DataFrame(
        [
            {
                "object": "returned_fit_counts",
                "status": "match" if np.array_equal(fit_counts, karospace_fit_counts) else "diff",
                "max_abs_numeric_diff": _array_max_abs_diff(fit_counts, karospace_fit_counts),
            },
            {
                "object": "returned_fit_meta",
                "status": "match" if fit_meta.equals(karospace_fit_meta) else "diff",
                "max_abs_numeric_diff": np.nan,
            },
            {
                "object": "fit_meta.attrs",
                "status": "match" if fit_meta.attrs == karospace_fit_meta.attrs else "diff",
                "max_abs_numeric_diff": np.nan,
            },
        ]
    )
    dds_report = compare_dds_objects(dds_notebook, dds_karospace)

    print("\nReturned argument comparison:")
    print(args_report.to_string(index=False))
    print("\nDDS comparison:")
    print(dds_report.to_string(index=False))
    print("\nDDS differences only:")
    differences = dds_report[dds_report["status"] != "match"]
    print(differences.to_string(index=False) if not differences.empty else "No differences.")

    print("\nRunning contrast from both dds objects...")
    raw_notebook = run_deseq2_pairwise(dds_notebook, args.source, args.reference)
    raw_karospace = run_deseq2_pairwise(dds_karospace, args.source, args.reference)
    result_report = compare_result_tables(raw_notebook, raw_karospace)

    print("\nDESeq2 contrast result delta summary:")
    print(result_report.to_string(index=False))

    selected_features = [
        feature.strip()
        for feature in str(args.show_features or "").split(",")
        if feature.strip() in raw_notebook.index and feature.strip() in raw_karospace.index
    ]
    common_columns = [
        column
        for column in ["baseMean", "log2FoldChange", "pvalue", "padj", "stat"]
        if column in raw_notebook.columns and column in raw_karospace.columns
    ]
    if selected_features:
        details = pd.concat(
            {
                "notebook_fit": raw_notebook.loc[selected_features, common_columns],
                "karospace_fit": raw_karospace.loc[selected_features, common_columns],
                "delta": raw_notebook.loc[selected_features, common_columns]
                - raw_karospace.loc[selected_features, common_columns],
            },
            axis=1,
        )
        print("\nSelected feature comparison:")
        print(details.to_string())

    if args.out_prefix is not None:
        prefix = args.out_prefix.expanduser()
        prefix.parent.mkdir(parents=True, exist_ok=True)
        args_report.to_csv(prefix.with_name(prefix.name + "_returned_args.csv"), index=False)
        dds_report.to_csv(prefix.with_name(prefix.name + "_dds.csv"), index=False)
        result_report.to_csv(prefix.with_name(prefix.name + "_contrast_delta.csv"), index=False)
        raw_notebook.to_csv(prefix.with_name(prefix.name + "_raw_notebook.csv"))
        raw_karospace.to_csv(prefix.with_name(prefix.name + "_raw_karospace.csv"))
        print(f"\nWrote CSV outputs with prefix: {prefix}")

    has_differences = bool((args_report["status"] != "match").any() or (dds_report["status"] != "match").any())
    has_result_differences = bool((result_report["max_abs_diff"].fillna(0) > 0).any())
    return 1 if has_differences or has_result_differences else 0


if __name__ == "__main__":
    raise SystemExit(main())
