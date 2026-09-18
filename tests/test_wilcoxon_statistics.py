import numpy as np
import pandas as pd
from anndata import AnnData

from karospace import wilcoxon as wilcoxon_module
from karospace.wilcoxon import (
    _fit_wilcoxon_runtime_model,
    apply_wilcoxon_runtime_calibration,
    calibrate_wilcoxon_runtime,
    compute_wilcoxon_group_de,
    decide_wilcoxon_runtime,
    estimate_wilcoxon_group_workload,
    parse_wilcoxon_runtime_limit,
    resolve_distribution_expression_matrix,
    resolve_wilcoxon_expression_matrix,
    wilcoxon_runtime_decision_log_lines,
)


def test_wilcoxon_runtime_limit_parser_accepts_hh_mm_ss():
    assert parse_wilcoxon_runtime_limit("00:30:00") == 1800.0
    assert parse_wilcoxon_runtime_limit("01:02:03") == 3723.0
    assert parse_wilcoxon_runtime_limit("45") == 45.0


def test_wilcoxon_runtime_decision_sums_rest_and_pairwise_estimates():
    items = [
        {"kind": "category_vs_rest", "annotation": "cell_type", "estimated_seconds": 18.0},
        {"kind": "category_vs_category", "annotation": "cell_type", "estimated_seconds": 15.0},
    ]

    auto_skip = decide_wilcoxon_runtime("auto", 30.0, items)
    assert auto_skip["should_run"] is False
    assert auto_skip["estimated_seconds"] == 33.0
    assert "exceeds runtime limit" in auto_skip["reason"]

    auto_run = decide_wilcoxon_runtime("auto", 40.0, items)
    assert auto_run["should_run"] is True
    assert auto_run["estimated_seconds"] == 33.0

    force = decide_wilcoxon_runtime("force", 1.0, items)
    assert force["should_run"] is True
    assert "force" in force["reason"]

    off = decide_wilcoxon_runtime("off", 40.0, items)
    assert off["should_run"] is False
    assert "disabled" in off["reason"]


def test_wilcoxon_group_preflight_counts_rest_and_requested_pairwise_work():
    adata = AnnData(
        np.ones((6, 3), dtype=float),
        obs=pd.DataFrame(
            {"cell_type": ["A", "A", "B", "B", "C", "C"]},
            index=[f"cell_{idx}" for idx in range(6)],
        ),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )

    items = estimate_wilcoxon_group_workload(
        adata,
        "cell_type",
        expression_matrix=adata.X,
        pairwise_categories=["A", "B"],
        min_cells=2,
    )

    rest = [item for item in items if item["kind"] == "category_vs_rest"][0]
    pairwise = [item for item in items if item["kind"] == "category_vs_category"][0]
    assert rest["comparison_count"] == 3
    assert np.isclose(rest["work_units"], 6 * np.log2(6) * 3 * 3)
    assert pairwise["comparison_count"] == 1
    assert np.isclose(pairwise["work_units"], 4 * np.log2(4) * 3)


def test_wilcoxon_runtime_calibration_separates_overhead_from_matrix_work():
    calibration = _fit_wilcoxon_runtime_model(
        [
            {"work_units": 1000.0, "elapsed_seconds": 1.10},
            {"work_units": 2000.0, "elapsed_seconds": 1.20},
            {"work_units": 3000.0, "elapsed_seconds": 1.30},
        ],
        safety_multiplier=1.25,
    )

    assert np.isclose(calibration["intercept_seconds"], 1.0)
    assert np.isclose(calibration["seconds_per_work_unit"], 0.0001)

    estimated = apply_wilcoxon_runtime_calibration(
        [
            {
                "kind": "category_vs_rest",
                "annotation": "cell_type",
                "comparison_count": 10,
                "work_units": 50000.0,
            }
        ],
        calibration,
        modality="rna",
    )[0]

    assert np.isclose(estimated["estimated_seconds"], (10.0 * 1.0 + 50000.0 * 0.0001) * 1.25)


def test_wilcoxon_runtime_estimate_accounts_for_parallel_workers():
    estimated = apply_wilcoxon_runtime_calibration(
        [
            {
                "kind": "category_vs_rest",
                "annotation": "cell_type",
                "comparison_count": 8,
                "work_units": 800.0,
            }
        ],
        {
            "intercept_seconds": 2.0,
            "seconds_per_work_unit": 0.1,
            "safety_multiplier": 1.0,
        },
        modality="rna",
        n_cpus=4,
    )[0]

    assert estimated["estimated_parallel_workers"] == 4
    assert np.isclose(estimated["estimated_serial_fixed_seconds"], 16.0)
    assert np.isclose(estimated["estimated_serial_variable_seconds"], 80.0)
    assert np.isclose(estimated["estimated_fixed_seconds"], 4.0)
    assert np.isclose(estimated["estimated_variable_seconds"], 20.0)
    assert np.isclose(estimated["estimated_seconds"], 24.0)

    decision = decide_wilcoxon_runtime("auto", 30.0, [estimated], n_cpus=4)
    assert decision["should_run"] is True

    lines = wilcoxon_runtime_decision_log_lines(decision)
    assert (
        "detail",
        1,
        "runtime_calculation=serial=(fixed=00:00:16 + variable=00:01:20) = 00:01:36; "
        "parallel_adjusted=(fixed=00:00:04 + variable=00:00:20) with n_cpus=4; "
        "* safety_multiplier=1 = 00:00:24; comparisons=8; work_units=800",
    ) in lines


def test_wilcoxon_runtime_calibration_discards_warmup_call():
    calls = []
    times = iter([0.0, 100.0, 100.0, 101.0, 101.0, 103.0, 103.0, 106.0])
    original_scanpy = wilcoxon_module._scanpy_wilcoxon_table
    original_timer = wilcoxon_module.time.perf_counter

    def fake_scanpy(matrix, labels, *, source, reference, feature_names):
        calls.append(tuple(matrix.shape))
        return pd.DataFrame(
            {
                "feature": list(feature_names),
                "score": np.zeros(len(feature_names)),
                "pvalue": np.ones(len(feature_names)),
            }
        )

    wilcoxon_module._scanpy_wilcoxon_table = fake_scanpy
    wilcoxon_module.time.perf_counter = lambda: next(times)
    try:
        calibration = calibrate_wilcoxon_runtime(
            np.ones((8, 8), dtype=float),
            max_cells=8,
            max_features=8,
        )
    finally:
        wilcoxon_module._scanpy_wilcoxon_table = original_scanpy
        wilcoxon_module.time.perf_counter = original_timer

    assert calls == [(2, 2), (2, 2), (4, 4), (8, 8)]
    assert len(calibration["calibration_samples"]) == 3
    assert calibration["work_units"] == 8 * 3 * 8
    assert calibration["calibration_total_elapsed_seconds"] == 6.0


def test_wilcoxon_runtime_calibration_measures_result_formatting():
    format_calls = []
    original_scanpy = wilcoxon_module._scanpy_wilcoxon_table
    original_format = wilcoxon_module._format_wilcoxon_result

    def fake_scanpy(matrix, labels, *, source, reference, feature_names):
        return pd.DataFrame(
            {
                "feature": list(feature_names),
                "score": np.zeros(len(feature_names)),
                "pvalue": np.ones(len(feature_names)),
            }
        )

    def fake_format(rank_table, **kwargs):
        format_calls.append(
            (
                len(rank_table),
                int(kwargs["source_mask"].sum()),
                int(kwargs["reference_mask"].sum()),
                len(kwargs["feature_names"]),
            )
        )
        return {"available": True}

    wilcoxon_module._scanpy_wilcoxon_table = fake_scanpy
    wilcoxon_module._format_wilcoxon_result = fake_format
    try:
        calibration = calibrate_wilcoxon_runtime(
            np.ones((8, 8), dtype=float),
            max_cells=8,
            max_features=8,
        )
    finally:
        wilcoxon_module._scanpy_wilcoxon_table = original_scanpy
        wilcoxon_module._format_wilcoxon_result = original_format

    assert format_calls == [(2, 1, 1, 2), (2, 1, 1, 2), (4, 2, 2, 4), (8, 4, 4, 8)]
    assert calibration["calibration_includes_formatting"] is True


def test_wilcoxon_runtime_decision_log_lines_explain_total_with_shallow_indent():
    decision = decide_wilcoxon_runtime(
        "auto",
        1800.0,
        [
            {
                "kind": "category_vs_rest",
                "comparison_count": 2,
                "work_units": 1000.0,
                "estimated_fixed_seconds": 4.0,
                "estimated_variable_seconds": 6.0,
                "estimated_seconds": 12.5,
            },
            {
                "kind": "category_vs_category",
                "comparison_count": 1,
                "work_units": 500.0,
                "estimated_fixed_seconds": 2.0,
                "estimated_variable_seconds": 3.0,
                "estimated_seconds": 6.25,
            },
        ],
    )

    lines = wilcoxon_runtime_decision_log_lines(decision)

    assert lines[0] == ("step", 0, "runtime_decision=run; estimated=00:00:19; limit=00:30:00")
    assert lines[1][0] == "step"
    assert lines[1][1] == 0
    assert lines[2] == (
        "detail",
        1,
        "runtime_calculation=(fixed=00:00:06 + variable=00:00:09) * safety_multiplier=1.25 = 00:00:19; comparisons=3; work_units=1,500",
    )
    assert (
        "detail",
        1,
        "estimated_by_kind=category_vs_rest; runtime=00:00:12; fixed=00:00:04; variable=00:00:06; comparisons=2; work_units=1,000",
    ) in lines
    assert (
        "detail",
        1,
        "estimated_by_kind=category_vs_category; runtime=00:00:06; fixed=00:00:02; variable=00:00:03; comparisons=1; work_units=500",
    ) in lines


def test_wilcoxon_expression_layer_prefers_normalized_then_log_normalized_raw_data():
    adata = AnnData(np.array([[1, 2], [3, 4]], dtype=float))
    adata.layers["counts"] = np.array([[10, 20], [30, 40]], dtype=float)
    adata.layers["normalized"] = np.array([[0.1, 0.2], [0.3, 0.4]], dtype=float)

    matrix, layer_name = resolve_wilcoxon_expression_matrix(adata, None)
    assert layer_name == "normalized"
    assert np.asarray(matrix)[0, 0] == 0.1

    del adata.layers["normalized"]
    matrix, layer_name = resolve_wilcoxon_expression_matrix(adata, None)
    assert layer_name == "counts_log1p_normalized"
    assert np.isclose(np.asarray(matrix)[0, 0], np.log1p(10 / 30 * 10000))

    del adata.layers["counts"]
    matrix, layer_name = resolve_wilcoxon_expression_matrix(adata, None)
    assert layer_name == "X_log1p_normalized"
    assert np.isclose(np.asarray(matrix)[0, 0], np.log1p(1 / 3 * 10000))


def test_wilcoxon_group_de_accepts_parallel_workers():
    adata = AnnData(
        np.array(
            [
                [4, 1, 0],
                [5, 1, 0],
                [1, 5, 0],
                [1, 4, 0],
                [0, 1, 5],
                [0, 1, 4],
            ],
            dtype=float,
        ),
        obs=pd.DataFrame(
            {"cell_type": ["A", "A", "B", "B", "C", "C"]},
            index=[f"cell_{idx}" for idx in range(6)],
        ),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )

    serial = compute_wilcoxon_group_de(
        adata,
        "cell_type",
        expression_matrix=adata.X,
        expression_layer_used="X",
        min_cells=2,
        n_cpus=1,
    )
    parallel = compute_wilcoxon_group_de(
        adata,
        "cell_type",
        expression_matrix=adata.X,
        expression_layer_used="X",
        min_cells=2,
        n_cpus=2,
    )

    assert parallel is not None
    assert sorted(key for key in parallel if not key.startswith("_")) == ["A", "B", "C"]
    assert parallel.keys() == serial.keys()
    assert parallel["A"].keys() == serial["A"].keys()


def test_distribution_expression_matrix_library_normalizes_without_log_transform():
    adata = AnnData(np.array([[1, 2], [3, 0]], dtype=float))
    adata.layers["counts"] = np.array([[10, 30], [5, 5]], dtype=float)
    adata.layers["normalized"] = np.array([[0.1, 0.2], [0.3, 0.4]], dtype=float)

    matrix, layer_name = resolve_distribution_expression_matrix(adata)

    assert layer_name == "counts_library_normalized"
    assert np.allclose(np.asarray(matrix), [[2500, 7500], [5000, 5000]])


def test_distribution_expression_matrix_supports_scale_log_and_direct_layer():
    adata = AnnData(np.array([[1, 2], [3, 0]], dtype=float))
    adata.layers["counts"] = np.array([[10, 30], [5, 5]], dtype=float)
    adata.layers["data"] = np.array([[0.1, 0.2], [0.3, 0.4]], dtype=float)

    matrix, layer_name = resolve_distribution_expression_matrix(
        adata,
        counts_layer="counts",
        normalization="RC",
        scale_factor=100,
    )
    assert layer_name == "counts_library_normalized"
    assert np.allclose(np.asarray(matrix), [[25, 75], [50, 50]])

    matrix, layer_name = resolve_distribution_expression_matrix(
        adata,
        counts_layer="counts",
        normalization="LogNormalize",
        scale_factor=100,
    )
    assert layer_name == "counts_log_normalized"
    assert np.allclose(np.asarray(matrix), np.log1p([[25, 75], [50, 50]]))

    matrix, layer_name = resolve_distribution_expression_matrix(
        adata,
        counts_layer="counts",
        normalization="RC",
        scale_factor=100,
        normalized_layer="data",
    )
    assert layer_name == "data_layer"
    assert np.allclose(np.asarray(matrix), adata.layers["data"])


def test_wilcoxon_missing_explicit_layer_falls_back_to_auto_selection():
    adata = AnnData(np.array([[1, 2], [3, 4]], dtype=float))

    matrix, layer_name = resolve_wilcoxon_expression_matrix(adata, "counts")

    assert layer_name == "X_log1p_normalized"
    assert np.isclose(np.asarray(matrix)[0, 0], np.log1p(1 / 3 * 10000))


def test_wilcoxon_group_de_emits_rest_and_pairwise_results():
    adata = AnnData(
        np.array(
            [
                [9, 0],
                [8, 0],
                [0, 7],
                [0, 6],
                [5, 1],
                [4, 1],
            ],
            dtype=float,
        ),
        obs=pd.DataFrame(
            {"cell_type": ["A", "A", "B", "B", "C", "C"]},
            index=[f"cell_{idx}" for idx in range(6)],
        ),
        var=pd.DataFrame(index=["marker_a", "marker_b"]),
    )
    result = compute_wilcoxon_group_de(
        adata,
        "cell_type",
        pairwise_categories=["A", "B"],
        expression_layer=None,
        min_cells=2,
        min_pct_expressed=0,
        p_adjust_method="fdr_bh",
        padj_cutoff=1,
        log2fc_cutoff=0,
        top_n_per_comparison=10,
    )

    assert result["_summary"]["source"] == "cell_wilcoxon"
    assert result["A"]["__rest__"]["method"] == "cell-wilcoxon-rest"
    assert result["A"]["B"]["method"] == "cell-wilcoxon-pairwise"
    assert result["B"]["A"]["method"] == "cell-wilcoxon-pairwise"
    assert "marker_a" in result["A"]["__rest__"]["features"]
    assert "marker_b" not in result["A"]["__rest__"]["features"]


def test_wilcoxon_group_de_applies_statistics_count_filters():
    adata = AnnData(
        np.array(
            [
                [50.0, 0.0, 2.0],
                [10.0, 0.0, 1.0],
                [9.0, 0.0, 1.0],
                [0.0, 50.0, 2.0],
                [0.0, 10.0, 1.0],
                [0.0, 9.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cell_type": ["A", "A", "A", "B", "B", "B"]},
            index=[f"cell_{idx}" for idx in range(6)],
        ),
        var=pd.DataFrame(index=["marker_a", "marker_b", "low_count"]),
    )
    adata.layers["normalized"] = adata.X.copy()
    filter_counts = np.array(
        [
            [1, 0, 0],
            [15, 0, 1],
            [15, 0, 1],
            [0, 1, 0],
            [0, 15, 1],
            [0, 15, 1],
        ],
        dtype=float,
    )
    calls = []
    original = wilcoxon_module._scanpy_wilcoxon_table

    def fake_wilcoxon_table(matrix, labels, *, source, reference, feature_names):
        calls.append((tuple(feature_names), int(matrix.shape[0]), int(matrix.shape[1])))
        return pd.DataFrame(
            {
                "feature": list(feature_names),
                "score": [2.0 for _ in feature_names],
                "pvalue": [0.01 for _ in feature_names],
            }
        )

    wilcoxon_module._scanpy_wilcoxon_table = fake_wilcoxon_table
    try:
        result = compute_wilcoxon_group_de(
            adata,
            "cell_type",
            pairwise_categories=["A", "B"],
            filter_expression_matrix=filter_counts,
            filter_expression_source="counts",
            min_cell_counts=10,
            min_feature_counts=10,
            min_cells=2,
            min_pct_expressed=0,
            p_adjust_method="none",
            padj_cutoff=1,
            log2fc_cutoff=0,
            top_n_per_comparison=10,
        )
    finally:
        wilcoxon_module._scanpy_wilcoxon_table = original

    summary = result["_summary"]
    assert summary["filter_expression_source"] == "counts"
    assert summary["n_cells_after_count_filter"] == 4
    assert summary["n_features_after_count_filter"] == 2
    assert summary["category_feature_means"]["features"] == ["marker_a", "marker_b"]
    assert summary["category_feature_means"]["n_cells"] == {"A": 2, "B": 2}
    assert all(call == (("marker_a", "marker_b"), 4, 2) for call in calls)
    assert "low_count" not in result["A"]["__rest__"]["features"]


def test_wilcoxon_pairwise_computes_each_category_pair_once():
    adata = AnnData(
        np.array(
            [
                [9, 1],
                [8, 2],
                [1, 9],
                [2, 8],
                [5, 5],
                [6, 4],
            ],
            dtype=float,
        ),
        obs=pd.DataFrame(
            {"cell_type": ["A", "A", "B", "B", "C", "C"]},
            index=[f"cell_{idx}" for idx in range(6)],
        ),
        var=pd.DataFrame(index=["marker_a", "marker_b"]),
    )
    adata.layers["normalized"] = adata.X.copy()
    calls = []
    original = wilcoxon_module._scanpy_wilcoxon_table

    def fake_wilcoxon_table(matrix, labels, *, source, reference, feature_names):
        calls.append((str(source), str(reference)))
        return pd.DataFrame(
            {
                "feature": list(feature_names),
                "score": [2.0, -1.0],
                "pvalue": [0.01, 0.02],
            }
        )

    wilcoxon_module._scanpy_wilcoxon_table = fake_wilcoxon_table
    try:
        result = compute_wilcoxon_group_de(
            adata,
            "cell_type",
            pairwise_categories=["A", "B", "C"],
            expression_layer=None,
            min_cells=2,
            min_pct_expressed=0,
            p_adjust_method="none",
            padj_cutoff=1,
            log2fc_cutoff=0,
            top_n_per_comparison=10,
        )
    finally:
        wilcoxon_module._scanpy_wilcoxon_table = original

    assert calls == [
        ("A", "rest"),
        ("B", "rest"),
        ("C", "rest"),
        ("A", "B"),
        ("A", "C"),
        ("B", "C"),
    ]
    forward = result["A"]["B"]
    reverse = result["B"]["A"]
    assert reverse["features"] == forward["features"]
    assert reverse["pvals"] == forward["pvals"]
    assert reverse["pvals_adj"] == forward["pvals_adj"]
    assert reverse["base_mean"] == forward["base_mean"]
    assert reverse["n_source"] == forward["n_reference"]
    assert reverse["n_reference"] == forward["n_source"]
    assert reverse["pct_source"] == forward["pct_reference"]
    assert reverse["pct_reference"] == forward["pct_source"]
    assert np.allclose(
        np.asarray(reverse["scores"], dtype=float),
        -np.asarray(forward["scores"], dtype=float),
    )
    assert np.allclose(
        np.asarray(reverse["log2foldchanges"], dtype=float),
        -np.asarray(forward["log2foldchanges"], dtype=float),
    )


def test_wilcoxon_rest_markers_use_log2fc_cutoff():
    adata = AnnData(
        np.array(
            [
                [10.0, 10.0],
                [10.0, 10.0],
                [0.0, 6.0],
                [0.0, 6.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cell_type": ["A", "A", "B", "B"]},
            index=[f"cell_{idx}" for idx in range(4)],
        ),
        var=pd.DataFrame(index=["strong_marker", "weak_marker"]),
    )
    adata.layers["normalized"] = adata.X.copy()
    result = compute_wilcoxon_group_de(
        adata,
        "cell_type",
        pairwise_categories=["A", "B"],
        expression_layer=None,
        min_cells=2,
        min_pct_expressed=0,
        p_adjust_method="none",
        padj_cutoff=1,
        log2fc_cutoff=1,
        top_n_per_comparison=10,
    )

    assert result["A"]["__rest__"]["features"] == ["strong_marker"]
    assert result["A"]["__rest__"]["min_pct_feature_count"] == 2
