import numpy as np
import pandas as pd
from anndata import AnnData

from karospace import wilcoxon as wilcoxon_module
from karospace.wilcoxon import (
    compute_wilcoxon_group_de,
    resolve_wilcoxon_expression_matrix,
)


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
