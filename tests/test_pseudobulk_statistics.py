import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from karospace.pseudobulk import (
    _aggregate_display_expression_matrix,
    _compute_category_feature_means_from_aggregate,
    _prepare_pseudobulk_distribution_aggregate,
    compute_pseudobulk_group_de,
)


def test_pseudobulk_category_means_can_use_display_scale_matrix():
    raw_counts = np.asarray([[220.0], [180.0]])
    display_values = np.asarray([[2.2], [1.8]])
    pb_meta = pd.DataFrame(
        {
            "_pb_replicate": ["S1", "S2"],
            "_pb_group": ["T", "B"],
            "n_cells": [2, 2],
        }
    )

    raw_summary = _compute_category_feature_means_from_aggregate(
        raw_counts,
        pb_meta,
        ["T", "B"],
        ["CD3"],
    )
    display_summary = _compute_category_feature_means_from_aggregate(
        display_values,
        pb_meta,
        ["T", "B"],
        ["CD3"],
        source="pseudobulk_display_aggregate",
    )

    assert raw_summary["means"]["T"][0] == 110.0
    assert display_summary["means"]["T"][0] == 1.1
    assert display_summary["source"] == "pseudobulk_display_aggregate"


def test_display_expression_aggregation_does_not_round_to_counts():
    incidence = sp.csr_matrix(
        np.asarray(
            [
                [1.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 1.0],
            ]
        )
    )
    display_matrix = np.asarray([[0.11], [0.12], [0.08], [0.10]])

    aggregate = _aggregate_display_expression_matrix(
        incidence,
        display_matrix,
        expected_shape=display_matrix.shape,
    )

    assert np.allclose(aggregate, [[0.23], [0.18]])


def test_pseudobulk_distribution_applies_feature_and_sample_filters_before_means():
    raw_aggregate = np.asarray(
        [
            [10.0, 1.0],
            [20.0, 2.0],
            [30.0, 3.0],
        ]
    )
    display_aggregate = np.asarray(
        [
            [10.0, 100.0],
            [20.0, 200.0],
            [30.0, 300.0],
        ]
    )
    pb_meta = pd.DataFrame(
        {
            "_pb_replicate": ["r1", "r2", "r1"],
            "_pb_group": ["A", "A", "B"],
            "n_cells": [2, 5, 3],
        }
    )
    pb_meta.attrs["feature_names"] = ["kept", "filtered"]

    filtered_counts, filtered_display, filtered_meta, summary = _prepare_pseudobulk_distribution_aggregate(
        raw_aggregate,
        display_aggregate,
        pb_meta,
        ["A", "B"],
        ["kept", "filtered"],
        min_feature_counts=20,
        min_cells=3,
    )

    assert filtered_counts.shape == (3, 1)
    assert filtered_display.shape == (3, 1)
    assert filtered_meta.attrs["feature_names"] == ["kept"]
    assert summary["features"] == ["kept"]
    assert summary["n_cells"] == {"A": 5, "B": 3}
    assert summary["means"]["A"] == [4.0]
    assert summary["means"]["B"] == [10.0]


def test_pseudobulk_distribution_survives_min_replicates_de_skip():
    adata = AnnData(
        X=np.asarray(
            [
                [10.0],
                [20.0],
                [30.0],
                [40.0],
            ]
        ),
        obs=pd.DataFrame(
            {
                "replicate": ["r1", "r1", "r2", "r2"],
                "cell_type": pd.Categorical(["A", "A", "B", "B"]),
            }
        ),
        var=pd.DataFrame(index=["gene_a"]),
    )

    result = compute_pseudobulk_group_de(
        adata,
        "cell_type",
        replicate="replicate",
        counts_layer=None,
        display_expression_matrix=np.asarray([[1.0], [3.0], [6.0], [10.0]]),
        min_cells=1,
        min_replicates=2,
    )

    assert list(result.keys()) == ["_summary"]
    means = result["_summary"]["category_feature_means"]
    assert means["features"] == ["gene_a"]
    assert means["n_cells"] == {"A": 2, "B": 2}
    assert means["means"] == {"A": [2.0], "B": [8.0]}
