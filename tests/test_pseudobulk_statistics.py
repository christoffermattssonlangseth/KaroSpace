import numpy as np
import pandas as pd
import scipy.sparse as sp

from karospace.pseudobulk import (
    _aggregate_display_expression_matrix,
    _compute_category_feature_means_from_aggregate,
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
