import numpy as np
import pandas as pd
from anndata import AnnData

from karospace.data_loader import load_spatial_data


def _make_section_metadata_adata():
    obs = pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s2", "s2"],
            "region": ["cortex", "cortex", "hippocampus", "hippocampus"],
            "cell_type": pd.Categorical(["A", "B", "A", "B"]),
        },
        index=[f"cell{i}" for i in range(4)],
    )
    adata = AnnData(
        X=np.ones((4, 2), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["gene_a", "gene_b"]),
    )
    adata.obsm["spatial"] = np.asarray(
        [
            [0, 0],
            [1, 0],
            [0, 1],
            [1, 1],
        ],
        dtype=np.float32,
    )
    return adata


def test_load_spatial_data_does_not_add_default_section_metadata():
    dataset = load_spatial_data(_make_section_metadata_adata(), section_key="sample_id")

    assert dataset.section_metadata == []
    assert dataset.section_metadata_extra == []
    assert all(section.metadata == {} for section in dataset.sections)

    data = dataset.to_json_data(annotation="cell_type", features=[])
    assert data["section_metadata"] == []
    assert data["section_metadata_extra"] == []
    assert data["metadata_filters"] == {}
    assert all(section["metadata"] == {} for section in data["sections"])


def test_load_spatial_data_keeps_only_existing_requested_section_metadata():
    dataset = load_spatial_data(
        _make_section_metadata_adata(),
        section_key="sample_id",
        section_metadata=["region", "missing_column"],
        section_metadata_extra=["also_missing"],
    )

    assert dataset.section_metadata == ["region"]
    assert dataset.section_metadata_extra == []
    assert [section.metadata for section in dataset.sections] == [
        {"region": "cortex"},
        {"region": "hippocampus"},
    ]

    data = dataset.to_json_data(annotation="cell_type", features=[])
    assert data["section_metadata"] == ["region"]
    assert data["section_metadata_extra"] == []
    assert data["metadata_filters"] == {"region": ["cortex", "hippocampus"]}
