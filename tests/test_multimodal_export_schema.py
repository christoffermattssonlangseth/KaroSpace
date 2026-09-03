import numpy as np
import pandas as pd
from anndata import AnnData

from karospace.data_loader import Modality, SectionData, SpatialDataset


def _make_multimodal_dataset():
    obs = pd.DataFrame(
        {
            "section": pd.Categorical(["s1", "s1", "s1", "s1"]),
            "replicate": pd.Categorical(["r1", "r2", "r1", "r2"]),
            "cell_type": pd.Categorical(["A", "A", "B", "B"]),
        },
        index=[f"cell{i}" for i in range(4)],
    )
    var = pd.DataFrame(index=["rna_a", "rna_b"])
    x = np.asarray(
        [
            [1, 0],
            [2, 1],
            [0, 3],
            [1, 4],
        ],
        dtype=np.float32,
    )
    adata = AnnData(X=x, obs=obs, var=var)
    adata.layers["normalized"] = x
    adata.obsm["spatial"] = np.asarray(
        [
            [0, 0],
            [1, 0],
            [0, 1],
            [1, 1],
        ],
        dtype=np.float32,
    )
    adata.obsp["spatial_connectivities"] = np.asarray(
        [
            [0, 1, 1, 0],
            [1, 0, 0, 1],
            [1, 0, 0, 1],
            [0, 1, 1, 0],
        ],
        dtype=np.float32,
    )

    protein_x = np.asarray(
        [
            [5, 1],
            [6, 1],
            [2, 8],
            [2, 9],
        ],
        dtype=np.float32,
    )
    modalities = {
        "rna": Modality(
            name="rna",
            matrix=x,
            var=var,
            layers={"normalized": x},
            value_kind="counts",
            label="RNA",
        ),
        "protein": Modality(
            name="protein",
            matrix=protein_x,
            var=pd.DataFrame(index=["protein_a", "protein_b"]),
            layers={"normalized": protein_x},
            value_kind="intensity",
            label="Protein",
        ),
    }
    return SpatialDataset(
        adata=adata,
        sections=[SectionData("s1", adata.obsm["spatial"])],
        section_key="section",
        obs_columns=["section", "replicate", "cell_type"],
        var_names=list(var.index),
        modalities=modalities,
        default_modality="rna",
    )


def test_multimodal_export_uses_only_by_modality_payloads():
    multimodal_dataset = _make_multimodal_dataset()
    data = multimodal_dataset.to_json_data(
        annotation="cell_type",
        features=["rna_a"],
        pseudobulk_de_annotations=[],
        interaction_marker_annotations=[],
        pseudobulk_modalities=["rna", "protein"],
    )

    assert data["default_modality"] == "rna"
    assert set(data["features_by_modality"]) == {"rna", "protein"}
    assert set(data["embedded_features_by_modality"]) == {"rna", "protein"}
    assert set(data["feature_state_by_modality"]) == {"rna", "protein"}
    assert set(data["pseudobulk_de_by_modality"]) == {"rna", "protein"}
    assert set(data["interaction_markers_by_modality"]) == {"rna", "protein"}

    rna_state = data["feature_state_by_modality"]["rna"]
    assert set(rna_state["features_meta"]) == {"rna_a"}
    assert rna_state["feature_encodings"]["rna_a"] in {"dense", "sparse"}
    assert data["embedded_features_by_modality"]["rna"] == ["rna_a"]
    assert data["embedded_features_by_modality"]["protein"] == []

    for removed_key in [
        "available_features",
        "embedded_features",
        "features_meta",
        "feature_encodings",
        "feature_value_encodings",
        "pseudobulk_de",
        "interaction_markers",
        "category_gene_means",
        "gene_correlations",
        "spatial_variable_genes",
        "pathway_settings",
    ]:
        assert removed_key not in data


def test_secondary_analytics_are_modality_scoped():
    multimodal_dataset = _make_multimodal_dataset()
    data = multimodal_dataset.to_json_data(
        annotation="cell_type",
        features=["rna_a", "protein_a"],
        pseudobulk_de_annotations=["cell_type"],
        pseudobulk_replicate_annotation="replicate",
        pseudobulk_modalities=["rna", "protein"],
        pseudobulk_min_replicates=1,
        pseudobulk_min_cells_per_pseudobulk=1,
        spatial_variable_genes_n=2,
        category_means_n_genes=2,
        gene_correlation_top_n=1,
    )

    assert set(data["category_feature_means_by_modality"]) == {"rna", "protein"}
    assert set(data["feature_correlations_by_modality"]) == {"rna", "protein"}
    assert set(data["spatial_variable_features_by_modality"]) == {"rna", "protein"}

    for modality, feature_name in [("rna", "rna_a"), ("protein", "protein_a")]:
        means = data["category_feature_means_by_modality"][modality]
        assert means is not None
        assert feature_name in means["genes"]

        correlations = data["feature_correlations_by_modality"][modality]
        assert feature_name in correlations

        spatial = data["spatial_variable_features_by_modality"][modality]
        assert any(row["gene"] == feature_name for row in spatial)
