import io
import json
import uuid
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData

from karospace.cli import _run_export_cli
from karospace.data_loader import Modality, SectionData, SpatialDataset, inspect_input_file
from karospace.exporter import _extract_embedded_viewer_data, export_to_html


def _make_multimodal_dataset():
    obs = pd.DataFrame(
        {
            "section": pd.Categorical(["s1", "s1", "s1", "s1"]),
            "replicate": pd.Categorical(["r1", "r2", "r1", "r2"]),
            "cell_type": pd.Categorical(["A", "A", "B", "B"]),
            "cell_state": pd.Categorical(["X", "Y", "X", "Y"]),
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
        obs_columns=["section", "replicate", "cell_type", "cell_state"],
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


def test_inspect_input_reports_feature_counts_by_modality(tmp_path=None):
    adata = _make_multimodal_dataset().adata.copy()
    adata.obsm["protein"] = np.asarray(
        [
            [5, 1],
            [6, 1],
            [2, 8],
            [2, 9],
        ],
        dtype=np.float32,
    )
    adata.uns["protein_var"] = pd.DataFrame(index=["protein_a", "protein_b"])

    report = inspect_input_file(adata)
    by_name = {entry["name"]: entry for entry in report["feature_modalities"]}
    assert by_name["rna"]["n_features"] == 2
    assert by_name["rna"]["is_default"] is True
    assert by_name["protein"]["n_features"] == 2


def test_cli_inspect_input_prints_feature_counts_by_modality(tmp_path=None):
    output_dir = Path(tmp_path) if tmp_path is not None else Path("/private/tmp")
    input_path = output_dir / "karospace-inspect-multimodal-test.h5ad"
    adata = _make_multimodal_dataset().adata.copy()
    adata.obsm["protein"] = np.asarray(
        [
            [5, 1],
            [6, 1],
            [2, 8],
            [2, 9],
        ],
        dtype=np.float32,
    )
    adata.uns["protein_var"] = pd.DataFrame(index=["protein_a", "protein_b"])
    adata.write_h5ad(input_path)

    stream = io.StringIO()
    with redirect_stdout(stream):
        _run_export_cli([str(input_path), "--inspect-input"])

    output = stream.getvalue()
    assert "Cells: 4\n" in output
    assert "Features by modality:" in output
    assert "  - rna [RNA] (default): 2 features" in output
    assert "  - protein: 2 features" in output
    assert "Genes:" not in output


def test_cli_help_prefers_feature_named_options():
    stream = io.StringIO()
    with redirect_stdout(stream):
        try:
            _run_export_cli(["--help"])
        except SystemExit as exc:
            assert exc.code == 0

    output = stream.getvalue()
    assert "Feature content and storage:" in output
    assert "--pseudobulk-min-feature-counts" in output
    assert "--interaction-markers-top-features" in output
    assert "--feature-correlation-top-n" in output
    assert "--category-means-n-features" in output
    assert "--spatial-variable-features-n" in output
    assert "--pseudobulk-min-gene-counts" not in output
    assert "--interaction-markers-top-genes" not in output
    assert "--gene-correlation-top-n" not in output
    assert "--category-means-n-genes" not in output
    assert "--spatial-variable-genes-n" not in output


def test_secondary_analytics_are_modality_scoped():
    multimodal_dataset = _make_multimodal_dataset()
    stream = io.StringIO()
    with redirect_stdout(stream):
        data = multimodal_dataset.to_json_data(
            annotation="cell_type",
            features=["rna_a"],
            pseudobulk_de_annotations=["cell_type"],
            pseudobulk_replicate_annotation="replicate",
            pseudobulk_modalities=["rna", "protein"],
            pseudobulk_min_replicates=1,
            pseudobulk_min_cells_per_pseudobulk=1,
            spatial_variable_genes_n=2,
            category_means_n_genes=2,
            gene_correlation_top_n=1,
        )

    log_text = stream.getvalue()
    assert "Stored pseudobulk result payload:" in log_text
    assert "Features selected for HTML embedding for modality rna:" in log_text
    assert "Features selected for HTML embedding for modality protein: 2." in log_text
    protein_log = log_text.split("Pseudobulk DE: modality=protein", 1)[1]
    protein_balanced_rest_log = protein_log.split("Preparing 1 pairwise", 1)[0]
    assert "log2FC >= 1" in protein_balanced_rest_log
    assert "|log2FC| >= 1" not in protein_balanced_rest_log
    assert data["embedded_features_by_modality"]["protein"] == ["protein_a", "protein_b"]
    assert set(data["feature_state_by_modality"]["protein"]["features_meta"]) == {"protein_a", "protein_b"}
    assert data["marker_features_by_modality"]["protein"]["cell_type"]["A"] == ["protein_a"]

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


def test_spatial_variable_features_warn_when_graph_missing():
    multimodal_dataset = _make_multimodal_dataset()
    multimodal_dataset.adata.obsp.clear()

    stream = io.StringIO()
    with redirect_stdout(stream):
        data = multimodal_dataset.to_json_data(
            annotation="cell_type",
            features=["rna_a"],
            pseudobulk_de_annotations=[],
            interaction_marker_annotations=[],
            spatial_variable_genes_n=2,
            category_means_n_genes=0,
            gene_correlation_top_n=0,
        )

    log_text = stream.getvalue()
    assert "Warning: Spatially variable feature calculation skipped for modality 'rna'" in log_text
    assert "no spatial graph found in adata.obsp" in log_text
    assert "spatial_connectivities, connectivities, neighbors, neighbor_graph" in log_text
    assert data["spatial_variable_features_by_modality"]["rna"] == []


def test_neighbor_and_interaction_annotations_include_cell_annotations():
    multimodal_dataset = _make_multimodal_dataset()
    data = multimodal_dataset.to_json_data(
        annotation="cell_type",
        cell_annotations=["cell_state"],
        features=["rna_a", "protein_a"],
        pseudobulk_de_annotations=[],
        interaction_marker_annotations=["cell_type", "cell_state"],
        neighbor_stats_annotations=["cell_type", "cell_state"],
        pseudobulk_replicate_annotation="replicate",
        pseudobulk_modalities=["rna", "protein"],
        pseudobulk_min_replicates=1,
        pseudobulk_min_cells_per_pseudobulk=1,
        interaction_markers_min_cells=1,
        interaction_markers_min_neighbors=1,
        interaction_markers_top_targets=2,
        interaction_markers_top_genes=2,
    )

    assert set(data["neighbor_stats"]) == {"cell_type", "cell_state"}
    assert set(data["interaction_markers_by_modality"]) == {"rna", "protein"}
    for modality_payload in data["interaction_markers_by_modality"].values():
        assert set(modality_payload) == {"cell_type", "cell_state"}


def test_export_defaults_analyze_neighbors_for_cell_annotations():
    multimodal_dataset = _make_multimodal_dataset()
    captured = {}

    class CapturedToJsonCall(RuntimeError):
        pass

    def capture_to_json_data(annotation, **kwargs):
        captured["annotation"] = annotation
        captured["kwargs"] = kwargs
        raise CapturedToJsonCall

    original_to_json_data = multimodal_dataset.to_json_data
    multimodal_dataset.to_json_data = capture_to_json_data
    try:
        try:
            export_to_html(
                multimodal_dataset,
                output_path="/private/tmp/karospace-export-defaults-capture.html",
                main_cell_annotation="cell_type",
                cell_annotations=["cell_state"],
                features=[],
                modalities=["rna"],
                pseudobulk=None,
                interaction_markers="auto",
                neighbor_stats_annotations=None,
                neighbor_stats_permutations=0,
                spatial_variable_genes_n=0,
                category_means_n_genes=0,
                gene_correlation_top_n=0,
                pathway_gsea_permutations=0,
                tutorial=False,
            )
        except CapturedToJsonCall:
            pass
    finally:
        multimodal_dataset.to_json_data = original_to_json_data

    assert captured["annotation"] == "cell_type"
    assert captured["kwargs"]["neighbor_stats_annotations"] == ["cell_type", "cell_state"]
    assert captured["kwargs"]["interaction_marker_annotations"] == ["cell_type", "cell_state"]


def test_export_neighbor_defaults_include_cell_annotations_without_interactions():
    multimodal_dataset = _make_multimodal_dataset()
    captured = {}

    class CapturedToJsonCall(RuntimeError):
        pass

    def capture_to_json_data(annotation, **kwargs):
        captured["annotation"] = annotation
        captured["kwargs"] = kwargs
        raise CapturedToJsonCall

    original_to_json_data = multimodal_dataset.to_json_data
    multimodal_dataset.to_json_data = capture_to_json_data
    try:
        try:
            export_to_html(
                multimodal_dataset,
                output_path="/private/tmp/karospace-neighbor-defaults-capture.html",
                main_cell_annotation="cell_type",
                cell_annotations=["cell_state"],
                features=[],
                modalities=["rna"],
                pseudobulk=None,
                interaction_markers=None,
                neighbor_stats_annotations=None,
                neighbor_stats_permutations=0,
                spatial_variable_genes_n=0,
                category_means_n_genes=0,
                gene_correlation_top_n=0,
                pathway_gsea_permutations=0,
                tutorial=False,
            )
        except CapturedToJsonCall:
            pass
    finally:
        multimodal_dataset.to_json_data = original_to_json_data

    assert captured["annotation"] == "cell_type"
    assert captured["kwargs"]["neighbor_stats_annotations"] == ["cell_type", "cell_state"]
    assert captured["kwargs"]["interaction_marker_annotations"] == []


def test_feature_sidecar_manifest_is_modality_only(tmp_path=None):
    output_dir = Path(tmp_path) if tmp_path is not None else Path("/private/tmp")
    output_path = output_dir / "karospace-sidecar-schema-test.html"
    manifest_path = output_dir / "karospace-sidecar-schema-test.features.json"

    export_to_html(
        _make_multimodal_dataset(),
        output_path=str(output_path),
        main_cell_annotation="cell_type",
        features=["rna_a"],
        modalities=["rna", "protein"],
        feature_storage="sidecar",
        feature_manifest_path=str(manifest_path),
        feature_sidecar_shard_size=1,
        pseudobulk=None,
        interaction_markers=None,
        spatial_variable_genes_n=0,
        category_means_n_genes=0,
        gene_correlation_top_n=0,
        pathway_gsea_permutations=0,
        tutorial=False,
    )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["format"] == "karospace-feature-sidecar-manifest-v4"
    assert set(manifest["modalities"]) == {"rna", "protein"}
    for removed_key in [
        "shards",
        "features_meta",
        "feature_encodings",
        "feature_value_encodings",
        "feature_to_shard",
        "section_order",
    ]:
        assert removed_key not in manifest

    for modality_entry in manifest["modalities"].values():
        assert "section_order" in modality_entry
        assert "features_meta" in modality_entry
        assert "feature_to_shard" in modality_entry


def test_embedded_storage_includes_explicit_extra_modalities_without_sidecar(tmp_path=None):
    output_dir = Path(tmp_path) if tmp_path is not None else Path("/private/tmp")
    stem = f"karospace-embedded-default-{uuid.uuid4().hex}"
    output_path = output_dir / f"{stem}.html"
    manifest_path = output_dir / f"{stem}.features.json"
    shard_dir = output_dir / f"{stem}.features"

    stream = io.StringIO()
    with redirect_stdout(stream):
        export_to_html(
            _make_multimodal_dataset(),
            output_path=str(output_path),
            main_cell_annotation="cell_type",
            features=["rna_a", "protein_a"],
            modalities=["rna", "protein"],
            feature_storage="embedded",
            pseudobulk=None,
            pseudobulk_modalities=["rna", "protein"],
            interaction_markers=None,
            spatial_variable_genes_n=0,
            category_means_n_genes=0,
            gene_correlation_top_n=0,
            pathway_gsea_permutations=0,
            tutorial=False,
        )

    data = _extract_embedded_viewer_data(output_path.read_text(encoding="utf-8"))
    log_text = stream.getvalue()
    assert data["feature_manifest_url"] is None
    assert set(data["features_by_modality"]) == {"rna", "protein"}
    assert [entry["name"] for entry in data["modalities"]] == ["rna", "protein"]
    assert data["requested_features_by_modality"] == {
        "rna": ["rna_a"],
        "protein": ["protein_a"],
    }
    assert data["embedded_features_by_modality"]["rna"] == ["rna_a"]
    assert data["embedded_features_by_modality"]["protein"] == ["protein_a"]
    assert set(data["feature_state_by_modality"]["protein"]["features_meta"]) == {"protein_a"}
    assert "Features embedded in HTML: rna=1, protein=1 (total 2)." in log_text
    assert "extra modalities" not in log_text
    assert "feature_storage='sidecar'" not in log_text
    assert "Skipping pathway enrichment for modality protein" not in log_text
    assert "pathway enrichment unavailable for modality protein" not in log_text
    assert not manifest_path.exists()
    assert not shard_dir.exists()


def test_embedded_storage_defaults_to_default_modality_without_sidecar(tmp_path=None):
    output_dir = Path(tmp_path) if tmp_path is not None else Path("/private/tmp")
    stem = f"karospace-embedded-default-only-{uuid.uuid4().hex}"
    output_path = output_dir / f"{stem}.html"
    manifest_path = output_dir / f"{stem}.features.json"
    shard_dir = output_dir / f"{stem}.features"

    stream = io.StringIO()
    with redirect_stdout(stream):
        export_to_html(
            _make_multimodal_dataset(),
            output_path=str(output_path),
            main_cell_annotation="cell_type",
            features=["rna_a", "protein_a"],
            feature_storage="embedded",
            pseudobulk=None,
            interaction_markers=None,
            spatial_variable_genes_n=0,
            category_means_n_genes=0,
            gene_correlation_top_n=0,
            pathway_gsea_permutations=0,
            tutorial=False,
        )

    data = _extract_embedded_viewer_data(output_path.read_text(encoding="utf-8"))
    log_text = stream.getvalue()
    assert data["feature_manifest_url"] is None
    assert set(data["features_by_modality"]) == {"rna"}
    assert data["requested_features_by_modality"] == {"rna": ["rna_a"]}
    assert data["embedded_features_by_modality"]["rna"] == ["rna_a"]
    assert data["embedded_features_by_modality"].get("protein", []) == []
    assert data["feature_state_by_modality"].get("protein", {}).get("features_meta", {}) == {}
    assert "Features embedded in HTML (rna): 1" in log_text
    assert "protein=1" not in log_text
    assert "extra modalities" not in log_text
    assert "feature_storage='sidecar'" not in log_text
    assert not manifest_path.exists()
    assert not shard_dir.exists()


def test_neighbor_interaction_dispersion_logs_follow_computation_tree():
    stream = io.StringIO()
    with redirect_stdout(stream):
        _make_multimodal_dataset().to_json_data(
            annotation="cell_type",
            cell_annotations=["cell_state"],
            features=["rna_a", "protein_a"],
            pseudobulk_de_annotations=[],
            interaction_marker_annotations=["cell_type", "cell_state"],
            neighbor_stats_annotations=["cell_type", "cell_state"],
            pseudobulk_replicate_annotation="replicate",
            pseudobulk_modalities=["rna", "protein"],
            pseudobulk_min_replicates=1,
            pseudobulk_min_cells_per_pseudobulk=1,
            interaction_markers_min_cells=1,
            interaction_markers_min_neighbors=1,
            interaction_markers_top_targets=2,
            interaction_markers_top_genes=2,
        )

    log_text = stream.getvalue()
    neighbor_idx = log_text.index("- Computing neighbor composition stats")
    interaction_idx = log_text.index("- Computing contact-conditioned pseudobulk interaction markers")
    dispersion_idx = log_text.index("- Computing full-cell spatial dispersion")
    assert neighbor_idx < interaction_idx < dispersion_idx

    lines = log_text.splitlines()
    assert "  - Neighbor stats: annotation column cell_type" in lines
    assert "    - A -> B: skipped, insufficient paired replicates (0; need >= 2)" in lines
    assert "  - A -> B: skipped, insufficient paired replicates (0; need >= 2)" not in lines

    dispersion_rows = [line for line in lines if "cell_type: stored dispersion rows" in line]
    assert dispersion_rows == ["  ↳ cell_type: stored dispersion rows for 2 categories."]
