import json
import re

import h5py
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import karospace.cli as cli_module
from karospace.data_loader import SectionData, SpatialDataset, load_spatial_data
from karospace.exporter import export_to_html


def _build_dataset() -> SpatialDataset:
    obs = pd.DataFrame(
        {
            "sample_id": pd.Categorical(["S1", "S1", "S2", "S2"]),
            "leiden": pd.Categorical(["A", "B", "A", "B"]),
            "condition": pd.Categorical(["ctrl", "ctrl", "treated", "treated"]),
        },
        index=[f"cell_{i}" for i in range(4)],
    )
    var = pd.DataFrame(index=["G1", "G2", "G3"])
    x = np.array(
        [
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 5.0],
            [3.0, 1.0, 0.0],
            [4.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    adata = AnnData(X=x, obs=obs, var=var)
    adata.layers["normalized"] = x.copy()
    adata.obsm["spatial"] = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ],
        dtype=float,
    )
    sections = [
        SectionData(section_id="S1", coordinates=adata.obsm["spatial"][:2], metadata={"condition": "ctrl"}),
        SectionData(section_id="S2", coordinates=adata.obsm["spatial"][2:], metadata={"condition": "treated"}),
    ]
    return SpatialDataset(
        adata=adata,
        sections=sections,
        groupby="sample_id",
        obs_columns=["sample_id", "leiden", "condition"],
        var_names=["G1", "G2", "G3"],
        metadata_columns=["condition"],
    )


def _extract_data_json(html_text: str) -> dict:
    match = re.search(
        r'<script id="karospace-data" type="application/json">(.*?)</script>',
        html_text,
        re.DOTALL,
    )
    assert match, "embedded data script not found"
    return json.loads(match.group(1).replace("<\\/", "</"))


def test_sidecar_export_writes_aux_and_updates_html_contract(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        color="leiden",
        genes=["G1"],
        use_hvgs=False,
        gene_encoding="sparse",
        gene_storage="sidecar",
    )

    aux_path = tmp_path / "viewer.genes.json"
    shard_dir = tmp_path / "viewer.genes"
    assert output_path.exists()
    assert aux_path.exists()
    assert shard_dir.exists()

    html_text = output_path.read_text(encoding="utf-8")
    embedded = _extract_data_json(html_text)
    manifest = json.loads(aux_path.read_text(encoding="utf-8"))
    shard_files = sorted(shard_dir.glob("*.json"))
    assert shard_files
    shard = json.loads(shard_files[0].read_text(encoding="utf-8"))

    assert embedded["available_genes"] == ["G1", "G2", "G3"]
    assert embedded["gene_aux_url"] == "viewer.genes.json"
    assert set(embedded["genes_meta"]) == {"G1"}
    assert set(manifest["genes_meta"]) == {"G2", "G3"}
    assert manifest["format"] == "karospace-gene-sidecar-manifest-v2"
    assert manifest["gene_to_shard"]["G2"].startswith("viewer.genes/")
    assert shard["format"] == "karospace-gene-sidecar-shard-v2"
    assert "G1" not in shard["genes"]
    assert "G2" in shard["genes"]
    assert "sparse" in shard["genes"]["G2"]["sections"]["S1"]
    assert "async function ensureGeneAvailable" in html_text
    assert "function requestModalBlendGene(gene)" in html_text
    assert "requestModalBlendGene(gene);" in html_text
    assert 'id="modal-controls-toggle"' in html_text
    assert ".modal-controls.hidden" in html_text
    assert 'data-modal-group="view"' in html_text
    assert "function updateModalToolbarState()" in html_text
    assert "function initModalControlsDragging()" in html_text
    assert 'id="modal-blend-loading"' in html_text
    assert 'id="gene-discovery-panel"' in html_text
    assert 'id="gene-panel-new"' in html_text
    assert 'id="genes-tab-compare"' in html_text
    assert "function getGeneSearchResults(query, limit = GENE_DISCOVERY_MAX_RESULTS)" in html_text
    assert "function renderGeneDiscoveryPanel()" in html_text
    assert "function renderClusterDE()" in html_text
    assert "function recordRecentGene(gene)" in html_text
    assert "function loadSavedGenePanels()" in html_text
    assert "const GENE_RECENTS_STORAGE_KEY = getViewerScopedStorageKey('gene-recents');" in html_text
    assert html_text.index("function getMarkerGenesForColorCategory(colorCol, category)") < html_text.index("function getGeneSuggestionGroups()")
    assert "was not pre-loaded" not in html_text


def test_embedded_mode_stays_single_file(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        color="leiden",
        genes=["G1"],
        use_hvgs=False,
        gene_storage="embedded",
    )

    html_text = output_path.read_text(encoding="utf-8")
    embedded = _extract_data_json(html_text)

    assert output_path.exists()
    assert not (tmp_path / "viewer.genes.json").exists()
    assert embedded["gene_aux_url"] is None
    assert embedded["available_genes"] == ["G1"]
    assert {section["id"]: section["rotation_deg"] for section in embedded["sections"]} == {
        "S1": 0.0,
        "S2": 0.0,
    }


def test_export_embeds_normalized_section_rotations(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        color="leiden",
        use_hvgs=False,
        section_rotations={"S1": 44.5, "S2": -50.25},
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    assert {section["id"]: section["rotation_deg"] for section in embedded["sections"]} == {
        "S1": 44.5,
        "S2": 309.75,
    }


def test_export_includes_initial_categorical_color_in_marker_genes(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        color="leiden",
        use_hvgs=False,
        marker_genes_groupby=["condition"],
        marker_genes_top_n=5,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    marker_genes = embedded.get("marker_genes") or {}
    assert "leiden" in marker_genes
    assert "condition" in marker_genes


def test_export_embeds_pairwise_cluster_de(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        color="leiden",
        use_hvgs=False,
        cluster_de_groupby=["leiden"],
        cluster_de_top_n=2,
        cluster_de_method="t-test",
        cluster_de_min_cells=1,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    cluster_de = embedded.get("cluster_de") or {}
    assert "leiden" in cluster_de
    assert set(cluster_de["leiden"]) == {"A", "B"}
    assert "A" not in cluster_de["leiden"]["A"]
    result = cluster_de["leiden"]["A"]["B"]
    assert result["available"] is True
    assert result["n_source"] == 2
    assert result["n_reference"] == 2
    assert 1 <= len(result["genes"]) <= 2
    assert len(result["genes"]) == len(result["logfoldchanges"]) == len(result["pvals_adj"]) == len(result["scores"])
    assert len(result["genes"]) == len(result["pct_source"]) == len(result["pct_reference"])


def test_export_marks_cluster_de_unavailable_for_small_clusters(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        color="leiden",
        use_hvgs=False,
        cluster_de_groupby=["leiden"],
        cluster_de_method="t-test",
        cluster_de_min_cells=3,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    result = embedded["cluster_de"]["leiden"]["A"]["B"]
    assert result["available"] is False
    assert result["reason"] == "insufficient_cells"
    assert result["min_cells_required"] == 3


def test_export_rejects_unknown_section_rotation_ids(tmp_path):
    dataset = _build_dataset()

    with pytest.raises(ValueError, match="unknown section_id"):
        export_to_html(
            dataset,
            output_path=str(tmp_path / "viewer.html"),
            color="leiden",
            use_hvgs=False,
            section_rotations={"missing": 45},
        )


def test_cli_passes_section_rotations_to_export(monkeypatch, tmp_path):
    input_path = tmp_path / "input.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    captured = {}

    def fake_load(path, groupby):
        captured["load"] = {"path": path, "groupby": groupby}
        return _build_dataset()

    def fake_export(dataset, **kwargs):
        captured["kwargs"] = kwargs
        return str(tmp_path / "viewer.html")

    monkeypatch.setattr("karospace.data_loader.load_spatial_data", fake_load)
    monkeypatch.setattr("karospace.exporter.export_to_html", fake_export)
    monkeypatch.setattr(
        cli_module.sys,
        "argv",
        [
            "karospace",
            str(input_path),
            "--section-rotations",
            "S1:45,S2:-90",
        ],
    )

    cli_module.main()

    assert captured["load"] == {"path": str(input_path), "groupby": "sample_id"}
    assert captured["kwargs"]["section_rotations"] == {"S1": 45.0, "S2": -90.0}


def test_cli_accepts_fractional_section_rotations(monkeypatch, tmp_path):
    input_path = tmp_path / "input.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    captured = {}

    def fake_load(path, groupby):
        return _build_dataset()

    def fake_export(dataset, **kwargs):
        captured["kwargs"] = kwargs
        return str(tmp_path / "viewer.html")

    monkeypatch.setattr("karospace.data_loader.load_spatial_data", fake_load)
    monkeypatch.setattr("karospace.exporter.export_to_html", fake_export)
    monkeypatch.setattr(
        cli_module.sys,
        "argv",
        [
            "karospace",
            str(input_path),
            "--section-rotations",
            "S1:37.5,S2:-12.25",
        ],
    )

    cli_module.main()

    assert captured["kwargs"]["section_rotations"] == {"S1": 37.5, "S2": -12.25}


def test_cli_passes_cluster_de_options_to_export(monkeypatch, tmp_path):
    input_path = tmp_path / "input.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    captured = {}

    def fake_load(path, groupby):
        return _build_dataset()

    def fake_export(dataset, **kwargs):
        captured["kwargs"] = kwargs
        return str(tmp_path / "viewer.html")

    monkeypatch.setattr("karospace.data_loader.load_spatial_data", fake_load)
    monkeypatch.setattr("karospace.exporter.export_to_html", fake_export)
    monkeypatch.setattr(
        cli_module.sys,
        "argv",
        [
            "karospace",
            str(input_path),
            "--cluster-de-groupby",
            "leiden,condition",
            "--cluster-de-top-n",
            "12",
            "--cluster-de-method",
            "t-test",
            "--cluster-de-layer",
            "normalized",
            "--cluster-de-min-cells",
            "7",
        ],
    )

    cli_module.main()

    assert captured["kwargs"]["cluster_de_groupby"] == ["leiden", "condition"]
    assert captured["kwargs"]["cluster_de_top_n"] == 12
    assert captured["kwargs"]["cluster_de_method"] == "t-test"
    assert captured["kwargs"]["cluster_de_layer"] == "normalized"
    assert captured["kwargs"]["cluster_de_min_cells"] == 7


def test_cli_rejects_invalid_section_rotations(monkeypatch, tmp_path, capsys):
    input_path = tmp_path / "input.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    monkeypatch.setattr(
        cli_module.sys,
        "argv",
        [
            "karospace",
            str(input_path),
            "--section-rotations",
            "bad-token",
        ],
    )

    with pytest.raises(SystemExit) as exc:
        cli_module.main()

    assert exc.value.code == 2
    assert "section_id:angle" in capsys.readouterr().err


def test_gene_correlations_embedded_in_export(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        color="leiden",
        genes=["G1", "G2", "G3"],
        use_hvgs=False,
        gene_storage="embedded",
        gene_correlation_top_n=5,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    gene_corr = embedded.get("gene_correlations")
    assert gene_corr is not None, "gene_correlations key missing from exported data"
    # All three embedded genes should have correlation entries
    assert set(gene_corr.keys()) == {"G1", "G2", "G3"}
    # Each entry is a list of dicts with 'gene' and 'r' keys
    for gene, corrs in gene_corr.items():
        assert isinstance(corrs, list)
        for entry in corrs:
            assert "gene" in entry
            assert "r" in entry
            assert entry["gene"] != gene  # no self-correlation
            assert isinstance(entry["r"], float)
            assert entry["r"] > 0  # only positive correlations included


def test_gene_correlations_disabled_when_top_n_zero(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        color="leiden",
        genes=["G1", "G2", "G3"],
        use_hvgs=False,
        gene_correlation_top_n=0,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    assert embedded.get("gene_correlations") == {}


def test_cli_passes_gene_correlation_top_n_to_export(monkeypatch, tmp_path):
    input_path = tmp_path / "input.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    captured = {}

    def fake_load(path, groupby):
        return _build_dataset()

    def fake_export(dataset, **kwargs):
        captured["kwargs"] = kwargs
        return str(tmp_path / "viewer.html")

    monkeypatch.setattr("karospace.data_loader.load_spatial_data", fake_load)
    monkeypatch.setattr("karospace.exporter.export_to_html", fake_export)
    monkeypatch.setattr(
        cli_module.sys,
        "argv",
        [
            "karospace",
            str(input_path),
            "--gene-correlation-top-n",
            "15",
        ],
    )

    cli_module.main()

    assert captured["kwargs"]["gene_correlation_top_n"] == 15


def test_load_spatial_data_handles_null_encoded_uns_entries(tmp_path):
    obs = pd.DataFrame(
        {"sample_id": pd.Categorical(["S1", "S1"])},
        index=["cell_0", "cell_1"],
    )
    var = pd.DataFrame(index=["G1", "G2"])
    adata = AnnData(X=np.array([[1.0, 0.0], [0.0, 1.0]], dtype=float), obs=obs, var=var)
    adata.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]], dtype=float)
    adata.uns["dummy"] = {"ok": 1}

    path = tmp_path / "broken.h5ad"
    adata.write_h5ad(path)

    with h5py.File(path, "r+") as handle:
        null_ds = handle["uns"].create_dataset("null_entry", data=np.float32(0))
        null_ds.attrs["encoding-type"] = "null"
        null_ds.attrs["encoding-version"] = "0.1.0"

    dataset = load_spatial_data(str(path), groupby="sample_id")

    assert dataset.n_sections == 1
    assert dataset.n_cells == 2
    assert dataset.var_names == ["G1", "G2"]
