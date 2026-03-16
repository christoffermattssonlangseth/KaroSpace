import json
import re

import numpy as np
import pandas as pd
from anndata import AnnData

from karospace.data_loader import SectionData, SpatialDataset
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
