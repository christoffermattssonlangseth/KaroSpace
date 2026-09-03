from pathlib import Path

from karospace.exporter import export_to_html
from tests.test_multimodal_export_schema import _make_multimodal_dataset


def _render_multimodal_html(tmp_path=None):
    output_path = (
        Path(tmp_path) / "multimodal-viewer.html"
        if tmp_path is not None
        else Path("/private/tmp/karospace-multimodal-runtime-test.html")
    )
    export_to_html(
        _make_multimodal_dataset(),
        output_path=str(output_path),
        main_cell_annotation="cell_type",
        features=["rna_a"],
        modalities=["rna"],
        pseudobulk=None,
        interaction_markers=None,
        spatial_variable_genes_n=0,
        category_means_n_genes=0,
        gene_correlation_top_n=0,
        pathway_gsea_permutations=0,
        tutorial=False,
    )
    return output_path.read_text(encoding="utf-8")


def test_generated_html_uses_modality_scoped_feature_helpers(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "function resolveCanonicalFeatureName(token, modality" in html
    assert "async function ensureFeatureAvailable(feature, options" in html
    assert "const FEATURE_INDEX_BY_MODALITY = new Map()" in html
    assert "function buildFeatureIndex(modality" in html
    assert "AVAILABLE_GENE_SET" not in html
    assert "resolveCanonicalGeneName" not in html


def test_generated_html_has_panel_scoped_modality_state(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "const PANEL_MODALITY_STATE" in html
    assert "visual:" in html
    assert "exploration:" in html
    assert "pseudobulk:" in html
    assert "interactions:" in html
    assert "let CURRENT_MODALITY" not in html
