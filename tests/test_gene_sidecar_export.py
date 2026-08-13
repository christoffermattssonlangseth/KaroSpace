import json
import py_compile
import re
import shutil
import struct
import subprocess
import tempfile
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import karospace.cli as cli_module
from karospace.annotations import integrate_polygon_annotations
from karospace.data_loader import SectionData, SpatialDataset, _read_h5ad_with_fallback, load_spatial_data
from karospace.exporter import export_to_html, package_sidecar_viewer


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
    adata.layers["counts"] = x.copy()
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
        metadata_section=["condition"],
    )


def _build_pseudobulk_dataset(cells_per_category_per_sample: int = 20) -> SpatialDataset:
    obs_rows = []
    x_rows = []
    coords = []
    for sample_idx, sample_id in enumerate(["S1", "S2"]):
        condition = "ctrl" if sample_id == "S1" else "treated"
        for category_idx, category in enumerate(["A", "B"]):
            for cell_idx in range(int(cells_per_category_per_sample)):
                obs_rows.append(
                    {
                        "sample_id": sample_id,
                        "leiden": category,
                        "condition": condition,
                    }
                )
                if category == "A":
                    x_rows.append([5.0, float(cell_idx % 3), 0.0])
                else:
                    x_rows.append([1.0, 4.0, float(cell_idx % 2)])
                coords.append(
                    [
                        float(sample_idx * 100 + cell_idx),
                        float(category_idx * 100 + cell_idx),
                    ]
                )

    obs = pd.DataFrame(
        obs_rows,
        index=[f"cell_{i}" for i in range(len(obs_rows))],
    )
    obs["sample_id"] = pd.Categorical(obs["sample_id"])
    obs["leiden"] = pd.Categorical(obs["leiden"])
    obs["condition"] = pd.Categorical(obs["condition"])
    var = pd.DataFrame(index=["G1", "G2", "G3"])
    x = np.asarray(x_rows, dtype=float)
    adata = AnnData(X=x, obs=obs, var=var)
    adata.layers["counts"] = x.copy()
    adata.layers["normalized"] = x.copy()
    adata.obsm["spatial"] = np.asarray(coords, dtype=float)
    sections = [
        SectionData(
            section_id=sample_id,
            coordinates=adata.obsm["spatial"][
                obs["sample_id"].astype(str).to_numpy() == sample_id
            ],
            metadata={"condition": condition},
        )
        for sample_id, condition in [("S1", "ctrl"), ("S2", "treated")]
    ]
    return SpatialDataset(
        adata=adata,
        sections=sections,
        groupby="sample_id",
        obs_columns=["sample_id", "leiden", "condition"],
        var_names=["G1", "G2", "G3"],
        metadata_section=["condition"],
    )


def _extract_data_json(html_text: str) -> dict:
    match = re.search(
        r'<script id="karospace-data" type="application/json">(.*?)</script>',
        html_text,
        re.DOTALL,
    )
    assert match, "embedded data script not found"
    return json.loads(match.group(1).replace("<\\/", "</"))


def _parse_binary_sidecar_index(path: Path) -> tuple[dict, bytes]:
    blob = path.read_bytes()
    assert blob[:4] == b"KSB1"
    version, _flags, gene_count, _reserved = struct.unpack_from("<HHII", blob, 4)
    assert version == 1
    offset = 16
    entries = {}
    for _ in range(gene_count):
        (name_len,) = struct.unpack_from("<H", blob, offset)
        offset += 2
        gene = blob[offset:offset + name_len].decode("utf-8")
        offset += name_len
        payload_kind, _reserved = struct.unpack_from("<BB", blob, offset)
        offset += 2
        payload_offset, payload_length = struct.unpack_from("<QQ", blob, offset)
        offset += 16
        entries[gene] = {
            "payload_kind": payload_kind,
            "payload_offset": payload_offset,
            "payload_length": payload_length,
        }
    return entries, blob


def _attach_companion_analytics(dataset: SpatialDataset, **payloads) -> None:
    companion_uns = {
        "analytics_storage": np.array("json-string-v1"),
        "analytics_columns": np.array([b"leiden"], dtype=object),
    }
    for key, value in payloads.items():
        companion_uns[f"{key}_json"] = json.dumps(value).encode("utf-8")
    dataset.adata.uns["karospace_companion"] = companion_uns


def _write_h5ad_with_nested_null_uns_entry(path: Path) -> None:
    obs = pd.DataFrame(
        {"sample_id": pd.Categorical(["S1", "S1"])},
        index=["cell_0", "cell_1"],
    )
    var = pd.DataFrame(index=["G1", "G2"])
    adata = AnnData(X=np.array([[1.0, 0.0], [0.0, 1.0]], dtype=float), obs=obs, var=var)
    adata.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]], dtype=float)
    adata.uns["dataset_source"] = {
        "chunksize": 256,
        "expression_path": "expr.h5",
        "metadata_path": "meta.csv",
    }
    adata.write_h5ad(path)

    with h5py.File(path, "r+") as handle:
        dataset_source = handle["uns"]["dataset_source"]
        null_ds = dataset_source.create_dataset("selected_fovs", data=np.float32(0))
        null_ds.attrs["encoding-type"] = "null"
        null_ds.attrs["encoding-version"] = "0.1.0"


def test_export_can_disable_pseudobulk_and_interaction_markers(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=["G1"],
        pseudobulk=None,
        interaction_markers=None,
        neighbor_stats_groupby=[],
        neighbor_stats_permutations=0,
        spatial_variable_genes_n=0,
        cluster_means_n_genes=0,
        gene_correlation_top_n=0,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    assert embedded["pseudobulk_de"] == {}
    assert embedded["interaction_markers"] == {}
    assert embedded["marker_genes"] == {}
    assert embedded["cluster_gene_means"] is None
    assert embedded["gene_correlations"] == {}


def test_sidecar_export_writes_aux_and_updates_html_contract(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=["G1"],
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
    shard_files = sorted(shard_dir.glob("*.bin"))
    assert shard_files
    entries, _blob = _parse_binary_sidecar_index(shard_files[0])

    assert embedded["available_genes"] == ["G1", "G2", "G3"]
    assert embedded["gene_aux_url"] == "viewer.genes.json"
    assert embedded["embedded_genes"] == []
    assert embedded["genes_meta"] == {}
    assert all(not (section.get("genes") or section.get("genes_sparse")) for section in embedded["sections"])
    assert set(manifest["genes_meta"]) == {"G1", "G2", "G3"}
    assert manifest["format"] == "karospace-gene-sidecar-manifest-v3"
    assert manifest["gene_sidecar_format"] == "binary"
    assert manifest["gene_to_shard"]["G1"].startswith("viewer.genes/")
    assert manifest["gene_to_shard"]["G2"].startswith("viewer.genes/")
    assert manifest["gene_to_shard"]["G2"].endswith(".bin")
    assert manifest["gene_value_encodings"]["G1"] == "uint16"
    assert manifest["gene_value_encodings"]["G2"] == "uint16"
    assert "G1" in entries
    assert "G2" in entries
    assert "async function ensureGeneAvailable" in html_text
    assert "function requestModalBlendGene(gene, modality = null)" in html_text
    assert "requestModalBlendGene(gene, modName);" in html_text
    assert html_text.index('id="zoom-out"') < html_text.index('id="zoom-in"') < html_text.index('id="zoom-info"')
    assert 'id="focused-modal-neighbor-toggle"' in html_text
    assert 'id="focused-modal-he-toggle"' in html_text
    assert 'id="focused-neighbor-panel"' in html_text
    assert 'id="focused-he-panel"' in html_text
    assert 'id="modal-controls-toggle"' not in html_text
    assert ".modal-controls.hidden" not in html_text
    assert 'data-modal-group="view"' not in html_text
    assert "function updateModalToolbarState()" in html_text
    assert "function initModalControlsDragging()" not in html_text
    assert 'id="modal-exit-subview-btn"' in html_text
    assert "function activateModalSubviewFromSelection()" in html_text
    assert "function exitModalSubview()" in html_text
    assert "function buildModalSubviewStateFromSelection(section = modalSection)" in html_text
    assert "Open focused view" in html_text
    assert "Focused view active. Use Back view to return to the full section." in html_text
    assert 'id="modal-gene-panel"' in html_text
    assert "Genes in selection" in html_text
    assert "function showModalGeneDiscoveryPanel(sectionId, selectedCellIndices)" in html_text
    assert "function hideModalGeneDiscoveryPanel()" in html_text
    assert "function scoreModalSelectionMarkerGenes(section, selectedCellIndices)" in html_text
    assert "function scoreLoadedGenesForModalSelection(section, selectedCellIndices)" in html_text
    assert "function createModalAnnotationFromSelection()" in html_text
    assert 'id="selection-create-annotation-btn"' in html_text
    assert 'id="selection-create-annotation-btn" type="button" title="Create annotation from this lasso selection" aria-label="Create annotation from this lasso selection" hidden' in html_text
    assert "function markPrimarySelectionChanged()" in html_text
    assert "function annotationAlreadyCreatedForCurrentSelection()" in html_text
    assert "function updateSelectionCreateAnnotationButtonState()" in html_text
    assert "function hasLassoDerivedAnnotationSelection()" in html_text
    assert "if (parentId === null || parentId === undefined || parentId === '') return null;" in html_text
    assert "Annotation already created for this selection" in html_text
    assert 'id="focused-modal-annotation-toggle"' not in html_text
    assert 'id="focused-modal-graph-toggle"' not in html_text
    assert 'id="umap-compare-toggle"' not in html_text
    assert 'id="umap-compare-wrap"' not in html_text
    assert "createModalAnnotationFromPath" not in html_text
    assert "modalAnnotationModeActive" not in html_text
    assert 'id="overview-gene-view-mode"' in html_text
    assert "function updateOverviewGeneViewState()" in html_text
    assert "function ensureSectionGeneDensityCache(section, transform, width, height, values, candidateIndices = null)" in html_text
    assert "function drawGeneDensityLayer(ctx, width, height, cache)" in html_text
    assert "class=\"modal-gene-panel-entry\"" in html_text
    assert "Search Google for this selection gene" in html_text
    assert "id=\"selection-show-genes-btn\"" in html_text
    assert "Genes in selection" in html_text
    assert 'id="genes-tab-spatial"' in html_text
    assert 'id="genes-tab-dotplot"' not in html_text
    assert 'id="compare-tab-cell-de"' in html_text
    assert 'id="dotplot-groupby"' not in html_text
    assert 'id="dotplot-genes"' not in html_text
    assert 'id="dotplot-grid"' not in html_text
    assert ">Dotplot<" not in html_text
    assert ">Spatial<" in html_text
    assert ">Cell DE<" in html_text
    assert "function renderSpatialVariableGenes()" in html_text
    assert "function getBarColor(" not in html_text
    assert 'data-marker-category="${escapeHtml(key)}"' in html_text
    assert 'title="Click to highlight this annotation in the viewer"' in html_text
    assert "Global Moran&apos;s I ranking precomputed at export." in html_text
    assert "Neighborhood Follow-up" not in html_text
    assert 'data-cluster-de-open="neighbors-enrichment"' not in html_text
    assert 'data-cluster-de-open="neighbors-interactions"' not in html_text
    assert "function renderComparisonCountSummary(" not in html_text
    assert "renderComparisonNeighborSummary(clusterDeGroupby, clusterDeSourceCategory, clusterDeReferenceCategory)" not in html_text
    assert "renderComparisonInteractionSummary(clusterDeGroupby, clusterDeSourceCategory, clusterDeReferenceCategory)" not in html_text
    assert "Pairwise DE for <strong>" not in html_text
    assert "let modalSpacePanActive = false;" in html_text
    assert "Pan inside the modal while Select is active." in html_text
    assert 'id="shortcuts-overlay"' in html_text
    assert "function initShortcutsOverlay()" in html_text
    assert "function toggleShortcutsOverlay()" in html_text
    assert "if (key === '?')" in html_text
    assert "Open the full keyboard shortcuts overlay." in html_text
    assert 'id="tutorial-overlay"' in html_text
    assert 'id="tutorial-trigger"' in html_text
    assert "function initTutorialOverlay()" in html_text
    assert "const tutorialSteps = [" in html_text
    assert "function buildDetailedTutorialSteps()" in html_text
    assert "Welcome the KaroSpace" in html_text
    assert "Recommended DE workflow" in html_text
    assert "I tried it" not in html_text
    assert 'id="tutorial-skip"' not in html_text
    assert "tutorial-task" in html_text
    assert html_text.count("step('") > 100
    embedded = _extract_data_json(html_text)
    assert embedded["tutorial"]["enabled"] is False
    assert embedded["tutorial"]["autostart"] is False
    assert "--selection-outline-color:" in html_text
    assert "function getSelectionOutlineColor()" in html_text
    assert "const metadataColor = getMetadataValueColor(OUTLINE_BY, value);" in html_text
    assert "function reorderModalAnnotation(sourceId, targetId)" in html_text
    assert "function moveModalAnnotation(sourceId, targetId, position = 'before')" in html_text
    assert "function moveModalAnnotationToRoot(sourceId)" in html_text
    assert "function createModalAnnotationGroup()" in html_text
    assert "parent_id: normalizeAnnotationParentId(annotation.parentId)" in html_text
    assert "type: annotation.isGroup ? 'group' : 'polygon'" in html_text
    assert 'id="modal-annotations-create-group"' in html_text
    assert "modal-annotation-children" not in html_text
    assert "modal-annotation-group-wrap" in html_text
    assert "data-annotation-group-id" in html_text
    assert "modal-annotation-group-placeholder" in html_text
    assert "annotation-group" in html_text
    assert "annotation-child" in html_text
    assert "drag-inside" in html_text
    assert "data-annotation-drag-handle" in html_text
    assert "event.target?.closest?.('input, button, select, textarea, a, [contenteditable=\"true\"]')" in html_text
    assert '.modal-annotation-row[draggable="true"]' in html_text
    assert 'id="focused-neighbor-panel"' in html_text
    assert 'id="focused-neighbor-panel-graph"' in html_text
    assert 'id="focused-neighbor-panel-neighbors"' in html_text
    assert 'id="focused-neighbor-hop-select"' in html_text
    assert "function updateFocusedNeighborPanelState()" in html_text
    assert "hopRow.hidden = !neighborHoverEnabled" in html_text
    assert ".focused-neighbor-param-row[hidden]" in html_text
    assert "panel.style.top = `${toggle.offsetTop}px`" in html_text
    assert '<button class="graph-toggle" id="graph-toggle"' not in html_text
    assert '<button class="graph-toggle" id="neighbor-hover-toggle"' not in html_text
    assert '<select id="neighbor-hop-select"' not in html_text
    assert ">Keyboard Shortcuts<" in html_text
    assert "<kbd>/</kbd>" in html_text

    assert "<kbd>I</kbd>" in html_text
    assert ".filter-group { display: flex; align-items: center; gap: 4px; flex-wrap: wrap; min-width: 0;" in html_text
    assert ".filter-group.outline-filter-group { margin-left: auto; }" in html_text
    assert "@media (max-width: 720px)" in html_text
    assert ".filter-group.outline-filter-group { margin-left: 0; width: 100%; }" in html_text
    assert "function toggleLegendPanel()" in html_text
    assert "function toggleInsightsPanel()" in html_text
    assert "function navigateModalSection(step)" in html_text
    assert "function initKeyboardShortcuts()" in html_text
    assert "function isEditableKeyboardTarget(target)" in html_text
    assert "function setModalButtonPriority(button, priority = 'default')" in html_text
    assert "function scheduleModalBlendRender()" in html_text
    assert "function updateModalBlendLabels()" in html_text
    assert "function ensureModalBlendRenderCache(section, transform, width, height, dpr, adjustedSpotSize, runtimes, candidateIndices = null)" in html_text
    assert "function renderModalBlendFromCache()" in html_text
    assert ".modal-controls button.modal-btn-primary" not in html_text
    assert ".modal-controls button.modal-btn-muted" not in html_text
    assert ".modal-controls button.modal-btn-danger" not in html_text
    assert ".focused-modal-tool-toggle" in html_text
    assert "buildColorPanel();" in html_text
    assert "colorToggle.addEventListener('click', toggleInsightsPanel);" in html_text
    assert "if (key === 'i' || key === 'I')" in html_text
    assert 'id="insights-footer-bar"' not in html_text
    assert '<div class="color-panel collapsed" id="color-panel"></div>' in html_text
    assert 'id="color-toggle"' in html_text
    assert "document.getElementById('legend-toggle').addEventListener('click', toggleLegendPanel);" in html_text
    assert "initKeyboardShortcuts();" in html_text
    assert ">Insights<" in html_text
    assert "function createViewerStorage()" in html_text
    assert "const VIEWER_STORAGE = createViewerStorage();" in html_text
    assert "VIEWER_STORAGE.getItem('spatial-viewer-theme')" in html_text
    assert "VIEWER_STORAGE.setItem('spatial-viewer-theme', currentTheme)" in html_text
    assert "--accent-text:" in html_text
    assert "--accent-fill:" in html_text
    assert "--panel-card-bg:" in html_text
    assert "--warning-bg:" in html_text
    assert "--split-a-text:" in html_text
    assert "button { color: var(--text-color); }" in html_text
    assert ":root.dark button" not in html_text
    assert "background: var(--panel-card-bg);" in html_text
    assert "function getAnnotationCountChipStyle(annotation)" in html_text
    assert "getTextColorForBackground(color)" in html_text
    assert "window.location.protocol === 'file:' && !window.__karospacePackageMode" in html_text
    assert 'id="modal-blend-loading"' in html_text
    assert 'id="gene-discovery-panel"' in html_text
    assert 'id="gene-panel-new"' in html_text
    assert 'id="color-tab-overview"' in html_text
    assert 'id="color-tab-genes"' in html_text
    assert 'id="color-tab-compare"' in html_text
    assert 'id="color-tab-neighbors"' in html_text
    assert ">Overview<" in html_text
    assert ">Summary<" in html_text
    assert ">Sections<" in html_text
    assert ">Neighbors<" in html_text
    assert ">Enrichment<" in html_text
    assert ">Interactions<" in html_text
    assert "function getGeneSearchResults(query, limit = GENE_DISCOVERY_MAX_RESULTS)" in html_text
    assert "function renderGeneDiscoveryPanel()" in html_text
    assert "function renderClusterDE()" in html_text
    assert "function isEmbeddedViewerGene(gene)" in html_text
    assert "Load marker gene into the viewer" in html_text
    assert "Click a marker gene to load it in the viewer." in html_text
    assert "Marker gene name only; this gene was not embedded in the viewer" in html_text
    assert "function applyModalInteractionPreview()" in html_text
    assert "function scheduleModalInteractionCommit(delayMs = 120)" in html_text
    assert "function ensureSectionSpatialIndex(section)" in html_text
    assert "function querySectionSpatialIndex(section, bbox)" in html_text
    assert "function getAvailableComparisonColors()" in html_text
    assert "function renderComparisonCountSummary(colorCol, sourceCategory, referenceCategory)" not in html_text
    assert "function renderComparisonNeighborSummary(colorCol, sourceCategory, referenceCategory)" in html_text
    assert "function renderComparisonInteractionSummary(colorCol, sourceCategory, referenceCategory)" in html_text
    assert "function renderClusterDEResultSection(colorCol, sourceCategory, referenceCategory)" in html_text
    assert "Find cells" in html_text
    assert 'data-selection-query-toggle' in html_text
    assert 'data-selection-query-input' in html_text
    assert "const bareTokenHasLetter = /[A-Za-z_]/.test(bareToken);" in html_text
    assert "tokens.push({ type: 'identifier', value: bareToken, pos: i });" in html_text
    assert 'data-selection-query-run' in html_text
    assert 'data-selection-query-add' in html_text
    assert 'data-selection-query-clear' in html_text
    assert "Click an example to fill the box, then edit it." in html_text
    assert "function tokenizeSelectionQuery(input)" in html_text
    assert "function parseSelectionQuery(input)" in html_text
    assert "function getSelectionQueryExamplePresets()" in html_text
    assert "function getSelectionQueryPlaceholder()" in html_text
    assert "Find cells by annotation, gene value, or section metadata." in html_text
    assert "function normalizeSelectionQueryAst(node)" in html_text
    assert "function evaluateSelectionQueryNode(node, section, cellIdx)" in html_text
    assert "async function runSelectionQuery(mode = 'replace')" in html_text
    assert "const MAX_SELECTION_QUERY_MATCHES = 150000;" in html_text
    assert "function getFilteredSections()" in html_text
    assert "function trimSelectionSetToFilteredSections(cells)" in html_text
    assert "function trimSelectionsToFilteredSections()" in html_text
    assert "function applyMetadataFilters()" in html_text
    assert "function isMissingDisplayValue(value)" in html_text
    assert "function getCategoricalValueInfo(config, value)" in html_text
    assert "const filteredSections = getFilteredSections();" in html_text
    assert "Most columns and genes work as bare names." in html_text
    assert "Use <code>obs(...)</code>, <code>gene(...)</code>," in html_text
    assert 'data-neighbor-view="table"' in html_text
    assert 'data-neighbor-view="bubble"' in html_text
    assert 'data-neighbor-view="chord"' in html_text
    assert "let neighborStatsView = 'table';" in html_text
    assert "function getNeighborStatsViewState()" in html_text
    assert "function renderNeighborNetworkView(viewState)" in html_text
    assert "function renderNeighborChordDiagram(viewState)" in html_text
    assert "function renderNeighborStatsTableView(viewState)" in html_text
    assert "function renderNeighborVisualizationToolbar(viewKey, label)" in html_text
    assert "function renderNeighborVisualizationControlPanel(viewKey, viewState)" in html_text
    assert "function bindNeighborVisualizationSettingControls(container, viewKey)" in html_text
    assert "function getNeighborMetricOptions(viewKey)" in html_text
    assert "function setNeighborMetric(viewKey, value)" in html_text
    assert "function getNeighborFocusedIndices(viewState, metric, controls)" in html_text
    assert "function sortNeighborIndicesForChord(indices, viewState, order)" in html_text
    assert "function applyNeighborVisualizationViewport(container, viewKey, options = {})" in html_text
    assert "function bindNeighborVisualizationControls(container, viewKey)" in html_text
    assert "function renderNeighborFocusDetail(focus, viewState)" in html_text
    assert "function bindNeighborVisualizationInteractions(container, viewState)" in html_text
    assert 'id="neighbor-focus-mode"' in html_text
    assert 'id="neighbor-metric-mode"' in html_text
    assert 'id="neighbor-threshold-range"' in html_text
    assert 'id="neighbor-max-categories"' in html_text
    assert 'id="neighbor-label-mode"' in html_text
    assert 'id="neighbor-chord-order"' in html_text
    assert 'data-neighbor-detail' in html_text
    assert 'data-neighbor-auto-sync="1"' in html_text
    assert 'data-neighbor-scroll' in html_text
    assert 'data-neighbor-svg' in html_text
    assert 'data-neighbor-zoom-action="in"' in html_text
    assert 'data-neighbor-zoom-action="expand"' in html_text
    assert "Neighbor enrichment network" in html_text
    assert "Neighbor connection chord diagram" in html_text
    assert "Selected + neighbors" in html_text
    assert ">Min strength<" in html_text
    assert ">Max categories<" in html_text
    assert ">Labels<" in html_text
    assert ">Order<" in html_text
    assert "function computeRegionAnnotationDE(annotationA, annotationB, options = {})" in html_text
    assert "function fetchGeneAuxShardForAnalysis(shardUrl)" in html_text
    assert "function decodeDenseSidecarSection(sectionEntry)" in html_text
    assert "function base64ToUint8Array(b64)" in html_text
    assert "|| getActiveFeatureList()[0]" not in html_text
    assert "function buildVolcanoPlot(genes, log2fc, pvalsAdj, padjCutoff = 0.05, log2fcCutoff = 0.5, colorCol = null, sourceCategory = null, referenceCategory = null)" in html_text
    assert "function buildMAPlot(genes, baseMean, log2fc, pvals, pvalsAdj, padjCutoff = 0.05, log2fcCutoff = 0.5)" in html_text
    assert "const pvals = Array.isArray(result.pvals) ? result.pvals : [];" in html_text
    assert "buildMAPlot(genes, baseMean, log2fc, pvals, pvalsAdj, padjCutoff, log2fcCutoff)" in html_text
    assert "baseMean" in html_text
    assert "baseMean + 1" not in html_text
    assert "Pseudobulk Analysis Diagnostics" in html_text
    assert "buildClusterDESummaryItem('p < 0.05', '#2f9e67', counts.green)" in html_text
    assert "buildClusterDESummaryItem('p >= 0.05', 'var(--border-color)', counts.grey)" in html_text
    assert "function buildPseudobulkPCAPlot(sampleInfo, colorCol)" in html_text
    assert "function buildPseudobulkDistanceHeatmap(sampleInfo, colorCol)" in html_text
    assert "function renderClusterPAResultSection(colorCol, sourceCategory, referenceCategory)" in html_text
    assert 'class="cluster-PA-result cluster-pa-result"' in html_text
    assert "Gene Enrichment" in html_text
    assert "data-pathway-annotation-select" in html_text
    assert 'data-pathway-annotation-panel="source"' in html_text
    assert 'data-pathway-annotation-panel="reference"' in html_text
    assert "function buildPathwayORAPlot(result, sourceCategory, referenceCategory, sourceColor, referenceColor)" in html_text
    assert "function buildPathwayORADotPlot(rows, sourceCategory, referenceCategory, sourceColor, referenceColor)" in html_text
    assert "GeneRatio" in html_text
    assert "-log10 adjusted p" in html_text
    assert "function buildPathwayGSEAPlot(result, sourceCategory, referenceCategory, sourceColor, referenceColor)" in html_text
    assert "function buildPathwayGSEAEnrichmentPlot(row, sourceCategory, referenceCategory, sourceColor, referenceColor)" in html_text
    assert "GSEA Enrichment" in html_text
    assert "data-pathway-gsea-select" in html_text
    assert "data-gsea-pathway-panel" in html_text
    assert "function buildPathwayBarPlot" not in html_text
    assert "GSEA Pathways" not in html_text
    assert "pathway_enrichment" in html_text
    assert "ORA Pathways" in html_text
    assert "Rank in ordered dataset" in html_text
    assert "const greyLabel = 'adj. p > ' + formatAdjustedPValue(padjThreshold)" in html_text
    assert "+ ' or |log2FC| < ' + formatScaleNumber(fcThr)" in html_text
    assert "function wrapClusterDEVolcanoPlot(volcanoHtml)" in html_text
    assert "buildVolcanoPlot(genes, log2fc, pvalsAdj, padjCutoff, log2fcCutoff, colorCol, sourceCategory, referenceCategory)" in html_text
    assert "function bindClusterDEPlotInteractions(container, rerenderFn)" in html_text
    assert "function buildGroupVolcanoPlot(entries)" in html_text
    assert "function bindVolcanoGroupInteraction(container, rerenderFn)" in html_text
    assert "function downloadCurrentClusterDEVolcanoSvg(container)" in html_text
    assert "function downloadCurrentClusterDEVolcanoPng(container)" in html_text
    assert "function downloadCurrentClusterDETable(format)" in html_text
    assert 'class="volcano-container"' in html_text
    assert 'class="cluster-de-plot-grid"' in html_text
    assert 'class="volcano-summary"' in html_text
    assert 'class="volcano-summary-chip"' in html_text
    assert 'class="volcano-summary-count"' in html_text
    assert 'class="volcano-tooltip"' in html_text
    assert "function formatAdjustedPValue(value)" in html_text
    assert "formatAdjustedPValue(pvalAdjValue)" in html_text
    assert 'data-tooltip-line2="adj. p: ' in html_text
    assert 'data-most-expressed="' in html_text
    assert "colorToRgbaCss(rowColor, 0.12)" in html_text
    assert "Save Volcano SVG" in html_text
    assert "Save Volcano PNG" in html_text
    assert "Save Table CSV" in html_text
    assert 'class="cluster-de-volcano-container"' in html_text
    assert "Save Table Excel" in html_text
    assert "Number.MIN_VALUE" in html_text
    assert "if (window.__karospacePackageMode) {" in html_text
    assert "function runFullRegionAnnotationDE(annotationA, annotationB)" in html_text
    assert "function exportAnnotationDEReport(annotationA, annotationB, exportState)" in html_text
    assert "function exportAnnotationDECsv(annotationA, annotationB, exportState)" in html_text
    assert "function computeCellSetDE(groupA, groupB, options = {})" in html_text
    assert "function runFullGroupDE(groupA, groupB)" in html_text
    assert "function renderGroupDE()" in html_text
    assert "function renderCompareInsights()" in html_text
    assert 'id="cluster-de-summary"' in html_text
    assert 'id="group-de-summary"' in html_text
    assert 'id="group-de-panel"' in html_text
    assert 'id="compare-tab-groups"' in html_text
    assert 'id="compare-tab-regions"' in html_text
    assert 'id="annotation-comparison"' in html_text
    assert 'id="group-de-source-spec"' in html_text
    assert 'id="group-de-source-value"' in html_text
    assert 'id="group-de-reference-value"' in html_text
    assert 'id="group-de-restrict-spec"' in html_text
    assert 'id="group-de-restrict-value"' in html_text
    assert 'id="group-de-scope"' in html_text
    assert "karospace-group-de-report-v1" in html_text
    assert 'if (/[",\\n\\r]/.test(text)) {' in html_text
    assert "function renderGeneGoogleSearchButton(gene, options = {})" in html_text
    assert "function renderAnnotationRegionDESection(annotations)" in html_text
    assert 'id="annotation-de-source"' in html_text
    assert 'id="annotation-de-reference"' in html_text
    assert 'id="annotation-de-topn"' in html_text
    assert 'id="annotation-de-run-full"' in html_text
    assert 'id="annotation-de-refresh-full"' in html_text
    assert 'id="annotation-de-cancel"' in html_text
    assert 'id="annotation-de-export-json"' in html_text
    assert 'id="annotation-de-export-csv"' in html_text
    assert 'data-gene-google-search="' in html_text
    assert "Region-to-Region DE" in html_text
    assert "Run Full Region DE" in html_text
    assert "Refresh Full DE" in html_text
    assert "Export JSON" in html_text
    assert "Export CSV" in html_text
    assert "Cached full sidecar DE" in html_text
    assert "karospace-region-de-report-v1" in html_text
    assert "Category A" in html_text
    assert "Category B" in html_text
    assert "function recordRecentGene(gene)" in html_text
    assert "function loadSavedGenePanels()" in html_text
    assert "const GENE_RECENTS_STORAGE_KEY = getViewerScopedStorageKey('gene-recents');" in html_text
    assert html_text.index("function getMarkerGenesForColorCategory(colorCol, category)") < html_text.index("function getGeneSuggestionGroups()")
    assert "was not pre-loaded" not in html_text

    if shutil.which("node"):
        script_matches = list(re.finditer(r"<script([^>]*)>(.*?)</script>", html_text, re.DOTALL | re.IGNORECASE))
        runtime_scripts = [
            match.group(2)
            for match in script_matches
            if 'application/json' not in (match.group(1) or '')
        ]
        assert runtime_scripts, "expected viewer runtime script"
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", suffix=".js", delete=False) as handle:
            handle.write(runtime_scripts[-1])
            script_path = handle.name
        try:
            subprocess.run(
                ["node", "--check", script_path],
                check=True,
                capture_output=True,
                text=True,
            )
        finally:
            Path(script_path).unlink(missing_ok=True)


def test_export_can_enable_guided_tutorial(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=["G1"],
        tutorial=True,
    )

    html_text = output_path.read_text(encoding="utf-8")
    embedded = _extract_data_json(html_text)
    assert embedded["tutorial"]["enabled"] is True
    assert embedded["tutorial"]["autostart"] is True
    assert embedded["tutorial"]["storage_key"].endswith(":viewer:v1")
    assert 'id="tutorial-trigger"' in html_text
    assert 'style=""' in html_text
    assert "Do not open automatically again" in html_text


def test_sidecar_export_with_no_embedded_genes_keeps_gene_catalog_and_warns(tmp_path, capsys):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=[],
        gene_storage="sidecar",
    )

    html_text = output_path.read_text(encoding="utf-8")
    embedded = _extract_data_json(html_text)
    manifest = json.loads((tmp_path / "viewer.genes.json").read_text(encoding="utf-8"))
    output = capsys.readouterr().out

    assert embedded["available_genes"] == ["G1", "G2", "G3"]
    assert embedded["genes_meta"] == {}
    assert all(not (section.get("genes") or section.get("genes_sparse")) for section in embedded["sections"])
    assert set(manifest["genes_meta"]) == {"G1", "G2", "G3"}
    assert "function getAvailableFeaturesForModality" in html_text
    assert "const availableGenes = DATA.available_genes || [];" in html_text
    assert "return getAvailableFeaturesForModality(CURRENT_MODALITY);" in html_text
    assert "Empty gene shard response from ${shardUrl}" in html_text
    assert "Invalid JSON in gene shard ${shardUrl}" in html_text
    assert "0 genes embedded in HTML; gene expression is sidecar-only" in output
    assert "sidecar viewers must be opened over HTTP(S)" in output


def test_sidecar_export_respects_custom_shard_size(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=["G1"],
        gene_storage="sidecar",
        gene_sidecar_shard_size=1,
    )

    aux_path = tmp_path / "viewer.genes.json"
    manifest = json.loads(aux_path.read_text(encoding="utf-8"))
    shard_files = sorted((tmp_path / "viewer.genes").glob("*.bin"))

    assert len(shard_files) == 2
    assert len(manifest["shards"]) == 2
    assert manifest["gene_to_shard"]["G2"] != manifest["gene_to_shard"]["G3"]


def test_sidecar_auto_encoding_prefers_compact_sparse_and_packed_dense():
    dataset = _build_dataset()

    sidecar = dataset.to_gene_sidecar_data(
        genes=["G1", "G3"],
        gene_encoding="auto",
    )

    assert sidecar["gene_encodings"]["G1"] == "dense"
    assert sidecar["gene_encodings"]["G3"] == "sparse"

    g1_s1 = sidecar["genes"]["G1"]["sections"]["S1"]
    g3_s1 = sidecar["genes"]["G3"]["sections"]["S1"]

    assert "dq16b64" in g1_s1
    assert "dense" not in g1_s1
    assert "sparse" in g3_s1


def test_sidecar_uint16_quantization_emits_compact_value_payloads(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=[],
        gene_storage="sidecar",
        gene_value_encoding="uint16",
    )

    manifest = json.loads((tmp_path / "viewer.genes.json").read_text(encoding="utf-8"))
    shard_path = sorted((tmp_path / "viewer.genes").glob("*.bin"))[0]
    entries, _blob = _parse_binary_sidecar_index(shard_path)

    assert manifest["gene_value_encodings"]["G1"] == "uint16"
    assert manifest["gene_value_encodings"]["G3"] == "uint16"
    assert set(entries) == {"G1", "G2", "G3"}


def test_sidecar_uint8_quantization_emits_compact_value_payloads(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=[],
        gene_storage="sidecar",
        gene_value_encoding="uint8",
    )

    manifest = json.loads((tmp_path / "viewer.genes.json").read_text(encoding="utf-8"))
    shard_path = sorted((tmp_path / "viewer.genes").glob("*.bin"))[0]
    entries, _blob = _parse_binary_sidecar_index(shard_path)

    assert manifest["gene_value_encodings"]["G1"] == "uint8"
    assert manifest["gene_value_encodings"]["G3"] == "uint8"
    assert set(entries) == {"G1", "G2", "G3"}


def test_sidecar_uint8_quantization_supports_sparse_payloads():
    dataset = _build_dataset()

    sidecar = dataset.to_gene_sidecar_data(
        genes=["G3"],
        gene_encoding="sparse",
        gene_value_encoding="uint8",
    )

    g3_s1 = sidecar["genes"]["G3"]["sections"]["S1"]["sparse"]

    assert sidecar["gene_value_encodings"]["G3"] == "uint8"
    assert "vb64" not in g3_s1
    assert "vq8b64" in g3_s1


def test_binary_sidecar_export_writes_indexed_bin_shards(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=["G1"],
        gene_storage="sidecar",
        gene_value_encoding="uint8",
    )

    html_text = output_path.read_text(encoding="utf-8")
    embedded = _extract_data_json(html_text)
    manifest = json.loads((tmp_path / "viewer.genes.json").read_text(encoding="utf-8"))
    shard_path = sorted((tmp_path / "viewer.genes").glob("*.bin"))[0]
    entries, blob = _parse_binary_sidecar_index(shard_path)

    assert embedded["gene_aux_url"] == "viewer.genes.json"
    assert manifest["format"] == "karospace-gene-sidecar-manifest-v3"
    assert manifest["gene_sidecar_format"] == "binary"
    assert manifest["section_order"] == ["S1", "S2"]
    assert manifest["gene_to_shard"]["G2"].endswith(".bin")
    assert set(entries) == {"G2", "G3"}
    assert entries["G2"]["payload_offset"] >= 16
    assert entries["G2"]["payload_offset"] + entries["G2"]["payload_length"] <= len(blob)
    assert "function parseBinaryShardIndexBuffer(buffer)" in html_text
    assert "async function readBinaryGenePayload(shardUrl, gene, indexInfo = null)" in html_text
    assert "function hydrateGeneFromBinary(gene, geneEntry, modality = null)" in html_text


def test_binary_karospace_package_stores_bin_shards_uncompressed(tmp_path):
    dataset = _build_dataset()
    package_path = tmp_path / "viewer.karospace"

    export_to_html(
        dataset,
        output_path=str(package_path),
        main_cells_annotation="leiden",
        genes=["G1"],
        gene_storage="sidecar",
        gene_value_encoding="uint8",
    )

    with zipfile.ZipFile(package_path) as zf:
        infos = {info.filename: info for info in zf.infolist()}
        assert "karospace-package.json" in infos
        bin_names = sorted(name for name in infos if name.endswith(".bin"))
        assert bin_names
        manifest = json.loads(zf.read("karospace-package.json").decode("utf-8"))
        gene_manifest_name = manifest["viewer"]["gene_manifest_path"]
        gene_manifest = json.loads(zf.read(gene_manifest_name).decode("utf-8"))

        assert gene_manifest["format"] == "karospace-gene-sidecar-manifest-v3"
        assert gene_manifest["gene_sidecar_format"] == "binary"
        assert infos[bin_names[0]].compress_type == zipfile.ZIP_STORED
        assert infos[gene_manifest_name].header_offset < infos[bin_names[0]].header_offset
        assert infos["index.html"].header_offset < infos[bin_names[0]].header_offset


def test_karospace_package_export_wraps_sidecar_assets(tmp_path):
    dataset = _build_dataset()
    package_path = tmp_path / "viewer.karospace"
    loader_path = tmp_path / "viewer.loader.html"

    returned = export_to_html(
        dataset,
        output_path=str(package_path),
        main_cells_annotation="leiden",
        genes=["G1"],
        gene_encoding="sparse",
        gene_storage="sidecar",
    )

    assert returned == str(package_path.resolve())
    assert package_path.exists()
    assert loader_path.exists()
    assert not (tmp_path / "index.html").exists()
    assert not (tmp_path / "viewer.genes.json").exists()
    assert not (tmp_path / "viewer.genes").exists()

    loader_html = loader_path.read_text(encoding="utf-8")
    assert "<title>KaroSpace</title>" in loader_html
    assert 'id="loading-label">KaroSpace</div>' in loader_html
    assert "drop a .karospace file to open" in loader_html
    assert "window.__karospacePackageSessions" in loader_html

    with zipfile.ZipFile(package_path) as zf:
        infos = {info.filename: info for info in zf.infolist()}
        names = set(infos)
        assert "karospace-package.json" in names
        assert "index.html" in names
        assert "viewer.genes.json" in names
        shard_names = sorted(name for name in names if name.startswith("viewer.genes/") and name.endswith(".bin"))
        assert shard_names

        package_manifest = json.loads(zf.read("karospace-package.json").decode("utf-8"))
        html_text = zf.read("index.html").decode("utf-8")
        embedded = _extract_data_json(html_text)
        manifest = json.loads(zf.read("viewer.genes.json").decode("utf-8"))

    assert package_manifest["format"] == "karospace-package-v1"
    assert package_manifest["entry_html"] == "index.html"
    assert package_manifest["viewer"]["mode"] == "sidecar-package"
    assert package_manifest["viewer"]["gene_storage"] == "sidecar"
    assert package_manifest["viewer"]["gene_manifest_path"] == "viewer.genes.json"
    assert package_manifest["viewer"]["gene_shard_dir"] == "viewer.genes"
    assert "index.html" in package_manifest["files"]
    assert infos[shard_names[0]].compress_type == zipfile.ZIP_STORED
    assert infos["viewer.genes.json"].header_offset < infos[shard_names[0]].header_offset
    assert infos["index.html"].header_offset < infos[shard_names[0]].header_offset


def test_package_sidecar_viewer_wraps_existing_sidecar_bundle(tmp_path, capsys):
    dataset = _build_dataset()
    html_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(html_path),
        main_cells_annotation="leiden",
        genes=["G1"],
        gene_storage="sidecar",
    )

    package_path = tmp_path / "viewer.karospace"
    returned = package_sidecar_viewer(html_path, output_path=package_path)
    output = capsys.readouterr().out

    assert returned == str(package_path.resolve())
    assert package_path.exists()
    assert (tmp_path / "viewer.loader.html").exists()
    assert "reading existing viewer HTML" in output
    assert "resolved gene manifest" in output
    assert "resolved shard directory" in output
    assert "staging package assets" in output
    assert "writing .karospace archive" in output
    assert "writing local package loader" in output

    with zipfile.ZipFile(package_path) as zf:
        names = set(zf.namelist())
        assert "karospace-package.json" in names
        assert "index.html" in names
        assert "viewer.genes.json" in names
        shard_names = sorted(name for name in names if name.startswith("viewer.genes/"))
        assert shard_names
        package_manifest = json.loads(zf.read("karospace-package.json").decode("utf-8"))
        embedded = _extract_data_json(zf.read("index.html").decode("utf-8"))
        manifest = json.loads(zf.read("viewer.genes.json").decode("utf-8"))
        assert shard_names[0].endswith(".bin")

    assert package_manifest["entry_html"] == "index.html"
    assert package_manifest["title"] == "Spatial Viewer"
    assert package_manifest["n_sections"] == dataset.n_sections
    assert package_manifest["total_cells"] == dataset.n_cells
    assert package_manifest["viewer"]["gene_manifest_path"] == "viewer.genes.json"
    assert package_manifest["viewer"]["gene_shard_dir"] == "viewer.genes"
    assert "viewer.genes.json" in package_manifest["files"]
    assert shard_names[0] in package_manifest["files"]

    assert embedded["gene_aux_url"] == "viewer.genes.json"
    assert embedded["available_genes"] == ["G1", "G2", "G3"]
    assert embedded["embedded_genes"] == []
    assert embedded["genes_meta"] == {}
    assert all(not (section.get("genes") or section.get("genes_sparse")) for section in embedded["sections"])
    assert set(manifest["genes_meta"]) == {"G1", "G2", "G3"}
    assert manifest["format"] == "karospace-gene-sidecar-manifest-v3"
    assert manifest["gene_sidecar_format"] == "binary"
    assert manifest["gene_to_shard"]["G2"].startswith("viewer.genes/")
    assert manifest["gene_to_shard"]["G2"].endswith(".bin")


def test_package_sidecar_viewer_rejects_unsupported_gene_manifest(tmp_path):
    dataset = _build_dataset()
    html_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(html_path),
        main_cells_annotation="leiden",
        genes=["G1"],
        gene_storage="sidecar",
    )

    manifest_path = tmp_path / "viewer.genes.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["format"] = "karospace-package-v1"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(ValueError, match="unsupported gene sidecar manifest format"):
        package_sidecar_viewer(html_path, output_path=tmp_path / "viewer.karospace")


def test_karospace_package_requires_sidecar_mode(tmp_path):
    dataset = _build_dataset()

    with pytest.raises(ValueError, match="gene_storage='sidecar'"):
        export_to_html(
            dataset,
            output_path=str(tmp_path / "viewer.karospace"),
            main_cells_annotation="leiden",
            genes=["G1"],
            gene_storage="embedded",
        )


def test_karospace_package_loader_page_contains_expected_runtime_hooks():
    loader_path = Path(__file__).resolve().parents[1] / "karospace-package-loader.html"
    loader_html = loader_path.read_text(encoding="utf-8")

    assert "<title>KaroSpace</title>" in loader_html
    assert 'id="loading-label">KaroSpace</div>' in loader_html
    assert "drop a .karospace file to open" in loader_html
    assert "karospace-package.json" in loader_html
    assert "window.__karospacePackageSessions" in loader_html
    assert "DecompressionStream" in loader_html
    assert "karospace-gene-sidecar-manifest-v2" in loader_html
    assert "karospace-gene-sidecar-manifest-v3" in loader_html
    assert "karospace-gene-sidecar-manifest-v4" in loader_html
    assert "Packaged gene manifest format is not supported at" in loader_html
    assert "JSON.stringify(packagedGeneManifest.format)" in loader_html
    assert "function injectBootstrap(htmlText, sessionId)" in loader_html
    assert "window.fetch = function packageAwareFetch(input, init)" in loader_html
    assert "window.__karospacePackageSession" in loader_html
    assert "document.write(hydratedHtml);" in loader_html
    assert "const LAZY_ARCHIVE_THRESHOLD_BYTES" in loader_html
    assert "const EAGER_ARCHIVE_COMPAT_MAX_BYTES" in loader_html
    assert "class LazyZipArchive" in loader_html
    assert "async function parseLazyZipArchive(file)" in loader_html
    assert "async function openZipArchive(file, options = {})" in loader_html
    assert "file.size > LAZY_ARCHIVE_THRESHOLD_BYTES" in loader_html
    assert "function shouldRetryWithEagerArchive(file, archive, error)" in loader_html
    assert "async ensureArchiveCompatibility(error)" in loader_html
    assert "async withArchiveRetry(operation)" in loader_html
    assert "retrying in compatibility mode" in loader_html
    assert "readBlobSlice(" in loader_html
    assert "centralDirectoryOffset + centralDirectorySize" in loader_html
    assert "supportsRange(path)" in loader_html
    assert "readRange(path, start, endExclusive)" in loader_html
    assert "'<scr' + 'ipt>'" in loader_html
    assert "'</scr' + 'ipt>'" in loader_html
    assert 'type="file" accept=".karospace,.zip,application/zip,application/x-karospace-package"' in loader_html


def test_embedded_mode_stays_single_file(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=["G1"],
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
        main_cells_annotation="leiden",
        section_rotations={"S1": 44.5, "S2": -50.25},
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    assert {section["id"]: section["rotation_deg"] for section in embedded["sections"]} == {
        "S1": 44.5,
        "S2": 309.75,
    }


def test_export_embeds_pairwise_pseudobulk_de(monkeypatch, tmp_path):
    dataset = _build_pseudobulk_dataset()
    output_path = tmp_path / "viewer.html"
    pathway_gmt = tmp_path / "pathways.gmt"
    pathway_gmt.write_text(
        "Source-favored pathway\tna\tG1\tG3\n"
        "Reference pathway\tna\tG2\n",
        encoding="utf-8",
    )

    class FakeDDS:
        def contrast(self, column, baseline, group_to_compare):
            return np.array([1.0 if str(group_to_compare) == "A" else -1.0])

    def fake_fit_shared(counts, metadata, categories, **kwargs):
        return FakeDDS(), counts, metadata

    def fake_shared_contrast(dds, contrast, n_cpus=1):
        return pd.DataFrame(
            {
                "baseMean": [4.0, 2.0, 1.0],
                "log2FoldChange": [2.0, -1.0, 0.5],
                "stat": [3.0, -1.0, 1.5],
                "pvalue": [0.001, 0.5, 0.02],
                "padj": [0.01, 0.6, 0.05],
            },
            index=["G1", "G2", "G3"],
        )

    monkeypatch.setattr("karospace.pseudobulk._fit_deseq2_shared_categories", fake_fit_shared)
    monkeypatch.setattr("karospace.pseudobulk._deseq2_shared_contrast", fake_shared_contrast)

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        pathway_gmt=str(pathway_gmt),
        pathway_organism="Human",
        pathway_min_overlap=1,
        pathway_gsea_permutations=0,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    pseudobulk_de = embedded.get("pseudobulk_de") or {}
    assert "leiden" in pseudobulk_de
    assert {"A", "B"}.issubset(set(pseudobulk_de["leiden"]))
    assert "A" not in pseudobulk_de["leiden"]["A"]
    result = pseudobulk_de["leiden"]["A"]["B"]
    assert result["available"] is True
    assert result["method"] == "pseudobulk-deseq2"
    assert result["n_source"] == 40
    assert result["n_reference"] == 40
    assert result["n_replicates"] == 2
    assert len(result["genes"]) == 3
    assert result["log2foldchanges"] == result["logfoldchanges"]
    assert result["p_adjust_method"] == "fdr_bh"
    assert result["padj_cutoff"] == 0.05
    assert result["log2fc_cutoff"] == 0.5
    assert len(result["genes"]) == len(result["log2foldchanges"]) == len(result["pvals_adj"]) == len(result["scores"])
    assert len(result["genes"]) == len(result["pct_source"]) == len(result["pct_reference"])
    rest_result = pseudobulk_de["leiden"]["A"]["__rest__"]
    assert rest_result["available"] is True
    assert rest_result["n_source"] == 40
    assert rest_result["n_reference"] == 40
    assert rest_result["n_replicates"] == 2
    summary = pseudobulk_de["leiden"]["_summary"]
    assert summary["diagnostics"] == "pairwise"
    assert "A" in summary["pair_diagnostics"]
    assert "B" in summary["pair_diagnostics"]["A"]
    diagnostic = summary["pair_diagnostics"]["A"]["B"]
    assert len(diagnostic["pca"]) == 4
    assert len(diagnostic["distance_matrix"]) == 4
    pathways = result["pathway_enrichment"]
    assert pathways["ora"]["up"][0]["term"] == "Source-favored pathway"
    gsea_row = pathways["gsea"]["positive"][0]
    assert gsea_row["running_profile"]
    assert gsea_row["hit_indices"]
    assert gsea_row["rank_metric_profile"]
    assert embedded["pathway_settings"]["source"] == "gmt"
    assert embedded["pathway_settings"]["enriched_comparisons"] >= 1
    assert embedded["marker_genes"]["leiden"]["A"] == ["G1"]


def test_pseudobulk_min_pct_expressed_does_not_truncate_de_table(monkeypatch, tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    def fake_deseq2_pair(counts, metadata, source, reference):
        return pd.DataFrame(
            {
                "baseMean": [4.0, 2.0, 1.0],
                "log2FoldChange": [2.0, -1.0, 0.5],
                "stat": [3.0, -1.0, 1.5],
                "pvalue": [0.001, 0.5, 0.02],
                "padj": [0.01, 0.6, 0.05],
            },
            index=["G1", "G2", "G3"],
        )

    monkeypatch.setattr("karospace.pseudobulk._fit_deseq2_pair", fake_deseq2_pair)

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        pseudobulk_min_pct_expressed=0.5,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    result = embedded["pseudobulk_de"]["leiden"]["A"]["B"]

    assert result["genes"] == ["G1", "G3", "G2"]
    assert result["pct_source"] == [1.0, 0.0, 0.5]
    assert result["pct_reference"] == [1.0, 0.5, 0.0]
    assert result["min_pct_expressed"] == 0.5


def test_export_marks_pseudobulk_de_unavailable_for_small_clusters(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    result = embedded["pseudobulk_de"]["leiden"]["A"]["B"]
    assert result["available"] is False
    assert result["reason"] == "insufficient_cells"
    assert result["min_cells_required"] == 20


def test_export_rejects_unknown_section_rotation_ids(tmp_path):
    dataset = _build_dataset()

    with pytest.raises(ValueError, match="unknown section_id"):
        export_to_html(
            dataset,
            output_path=str(tmp_path / "viewer.html"),
            main_cells_annotation="leiden",
            section_rotations={"missing": 45},
        )


def test_cli_passes_section_rotations_to_export(monkeypatch, tmp_path):
    input_path = tmp_path / "input.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    captured = {}

    def fake_load(path, groupby, spatial_key="spatial"):
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

    def fake_load(path, groupby, spatial_key="spatial"):
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


def test_cli_omits_removed_cluster_de_options(monkeypatch, tmp_path):
    input_path = tmp_path / "input.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    captured = {}

    def fake_load(path, groupby, spatial_key="spatial"):
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
            "--pseudobulk",
            "None",
            "--pseudobulk-additional-annotations",
            "subclass,region",
            "--pseudobulk-simple-constrast-categories",
            '{"leiden":["A","B"],"subclass":["T cell"],"region":["cortex","white matter"]}',
            "--pseudobulk-counts-layer",
            "counts",
            "--pseudobulk-min-replicates",
            "3",
            "--pseudobulk-min-pct-expressed",
            "0.1",
            "--pseudobulk-p-adjust-method",
            "holm",
            "--pseudobulk-padj-cutoff",
            "0.01",
            "--pseudobulk-log2fc-cutoff",
            "1.2",
            "--pseudobulk-deseq2-fit-type",
            "mean",
            "--pseudobulk-embed-top-n-per-comparison",
            "7",
            "--pathway-gmt",
            "reactome.gmt,go_bp.gmt",
            "--pathway-organism",
            "Mouse",
            "--pathway-top-n",
            "9",
            "--pathway-min-overlap",
            "2",
            "--pathway-gsea-permutations",
            "0",
            "--interaction-markers",
            "None",
            "--tutorial",
        ],
    )

    cli_module.main()

    assert "cluster_de_groupby" not in captured["kwargs"]
    assert "cluster_de_top_n" not in captured["kwargs"]
    assert "cluster_de_method" not in captured["kwargs"]
    assert "cluster_de_layer" not in captured["kwargs"]
    assert "cluster_de_min_cells" not in captured["kwargs"]
    assert captured["kwargs"]["pseudobulk_additional_annotations"] == ["subclass", "region"]
    assert captured["kwargs"]["pseudobulk_simple_constrast_categories"] == {
        "leiden": ["A", "B"],
        "subclass": ["T cell"],
        "region": ["cortex", "white matter"],
    }
    assert captured["kwargs"]["pseudobulk"] is None
    assert captured["kwargs"]["interaction_markers"] is None
    assert "pseudobulk_groupby" not in captured["kwargs"]
    assert "pseudobulk_replicate" not in captured["kwargs"]
    assert "interaction_markers_groupby" not in captured["kwargs"]
    assert captured["kwargs"]["pseudobulk_counts_layer"] == "counts"
    assert captured["kwargs"]["pseudobulk_min_replicates"] == 3
    assert captured["kwargs"]["pseudobulk_min_pct_expressed"] == 0.1
    assert captured["kwargs"]["pseudobulk_p_adjust_method"] == "holm"
    assert captured["kwargs"]["pseudobulk_padj_cutoff"] == 0.01
    assert captured["kwargs"]["pseudobulk_log2fc_cutoff"] == 1.2
    assert captured["kwargs"]["pseudobulk_deseq2_fit_type"] == "mean"
    assert captured["kwargs"]["pseudobulk_embed_top_n_per_comparison"] == 7
    assert captured["kwargs"]["pathway_gmt"] == ["reactome.gmt", "go_bp.gmt"]
    assert captured["kwargs"]["pathway_organism"] == "Mouse"
    assert captured["kwargs"]["pathway_top_n"] == 9
    assert captured["kwargs"]["pathway_min_overlap"] == 2
    assert captured["kwargs"]["pathway_gsea_permutations"] == 0
    assert captured["kwargs"]["tutorial"] is True


def test_cli_rejects_flat_pseudobulk_categories_with_multiple_annotations(monkeypatch, tmp_path, capsys):
    input_path = tmp_path / "input.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    monkeypatch.setattr(
        cli_module.sys,
        "argv",
        [
            "karospace",
            str(input_path),
            "--pseudobulk-additional-annotations",
            "subclass,region",
            "--pseudobulk-simple-constrast-categories",
            "A,B",
        ],
    )

    with pytest.raises(SystemExit) as excinfo:
        cli_module.main()

    assert excinfo.value.code == 2
    captured = capsys.readouterr()
    assert "ambiguous with multiple pseudobulk annotations" in captured.err


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
        main_cells_annotation="leiden",
        genes=["G1", "G2", "G3"],
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
        main_cells_annotation="leiden",
        genes=["G1", "G2", "G3"],
        gene_correlation_top_n=0,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    assert embedded.get("gene_correlations") == {}


def test_cli_passes_gene_correlation_top_n_to_export(monkeypatch, tmp_path):
    input_path = tmp_path / "input.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    captured = {}

    def fake_load(path, groupby, spatial_key="spatial"):
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


def test_cli_sidecar_done_message_points_to_http_server(monkeypatch, tmp_path, capsys):
    input_path = tmp_path / "input.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    output_path = tmp_path / "viewer.html"

    def fake_load(path, groupby, spatial_key="spatial"):
        return _build_dataset()

    def fake_export(dataset, **kwargs):
        return str(output_path)

    monkeypatch.setattr("karospace.data_loader.load_spatial_data", fake_load)
    monkeypatch.setattr("karospace.exporter.export_to_html", fake_export)
    monkeypatch.setattr(
        cli_module.sys,
        "argv",
        [
            "karospace",
            str(input_path),
            "-o",
            str(output_path),
            "--gene-storage",
            "sidecar",
        ],
    )

    cli_module.main()

    output = capsys.readouterr().out
    assert "Sidecar gene loading requires HTTP(S)" in output
    assert f"python -m http.server --directory {output_path.parent}" in output
    assert f"http://localhost:8000/{output_path.name}" in output


def test_spatial_variable_genes_empty_without_spatial_graph(tmp_path):
    """spatial_variable_genes is [] when no obsp spatial graph exists."""
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        spatial_variable_genes_n=100,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    assert "spatial_variable_genes" in embedded
    assert embedded["spatial_variable_genes"] == []


def _build_dataset_with_spatial_graph() -> "SpatialDataset":
    """Build a dataset that has a spatial connectivity matrix in obsp."""
    import scipy.sparse as sp
    from karospace.data_loader import SectionData, SpatialDataset

    obs = pd.DataFrame(
        {
            "sample_id": pd.Categorical(["S1", "S1", "S2", "S2"]),
            "leiden": pd.Categorical(["A", "B", "A", "B"]),
        },
        index=[f"cell_{i}" for i in range(4)],
    )
    var = pd.DataFrame(index=["G1", "G2", "G3"])
    x = np.array(
        [[1.0, 0.0, 5.0], [2.0, 0.5, 0.0], [3.0, 1.0, 4.0], [4.0, 0.0, 1.0]],
        dtype=float,
    )
    adata = AnnData(X=x, obs=obs, var=var)
    adata.obsm["spatial"] = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]], dtype=float
    )
    # Add a simple spatial connectivity matrix (ring graph)
    W = sp.csr_matrix(
        np.array(
            [[0, 1, 1, 0], [1, 0, 0, 1], [1, 0, 0, 1], [0, 1, 1, 0]], dtype=float
        )
    )
    adata.obsp["spatial_connectivities"] = W
    sections = [
        SectionData(section_id="S1", coordinates=adata.obsm["spatial"][:2]),
        SectionData(section_id="S2", coordinates=adata.obsm["spatial"][2:]),
    ]
    return SpatialDataset(
        adata=adata,
        sections=sections,
        groupby="sample_id",
        obs_columns=["sample_id", "leiden"],
        var_names=["G1", "G2", "G3"],
        metadata_section=[],
    )


def test_spatial_variable_genes_computed_with_spatial_graph(tmp_path):
    """spatial_variable_genes is a sorted list of {gene, I} when obsp graph is present."""
    dataset = _build_dataset_with_spatial_graph()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        spatial_variable_genes_n=10,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    svg = embedded.get("spatial_variable_genes")
    assert svg is not None
    assert isinstance(svg, list)
    assert len(svg) > 0
    for entry in svg:
        assert "gene" in entry
        assert "I" in entry
        assert isinstance(entry["I"], float)
        assert -1.0 <= entry["I"] <= 1.0
    # Should be sorted descending
    scores = [e["I"] for e in svg]
    assert scores == sorted(scores, reverse=True)


def test_spatial_variable_genes_disabled_when_n_zero(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        spatial_variable_genes_n=0,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    assert embedded.get("spatial_variable_genes") == []


def test_cli_passes_spatial_variable_genes_n_to_export(monkeypatch, tmp_path):
    import karospace.cli as cli_module

    input_path = tmp_path / "input.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    captured = {}

    def fake_load(path, groupby, spatial_key="spatial"):
        return _build_dataset()

    def fake_export(dataset, **kwargs):
        captured["kwargs"] = kwargs
        return str(tmp_path / "viewer.html")

    monkeypatch.setattr("karospace.data_loader.load_spatial_data", fake_load)
    monkeypatch.setattr("karospace.exporter.export_to_html", fake_export)
    monkeypatch.setattr(
        cli_module.sys,
        "argv",
        ["karospace", str(input_path), "--spatial-variable-genes-n", "50"],
    )

    cli_module.main()

    assert captured["kwargs"]["spatial_variable_genes_n"] == 50


def test_cli_package_sidecar_dispatches_to_packager(monkeypatch, tmp_path):
    html_path = tmp_path / "viewer.html"
    html_path.write_text("<html><head><title>KaroSpace</title></head><body></body></html>", encoding="utf-8")
    package_path = tmp_path / "viewer.karospace"
    captured = {}

    def fake_package(html, **kwargs):
        captured["html"] = html
        captured["kwargs"] = kwargs
        return str(package_path)

    monkeypatch.setattr("karospace.exporter.package_sidecar_viewer", fake_package)
    monkeypatch.setattr(
        cli_module.sys,
        "argv",
        [
            "karospace",
            "package-sidecar",
            str(html_path),
            "--output",
            str(package_path),
        ],
    )

    cli_module.main()

    assert Path(captured["html"]) == html_path
    assert captured["kwargs"]["output_path"] == str(package_path)
    assert captured["kwargs"]["gene_manifest_path"] is None
    assert captured["kwargs"]["gene_shard_dir"] is None


def test_export_uses_companion_analytics_when_present(monkeypatch, tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"
    pseudobulk_de = {
        "leiden": {
            "A": {
                "B": {
                    "available": True,
                    "genes": ["G1"],
                    "logfoldchanges": [1.2],
                    "pvals_adj": [0.01],
                    "scores": [2.1],
                    "pct_source": [1.0],
                    "pct_reference": [0.5],
                    "n_source": 2,
                    "n_reference": 2,
                    "pseudobulk_samples": {
                        "pca": [[0.0, 0.0], [1.0, 1.0]],
                        "pca_variance": [0.9, 0.1],
                        "distance_matrix": [[0.0, 1.0], [1.0, 0.0]],
                    },
                }
            }
        }
    }
    neighbor_stats = {
        "leiden": {
            "categories": ["A", "B"],
            "counts": [[1.0, 2.0], [2.0, 1.0]],
            "n_cells": [2, 2],
            "mean_degree": [1.5, 1.5],
        }
    }
    interaction_markers = {
        "leiden": {
            "A": {
                "B": {
                    "available": True,
                    "genes": ["G1"],
                    "logfoldchanges": [0.9],
                    "pvals_adj": [0.02],
                    "n_contact": 1,
                    "n_non_contact": 1,
                    "pct_contact": 50.0,
                    "mean_target_neighbors_contact": 1.0,
                    "mean_target_neighbors_non_contact": 0.0,
                    "target_edge_count": 2.0,
                    "target_zscore": None,
                }
            }
        }
    }
    companion_gene_correlations = {"G1": [{"gene": "G2", "r": 0.75}]}
    companion_spatial_variable_genes = [{"gene": "G1", "I": 0.42}]
    companion_cluster_gene_means = {
        "genes": ["G1"],
        "columns": {
            "leiden": {
                "categories": ["A", "B"],
                "means": {"A": [1.0], "B": [2.0]},
                "background": [1.5],
            }
        },
    }
    recomputed_gene_correlations = {"G1": [{"gene": "G3", "r": 0.66}]}
    recomputed_spatial_variable_genes = [{"gene": "G2", "I": 0.31}]
    recomputed_cluster_gene_means = {
        "genes": ["G1"],
        "columns": {
            "leiden": {
                "categories": ["A", "B"],
                "means": {"A": [3.0], "B": [4.0]},
                "background": [3.5],
            }
        },
    }
    _attach_companion_analytics(
        dataset,
        pseudobulk_de=pseudobulk_de,
        neighbor_stats=neighbor_stats,
        interaction_markers=interaction_markers,
        gene_correlations=companion_gene_correlations,
        spatial_variable_genes=companion_spatial_variable_genes,
        cluster_gene_means=companion_cluster_gene_means,
    )

    morans_calls = []

    def fake_morans_i(adata, genes, n_genes):
        morans_calls.append((adata.n_obs, list(genes), n_genes))
        return recomputed_spatial_variable_genes

    monkeypatch.setattr(
        "karospace.exporter._compute_gene_correlations_from_category_means",
        lambda *args, **kwargs: recomputed_gene_correlations,
    )
    monkeypatch.setattr("karospace.exporter._compute_morans_i", fake_morans_i)
    monkeypatch.setattr(
        "karospace.exporter._cluster_gene_means_from_pseudobulk_de",
        lambda *args, **kwargs: recomputed_cluster_gene_means,
    )

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=["G1"],
        neighbor_stats_groupby=["leiden"],
        gene_correlation_top_n=5,
        spatial_variable_genes_n=25,
        cluster_means_n_genes=10,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    expected_pseudobulk_de = json.loads(json.dumps(pseudobulk_de))
    del expected_pseudobulk_de["leiden"]["A"]["B"]["pseudobulk_samples"]
    assert embedded["pseudobulk_de"] == expected_pseudobulk_de
    assert "pseudobulk_samples" in pseudobulk_de["leiden"]["A"]["B"]
    assert embedded["neighbor_stats"] == neighbor_stats
    assert embedded["interaction_markers"] == interaction_markers
    assert embedded["gene_correlations"] == recomputed_gene_correlations
    assert embedded["spatial_variable_genes"] == recomputed_spatial_variable_genes
    assert embedded["cluster_gene_means"] == recomputed_cluster_gene_means
    assert morans_calls == [(dataset.adata.n_obs, list(dataset.var_names), 25)]


def test_pseudobulk_embed_top_n_per_comparison_limits_auto_embedded_genes(tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"
    pseudobulk_de = {
        "leiden": {
            "A": {
                "B": {
                    "available": True,
                    "genes": ["G2", "G1", "G3"],
                    "logfoldchanges": [2.0, 1.5, 1.2],
                    "log2foldchanges": [2.0, 1.5, 1.2],
                    "pvals_adj": [0.003, 0.001, 0.002],
                    "scores": [3.0, 5.0, 4.0],
                    "pct_source": [1.0, 1.0, 1.0],
                    "pct_reference": [0.0, 0.0, 0.0],
                    "n_source": 2,
                    "n_reference": 2,
                    "padj_cutoff": 0.05,
                    "log2fc_cutoff": 0.5,
                }
            }
        }
    }
    _attach_companion_analytics(dataset, pseudobulk_de=pseudobulk_de)

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=["G3"],
        pseudobulk_embed_top_n_per_comparison=1,
        interaction_markers=None,
        neighbor_stats_groupby=[],
        spatial_variable_genes_n=0,
        cluster_means_n_genes=0,
        gene_correlation_top_n=0,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    assert embedded["pseudobulk_de"]["leiden"]["A"]["B"]["genes"] == ["G2", "G1", "G3"]
    assert set(embedded["available_genes"]) == {"G1", "G3"}
    assert set(embedded["embedded_genes"]) == {"G1", "G3"}
    assert set(embedded["genes_meta"]) == {"G1", "G3"}
    assert embedded["pseudobulk_settings"]["embed_top_n_per_comparison"] == 1


def test_export_falls_back_when_companion_analytics_absent(monkeypatch, tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"
    expected_gene_correlations = {"G1": [{"gene": "G2", "r": 0.33}]}

    cluster_gene_means = {
        "genes": ["G1", "G2"],
        "columns": {"leiden": {"means": {"A": [1.0, 2.0], "B": [2.0, 1.0]}}},
    }
    monkeypatch.setattr(
        "karospace.exporter._cluster_gene_means_from_pseudobulk_de",
        lambda *args, **kwargs: cluster_gene_means,
    )
    monkeypatch.setattr(
        "karospace.exporter._compute_gene_correlations_from_category_means",
        lambda *args, **kwargs: expected_gene_correlations,
    )

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=["G1"],
        gene_correlation_top_n=5,
        spatial_variable_genes_n=0,
        cluster_means_n_genes=10,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    assert embedded["gene_correlations"] == expected_gene_correlations


def test_export_computes_correlation_means_without_exporting_cluster_means(monkeypatch, tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"
    expected_gene_correlations = {"G1": [{"gene": "G2", "r": 0.44}], "G2": [{"gene": "G1", "r": 0.44}]}
    cluster_mean_calls = []
    cluster_gene_means = {
        "genes": ["G1", "G2"],
        "columns": {"leiden": {"means": {"A": [1.0, 2.0], "B": [2.0, 1.0]}}},
    }

    def fake_cluster_means(pseudobulk_de, genes, max_genes):
        cluster_mean_calls.append((list(genes), max_genes))
        return cluster_gene_means

    monkeypatch.setattr(
        "karospace.exporter._cluster_gene_means_from_pseudobulk_de",
        fake_cluster_means,
    )
    monkeypatch.setattr(
        "karospace.exporter._compute_gene_correlations_from_category_means",
        lambda *args, **kwargs: expected_gene_correlations,
    )

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=["G1", "G2"],
        pseudobulk=None,
        interaction_markers=None,
        neighbor_stats_groupby=[],
        gene_correlation_top_n=5,
        spatial_variable_genes_n=0,
        cluster_means_n_genes=0,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    assert embedded["cluster_gene_means"] is None
    assert embedded["gene_correlations"] == expected_gene_correlations
    assert cluster_mean_calls == [(["G1", "G2"], 2)]


def test_export_recomputes_only_missing_companion_analytics(monkeypatch, tmp_path):
    dataset = _build_dataset()
    output_path = tmp_path / "viewer.html"
    expected_gene_correlations = {"G1": [{"gene": "G2", "r": 0.61}]}
    cluster_gene_means = {
        "genes": ["G1", "G2"],
        "columns": {"leiden": {"means": {"A": [1.0, 2.0], "B": [2.0, 1.0]}}},
    }
    monkeypatch.setattr(
        "karospace.exporter._cluster_gene_means_from_pseudobulk_de",
        lambda *args, **kwargs: cluster_gene_means,
    )
    monkeypatch.setattr(
        "karospace.exporter._compute_gene_correlations_from_category_means",
        lambda *args, **kwargs: expected_gene_correlations,
    )

    export_to_html(
        dataset,
        output_path=str(output_path),
        main_cells_annotation="leiden",
        genes=["G1"],
        gene_correlation_top_n=5,
        spatial_variable_genes_n=0,
        cluster_means_n_genes=10,
    )

    embedded = _extract_data_json(output_path.read_text(encoding="utf-8"))
    assert embedded["marker_genes"] == {}
    assert embedded["gene_correlations"] == expected_gene_correlations


def test_load_spatial_data_handles_null_encoded_uns_entries(tmp_path):
    path = tmp_path / "broken.h5ad"
    _write_h5ad_with_nested_null_uns_entry(path)

    dataset = load_spatial_data(str(path), groupby="sample_id")

    assert dataset.n_sections == 1
    assert dataset.n_cells == 2
    assert dataset.var_names == ["G1", "G2"]
    assert dataset.adata.uns["dataset_source"]["chunksize"] == 256
    assert dataset.adata.uns["dataset_source"]["expression_path"] == "expr.h5"
    # Older anndata builds crash on null-encoded fields, so we strip them (see
    # _strip_null_encoded_h5ad_entries, covered by the sibling fallback test);
    # newer builds read them natively as None. Either way the entry must not
    # carry a garbage value and the surrounding fields must stay intact.
    assert dataset.adata.uns["dataset_source"].get("selected_fovs") is None


def test_read_h5ad_with_fallback_uses_sanitized_copy_and_cleans_it_up(tmp_path, monkeypatch):
    path = tmp_path / "broken.h5ad"
    _write_h5ad_with_nested_null_uns_entry(path)

    expected = AnnData(
        X=np.array([[1.0]], dtype=float),
        obs=pd.DataFrame(index=["cell_0"]),
        var=pd.DataFrame(index=["G1"]),
    )
    read_calls = []

    def fake_read_h5ad(target_path):
        read_calls.append(Path(target_path))
        if len(read_calls) == 1:
            raise RuntimeError(
                "No read method registered for IOSpec(encoding_type='null', encoding_version='0.1.0')"
            )
        assert Path(target_path) != path
        assert Path(target_path).exists()
        with h5py.File(target_path, "r") as handle:
            assert "selected_fovs" not in handle["uns"]["dataset_source"]
            expression_path = handle["uns"]["dataset_source"]["expression_path"][()]
            if isinstance(expression_path, bytes):
                expression_path = expression_path.decode("utf-8")
            assert expression_path == "expr.h5"
        return expected

    monkeypatch.setattr("karospace.data_loader.sc.read_h5ad", fake_read_h5ad)

    loaded = _read_h5ad_with_fallback(str(path))

    assert loaded is expected
    assert read_calls[0] == path
    assert len(read_calls) == 2
    assert not read_calls[1].exists()


def test_read_h5ad_with_fallback_clears_legacy_uns_metadata(tmp_path, monkeypatch):
    path = tmp_path / "legacy_uns.h5ad"
    obs = pd.DataFrame(index=["cell_0"])
    var = pd.DataFrame(index=["G1"])
    adata = AnnData(X=np.array([[1.0]], dtype=float), obs=obs, var=var)
    adata.uns["legacy_tool_state"] = {"payload": "old"}
    adata.write_h5ad(path)

    expected = AnnData(X=np.array([[1.0]], dtype=float), obs=obs, var=var)
    read_calls = []

    def fake_read_h5ad(target_path):
        read_calls.append(Path(target_path))
        if len(read_calls) == 1:
            raise RuntimeError(
                "No read method registered for IOSpec("
                "encoding_type='legacy-json', encoding_version='0.0.1')"
            )
        assert Path(target_path) != path
        assert Path(target_path).exists()
        with h5py.File(target_path, "r") as handle:
            assert list(handle["uns"].keys()) == []
            assert handle["uns"].attrs["encoding-type"] == "dict"
            assert handle["X"].shape == (1, 1)
            assert list(handle["obs"].keys())
            assert list(handle["var"].keys())
        return expected

    monkeypatch.setattr("karospace.data_loader.sc.read_h5ad", fake_read_h5ad)

    loaded = _read_h5ad_with_fallback(str(path))

    assert loaded is expected
    assert read_calls[0] == path
    assert len(read_calls) == 2
    assert not read_calls[1].exists()


def test_read_h5ad_with_fallback_raises_external_volume_guidance(monkeypatch):
    blocked_path = "/Volumes/processing2/example.companion.ready.h5ad"

    def fake_read_h5ad(_target_path):
        raise PermissionError(
            "[Errno 1] Unable to synchronously open file (error message = 'Operation not permitted')"
        )

    monkeypatch.setattr("karospace.data_loader.sc.read_h5ad", fake_read_h5ad)

    with pytest.raises(PermissionError, match="mounted volume"):
        _read_h5ad_with_fallback(blocked_path)


def test_integrate_polygon_annotations_produces_writeable_uns(tmp_path):
    adata = AnnData(
        X=np.array([[1.0], [2.0], [3.0]], dtype=float),
        obs=pd.DataFrame(
            {"sample_id": pd.Categorical(["S1", "S1", "S2"])},
            index=["cell_0", "cell_1", "cell_2"],
        ),
        var=pd.DataFrame(index=["G1"]),
    )
    adata.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float)

    payload = {
        "format": "karospace-polygon-annotations-v1",
        "created_at": "2026-04-08T19:51:08.859Z",
        "groupby": "sample_id",
        "polygons": [
            {
                "id": 1,
                "label": "Region A",
                "section_id": "S1",
                "cell_local_indices": [0, 1],
                "vertices": [{"x": 0.0, "y": 0.0}, {"x": 1.0, "y": 0.0}, {"x": 1.0, "y": 1.0}],
                "color": "#ff7f50",
            }
        ],
    }

    adata = integrate_polygon_annotations(adata, payload)

    assert adata.obs["karospace_polygon_labels"].tolist() == ["Region A", "Region A", None]
    assert adata.obs["karospace_polygon_count"].tolist() == [1, 1, 0]
    assert adata.uns["karospace_polygon_annotations"]["polygons_storage"] == "columnar-json-v1"
    poly_table = adata.uns["karospace_polygon_annotations"]["polygons"]
    assert poly_table["label"] == ["Region A"]
    assert poly_table["section_id"] == ["S1"]
    assert json.loads(poly_table["cell_global_indices"][0]) == [0, 1]
    assert json.loads(poly_table["vertices"][0])[0] == {"x": 0.0, "y": 0.0}

    output_path = tmp_path / "with_polygons.h5ad"
    adata.write_h5ad(output_path)
    assert output_path.exists()


def test_all_example_scripts_compile():
    examples_dir = Path(__file__).resolve().parents[1] / "examples"
    example_paths = sorted(examples_dir.glob("*.py"))

    assert example_paths, "expected at least one example script"

    failures = []
    for path in example_paths:
        try:
            py_compile.compile(str(path), doraise=True)
        except Exception as exc:  # pragma: no cover - exercised only on failure
            failures.append((path.name, str(exc)))

    assert failures == []
