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
        spatial_variable_features_n=0,
        feature_correlation_top_n=0,
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
    assert ("AVAILABLE_" + "G" + "ENE_SET") not in html
    assert "DATA.available_features" not in html
    assert ("resolveCanonical" + "G" + "eneName") not in html


def test_generated_html_drops_removed_old_paths(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "renderOldGroupDE" not in html
    assert "computeCellSetDEAsync" not in html
    assert "function computeCellSetDE(" not in html
    assert "section.he_image)" not in html
    assert "section.he_image " not in html
    assert "logfoldchanges" not in html
    assert "karospace-feature-sidecar-manifest-v3" not in html
    assert "targetModality === DEFAULT_MODALITY_NAME ? manifest : null" not in html
    assert "manifest.section_order" not in html


def test_generated_html_has_panel_scoped_modality_state(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "const PANEL_MODALITY_STATE" in html
    assert "visual:" in html
    assert "exploration:" in html
    assert "module:" in html
    assert "interactions:" in html
    assert "let CURRENT_MODALITY" not in html


def test_visual_controls_have_feature_namespace_select(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert 'id="visual-feature-namespace-select"' in html
    assert 'id="feature-input"' in html
    assert 'id="feature-list"' in html
    assert 'id="feature-discovery-panel"' in html
    assert ('id="' + 'g' + 'ene-input"') not in html
    assert 'id="modality-select"' not in html


def test_feature_discovery_uses_visual_namespace(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "const discoveryModality = getVisualModality();" in html
    assert "getFeatureSuggestionGroups(discoveryModality)" in html
    assert "getFeatureTokensForModality(recentFeatures, discoveryModality)" in html
    assert "getFeatureTokensForModality(panel?.features || [], discoveryModality)" in html


def test_exploration_feature_controls_accept_manual_input(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert '<input type="text" id="feature-module-feature-picker"' in html
    assert 'list="feature-module-feature-list"' in html
    assert '<select id="feature-module-feature-picker"' not in html
    assert 'id="marker-feature-search" type="text" list="marker-feature-search-list"' in html
    assert '<select class="marker-search" id="marker-feature-search"' not in html
    assert "moduleFeaturePicker?.addEventListener('keydown'" in html
    assert "markerSearch?.addEventListener('keydown'" in html
    assert "resolveCanonicalFeatureName(moduleFeaturePicker.value, moduleFocusedModality)" in html


def test_feature_modules_are_scoped_to_focused_modality(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert 'id="feature-module-modality-select"' in html
    assert "function getModuleBuilderModality()" in html
    assert "function setModuleBuilderModality(modality)" in html
    assert "function getFeatureModulesForModality(modality = getModuleBuilderModality())" in html
    assert "const activeModules = getFeatureModulesForModality(moduleFocusedModality);" in html
    assert "getFeatureDatalistValuesForModality(moduleFocusedModality)" in html
    assert "createFeatureModule('', features, moduleFocusedModality)" in html
    assert "format: 'karospace-feature-modules-v2'" in html
    assert "modality: getFeatureModuleModality(module)" in html
    assert "features: module.features.slice()" in html
    assert "getSectionFeatureValues(section, moduleFeature, sourceModality)" in html
    assert "getFeatureScaleRange(moduleFeature, sourceModality)" in html


def test_feature_module_changes_refresh_visual_namespace_controls(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "function syncVisualFeatureNamespaceSelect()" in html
    assert "if (typeof syncVisualFeatureNamespaceSelect === 'function') syncVisualFeatureNamespaceSelect();" in html
    assert "const syncFeatureNamespaceSelect = () =>" in html
    assert "syncVisualFeatureNamespaceSelect();" in html
    assert "setSelectOptions(select, options, selected);" in html


def test_feature_dropdowns_use_full_catalog_only_with_sidecar(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "const base = DATA.feature_manifest_url" in html
    assert "? getFeatureCatalog(modality)" in html
    assert ": getLoadedFeaturesForModality(modality);" in html
    assert '"rna_b"' in html


def test_marker_search_datalist_uses_embedded_features_only(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "function getEmbeddedFeatureDatalistValuesForModality" in html
    assert "if (subtab === 'distribution') return embedded;" in html
    assert "subtab === 'means'" in html
    assert "getPseudobulkMeanFeatureNames(explorationColorCol || currentAnnotation || '', modality)" in html


def test_marker_search_does_not_mutate_visual_feature_controls(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "Marker search filters Insights panels only." in html
    assert "if (isViewerFeatureLoadable(feature, modality))" not in html


def test_split_controls_have_independent_feature_namespaces(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert 'id="overview-blend-a-namespace"' in html
    assert 'id="overview-blend-b-namespace"' in html
    assert "overviewBlendSpec.a.modality" in html
    assert "overviewBlendSpec.b.modality" in html


def test_exploration_uses_features_and_modality_payloads(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert 'data-insights-tree-parent="features"' in html
    assert 'id="exploration-feature-modality-select"' in html
    assert "Focused modality" in html
    assert html.index('id="exploration-annotation-select"') < html.index('id="exploration-feature-modality-select"')
    assert html.index('id="exploration-feature-modality-select"') < html.index('id="visualization-menu-label"')
    assert "DATA.pseudobulk_de_by_modality" in html
    assert "DATA.category_feature_means_by_modality" in html
    assert ("DATA." + "g" + "ene_correlations") not in html
    assert ("DATA.spatial_variable_" + "g" + "enes") not in html


def test_compare_group_de_uses_exploration_modality(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "function shouldRunFullSidecarDE(modality = getExplorationModality())" in html
    assert "const activeModality = getExplorationModality();" in html
    assert "getSectionFeatureValues(section, feature, modality)" in html
    assert "getFeatureSidecarManifestEntryForModality(manifest, targetModality)" in html
    assert "return `${getExplorationModality()}::${groupA.key}::${groupB.key}`;" in html
    assert "return `${getExplorationModality()}::${groupAKey}::${groupBKey}`;" in html
    assert "async function runAnnotationGroupDE(sourceSpec, valueA, valueB, restrictSpec, restrictValue)" in html
    assert "runAnnotationGroupDE(sourceSpec, groupDeSourceValue, groupDeReferenceValue, restrictSpec, groupDeRestrictValue)" in html
    assert "if (shouldRunFullSidecarDE(getExplorationModality()))" in html
    assert "runFullRegionAnnotationDE(source, reference);" in html
    assert "const resultModality = result?.modality || getExplorationModality();" in html
    assert "allowUnknown: true" in html
    assert "modality: resultModality" in html


def test_selection_comparison_reruns_when_focused_modality_changes(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "return `${getExplorationModality()}:${groupB === null ? 'all' : 'region-b'}:${selectionRevision}:${selectionWelchRevision}`;" in html
    assert "const shouldRerunSelectionComparison = selectionWelchRunRequested && selectedCells.size > 0;" in html
    assert "resetSelectionWelchState({ keepRequested: shouldRerunSelectionComparison });" in html
    assert "updateSelectionInfo?.();" in html
    assert "await runSelectionWelchComparison();" in html


def test_selection_comparison_uses_full_sidecar_features(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "async function runSelectionWelchComparison()" in html
    assert "const fullResult = await runFullCellSetDE(cellSetA, cellSetB, {" in html
    assert "if (shouldRunFullSidecarDE(targetModality))" in html
    assert "selectionWelchCache.set(key, normalized);" in html
    assert "const selectionWelchResult = getCachedSelectionWelchResult(selectedCells, compareAllCells ? null : selectedCellsB);" in html
    assert "const cachedResult = getCachedSelectionWelchResult(selectedCells, compareAllCells ? null : selectedCellsB);" in html
    assert "Scanning all ${getModalityDisplayLabel(resultModality)} features from the sidecar." in html
    assert "Full sidecar comparison across ${Number(selectionWelchResult.totalFeatureCount || 0).toLocaleString()} ${getModalityDisplayLabel(resultModality)} features." in html


def test_compare_pseudobulk_follows_exploration_controls(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert 'id="pseudobulk-de-modality-select"' not in html
    assert 'id="pseudobulk-de-annotation"' not in html
    assert "return getExplorationModality();" in html
    assert "function setPseudobulkPanelModality" not in html
    assert "getPseudobulkDEPayloadForModality" in html
    assert "(DATA.pseudobulk_de || {})" not in html


def test_interaction_markers_use_modality_payload(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert 'id="interaction-marker-modality-select"' in html
    assert "DATA.interaction_markers_by_modality" in html
    assert "(DATA.interaction_markers || {})" not in html


def test_html_copy_uses_feature_labels(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "Feature discovery" in html
    assert "Marker features" in html
    assert "Spatial features" in html
    assert "Pseudobulk feature differential analysis" in html
    assert "Features in selection" in html
    assert "Feature values - annotation A vs annotation B" in html
    assert "No features matched" in html
    assert ("No " + "g" + "enes matched") not in html
    assert ("G" + "ene symbol") not in html
    assert ("G" + "enes in selection") not in html
    assert ("G" + "ene expression") not in html
    assert ("Sidecar " + "g" + "ene loading") not in html
    assert ("No " + "g" + "enes are currently loaded") not in html


def test_download_filenames_include_modality(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "'karospace-pseudobulk-de'," in html
    assert "sanitizeFilenamePart(getPseudobulkPanelModality()" in html
    assert "karospace-pseudobulk-de-features-${modName}-" in html
    assert "['modality', 'feature', 'base_mean'" in html
    assert "['modality', 'annotation_column', 'category', 'reference', 'rank', 'feature'" in html
