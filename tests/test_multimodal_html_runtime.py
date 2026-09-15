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
    assert "magma(" not in html
    assert "magmaRgb(" not in html


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
    assert "function getExplorationDistributionFeatureNames(modality = getExplorationModality())" in html
    assert "return getFeatureDatalistValuesForModality(modality);" in html
    assert "function getStatisticsDistributionFeatureNames(annotationCol = explorationColorCol || currentAnnotation || '', modality = getExplorationModality())" in html
    assert "if (subtab === 'means') return getStatisticsDistributionFeatureNames" in html
    assert "await ensureFeatureAvailable(feature, {" in html


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


def test_feature_google_search_uses_modality_keyword(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "const modality = options.modality || getExplorationModality?.() || getVisualModality?.() || '';" in html
    assert "const modalityLabel = String(getModalityDisplayLabel?.(modality) || modality || '').trim();" in html
    assert "const query = options.query || [label, modalityLabel].filter(Boolean).join(' ');" in html
    assert "`${label} feature`" not in html


def test_marker_search_datalist_uses_embedded_features_only(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "function getEmbeddedFeatureDatalistValuesForModality" in html
    assert "if (subtab === 'distribution') return getExplorationDistributionFeatureNames(modality);" in html
    assert "if (subtab === 'means') return getStatisticsDistributionFeatureNames" in html
    assert "return getPseudobulkMeanFeatureNames(annotationCol, modality);" in html


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
    assert "DATA.wilcoxon_de_by_modality" in html
    assert "DATA.marker_features_by_method_by_modality" in html
    assert "DATA.category_feature_means_by_method_by_modality" in html
    assert ("DATA." + "g" + "ene_correlations") not in html
    assert ("DATA.spatial_variable_" + "g" + "enes") not in html


def test_statistics_method_selector_runtime_is_available(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "let statisticsMethodByPanel" in html
    assert "function getStatisticsMethodOptions" in html
    assert "function getActiveStatisticsMethod" in html
    assert "statistics-method-select" in html
    assert "const labelFor = (method) => method === 'pseudobulk' ? 'Pseudobulk' : 'Wilcoxon';" in html
    assert "activeMethod === 'pseudobulk' ? 'Pseudobulk' : 'Wilcoxon'" in html


def test_simple_design_marker_chips_can_spotlight_categories(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "function toggleAnnotationCategorySpotlight(annotationCol, category)" in html
    assert "data-pseudobulk-de-marker-annotation" in html
    assert "data-pseudobulk-de-marker-category" in html
    assert "toggleAnnotationCategorySpotlight(annotationCol, cat);" in html
    assert "toggleAnnotationCategorySpotlight(markerColorCol, cat);" in html


def test_simple_design_calc_info_matches_active_analysis_method(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "wilcoxon_marker_features" in html
    assert "wilcoxon_simple_de_section" in html
    assert "wilcoxon_simple_de_table" in html
    assert "KaroSpace uses the normalized layer when available" in html
    assert "Pseudobulk DESeq2 marker features" in html
    assert "pseudobulk sample = sum raw counts for replicate x category" in html
    assert "Pseudobulk DESeq2 feature table" in html
    assert "function getAnalysisCalcInfoKey(method, subject)" in html
    assert "getAnalysisCalcInfoKey(activeMethod, 'simple_section')" in html
    assert "getAnalysisCalcInfoKey(activeMethod, 'simple_table')" in html
    assert "getAnalysisCalcInfoKey(activeMethod, 'volcano_plot')" in html


def test_statistics_feature_calc_info_describes_wilcoxon_and_pseudobulk_means(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "Values are grouped by the selected Annotation dropdown" in html
    assert "display values are normalized Relative Counts (RC)" in html
    assert "The export can instead use LogNormalize or a pre-normalized layer" in html
    assert "Optional restriction filters limit which cells enter the summary" in html
    assert "Tiles use category means from the active Statistics method" in html
    assert "RC value = raw count * scale_factor / cell library size" in html
    assert "LogNormalize value = log1p(raw count * scale_factor / cell library size)" in html
    assert "mean = sum(values) / cells in group" in html
    assert "% Expr = 100 * cells with value > 0 / cells in group" in html
    assert "Wilcoxon/Pseudobulk category means" in html
    assert "Low-count cells and features are first filtered out" in html
    assert "Wilcoxon and Pseudobulk use normalized Relative Counts (RC) values" in html
    assert "replicate-category samples with too few cells are removed prior average calculation" in html
    assert "DESeq2 applies the minimum replicate filter afterward" in html
    assert "LogNormalize value = log1p(raw count * scale_factor / cell library size); Wilcoxon mean = mean per-cell value in category; pseudobulk mean = mean over retained replicate-category samples of mean value; delta = category mean - background" in html
    assert "LogNormalize value = log1p(raw count * 10000 / cell library size); Wilcoxon mean = mean per-cell display value in category; pseudobulk display mean = mean over retained replicate-category samples of mean display value; delta = category mean - background" not in html


def test_feature_distribution_panels_show_distribution_settings(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "function renderFeatureDistributionSettingsInfo(modality, method = null)" in html
    assert "function renderFeatureDistributionCountFilterWarning(method = null)" in html
    assert "feature-distribution-info" in html
    assert "feature-distribution-warning" in html
    assert "<strong>Assay:</strong>" not in html
    assert "<strong>Count filters:</strong> cells >=" in html
    assert "renderWarningDiv('feature-distribution-warning', content)" in html
    assert "<strong>Matrix:</strong> counts layer" in html
    assert "<strong>Normalization:</strong> RC library-size normalization without log1p" in html
    assert "<strong>Scale factor:</strong>" in html
    assert "<strong>Log transform:</strong> no" in html
    assert "renderFeatureDistributionSettingsInfo(modality, methodLabel)" in html
    assert "${countFilterWarningHtml}\n            ${settingsInfoHtml}" in html


def test_exploration_distribution_cells_column_is_last(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert '<th data-feature-dist-sort="cat">Group${arrow(\'cat\')}</th>\n                        <th data-feature-dist-sort="mean">Mean${arrow(\'mean\')}</th>\n                        <th data-feature-dist-sort="median">Median${arrow(\'median\')}</th>\n                        <th data-feature-dist-sort="pctExpr">% Expr${arrow(\'pctExpr\')}</th>\n                        <th data-feature-dist-sort="n">Cells${arrow(\'n\')}</th>' in html
    assert '<th data-feature-dist-sort="n">n${arrow(\'n\')}</th>' not in html
    assert '<td>${fmtP(s.pctExpr)}</td>\n                <td>${fmtN(s.n)}</td>' in html


def test_warning_blocks_use_icon_without_warning_prefix(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "function renderWarningDiv(className, contentHtml)" in html
    assert 'class="exploration-embedded-warning-icon"' in html
    assert "renderWarningDiv('comparison-info-warning'" in html
    assert "renderWarningDiv('pseudobulk-comparison-warning'" in html
    assert "renderWarningDiv('neighbor-warning'" in html
    assert "renderWarningDiv('features-warning'" in html
    assert "Warning.</strong>" not in html
    assert "warning.</strong>" not in html


def test_insights_has_separate_statistics_menu(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert 'id="insights-mode-statistics" data-insights-mode="statistics"' in html
    assert 'id="exploration-embedded-warning"' in html
    assert "embedded cells and features inside this HTML file" in html
    assert "classList.toggle('hidden', mode !== 'exploration')" in html
    assert "exploration: {" in html
    assert "features: ['distribution']" in html
    assert "compare: ['groups', 'regions', 'selection', 'river']" in html
    assert "statistics: {" in html
    assert "features: ['means', 'de-features', 'spatial']" in html
    assert "compare: ['cell-de', 'complex-contrast']" in html
    assert "neighbors: ['enrichment', 'interactions', 'dispersion']" in html
    assert "getInsightsModeForLeaf(topLevel, subtab)" in html
    assert ">Per cell<" not in html
    assert ">Per sample<" not in html


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


def test_selection_comparison_uses_selection_labels_and_outside_bar_labels(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "Cell Composition — SELECTION A vs SELECTION B" in html
    assert "const compositionTitle = compareAllCells" in html
    assert "Feature values - selection a vs selection b" in html
    assert "const labelA = compareAllCells ? 'Selected cells' : 'Selection A';" in html
    assert "const labelB = compareAllCells ? 'All cells' : 'Selection B';" in html
    assert "Selection comparison" in html
    assert "Selection A and B category percentages" in html
    assert "const barLabelA = `${formatCompactNumber(meanA)} (${pctA.toFixed(0)}%)`;" in html
    assert "const barLabelB = `${formatCompactNumber(meanB)} (${pctB.toFixed(0)}%)`;" in html
    assert "style=\"width:${clampPercent(100 * meanA / vmax)}%;color:#fff;\"" in html
    assert "style=\"width:${clampPercent(100 * meanB / vmax)}%;color:#fff;\"" in html
    assert "barLabel(barLabelA, outsideA)" in html
    assert "barLabel(barLabelB, outsideB)" in html


def test_annotation_comparison_moves_tiny_bar_labels_outside(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert ".selection-summary-expr-bars.has-outside-label" in html
    assert ".selection-summary-expr-bar-label-outside" in html
    assert "const outsideA = Number.isFinite(ratio) && ratio < 0.2;" in html
    assert "const outsideB = ratio > 5;" in html
    assert "barLabel(labelA, outsideA)" in html
    assert "barLabel(labelB, outsideB)" in html


def test_region_comparison_moves_tiny_bar_labels_outside(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "Feature Values by Region" in html
    assert "const barsClass = `selection-summary-expr-bars${outsideA || outsideB ? ' has-outside-label' : ''}`;" in html
    assert 'Region A mean: ${formatCompactNumber(meanA)}">${barLabel(labelA, outsideA)}</div>' in html
    assert 'Region B mean: ${formatCompactNumber(meanB)}">${barLabel(labelB, outsideB)}</div>' in html


def test_compare_pseudobulk_follows_exploration_controls(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert 'id="pseudobulk-de-modality-select"' not in html
    assert 'id="pseudobulk-de-annotation"' not in html
    assert "return getExplorationModality();" in html
    assert "function setPseudobulkPanelModality" not in html
    assert "getPseudobulkDEPayloadForModality" in html


def test_interaction_markers_use_modality_payload(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert 'id="interaction-marker-modality-select"' in html
    assert "DATA.interaction_markers_by_modality" in html


def test_html_copy_uses_feature_labels(tmp_path=None):
    html = _render_multimodal_html(tmp_path)

    assert "Feature discovery" in html
    assert "Marker features" in html
    assert "Spatial features" in html
    assert "activeMethod === 'pseudobulk' ? 'Pseudobulk' : 'Wilcoxon'" in html
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
    assert "['method', 'modality', 'annotation_column', 'category', 'reference', 'rank', 'feature'" in html
