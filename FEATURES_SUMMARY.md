# KaroSpace Feature Summary

This document summarizes what the generated KaroSpace HTML viewer currently displays and how the user is expected to move through it.

## 1. Product Surfaces

- Python API for inspecting input metadata, loading `.h5ad`, `AnnData`, or SpatialData input, and exporting an HTML viewer.
- Command-line interface (`karospace`) for scriptable exports and inspect-only metadata checks.
- Desktop GUI (`karospace.gui`) for non-code export configuration.
- Standalone HTML output that can be opened in a browser without a Python server when all required gene data is embedded.
- Optional sidecar and `.karospace` package outputs for large gene payloads.

## 2. Input and Export Payload

- Input can be an `.h5ad` path, an in-memory `AnnData`, a SpatialData `.zarr` store, or an in-memory SpatialData object.
- SpatialData input is normalized to one AnnData table via `spatialdata_table` / `--spatialdata-table`.
- If full SpatialData construction fails because unrelated image or label elements are invalid, KaroSpace can fall back to reading the selected AnnData table directly from `tables/<table>`.
- `inspect_input_file(...)` / `--inspect-input` reports available `adata.obs` metadata, value types, examples, and missing-value counts without building sections, downsampling, exporting HTML, or running analytics.
- Section grouping uses `section_key`; SpatialData region metadata can be used automatically when `sample_id` is absent.
- Section metadata is split into `section_metadata` for visual filter chips and `section_metadata_extra` for stored section metadata that is not displayed as filter chips.
- Cell annotation dropdowns use `main_cell_annotation` plus `cell_annotations`.
- Optional payloads include UMAP coordinates, neighborhood graphs, image overlays, deconvolution keys, modules, sidecar genes, pseudobulk DE, pathways, spatial genes, category means, and gene correlations.

## 3. Introduction

- The central grid is the primary spatial browsing surface.
- Each visible section appears as a card with metadata, rotation controls, scalebar, outline color, and a canvas-rendered cell layer.
- The section filter bar displays exported section-level metadata chips and can reduce the grid to specific stages, samples, conditions, regions, or batches.
- The status area reports visible sections, exported cells, and exported genes.
- Downsampling status is displayed when the export contains fewer cells than the original input.

## 4. Helpers

- The Info button opens dataset notes, viewer context, and keyboard shortcut reminders.
- The message-circle-question-mark button opens the button guide.
- The button guide explains icon-only controls using the actual icons shown in the viewer.
- Keyboard shortcuts include `I` for toggling the Insights panel.

## 5. Session Tools

- The theme button switches light and dark modes.
- The screenshot menu exports the current grid view only, with size and transparent-background options.
- Session export downloads a JSON state file containing annotations, hidden categories, and current views.
- Session import restores a previously exported session JSON.
- The annotation export menu contains broad data, category, and annotation export actions.

## 6. Visual Setup

- The visual parameters button opens cell display controls such as size and opacity.
- Default mode displays one layer across the grid.
- Split mode compares two aligned visual layers in the same spatial panels.
- Default source can be a cell annotation layer or a gene/expression layer.
- The annotation selector chooses which cell-level annotation colors the grid.
- Split layer A and split layer B can represent annotations, genes, or exported modalities.
- The split boundary slider controls where layer A ends and layer B begins within each spatial panel.

## 7. Gene Expression

- Gene source changes the grid from annotation colors to gene or feature expression.
- The gene discovery panel supports fuzzy search, keyboard navigation, recent genes, marker suggestions, spatial genes, correlated-gene suggestions, and module-related suggestions when those payloads exist.
- The gene input accepts embedded genes and sidecar-loadable genes.
- Gene expression scale controls adjust the visible color range and can be propagated across split gene comparison.
- Sidecar mode keeps all gene expression vectors outside the HTML and fetches them from nearby sidecar files.
- Multiple modalities can be selected when exported, for example RNA genes and protein features.
- Genes not embedded in the HTML can still appear in DE tables or marker lists; sidecar-loadable genes can be fetched on demand.

## 8. Spatial Selection

- Pan mode is the default movement mode for navigating section panels.
- Lasso mode selects cells directly in spatial panels.
- The compare-selection workflow creates Region A and Region B cell sets for side-by-side comparison.
- Region B can be cleared with the cross-format compare button.
- Selected cells can be saved as a region annotation.
- Selected cells can be cleared from the lasso/cross control or from the selection chip.
- The query tool can select cells from annotation values, genes, or section metadata.

## 9. Modal View

- Clicking a section opens a high-detail modal view.
- Modal entry hides stretched intermediate cells and renders the section after the modal frame is open.
- Modal zoom controls change magnification.
- Modal rotation controls adjust visual section orientation.
- Neighbor controls appear when a neighbor graph is available and can display graph edges or hover-neighbor rings.
- H&E/image overlay controls can load and adjust image overlays with opacity, scale, rotation, flip, alignment, and related parameters.
- Modal display state shares section rotation and selection behavior with the overview grid.

## 10. Legend

- The Legend button shows or hides category rows.
- Legend rows show category colors and labels for the active annotation.
- Clicking a category color dot hides or shows that specific category.
- The global visibility button hides or shows all categories.
- Spotlight mode focuses one category while muting the rest.
- Legend export downloads palettes or category labels.
- Legend import restores previously exported palette or category-label files.

## 11. UMAP

- If `adata.obsm["X_umap"]` was exported, the UMAP button opens a linked embedding panel.
- The UMAP panel can be docked to corners, resized, and kept pinned while scrolling the grid.
- UMAP supports pan, zoom, point-size control, and lasso selection.
- UMAP lasso selections synchronize with the spatial grid.
- Selected UMAP cells can be saved as a cell-set region, not as a polygon region tied to one section.

## 12. Insights Modes

- Insights is the right-side workspace for selected cells, regions, gene modules, and built-in analysis panels.
- Selection mode shows compact summaries for active lasso, UMAP, or query selections.
- Region mode stores, selects, groups, imports, exports, recolors, and deletes user-created spatial regions or cell sets.
- Module mode creates gene modules and displays module scores like expression layers.
- Exploration mode contains the Visualization menu tree and the built-in analysis panels.
- The Exploration annotation selector is independent from the main grid annotation selector.

## 13. Selection Workflow

- The Selection summary reports active selected cells by section and main annotation.
- Find More opens the full `Compare > Per cell > Selections` workflow.
- Selection marker genes are calculated from selected cells using a two-sided Welch test.
- The selection comparison panel displays composition and expression summaries.
- Expression summaries show mean expression and percent expressed for displayed genes.

## 14. Region Workflow

- Region creation starts from selected cells.
- Regions can be created from modal or grid lasso selections.
- Each region row represents a saved spatial region or cell group.
- Clicking a region row selects its cells.
- Clicking a group row selects the union of cells from nested regions.
- Regions can be recolored, grouped, exported, imported, or deleted.
- Imported region JSON restores saved regions and associated selected cells in the viewer.
- Region comparison in Exploration uses saved regions.

## 15. Module Workflow

- The module panel has a gene picker and a draft gene list.
- Creating a module scales each gene and computes an average module score.
- Module scores can be loaded into the spatial viewer like gene expression.
- Module definitions can be imported or exported as JSON.

## 16. Exploration Menu

- The Visualization menu tree opens Overview, Genes, Compare, and Neighbors panels.
- Overview summarizes section composition and metadata trends.
- Genes focuses marker genes, spatial genes, per-cell distributions, and per-sample/category mean expression.
- Compare contains selection, region, annotation, pseudobulk, pathway, and relationship comparisons.
- Neighbors contains spatial adjacency enrichment, interaction markers, and dispersion analysis.

## 17. Exploration > Overview

- `Overview > Summary` aggregates cells across section metadata.
- Summary calculations use cells embedded in the HTML file.
- The per-annotation trend selector focuses one category across section metadata for abundance checks.
- `Overview > Sections` compares section-level composition.
- Section composition can be shown as stacked bars or a heatmap.

## 18. Exploration > Genes

- `Genes > Markers` lists pseudobulk-derived marker genes by category when available.
- Marker genes can be displayed as compact lists or heatmaps.
- `Genes > Spatial` shows Moran Index rankings computed at export on the full input cell set.
- Spatial genes can be displayed as a ranked list or graph.
- `Genes > Distribution > Per cell` summarizes expression distributions across categories for a selected gene.
- Per-cell distribution calculations use cells embedded in the HTML.
- Per-cell distributions can be shown as a table or violin/boxplot.
- `Genes > Distribution > Per sample` uses pseudobulk/category mean summaries for selected genes.
- Per-sample/category means can be shown as a table or barplot.

## 19. Exploration > Compare

- `Compare > Per cell > Selections` analyzes active lasso/query selections and is the detailed view behind Selection Find More.
- `Compare > Per cell > Regions` compares saved region annotations and can run region comparisons similar to selection comparisons.
- `Compare > Per cell > Annotations` compares categories within the selected annotation using per-cell summaries.
- `Compare > Per sample > Simple design` displays category-versus-category pseudobulk DE when exported.
- Simple design includes explanatory calculation popovers, double-dipping warnings, expression-percent filtering warnings, raw tables, markers, MA plots, volcano plots, PCA, distance matrix diagnostics, and pathway enrichment.
- Simple design has a full-width switch with Raw table, Genes, and Samples views.
- DE genes are filtered with `padj < cutoff` and `abs(log2FC) >= cutoff`.
- Genes below the minimum percent-expressed threshold in both compared groups are removed before `DeseqStats`, so they do not enter contrast-level multiple-testing correction.
- The pseudobulk marker lists are ordered by adjusted p-value then log2FC and can expand from the first displayed rows.
- MA and volcano plots show non-significant genes, minimum-percent-expression filtered genes, and significant genes with updated legends.
- Pathway Enrichment displays ORA pathways and GSEA enrichment in a separate panel.
- ORA uses significant DE genes favoring the selected annotation.
- ORA is shown as a dot plot with GeneRatio on x, pathway labels on y, dot size as gene count, fill color as `-log10(adjusted p-value)`, and border color indicating adjusted p-value threshold.
- ORA tables chip adjusted p-values and `-log10` adjusted p-values by significance.
- GSEA uses the retained ranked gene list after model and expression-percent filtering.
- GSEA profiles can be selected from a pathway dropdown and are displayed as enrichment-profile plots with hits and ranked-metric context.
- ORA and GSEA sections can switch between plot and raw table and include download buttons.
- `Compare > Relationships` shows how categories from two annotations map to each other, with controls to select the second annotation, swap direction, and export the correspondence table.

## 20. Exploration > Neighbors

- `Neighbors > Enrichment` summarizes which annotation categories are spatially adjacent more or less often than expected.
- Enrichment can be displayed as a table, network, or chord view.
- Optional neighbor permutations add enrichment z-scores; with zero permutations, observed counts, shares, cell counts, and mean degree still remain available.
- Yellow warnings appear when no neighbor graph exists or when the selected Exploration annotation has no neighbor stats; warnings list annotations with available stats.
- `Neighbors > Interactions` compares source cells based on which target categories they touch.
- Interaction markers are contact-conditioned pseudobulk marker results when exported.
- Interaction controls choose the source category and filter target names.
- `Neighbors > Dispersion` summarizes whether categories are clustered, dispersed, or close to random across all cells before HTML downsampling.
- Dispersion complements immediate neighbor enrichment by describing whole-section spatial arrangement.

## 21. Exported Analytics

- Pseudobulk DE uses raw counts grouped by replicate and annotation.
- Pseudobulk samples require at least `pseudobulk_min_cells_per_sample` cells before entering the shared DESeq2 fit.
- Category-versus-category contrasts are extracted from a shared fit per annotation column.
- Balanced-rest contrasts compare one category against the equally weighted mean of retained other categories.
- Pairwise PCA and distance diagnostics are generated for selected pseudobulk comparisons.
- ORA and GSEA are computed after Simple design pseudobulk DE and feed the Pathway Enrichment panel.
- Spatially variable genes are computed with Moran's I for up to `spatial_variable_genes_n` variable genes on the full input cell set.
- Category gene means are derived from embedded pseudobulk DE genes and feed per-sample/category distribution panels.
- Gene correlations are computed from category means and feed related-gene suggestions.
- Full-cell spatial dispersion is computed before HTML downsampling for the main cell annotation and requested `cell_annotations`.

## 22. Sharing and Storage

- Embedded HTML stores viewer data directly in the HTML file.
- Sidecar mode writes all gene expression vectors to a manifest and binary shard directory for lazy loading over HTTP(S).
- `.karospace` packages wrap sidecar exports into one shareable ZIP-based package with a local or hosted loader.
- Generated HTML is static; re-export is required to pick up new viewer code or new analytics.
- Old exported HTML files do not gain new behavior automatically.

## 23. API, CLI, and GUI Configuration

- API entry points include:
  - `inspect_input_file(...)`
  - `load_spatial_data(...)`
  - `export_to_html(...)`
  - `integrate_polygon_annotations(...)`
- CLI supports:
  - H5AD and SpatialData `.zarr` inputs
  - `--inspect-input`
  - `--main-cell-annotation`
  - `--cell-annotations`
  - `--section-metadata`
  - `--section-metadata-extra`
  - `--spatialdata-table`
  - gene encoding and sidecar storage controls
  - pseudobulk DE controls
  - pathway enrichment controls
  - neighbor stats and interaction-marker controls
  - section rotations and image overlays
- GUI supports inspect/validate/export workflows, searchable annotation/gene editors, presets, and advanced settings for feature storage, pseudobulk DE, neighbor stats, and interaction markers.

## 24. Practical Notes

- Large datasets may produce large HTML files; downsampling, sparse encoding, sidecar storage, and `.karospace` packages are the main size/performance levers.
- Sidecar HTML should be served over HTTP(S); direct `file://` loading can block sidecar fetches.
- Many Exploration panels depend on optional exported payloads. When data is missing, the HTML shows warnings or hides unavailable controls instead of crashing.
- Most plot and table download buttons export the currently visible analysis panel, while broad session and annotation exports live in the top toolbar.
