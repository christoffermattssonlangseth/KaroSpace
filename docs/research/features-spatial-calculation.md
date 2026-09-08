# Features Spatial Calculation

## Short Answer

The concept named approximately "Features Spatial" is the HTML viewer's `Features > Spatial` panel. It is not calculated in the browser; it is precomputed during export as `spatial_variable_features_by_modality` and rendered later as Moran Index rankings. The generated UI has a `Features` branch with a `Spatial` leaf, and the tutorial calls this the "Features Spatial panel" that shows Moran Index rankings computed at export. Citations: `karospace/exporter.py:25615-25627`, `karospace/exporter.py:9664-9672`, `karospace/data_loader.py:3342-3349`, `karospace/data_loader.py:3411-3417`.

## Calculation Site

The actual calculation is in `karospace/data_loader.py` inside `_compute_morans_i_for_features`. The function looks for the first available graph in `adata.obsp` under `spatial_connectivities`, `connectivities`, `neighbors`, or `neighbor_graph`; if no graph is found, it returns an empty list. Citations: `karospace/data_loader.py:853-861`.

The function converts the graph to CSR, row-normalizes it, rejects an all-zero graph, selects the top variable features from `adata.X`, centers the selected expression matrix, computes Moran's I as `n_obs * (Z * WZ).sum(axis=0) / (s0 * (Z * Z).sum(axis=0))`, clips results to `[-1, 1]`, rounds to four decimals, and sorts descending by `I`. Citations: `karospace/data_loader.py:863-899`.

The feature candidate list comes from `_select_top_variable_features`, which ranks features by variance across cells using sparse or dense matrix handling. Citations: `karospace/data_loader.py:840-850`.

## Export-Time Call Path

`export_to_html` currently keeps the legacy `spatial_variable_genes_n` parameter with default `20` for compatibility; its docstring describes the value as the number of top variable features scored with Moran's I, requires a spatial weight matrix in `adata.obsp`, and can be disabled with `0`. Citations: `karospace/exporter.py:34709-34765`, `karospace/exporter.py:34830-34838`.

During export, `export_to_html` calls `dataset.to_json_data(...)` while building the viewer payload and passes `analytics_modalities=selected_modalities` plus `spatial_variable_genes_n=int(spatial_variable_genes_n)`. Citations: `karospace/exporter.py:35199-35205`, `karospace/exporter.py:35257-35295`.

Inside `SpatialDataset.to_json_data`, all available features are recorded per modality in `features_by_modality`; requested/embedded features are handled separately, but the spatial calculation is passed `list(features_by_modality.get(modality_name) or [])`, so the spatial ranking uses the modality's available feature universe, not only the manually requested features. Citations: `karospace/data_loader.py:2422-2459`, `karospace/data_loader.py:3342-3349`.

For each selected analytics modality, `to_json_data` obtains a modality-specific AnnData object via `_adata_for_pseudobulk_modality`. For non-default modalities, that object is built from the modality matrix and var table while copying compatible `.obsm` and `.obsp` entries from the base AnnData, preserving the spatial graph used by Moran's I. Citations: `karospace/data_loader.py:2390-2412`, `karospace/data_loader.py:3342-3349`.

The computed lists are stored in the JSON payload under `spatial_variable_features_by_modality`; every modality is initialized with an empty list, then selected analytics modalities are filled when `spatial_variable_genes_n > 0`. Citations: `karospace/data_loader.py:3311-3313`, `karospace/data_loader.py:3342-3349`, `karospace/data_loader.py:3411-3417`.

## User-Facing Entry Points

CLI users control the calculation with `--spatial-variable-features-n`, whose help text says it scores top variable features with Moran's I, requires a spatial graph in `obsp`, and can be disabled with `0`. The legacy `--spatial-variable-genes-n` alias remains accepted and both forms are passed to `export_to_html` as `spatial_variable_genes_n=args.spatial_variable_genes_n`. Citations: `karospace/cli.py:550-555`, `karospace/cli.py:810-820`, `karospace/cli.py:830-868`.

GUI users control the same value through the "Moran features N" field, initialized to `20`, parsed as a non-negative integer, and added to `export_kwargs` as `spatial_variable_genes_n`; the GUI worker then calls `load_spatial_data(...)` followed by `export_to_html(dataset, **export_kwargs)`. Citations: `karospace/gui.py:277-292`, `karospace/gui.py:641-646`, `karospace/gui.py:1256-1283`, `karospace/gui.py:1328-1335`.

Python API users reach the same code through `karospace.export_to_html`; the package exposes `export_to_html`, `load_spatial_data`, and `SpatialDataset` from `karospace.__init__`, and `export_to_html` carries the same `spatial_variable_genes_n` parameter. Citations: `karospace/__init__.py:40-57`, `karospace/exporter.py:34709-34765`.

## Viewer Consumption

The generated JavaScript reads `DATA.spatial_variable_features_by_modality` via `getModalityPayload(...)`, and `getExplorationSpatialVariableFeaturesPayload(...)` returns the modality-specific list or `[]`. Citations: `karospace/exporter.py:7761-7770`, `karospace/exporter.py:7838-7841`.

When the active Insights tab is `features` and its subtab is `spatial`, the viewer calls `renderSpatialVariableGenes()`. Citations: `karospace/exporter.py:24765-24773`.

`renderSpatialVariableGenes()` renders an empty message when no spatial rankings were precomputed, otherwise it reads each row's `gene` and `I` value, presents list rows with rank and `I=...`, and can switch to a graph view. Citations: `karospace/exporter.py:28013-28024`, `karospace/exporter.py:28028-28049`, `karospace/exporter.py:28055-28073`.

The graph view plots the top entries on a "Moran Index" axis, and the calculation-info popover describes spatial features as ranked by Moran index, measuring whether nearby cells have similar feature values. Citations: `karospace/exporter.py:27137-27165`, `karospace/exporter.py:8433-8437`.

The feature discovery sidebar also surfaces the top 12 entries from `getExplorationSpatialVariableFeaturesPayload(...)` under "Spatially variable features" with `I` metadata. Citations: `karospace/exporter.py:14649-14666`.

## Documentation and Tests

The README documents `--spatial-variable-features-n` as "Top variable features scored with Moran's I" with default `20`, and documents that `Insights -> Features` includes spatial features scoped by the feature namespace selector. Citations: `README.md:329-331`, `README.md:408-414`.

`FEATURES_SUMMARY.md` documents that `Features > Spatial` shows Moran Index rankings computed at export for each selected modality and that spatially variable features are computed with Moran's I per selected modality for up to `spatial_variable_features_n` variable features on the full input cell set. Citations: `FEATURES_SUMMARY.md:155-166`, `FEATURES_SUMMARY.md:202-215`.

The multimodal export schema test builds an AnnData object with `obsm["spatial"]` and `obsp["spatial_connectivities"]`, calls `to_json_data(..., spatial_variable_genes_n=2, ...)`, and asserts that `spatial_variable_features_by_modality` has both `rna` and `protein` keys with expected features in each modality's spatial list. Citations: `tests/test_multimodal_export_schema.py:38-55`, `tests/test_multimodal_export_schema.py:182-224`.

The HTML runtime test checks that generated copy includes "Spatial features", confirming the browser output retains the user-facing feature label. Citations: `tests/test_multimodal_html_runtime.py:235-243`.

## Practical Interpretation

If "Features Spatial" is missing or empty in a generated viewer, the primary causes visible in the code are: `spatial_variable_genes_n` was set to `0`; the selected analytics modality did not run; or the AnnData object lacked one of the accepted `.obsp` graph keys. Citations: `karospace/data_loader.py:853-861`, `karospace/data_loader.py:3342-3349`, `karospace/exporter.py:35032-35039`, `karospace/exporter.py:35288-35292`.
