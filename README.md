# KaroSpace

**KaroSpace** is an interactive HTML viewer for exploring spatial transcriptomics data. It generates standalone HTML files from AnnData/H5AD or SpatialData inputs that can be shared and viewed in any web browser — no server or Python installation required.

Originally developed at Karolinska Institutet for visualizing Xenium spatial transcriptomics data across multiple tissue sections.

Visit [KaroSpace Website](https://karospace.se/).

![KaroSpace viewer showing a multi-section spatial transcriptomics grid](assets/karospace-readme-visual.png)

> [!NOTE]
> Live Demo ! **Pancreas viewer**: [Open hosted demo](https://christoffermattssonlangseth.github.io/KaroSpace/pancreas.html).

## Features

- [x] **Modal interaction** — Browse many sections in a responsive grid, then zoom and pan any section in detail
- [x] **Section filters** — Filter visible sections by exported metadata such as stage, condition, region, sex, model, sample, or batch
- [x] **UMAP to sections connection** — Lasso selection works in UMAP and modal view with synced highlights
- [x] **Cells selection composition** — Selected-cell totals and per-type counts with expandable scrollable lists
- [x] **Polygon regions** — Save lasso selections as persistent regions, reorder labels, and export JSON for downstream integration
- [x] **Region-to-region comparison** — Compare saved regions directly in the viewer, export JSON/CSV reports, and search top hits
- [x] **Cell search** — Select cells with query syntax based on annotations, features, or section metadata, then reuse the selection in summaries and comparisons
- [x] **Split screen** — Compare two variables side-by-side in the modal (`Annotation`, a selected feature modality, or `Module`)
- [x] **Feature modules** — Build custom feature sets, compute averaged module scores, display them like feature layers, and import/export module definitions
- [x] **Legend controls** — Toggle/hide categories and spotlight one class across grid and UMAP
- [x] **Feature exploration** — Search within a selected modality, inspect value distributions, review marker features, spatial features, category means, and related-feature suggestions
- [x] **Per cell comparison** — Live comparison of cell selections or regions with table and graph visualization (Welch test scores, log2FC, mean, percent detected)
- [x] **Per sample comparison** — Precomputed pseudobulk differential feature analysis using DESeq2 (PCA, distance matrices, volcano plots) with pathway enrichment for feature-supported modalities.
- [x] **Neighbor graph tools** — Graph overlay, hover rings (1–3 hops), enrichment, interactions, and dispersion summaries when `adata.obsp` contains a spatial graph
- [x] **Quality-of-life controls** — Hideable toolbar, screenshots, light/dark theme toggle, buttons explanation, keyboard shortcuts, and adjustable spot size
- [x] **Standalone export** — One self-contained HTML file, no backend required
- [x] **Compact sidecar** — Keep large feature matrices outside the HTML with lazy-loaded sidecar manifests and binary shards for lighter initial viewer files
- [x] **Shareable packages** — Export as `.karospace` bundles (ZIP + viewer HTML)

## Quick Start

A GUI version of KaroSpace (`KaroSpaceBuilder`) has been developed to allow researchers with moderate computational skills to create HTML file.

Prebuilt executables are available from the
[KaroSpaceBuilder Releases](https://github.com/christoffermattssonlangseth/KaroSpaceBuilder/releases) page:
- Apple Silicon: `KaroSpaceBuilder-macos-arm64.zip`
- Windows: `KaroSpaceBuilder-windows.zip`
- Linux: `KaroSpaceBuilder-linux.zip`

Download, unzip, and run — no Python required.

To install from source:

```bash
git clone https://github.com/christoffermattssonlangseth/KaroSpaceBuilder
cd KaroSpaceBuilder
python -m pip install "git+https://github.com/christoffermattssonlangseth/KaroSpace.git"
python -m pip install -e .
KaroSpaceBuilder
```

If KaroSpaceBuilder is already installed, launch with:

```bash
karospacebuilder   # or: karospace-gui
```

## Installation

```bash
git clone https://github.com/christoffermattssonlangseth/karospace.git
cd karospace
pip install -e .
```

> [!WARNING]
> Dependencies :
> ```bash
> - Python >= 3.9
> - scanpy >= 1.9.0
> - anndata >= 0.8.0
> - numpy >= 1.20.0
> - pandas >= 1.3.0
> - pydeseq2 >= 0.5.0
> - scipy >= 1.7.0
> - gseapy >= 1.1.0
> - tqdm >= 4.66.0
> ```

SpatialData input is optional. Install only when you want to load SpatialData `.zarr` :

```bash
pip install -e ".[spatialdata]"
```

If KaroSpace is already installed and you want to reinstall the local checkout after editing the source:

```bash
python -m pip uninstall karospace -y
python -m pip install -e .
```

Install the SpatialData extra in the same environment if `import spatialdata` fails:

```bash
python -m pip install -e ".[spatialdata]"
```

## Usage

### Python API

```python
from karospace import inspect_input_file, load_spatial_data, export_to_html

# Optional: inspect obs metadata first without exporting HTML or running analytics.
report = inspect_input_file("your_data.h5ad")
for column in report["metadata"]:
    print(column["name"], column["type"], column["examples"][:3])

dataset = load_spatial_data(
    "your_data.h5ad",
    section_key="sample_id",  # Column identifying each section
    section_metadata=["course", "region", "condition"],  # Section metadata shown in the visual params bar/filter chips
    section_metadata_extra=["patient_id", "slide_id"],  # Section metadata stored without visual params bar chips
    metadata_value_order={
        "course": ["naive", "peak_I", "peak_II", "peak_III"],
    },
)

export_to_html(
    dataset,
    output_path="viewer.html",
    main_cell_annotation="cell_type",    # Main cell-annotation column shown first
    title="KaroSpace",
    min_panel_size=150,          # Min panel width (responsive autoscaling)
    spot_size="auto",            # Adaptive by section density (or set a fixed number)
    downsample=30000,            # Max cells per section
    cell_annotations=[          # Extra cell obs annotation columns for annotation dropdowns
        "leiden",
        "niche",
    ],
    features=[                      # Pre-load features for visualization; matched across selected modalities
        "Cd4",
        "Cd8a",
        "Gfap",
    ],
    features_list=None,             # Optional text file with one feature per line
    feature_encoding="auto",        # "auto" | "dense" | "sparse"
    feature_storage="embedded",     # "embedded" | "sidecar"
    feature_manifest_path=None,          # Optional manifest path; defaults to viewer.features.json
    feature_sparse_zero_threshold=0.8,
    modalities=["rna", "protein"],  # Feature namespaces to export into the viewer
    neighbor_stats_annotations=["cell_type"],
    neighbor_stats_permutations=20,
    pseudobulk="auto",           # Use None to disable category pseudobulk DE in Python
    pseudobulk_additional_annotations=["niche"],
    pseudobulk_counts_layer="counts",
    pseudobulk_modalities=["rna"],  # Use ["all"] or e.g. ["rna", "protein"] to run DE on multiple modalities
    pseudobulk_min_cells_per_pseudobulk=20,
    pseudobulk_min_pct_expressed=0.2,
    pseudobulk_p_adjust_method="fdr_bh",
    pseudobulk_padj_cutoff=0.05,
    pseudobulk_log2fc_cutoff=1,
    pseudobulk_deseq2_fit_type="parametric",
    pseudobulk_n_cpus=1,
    pseudobulk_embed_top_n_per_comparison=2,
    pathway="auto",             # Use None to disable pathway enrichment in Python
    pathway_gmt=None,            # default cached Reactome; or pass "reactome.gmt"
    pathway_organism="Mouse",
    pathway_top_n=10,
    pathway_min_overlap=3,
    pathway_gsea_permutations=100,
    interaction_markers="auto",  # Use None to disable contact-conditioned marker DE in Python
    embed_reproducibility_info=True,  # Embed export arguments/settings in the HTML header popover
    source_input_path="your_data.h5ad",  # Optional provenance path shown in the reproducibility popover
    section_rotations={
        "sample_a": 37.5,
        "sample_b": -90,
    },
)
```

You can also pass an already-loaded `AnnData` object:

```python
dataset = load_spatial_data(
    adata,
    section_key="sample_id",
)
```

For SpatialData, pass either the object or a `.zarr` store. If the SpatialData object has multiple AnnData tables, select the table explicitly:

```python
import spatialdata as sd
from karospace import inspect_input_file, load_spatial_data, export_to_html

sdata = sd.read_zarr("your_spatialdata.zarr")

report = inspect_input_file(sdata, spatialdata_table="table")

dataset = load_spatial_data(
    sdata,
    spatialdata_table="table",  # or "cells", depending on your object
    section_key="sample_id",
)

export_to_html(dataset, "viewer.html", main_cell_annotation="cell_type")
```

### Command Line

Inspect an `.h5ad` without generating HTML or running analytics:

```bash
karospace your_data.h5ad --inspect-input
```

SpatialData `.zarr` stores are also supported:

```bash
karospace your_spatialdata.zarr -o viewer.html --main-cell-annotation cell_type --spatialdata-table table
```

```bash
karospace your_data.h5ad \
  -o viewer.html \
  --section-key sample_id \
  --section-metadata course,region,condition \
  --section-metadata-extra patient_id,slide_id \
  --metadata-value-order '{"course":["naive","peak_I","peak_II","peak_III"]}' \
  --main-cell-annotation cell_type \
  --title "KaroSpace" \
  --min-panel-size 150 \
  --spot-size auto \
  --downsample 30000 \
  --cell-annotations leiden, niche \
  --features Cd4,Cd8a,Gfap \
  --features-list features.txt \
  --feature-encoding auto \
  --feature-storage embedded \
  --feature-sparse-zero-threshold 0.8 \
  --modalities rna,protein \
  --neighbor-stats-annotations cell_type \
  --neighbor-permutations 20 \
  --pseudobulk auto \
  --pseudobulk-additional-annotations niche \
  --pseudobulk-counts-layer counts \
  --pseudobulk-modalities rna \
  --pseudobulk-min-cells-per-pseudobulk 20 \
  --pseudobulk-min-pct-expressed 0.2 \
  --pseudobulk-p-adjust-method fdr_bh \
  --pseudobulk-padj-cutoff 0.05 \
  --pseudobulk-log2fc-cutoff 1 \
  --pseudobulk-deseq2-fit-type parametric \
  --pseudobulk-n-cpus 1 \
  --pseudobulk-embed-top-n-per-comparison 2 \
  --pathway auto \
  --pathway-organism Mouse \
  --pathway-top-n 10 \
  --pathway-min-overlap 3 \
  --pathway-gsea-permutations 100 \
  --interaction-markers auto \
  --section-rotations sample_a:37.5,sample_b:-90
```

#### CLI Options

CLI value conventions:
- Use `auto` when KaroSpace should choose behavior automatically.
- Use `off` to disable analysis modes such as `--pseudobulk`, `--pathway`, and `--interaction-markers`.
- Use `none` where an option documents a nullable string value, such as `--outlineby` and `--pseudobulk-counts-layer`.
- Omit comma-separated/JSON options, or pass `""`, for an empty list or object.
- Use `0` for numeric disable switches.

##### Required

| Option | Description | Default |
|--------|-------------|---------|
| `input` | Path to input `.h5ad` file or SpatialData `.zarr` store | required |
| `-o, --output` | Output HTML file path | `karospace.html` |
| `--section-key` | Column to identify sections | `sample_id` |
| `--section-order` | Comma-separated section IDs to control section order | empty string |
| `--spatial-key` | Key in `adata.obsm` containing spatial coordinates, or target key created from `--spatial-x/--spatial-y` | `spatial` |
| `--main-cell-annotation` | Main cell-annotation column shown first in the viewer | `leiden` |
| `--section-metadata` | Comma-separated obs columns to use as section metadata shown in the visual params bar/filter chips | empty string |
| `--modalities` | Comma-separated modalities to export | all detected |

##### Inspection

| Option | Description | Default |
|--------|-------------|---------|
| `--inspect-input` | Read input metadata and exit without building sections, downsampling, exporting HTML, or running analytics | flag off |

##### Coordinates

| Option | Description | Default |
|--------|-------------|---------|
| `--spatial-x` | Obs/metadata column to use as X coordinates; requires `--spatial-y` | not set |
| `--spatial-y` | Obs/metadata column to use as Y coordinates; requires `--spatial-x` | not set |
| `--spatialdata-table` | AnnData table key to use when the input is a SpatialData object/store; required when multiple tables are present and no table named `table` exists | not set |

##### Annotations

| Option | Description | Default |
|--------|-------------|---------|
| `--cell-annotations` | Comma-separated extra cell obs annotation columns to embed as selectable cell annotations | empty string |
| `--section-metadata-extra` | Comma-separated obs columns to store as section metadata without visual params bar/filter chips | empty string |
| `--metadata-value-order` | JSON object mapping metadata columns to ordered value lists | empty string |
| `--metadata-max-columns` | Limit metadata columns used, preserving order | not set |
| `--metadata-labels` | JSON object mapping metadata/obs column keys to display labels in the viewer UI | empty string |

##### Viewer layout

| Option | Description | Default |
|--------|-------------|---------|
| `--outlineby` | Metadata column used to paint panel outlines; use `none` to disable | `none` |
| `--downsample` | Max cells per section | not set |
| `--title` | Page title | `KaroSpace` |
| `--tutorial` | (in development) Embed the static Story Mode HTML tutorial; users start it from the graduation-cap control and move through prepared viewer states with Next/Back | flag off |
| `--no-reproducibility-info` | Do not embed export arguments, thresholds, cutoffs, inputs, and resolved settings in the HTML reproducibility popover | flag off |
| `--min-panel-size` | Minimum panel width in pixels | `150` |
| `--spot-size` | Cell/spot size (`auto` or positive number) | `auto` |
| `--deconvolutions` | JSON object mapping deconvolution labels to obs/obsm keys | empty string |
| `--scalebar-unit` | Unit label for the scalebar | `μm` |
| `--viewer-info-html` | HTML string shown in the viewer Info tab | default info |
| `--viewer-info-html-file` | Path to an HTML fragment shown in the viewer Info tab | not set |

##### Feature content and storage

| Option | Description | Default |
|--------|-------------|---------|
| `--features` | Comma-separated features to preload. Requested names are matched against every exported `--modalities` namespace, so the same feature name is included in each selected modality where it exists. Significant pseudobulk DE features are embedded automatically up to the per-comparison cap | empty string |
| `--features-list` | Text file with one feature per line; combined with `--features`, deduplicated, and resolved across selected modalities | not set |
| `--feature-storage` | Feature storage mode: `embedded` stores requested/top DE feature vectors in the HTML; `sidecar` stores all feature vectors outside the HTML | `embedded` |
| `--feature-encoding` | Feature vector encoding (`auto`, `dense`, `sparse`) | `auto` |
| `--feature-value-encoding` | Sidecar/package feature value encoding for binary shards (`uint16`, `uint8`) | `uint16` |
| `--feature-manifest-path` | Path for the feature sidecar manifest JSON | derived from output path |
| `--feature-sidecar-shard-size` | Features per sidecar shard | `256` |
| `--feature-sparse-zero-threshold` | Zero fraction threshold for `auto` sparse encoding | `0.8` |

##### Neighborhoods

| Option | Description | Default |
|--------|-------------|---------|
| `--neighbor-permutations` | Permutations for neighbor enrichment z-scores | `auto` |
| `--neighbor-stats-annotations` | Obs columns for neighbor composition stats (`auto` or comma-separated; pass `""` to disable standalone neighbor enrichment) | `auto` |
| `--neighbor-stats-seed` | Random seed for neighbor enrichment permutations | `0` |

##### Interactions

| Option | Description | Default |
|--------|-------------|---------|
| `--interaction-markers` | Contact-conditioned pseudobulk marker mode (`auto`, `off`) | `auto` |
| `--interaction-markers-top-targets` | Target categories evaluated per source for contact-conditioned markers | `5` |
| `--interaction-markers-top-features` | Top DE features kept per source-target interaction | `20` |
| `--interaction-markers-min-cells` | Minimum cells per replicate contact+ and contact- pseudobulk sample | `30` |
| `--interaction-markers-min-neighbors` | Minimum target neighbors to classify contact+ source cells | `1` |

##### Connections

| Option | Description | Default |
|--------|-------------|---------|
| `--feature-correlation-top-n` | Correlated features shown per embedded feature in discovery panel | `5` |
| `--spatial-variable-features-n` | Top variable features scored with Moran's I; use `0` to disable | `20` |

##### Pseudobulk DE

| Option | Description | Default |
|--------|-------------|---------|
| `--pseudobulk` | Category pseudobulk DE mode (`auto`, `off`) | `auto` |
| `--pseudobulk-additional-annotations` | Additional annotation columns to analyze when pseudobulk or interaction markers are enabled. `--main-cell-annotation` is included automatically | empty string |
| `--pseudobulk-replicate-annotation` | Obs annotation to use as the biological replicate for pseudobulk analyses; defaults to `--section-key` | `--section-key` |
| `--pseudobulk-counts-layer` | Raw-count AnnData layer for pseudobulk aggregation; use `none` for `adata.X` | `counts` |
| `--pseudobulk-modalities` | Comma-separated modalities to run category pseudobulk DE and contact-conditioned interaction markers on. Use `all` for all detected modalities. This is independent of `--modalities`, which controls feature export for the viewer | dataset default modality |
| `--pseudobulk-simple-constrast-categories` | Categories to report in category-versus-category contrasts. With `--pseudobulk-additional-annotations`, use annotation-specific JSON wrapped in single quotes, such as `'{"cell_type":["Astrocyte","B cell"],"region":["Cortex"]}'`, or a nested list matching `[main-cell-annotation, additional...]` | empty string |
| `--pseudobulk-min-cell-counts` | Exclude cells with fewer than this many total raw counts before pseudobulk aggregation; use `0` to disable | `0` |
| `--pseudobulk-min-feature-counts` | Exclude features with fewer than this many total raw pseudobulk counts in the shared DESeq2 fit; use `0` to disable | `0` |
| `--pseudobulk-min-cells-per-pseudobulk` | Minimum cells required in each replicate × annotation pseudobulk sample before it can enter the shared DESeq2 fit | `20` |
| `--pseudobulk-min-replicates` | Minimum paired replicates required for each reported contrast | `2` |
| `--pseudobulk-min-pct-expressed` | Minimum fraction of cells with nonzero feature values required in at least one compared group before DE results are reported; values >1 are interpreted as percentages | `0` |
| `--pseudobulk-p-adjust-method` | Multiple-testing correction method (`fdr_bh`, `bonferroni`, `holm`, `none`) | `fdr_bh` |
| `--pseudobulk-padj-cutoff` | Adjusted p-value threshold for DE calls and plot coloring; DE features must pass `padj < cutoff` | `0.05` |
| `--pseudobulk-log2fc-cutoff` | Absolute log2FC cutoff for volcano highlighting and DE table inclusion | `1` |
| `--pseudobulk-deseq2-fit-type` | PyDESeq2 dispersion trend fit type; use `mean` to avoid parametric trend fallback warnings | `parametric` |
| `--pseudobulk-n-cpus` | CPU workers for the shared DESeq2 fit and maximum parallel shared-fit contrasts | `1` |
| `--pseudobulk-embed-top-n-per-comparison` | Significant DE features to auto-embed per category/contact comparison in embedded mode; ignored by sidecar mode because feature vectors are sidecar-loaded | `2` |

##### Pathway enrichment

| Option | Description | Default |
|--------|-------------|---------|
| `--pathway` | Pathway enrichment mode (`auto`, `off`) | `auto` |
| `--pathway-gmt` | GMT pathway file(s) for ORA/GSEA after Simple design DE; omitted uses cached/default Reactome when available, then falls back to GSEApy/Enrichr | Reactome |
| `--pathway-organism` | Organism used for default Reactome loading, e.g. `Human` or `Mouse` | `Mouse` |
| `--pathway-top-n` | Maximum ORA/GSEA pathways stored per direction and comparison | `10` |
| `--pathway-min-overlap` | Minimum pathway/query feature overlap for ORA/GSEA reporting | `3` |
| `--pathway-gsea-permutations` | Permutations for compact preranked GSEA p-values | `100` |

##### Images, overlays, and utilities

| Option | Description | Default |
|--------|-------------|---------|
| `--section-rotations` | Comma-separated `section_id:angle` pairs | empty string |
| `--section-images` | JSON object mapping section IDs to image paths/specs | empty string |
| `--section-images-max-px` | Maximum image dimension when embedding section images | `4096` |

> [!NOTE]
> See [FEATURES_SUMMARY.md](FEATURES_SUMMARY.md) for a guided overview of the main HTML viewer features.

## Data Requirements

KaroSpace accepts:

- An `.h5ad` file path
- An `AnnData` object passed through the Python API
- A SpatialData `.zarr` store path
- A SpatialData object passed through the Python API

Use `inspect_input_file(...)` in Python or `--inspect-input` on the CLI to list available `adata.obs` metadata, value types, example values, and missing-value counts before choosing annotation, section metadata, and pseudobulk arguments. This inspect path only reads the AnnData table and does not validate coordinates or run the export pipeline.

Internally, SpatialData input is normalized to one AnnData table before export. The selected table must satisfy the same requirements as a regular AnnData input:

- **`adata.obsm['spatial']`** — 2D coordinates for each cell (x, y)
- **`adata.obs[section_key]`** — Column identifying which section each cell belongs to
- **Categorical or numeric columns in `adata.obs`** — For assigning cell annotations and visualizing cells

For SpatialData tables, use `spatialdata_table="..."` / `--spatialdata-table ...` when the object contains more than one table. If the default `section_key="sample_id"` is missing, KaroSpace uses the table's SpatialData `region_key` automatically when available. If no per-cell region key exists, the table is exported as one section.

If a SpatialData `.zarr` store contains invalid non-table elements, such as images with missing transformations, KaroSpace falls back to reading the selected AnnData table directly from `tables/<table>` instead of failing the whole export. The image/label elements are ignored by this fallback; pass separate section images with `--section-images` if you want image overlays in the viewer.

If coordinates are stored as separate obs columns instead of an `obsm` matrix,
pass them on the CLI:

```bash
karospace your_data.h5ad -o viewer.html --spatial-x x_centroid --spatial-y y_centroid
```

This creates `adata.obsm["spatial"]` during loading. Use `--spatial-key` to pick
a different target key.

> [!TIP]
> See [`examples/`](examples/) for complete dataset-specific export scripts.

### Optional metadata

Use `section_metadata=[...]` / `--section-metadata ...` for section-level obs columns that should appear in the visual params bar and filter chips. If omitted, no section metadata columns are added by default. Use `section_metadata_extra=[...]` / `--section-metadata-extra ...` for section-level metadata that should be stored in the viewer payload but not shown as filter chips.

- `course` — Experimental phase (e.g., `"naive"`, `"peak_I"`); sections are outlined by this column when `outlineby="course"` / `--outlineby course` is used
- `region`, `condition`, `timepoint` — Typical section metadata shown as filter chips

Use `cell_annotations=[...]` / `--cell-annotations ...` for additional cell-level annotation columns that should be available in annotation dropdowns and comparison panels.

Control display order of metadata values and section ordering via `metadata_value_order`:

```python
dataset = load_spatial_data(
    "your_data.h5ad",
    section_key="sample_id",
    section_metadata="course,
    metadata_value_order={
        "course": ["naive", "peak_I", "peak_II", "peak_III"],
    },
)
```

### Optional category palettes

If `adata.uns["{col}_colors"]` exists (scanpy convention — list of hex aligned to `adata.obs[col].cat.categories`), KaroSpace uses it for that column everywhere (legend, spots, neighbor views, samples panel). Length mismatch or missing key falls back to the default palette.

### Optional neighborhood graph

If `adata.obsp` contains a neighbor graph (`spatial_connectivities`, `connectivities`, `neighbors`, or `neighbor_graph`), KaroSpace exposes graph overlay and neighbor-hover controls.

If the export is downsampled, the visible graph overlay and neighbor-hover controls are downsampled with the displayed cells, so edges to cells that were not exported are not shown.

`Insights → Neighbors → Enrichment` and `Interactions` use neighbor composition statistics for the selected Exploration annotation. If no graph or no stats exist for that annotation, the viewer shows a yellow warning and lists the annotations that do have neighbor stats. `Insights → Neighbors → Dispersion` is computed from all cells before HTML downsampling for the main cells annotation and any requested `cell_annotations`, then summarizes whether each category is clustered, random, or dispersed relative to the observed all-cell layout.

### Multimodal feature viewer

When multiple modalities are exported with `modalities=["rna", "protein"]` or `--modalities rna,protein`, the viewer treats each modality as its own feature namespace. In Default mode, switch from `Annotation` to `Feature`, then use the feature-namespace dropdown beside the source switch to choose RNA, protein, or `Module` when modules exist. The search box is scoped to that namespace, so a feature name shared by RNA and protein is loaded from the active namespace rather than shown with modality badges.

Split view keeps an independent source and feature namespace for layer A and layer B. This allows comparisons such as RNA feature versus protein feature, annotation versus protein, or module score versus RNA without changing the main visual namespace.

`Insights → Features` has its own feature-namespace selector for marker features, spatial features, per-cell distributions, per-sample/category means, and related-feature suggestions. `Insights → Compare → Per sample → Simple design` has a pseudobulk modality selector, and `Insights → Neighbors → Interactions` has an interaction-marker modality selector. Exported CSV/SVG filenames include the selected modality where a result is modality-specific.

Pathway enrichment is feature-supported only. It is computed for RNA-like modalities and reported as unavailable for unsupported modalities unless a future export provides an explicit feature-to-pathway mapping.

### Optional pseudobulk category selection

Pseudobulk category DE is precomputed automatically for the initial `main cells annotation` column unless `pseudobulk=None` in Python or `--pseudobulk off` on the CLI is used, and shown in `Insights → Compare → Per sample → Simple design`. KaroSpace aggregates raw counts by replicate and annotation, keeps replicate × annotation pseudobulk samples with at least `pseudobulk_min_cells_per_pseudobulk` / `--pseudobulk-min-cells-per-pseudobulk` cells, fits one shared `~ replicate + annotation` DESeq2 model per annotation column, then extracts category-versus-category contrasts. It also extracts a balanced-rest contrast for every category: the category minus the equally weighted mean of all other retained annotation categories. Features that do not reach `pseudobulk_min_pct_expressed` / `--pseudobulk-min-pct-expressed` in at least one compared group are removed from reported DE results, so they do not enter the contrast-level multiple-testing correction applied by KaroSpace. Pairwise PCA/distance diagnostics are generated automatically for selected category pairs.

By default, pseudobulk DE and contact-conditioned interaction markers run on the dataset default modality, usually `rna`. Use `pseudobulk_modalities=["rna", "protein"]` in Python or `--pseudobulk-modalities rna,protein` on the CLI to run those analyses on selected modalities, or use `all` for every detected modality. Results are stored only in modality-keyed payloads such as `pseudobulk_de_by_modality`, `interaction_markers_by_modality`, `category_feature_means_by_modality`, `feature_correlations_by_modality`, `spatial_variable_features_by_modality`, and `pathway_settings_by_modality`.

When selecting specific pairwise categories from the command line, wrap listed values in single quotes:

```bash
--main-cell-annotation cell_type \
--pseudobulk-additional-annotations region \
--pseudobulk-simple-constrast-categories '{"cell_type":["Astrocyte","B cell"],"region":["Cortex"]}'
```

## Deployment and Sharing

KaroSpace has three practical export modes:

| Mode | Output | Best for | How to open |
| --- | --- | --- | --- |
| Embedded HTML | `viewer.html` | Small to medium feature payloads, easiest sharing | Double-click or open the file in a browser |
| Sidecar viewer | `viewer.html` + `viewer.features.json` + `viewer.features/` | Large feature payloads with lazy loading | Serve the directory over HTTP(S), then open the HTML URL |
| `.karospace` package | `viewer.karospace` + optional `viewer.loader.html` | One-file sharing of a sidecar viewer | Drop the package into the hosted loader or the generated local loader |

Sidecar mode keeps the initial HTML smaller by moving feature vectors into a manifest and binary shard files. The viewer fetches those shards only when a feature is needed. This is useful when many features or modalities would make a single HTML file too large.

### Create a sidecar viewer

CLI:

```bash
karospace your_data.h5ad \
  -o viewer.html \
  --main-cell-annotation cell_type \
  --features Cd4,Cd8a,Gfap \
  --feature-storage sidecar \
  --feature-manifest-path viewer.features.json
```

Python API:

```python
from karospace import load_spatial_data, export_to_html

dataset = load_spatial_data("your_data.h5ad", section_key="sample_id")

export_to_html(
    dataset,
    output_path="viewer.html",
    main_cell_annotation="cell_type",
    features=["Cd4", "Cd8a", "Gfap"],
    features_list="features.txt",
    feature_storage="sidecar",
    feature_manifest_path="viewer.features.json"
)
```

This writes:

```text
viewer.html
viewer.features.json
viewer.features/
  000.bin
  001.bin
  ...
```

> [!IMPORTANT]
> Keep all three elements together. The HTML contains the viewer and embedded summary data, but no feature vectors in sidecar mode; `viewer.features.json` is the sidecar manifest; `viewer.features/` contains the binary feature shards.

### Open a sidecar viewer

> [!CAUTION]
> Do not open sidecar HTML directly with `file://`; browsers block local shard loading.

Serve the output directory instead:

```bash
python -m http.server --directory /path/to/output-dir 8000
```

Then open:

```text
http://localhost:8000/viewer.html
```

> [!TIP]
> For deployment, upload the HTML, manifest, and shard directory with the same relative paths to GitHub Pages, S3, an institutional web server, or a lab intranet. If `viewer.html` references `viewer.features.json`, then `viewer.features.json` must be next to the HTML and its `viewer.features/` shard directory must also be next to the HTML unless you intentionally used matching custom paths.

### Create a `.karospace` package directly

Use `.karospace` output when you want sidecar loading with a more compact shareable file:

CLI :

```bash
karospace your_data.h5ad \
  -o viewer.karospace \
  --main-cell-annotation cell_type \
  --feature-storage sidecar \
  --feature-manifest-path viewer.features.json
```

Python API:

```python
export_to_html(
    dataset,
    output_path="viewer.karospace",
    main_cell_annotation="cell_type",
    feature_storage="sidecar",
    feature_manifest_path="viewer.features.json",
)
```

Direct package export writes:

```text
viewer.karospace
viewer.loader.html
```

The `.karospace` file is a ZIP-based package containing `index.html`, `karospace-package.json`, the sidecar manifest, and the binary shard directory. The sibling `viewer.loader.html` is a local opener; it is not part of the package itself.

> [!NOTE]
> Open the package by visiting the hosted loader at [karospace.se/open](https://karospace.se/open) or by opening `viewer.loader.html` and dropping/selecting `viewer.karospace`. Package loading happens in the browser; the package is read locally by the browser and is not uploaded by the local loader.

### Package an existing sidecar into `.karospace`

If you already have an unpacked sidecar viewer, package it without recomputing analytics:

```bash
# Short form: auto-detect sidecar paths from the HTML.
karospace package-sidecar viewer.html --output viewer.karospace

# Explicit form: use this when the manifest or shard directory is not next to the HTML.
karospace package-sidecar viewer.html \
  --output viewer.karospace \
  --feature-manifest-path viewer.features.json \
  --feature-shard-dir viewer.features \
  --loader-output viewer.loader.html
```

### Sidecar troubleshooting

> [!TIP]
> - **The viewer opens but features do not load**: check that the HTML is served over HTTP(S), not opened with `file://`.
> - **404 for `viewer.features.json` or `.bin` shards**: keep `viewer.html`, `viewer.features.json`, and `viewer.features/` in the same relative layout used at export time.
> - **Custom `--feature-manifest-path`**: for normal sidecar HTML it may be a path; for direct `.karospace` export it must be a filename inside the package.
> - **Large packages**: increase `--feature-sidecar-shard-size` for fewer shard files, decrease it for smaller individual requests.

## Integrate polygon regions back into AnnData

```python
import scanpy as sc
from karospace import integrate_polygon_annotations

adata = sc.read_h5ad("your_data.h5ad")
integrate_polygon_annotations(
    adata,
    "karospace-annotations-2026-02-12T12-00-00-000Z.json",
    label_key="lesion_labels",
    count_key="lesion_label_count",
    uns_key="lesion_polygons",
)
adata.write_h5ad("your_data_with_polygons.h5ad")
```

## Browser Considerations

> [!WARNING]
> KaroSpace is canvas-heavy. Chrome/Chromium is generally fastest; Safari can be noticeably slower on large datasets. The viewer caps canvas DPR at `1.0` in Safari by default to reduce pixel work on Retina displays.
> For large datasets:
> - Use `downsample` to limit cells per section
> - Lower `min_panel_size` to reduce pixels drawn per thumbnail
> - Keep the neighbor graph toggle off unless needed

## License

MIT License

## Author

Christoffer Mattsson Langseth — Karolinska Institutet
