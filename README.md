<p align="center">
  <img src="assets/logo.png" alt="KaroSpace Logo" width="800">
</p>

**KaroSpace** is an interactive HTML viewer for exploring spatial transcriptomics data. It generates standalone HTML files from h5ad files that can be shared and viewed in any web browser without requiring a server or Python installation.

Originally developed at Karolinska Institutet for visualizing Xenium spatial transcriptomics data across multiple tissue sections.

## Live Demo
- [KaroSpace website](https://karospace.se/)
- **Pancreas viewer (GitHub Pages)**: [Open hosted demo](https://christoffermattssonlangseth.github.io/KaroSpace/pancreas.html)

## Features

For a full feature inventory, see [FEATURES_SUMMARY.md](FEATURES_SUMMARY.md).
For the shareable package format, see [KAROSPACE_PACKAGE_FORMAT_SPEC.md](KAROSPACE_PACKAGE_FORMAT_SPEC.md).
The reference browser opener for `.karospace` packages is [karospace-package-loader.html](karospace-package-loader.html).
When exporting `something.karospace`, KaroSpace also writes `something.loader.html` as a local no-install opener for that package workflow.

- **Grid + modal exploration** - Browse many sections in a responsive grid, then zoom/pan any section in detail
- **Per-section rotation** - Set exact initial section angles at export time and adjust them interactively in the viewer
- **Linked UMAP + section selection** - Magic Wand lasso works in UMAP and modal view with synced highlights
- **Selection summaries** - Selected-cell totals and per-type counts with expandable, scrollable lists
- **Polygon annotation workflow** - Draw multiple persistent polygons, manage labels, and export JSON for downstream `adata` integration
- **Region-to-region DE** - Compare drawn annotations directly in the viewer, export JSON/CSV reports, and Google-search top hits
- **Split compare slider** - Compare A/B variables in modal view (`Cell type` or `Gene`, including `All categories`)
- **Legend controls + spotlight** - Toggle/hide categories and spotlight one class across grid and UMAP
- **Flexible coloring + gene discovery** - Switch between multiple annotation columns, fuzzy-search genes, and reuse recent or saved genes
- **Insights panel** - `Summary`, `Compare`, `Genes`, `Neighborhood`, and `Regions` tabs with marker genes, pairwise cluster DE, neighbor composition, interaction markers, and region comparison tools
- **Modal selection workflow** - Lasso in the sample view, open a focused subview, browse `Genes in selection`, and use `Space` + drag to pan while Select or Annotate is active
- **Shareable sidecar packaging** - Export sidecar viewers as `.karospace` bundles with matching local loader HTML, or package an existing sidecar bundle later via the CLI
- **Compact sidecar options** - JSON and binary sidecar formats, sparse-first encoding, and optional `uint16` / `uint8` quantization for large datasets
- **Metadata-aware browsing** - Filter by metadata and optionally outline sections by `course` (or another column)
- **Optional neighbor graph tools** - Graph overlay + hover rings (1–3 hops) when `adata.obsp` graph exists
- **Quality-of-life controls** - Compact hideable sample-view toolbar, screenshots, theme toggle, keyboard shortcuts, and adjustable spot size
- **Standalone export** - One self-contained HTML file, no backend required

## Browser Considerations

KaroSpace is canvas-heavy. Chrome/Chromium browsers are generally fastest; Safari can be noticeably slower on large datasets.
To mitigate this, the viewer caps canvas device pixel ratio (DPR) in Safari at `1.0` by default, reducing pixel work on Retina displays.

Performance tips for large datasets:
- Reduce `downsample` to limit cells per section.
- Lower `min_panel_size` to reduce the number of pixels drawn.
- Keep the neighbor graph toggle off unless needed.

## Installation

### From source

```bash
git clone https://github.com/christoffermattssonlangseth/karospace.git
cd karospace
pip install -e .
```

### Dependencies

- Python >= 3.9
- scanpy >= 1.9.0
- anndata >= 0.8.0
- numpy >= 1.20.0
- pandas >= 1.3.0
- scipy >= 1.7.0

## Quick Start

### Desktop GUI (KaroSpaceBuilder)

#### Generate HTML via the KaroSpaceBuilder repository (straightforward setup)

If you prefer a guided GUI workflow, use the dedicated
[`KaroSpaceBuilder`](https://github.com/christoffermattssonlangseth/KaroSpaceBuilder)
repository:

Prebuilt executables are available from
[`KaroSpaceBuilder` Releases](https://github.com/christoffermattssonlangseth/KaroSpaceBuilder/releases):
- Apple Silicon: `KaroSpaceBuilder-macos-arm64.zip`
- Windows: `KaroSpaceBuilder-windows.zip`
- Linux: `KaroSpaceBuilder-linux.zip`

If you want the no-code path, download the correct zip for your machine from the
Releases page, unzip it, and run the app.

```bash
git clone https://github.com/christoffermattssonlangseth/KaroSpaceBuilder
cd KaroSpaceBuilder
python -m pip install "git+https://github.com/christoffermattssonlangseth/KaroSpace.git"
python -m pip install -e .
KaroSpaceBuilder
```

If you are developing both repos locally, install KaroSpace from source instead:

```bash
python -m pip install -e /path/to/spatial-viewer
```

In the app:
1. Set input `.h5ad` and output directory
2. Click `Inspect H5AD`
3. Choose groupby/color settings (or a preset)
4. Click `Export`

Output is written to:

```text
<output_dir>/index.html
```

If KaroSpaceBuilder is already installed in your current environment, launch directly with:

```bash
karospacebuilder
```

Compatibility launcher (same GUI):

```bash
karospace-gui
```

#### Quick Start (GUI)

1. Launch `karospacebuilder`
2. Pick a preset (`Default`, `Pancreas`, or `Lightweight`)
3. In `Basic`, choose input `.h5ad` and click `Inspect`
4. In `Colors & Genes`, build `additional_colors` and `genes` using searchable `+ Add` / `Remove` editors
5. Open `Advanced` only if needed (collapsed by default), then export

#### Tabs

- `Basic`: input/output paths, grouping/color setup, metadata loading options
- `Colors & Genes`: searchable list builders for color columns and hand-picked genes
- `Advanced`: expandable panel for encoding, neighbor stats, marker genes, interaction markers
- `Help`: in-app field guide with expected input formats and examples

#### Searchable List Builders

- `additional_colors`: search `adata.obs` columns, add/remove entries interactively
- `genes`: search `adata.var_names`, add/remove hand-picked genes interactively
- Advanced groupby lists:
  - `neighbor_stats_groupby`
  - `marker_genes_groupby`
  - `interaction_markers_groupby`

These choices are populated after `Inspect` reads the `.h5ad` object.

#### Presets

- `Default`: balanced settings and starter examples
- `Pancreas`: prefilled to mirror `examples/pancreas.py`-style workflow
- `Lightweight`: smaller/faster export defaults for quick iteration

The GUI also includes:
- Built-in field validation and error dialogs
- Scrollable layout for smaller windows
- Export log panel with progress and failures

### Python API

```python
from karospace import load_spatial_data, export_to_html

# Load your h5ad file
dataset = load_spatial_data(
    "your_data.h5ad",
    groupby="sample_id",  # Column identifying each section
    # Optional: custom ordering for metadata values (affects filter chips + outlines)
    # If group_order isn't set, the first key here is also used to order sections.
    metadata_value_order={
        "course": ["naive", "peak_I", "peak_II", "peak_III"],
    },
)

# Export to HTML
export_to_html(
    dataset,
    output_path="viewer.html",
    color="cell_type",           # Initial color column
    title="KaroSpace",
    min_panel_size=150,          # Min panel width (responsive autoscaling)
    spot_size="auto",            # Adaptive by section density (or set a fixed number)
    downsample=30000,            # Max cells per section (for large datasets)
    theme="light",               # "light" or "dark"
    additional_colors=[          # Extra columns for color dropdown
        "leiden",
        "condition",
    ],
    genes=[                      # Pre-load genes for expression view
        "Cd4",
        "Cd8a",
        "Gfap",
    ],
    # For zero-inflated expression matrices, store gene vectors sparsely (smaller HTML)
    gene_encoding="auto",        # "auto" | "dense" | "sparse"
    gene_storage="embedded",     # "embedded" | "sidecar" (manifest + shards for lazy loading)
    gene_aux_path=None,          # Optional manifest path; defaults to viewer.genes.json
    gene_sparse_zero_threshold=0.8,
    pack_arrays=True,            # Pack coords/colors/UMAP as base64 typed arrays (smaller + faster load)
    pack_arrays_min_len=1024,
    use_hvgs=True,               # Use adata.var['highly_variable'] when available
    hvg_limit=20,                # Max number of HVGs to include
    marker_genes_groupby=None,   # Marker genes off by default (enable with a categorical column)
    marker_genes_top_n=30,       # Top N markers per group
    neighbor_stats_groupby=[         # Compute neighbor composition stats for these categorical obs columns
        "cell_type",
    ],
    neighbor_stats_permutations=20,   # Permutations for neighbor enrichment z-scores (0 disables)
    neighbor_stats_seed=0,       # Random seed for permutations
    interaction_markers_groupby=None,           # Interaction markers off by default
    interaction_markers_top_targets=8,          # Targets evaluated per source type
    interaction_markers_top_genes=20,           # Genes shown per source-target interaction
    interaction_markers_min_cells=30,           # Minimum cells in contact+ and contact- groups
    interaction_markers_min_neighbors=1,        # Source cell needs >= this many target neighbors to be contact+
    interaction_markers_method="wilcoxon",      # DE method for contact+ vs contact-
    interaction_markers_layer="normalized",     # Layer used for DE (falls back to adata.X if missing)
    cluster_de_groupby=["cell_type"],           # Optional pairwise cluster-vs-cluster DE for Insights -> Genes -> Compare
    cluster_de_top_n=20,                        # Top genes kept per source/reference pair
    cluster_de_method="wilcoxon",               # scanpy.tl.rank_genes_groups method
    cluster_de_layer="normalized",              # Layer used for pairwise DE (falls back to adata.X if missing)
    cluster_de_min_cells=20,                    # Minimum cells required in both clusters
    section_rotations={                         # Optional exact initial per-section rotation angles
        "sample_a": 37.5,
        "sample_b": -90,
    },
)
```

### Command Line

```bash
karospace your_data.h5ad -o viewer.html --color leiden
```

#### CLI Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --output` | Output HTML file path | `karospace.html` |
| `-c, --color` | Initial color column | `leiden` |
| `-g, --groupby` | Column to group sections by | `sample_id` |
| `--min-panel-size` | Minimum panel width in pixels (responsive autoscaling) | `150` |
| `--spot-size` | Cell/spot size (`auto` or positive number) | `auto` |
| `--downsample` | Max cells per section | None |
| `--theme` | Color theme (`light` or `dark`) | `light` |
| `--title` | Page title | `KaroSpace` |
| `--gene-encoding` | Gene vector encoding (`auto`, `dense`, `sparse`) | `auto` |
| `--gene-storage` | Gene storage mode (`embedded`, `sidecar`) | `embedded` |
| `--gene-aux-path` | Optional path for the gene sidecar manifest JSON | auto |
| `--gene-sparse-zero-threshold` | Zero fraction threshold for `auto` sparse encoding | `0.8` |
| `--no-pack-arrays` | Disable base64 packing of large per-section arrays | off |
| `--pack-arrays-min-len` | Only pack arrays when section cell count ≥ this value | `1024` |
| `--neighbor-permutations` | Permutations for neighbor enrichment z-scores (`auto`, `0`, `100`, …) | `auto` |
| `--neighbor-stats-groupby` | Comma-separated obs columns to compute neighbor composition stats for (`auto`, empty disables) | `auto` |
| `--marker-genes-groupby` | Comma-separated obs columns to compute marker genes for (empty disables) | empty |
| `--interaction-markers-groupby` | Comma-separated obs columns to compute contact-conditioned markers for (empty disables) | empty |
| `--cluster-de-groupby` | Comma-separated categorical obs columns to precompute pairwise cluster DE for (`Insights -> Genes -> Compare`) | empty |
| `--cluster-de-top-n` | Top genes kept per source/reference cluster pair | `20` |
| `--cluster-de-method` | `scanpy.tl.rank_genes_groups` method used for pairwise cluster DE | `wilcoxon` |
| `--cluster-de-layer` | AnnData layer used for pairwise cluster DE when present | `normalized` |
| `--cluster-de-min-cells` | Minimum cells required in both clusters to report pairwise DE | `20` |
| `--section-rotations` | Comma-separated `section_id:angle` pairs for exact initial section rotations | empty |
| `--gene-correlation-top-n` | Number of correlated genes shown per embedded gene in the discovery panel (`0` disables) | `10` |

## Data Requirements

Your h5ad file should have:

- **`adata.obsm['spatial']`** - 2D coordinates for each cell (x, y)
- **`adata.obs[groupby]`** - Column identifying which section each cell belongs to
- **Categorical or numeric columns in `adata.obs`** - For coloring cells (e.g., cell types, clusters)

### Optional metadata columns

For filtering functionality, include these columns in `adata.obs`:
- `course` - Experimental course/phase (e.g., "peak_I", "peak_II", "naive")
- `region` - Tissue region
- `condition` - Experimental condition
- `timepoint` - Time point

Sections will be outlined with colors based on their `course` metadata if present.

### Optional metadata ordering

You can control the display order of metadata values (filter chips + outline legend)
and section ordering by passing `metadata_value_order` to `load_spatial_data`.
If `group_order` is not provided, sections are ordered by the first key in
`metadata_value_order` (unknown values appear after the custom list).

```python
dataset = load_spatial_data(
    "your_data.h5ad",
    groupby="sample_id",
    metadata_value_order={
        "course": ["naive", "peak_I", "peak_II", "peak_III"],
    },
)
```

### Optional neighborhood graph

If `adata.obsp` contains a neighbor graph (e.g., `spatial_connectivities`, `connectivities`,
`neighbors`, or `neighbor_graph`), KaroSpace will expose graph/neighbor-hover controls.
Neighbor composition stats are **enabled by default** for the initial color. To customize,
pass `neighbor_stats_groupby=[...]` to `export_to_html` (or `--neighbor-stats-groupby`
in the CLI).

To enable contact-conditioned interaction markers in the interaction browser, pass
`interaction_markers_groupby=[...]` to `export_to_html` (typically the same categorical
column used for coloring, e.g. `cell_type`).

Interaction markers are defined per source-target pair as:
- **contact+**: source cells with at least `interaction_markers_min_neighbors` neighbors of the target type
- **contact-**: source cells of the same source type with zero neighbors of that target type

To enable pairwise cluster-vs-cluster differential expression in `Insights -> Genes -> Compare`,
pass `cluster_de_groupby=[...]` to `export_to_html` (or `--cluster-de-groupby` in the CLI).
Each requested categorical column is precomputed at export time, and the viewer lets you compare
any source/reference cluster pair without a backend.

To expose region-to-region differential expression in `Insights -> Regions`, keep genes available
either embedded in the viewer or via sidecar mode. Sidecar viewers can also run full region DE in
the browser by streaming the existing sidecar data without a backend.

## Examples

See the scripts in [`examples/`](examples/) for complete dataset-specific exports.

## Viewer Controls

### Grid View
- **Click a section** - Open detailed modal view
- **Color dropdown** - Switch between annotation columns
- **Gene input** - Opens a discovery panel with fuzzy search, keyboard navigation, recent genes, saved panels, and marker-based suggestions
- **Size slider** - Adjust spot size (drag or +/- buttons)
- **Filter chips** - Click to filter sections by metadata
- **Legend items + Spotlight** - Toggle categories and optionally spotlight one across grid + UMAP
- **Legend button** - Show/hide the legend panel
- **Insights button** - Toggle the insights panel
- **Insights tabs** - `Summary`, `Compare`, `Genes`, `Neighborhood`, and `Regions`
- **Genes subtabs** - `Dotplot` and `Markers`
- **Screenshot button** - Download a snapshot of the current view
- **Theme button** - Toggle dark/light mode
- **Footer badge** - “KaroSpace” label with a GitHub link

### Modal View (Detailed Section)
- **Mouse wheel** - Zoom in/out
- **Click and drag** - Pan around
- **Space + drag** - Temporarily pan even when `Select` or `Annotate` is active
- **Zoom buttons** - +/- zoom controls
- **Rotate buttons** - Rotate each section view with exact stored angles and quick step buttons
- **Reset button** - Return to default zoom/pan
- **Split button** - Open the A/B comparison panel in modal view
- **A/B variable selectors** - Choose per side: `Cell type` (single category or `All categories`) or `Gene`
- **Split slider** - Controls left/right display share (e.g. 10% = left 10% uses A, right 90% uses B)
- **Magic Wand** - Draw lasso selection directly on the section
- **Genes in selection** - Open an on-demand ranked gene panel from the lasso selection summary
- **Focused view** - Open the current single-section selection as a temporary filtered subview and return with `Back view`
- **Annotate + Annotation panel** - Draw persistent polygons, rename/select/delete, clear section/all
- **Hide tools** - Collapse the bottom sample-view toolbar when you want an unobstructed canvas
- **Draggable grouped toolbar** - View, Tools, Graph, and Actions controls live in a compact dock that can be repositioned
- **Export JSON** - Download all polygon annotations including vertices + mapped cell indices
- **Screenshot button** - Download a PNG of the current sample (modal) view
- **Clear selection** - Remove current selected cells
- **Graph / Neighbors / Hop selector** - Available when a neighbor graph exists
- **Size slider** - Adjust spot size for this view
- **Selection summary** - Selected-cell count + per-type counts; expandable and scrollable
- **Hover over cells** - See cell type or expression value
- **Escape or click outside** - Close modal

### UMAP View (if available)
If your h5ad file contains UMAP coordinates (`adata.obsm['X_umap']`), a UMAP toggle button appears:
- **UMAP button** - Toggle the UMAP panel on/off
- **Dock button (TR/TL/BR/BL)** - Cycle the panel corner placement
- **Panel +/-** - Resize the UMAP panel
- **Magic Wand** - Activate lasso selection mode
- **Draw selection** - Click and drag to draw a selection area
- **Clear** - Clear the current selection
- **Size slider** - Adjust point size in the UMAP view
- **Mouse wheel** - Zoom the UMAP view
- **Selection summary** - Same synced summary behavior as modal view (count + expandable/scrollable category list)

## Deployment and Sharing

- **Single-file HTML**: Default `export_to_html(..., gene_storage="embedded")` writes one standalone HTML file with embedded data + viewer logic
- **Sidecar mode**: `gene_storage="sidecar"` writes the HTML plus `<name>.genes.json` and a `<name>.genes/` shard directory, allowing downstream lazy loading of non-embedded genes
- **Local use**: Embedded mode opens directly in any modern browser. Sidecar mode requires serving the files over HTTP(S) because browsers block `fetch()` from `file://`.
- **Static hosting**: Both modes work on static hosting (e.g., GitHub Pages, S3, lab intranet). Sidecar mode needs the HTML, manifest JSON, and shard directory deployed together.
- **GitHub Pages (this repo)**: Use **Source = GitHub Actions** (static HTML artifact deploy). Do **not** use Jekyll/deploy-from-branch mode for this workflow.
- **Pancreas deployment workflow**: `.github/workflows/pages-balo.yml` publishes `pancreas.html` and creates an `index.html` redirect.
- **Pancreas test URL**: `https://christoffermattssonlangseth.github.io/KaroSpace/pancreas.html` (available after a successful Pages workflow run).
- **Sharing**: In embedded mode, send the HTML file directly and it will work offline. In sidecar mode, share the HTML, `.genes.json`, and matching `.genes/` directory together via a served/static-hosted location.
- **Package existing sidecar**: If you already have a sidecar export and want a shareable `.karospace` bundle without recomputing the viewer, package it afterwards with the CLI. This writes both `something.karospace` and `something.loader.html`.
- **Updates**: Re-run `export_to_html` to refresh the file when annotations or metadata change.
- **File size note**: Large datasets create large HTML files. Consider `downsample` and limiting `genes`/`additional_colors`.

### Package Existing Sidecar Into `.karospace`

Short form, when the HTML already points at the correct sidecar files in the same directory:

```bash
karospace package-sidecar BALO.html --output BALO.karospace
```

Explicit form, when you want to pass the existing sidecar paths directly:

```bash
karospace package-sidecar BALO.html \
  --output BALO.karospace \
  --gene-aux-path BALO.genes.json \
  --gene-shard-dir BALO.genes \
  --loader-output BALO.loader.html
```

You can also use the shorter alias:

```bash
karospace package BALO.html --output BALO.karospace
```

The packaging command is intended for existing sidecar bundles and does not recompute viewer analytics or rewrite the HTML viewer. It simply wraps the current sidecar assets into a `.karospace` archive and writes the matching local loader HTML.

Selected cells are highlighted with a yellow/gold outline in both UMAP and spatial views.

### Integrate Polygon Annotations Back Into AnnData

Polygon exports use global cell indices when available, so you can map annotations directly back to the original `adata`.

```python
import scanpy as sc
from karospace import integrate_polygon_annotations

adata = sc.read_h5ad("your_data.h5ad")
integrate_polygon_annotations(
    adata,
    "karospace-annotations-2026-02-12T12-00-00-000Z.json",
    label_key="lesion_labels",           # per-cell joined labels
    count_key="lesion_label_count",      # number of polygons covering each cell
    uns_key="lesion_polygons",           # full polygon metadata
)

adata.write_h5ad("your_data_with_polygons.h5ad")
```

## Performance Tips

- Use `downsample` parameter for datasets with >50,000 cells per section
- Limit `genes` list to only essential genes (each adds to file size)
- If you enable `use_hvgs`, the viewer preloads up to 20 HVGs to limit file size
- Consider splitting very large datasets into multiple viewers

## License

MIT License

## Author

Christoffer Mattsson Langseth - Karolinska Institutet
