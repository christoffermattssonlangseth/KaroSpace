# KaroSpace

**KaroSpace** is an interactive HTML viewer for exploring spatial transcriptomics data. It generates standalone HTML files from h5ad files that can be shared and viewed in any web browser — no server or Python installation required.

Originally developed at Karolinska Institutet for visualizing Xenium spatial transcriptomics data across multiple tissue sections.

## Live Demo
- [KaroSpace website](https://karospace.se/)
- **Pancreas viewer**: [Open hosted demo](https://christoffermattssonlangseth.github.io/KaroSpace/pancreas.html)

## Features

- **Grid + modal exploration** — Browse many sections in a responsive grid, then zoom and pan any section in detail
- **Per-section rotation** — Set exact initial section angles at export time and adjust them interactively
- **Linked UMAP + section selection** — Magic Wand lasso works in UMAP and modal view with synced highlights
- **Selection summaries** — Selected-cell totals and per-type counts with expandable scrollable lists
- **Polygon annotation workflow** — Draw persistent polygons, manage labels, and export JSON for downstream `adata` integration
- **Region-to-region DE** — Compare drawn annotations directly in the viewer, export JSON/CSV reports, and search top hits
- **Split compare slider** — Compare two variables side-by-side in the modal (`Cell type` or `Gene`, including `All categories`), draggable directly on the canvas
- **Legend controls + spotlight** — Toggle/hide categories and spotlight one class across grid and UMAP
- **Flexible coloring + gene discovery** — Switch between annotation columns, fuzzy-search genes, and reuse recent or saved gene panels
- **Insights panel** — `Summary`, `Compare`, `Genes`, `Neighborhood`, and `Regions` tabs with marker genes, pairwise cluster DE, neighbor composition, interaction markers, and region comparison
- **Modal selection workflow** — Lasso in the sample view, open a focused subview, browse `Genes in selection`, and use `Space` + drag to pan while Select or Annotate is active
- **Shareable packages** — Export as `.karospace` bundles (ZIP + viewer HTML) for offline sharing; open via the hosted loader at [karospace.se/open](https://karospace.se/open) or a local `loader.html`. See [KAROSPACE_PACKAGE_FORMAT_SPEC.md](KAROSPACE_PACKAGE_FORMAT_SPEC.md)
- **Compact sidecar options** — JSON and binary shard formats, sparse-first encoding, and optional `uint16`/`uint8` quantization for large datasets
- **Metadata-aware browsing** — Filter sections by metadata and outline by course or another column
- **Neighbor graph tools** — Graph overlay and hover rings (1–3 hops) when `adata.obsp` contains a spatial graph
- **Quality-of-life controls** — Hideable toolbar, screenshots, theme toggle, keyboard shortcuts, and adjustable spot size
- **Standalone export** — One self-contained HTML file, no backend required

## Browser Considerations

KaroSpace is canvas-heavy. Chrome/Chromium is generally fastest; Safari can be noticeably slower on large datasets. The viewer caps canvas DPR at `1.0` in Safari by default to reduce pixel work on Retina displays.

For large datasets:
- Use `downsample` to limit cells per section
- Lower `min_panel_size` to reduce pixels drawn per thumbnail
- Keep the neighbor graph toggle off unless needed

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

Prebuilt executables are available from the
[`KaroSpaceBuilder` Releases](https://github.com/christoffermattssonlangseth/KaroSpaceBuilder/releases) page:
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

**Workflow:**
1. Launch the app and pick a preset (`Default`, `Pancreas`, or `Lightweight`)
2. In `Basic`, set input `.h5ad` and output path, then click `Inspect`
3. In `Colors & Genes`, build color columns and gene lists with the searchable editors
4. Open `Advanced` only if needed, then export

Output is written to `<output_dir>/index.html`.

### Python API

```python
from karospace import load_spatial_data, export_to_html

dataset = load_spatial_data(
    "your_data.h5ad",
    groupby="sample_id",  # Column identifying each section
    metadata_value_order={
        "course": ["naive", "peak_I", "peak_II", "peak_III"],
    },
)

export_to_html(
    dataset,
    output_path="viewer.html",
    color="cell_type",           # Initial color column
    title="KaroSpace",
    min_panel_size=150,          # Min panel width (responsive autoscaling)
    spot_size="auto",            # Adaptive by section density (or set a fixed number)
    downsample=30000,            # Max cells per section
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
    gene_encoding="auto",        # "auto" | "dense" | "sparse"
    gene_storage="embedded",     # "embedded" | "sidecar"
    gene_aux_path=None,          # Optional manifest path; defaults to viewer.genes.json
    gene_sparse_zero_threshold=0.8,
    pack_arrays=True,
    use_hvgs=True,
    hvg_limit=20,
    marker_genes_groupby=None,
    marker_genes_top_n=30,
    neighbor_stats_groupby=["cell_type"],
    neighbor_stats_permutations=20,
    interaction_markers_groupby=None,
    cluster_de_groupby=["cell_type"],
    cluster_de_top_n=20,
    cluster_de_method="wilcoxon",
    section_rotations={
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
| `--min-panel-size` | Minimum panel width in pixels | `150` |
| `--spot-size` | Cell/spot size (`auto` or positive number) | `auto` |
| `--downsample` | Max cells per section | None |
| `--theme` | Color theme (`light` or `dark`) | `light` |
| `--title` | Page title | `KaroSpace` |
| `--gene-encoding` | Gene vector encoding (`auto`, `dense`, `sparse`) | `auto` |
| `--gene-storage` | Gene storage mode (`embedded`, `sidecar`) | `embedded` |
| `--gene-aux-path` | Path for the gene sidecar manifest JSON | auto |
| `--gene-sparse-zero-threshold` | Zero fraction threshold for `auto` sparse encoding | `0.8` |
| `--no-pack-arrays` | Disable base64 packing of large per-section arrays | off |
| `--pack-arrays-min-len` | Only pack arrays when section cell count ≥ this value | `1024` |
| `--neighbor-permutations` | Permutations for neighbor enrichment z-scores | `auto` |
| `--neighbor-stats-groupby` | Obs columns for neighbor composition stats (`auto` or comma-separated) | `auto` |
| `--marker-genes-groupby` | Obs columns to compute marker genes for | empty |
| `--interaction-markers-groupby` | Obs columns to compute contact-conditioned markers for | empty |
| `--cluster-de-groupby` | Categorical obs columns for pairwise cluster DE (`Insights → Compare`) | empty |
| `--cluster-de-top-n` | Top genes kept per cluster pair | `20` |
| `--cluster-de-method` | `scanpy.tl.rank_genes_groups` method | `wilcoxon` |
| `--cluster-de-layer` | AnnData layer for pairwise cluster DE | `normalized` |
| `--cluster-de-min-cells` | Minimum cells required in both clusters | `20` |
| `--section-rotations` | Comma-separated `section_id:angle` pairs | empty |
| `--gene-correlation-top-n` | Correlated genes shown per embedded gene in discovery panel | `10` |

## Data Requirements

- **`adata.obsm['spatial']`** — 2D coordinates for each cell (x, y)
- **`adata.obs[groupby]`** — Column identifying which section each cell belongs to
- **Categorical or numeric columns in `adata.obs`** — For coloring cells

### Optional metadata

- `course` — Experimental phase (e.g., `"naive"`, `"peak_I"`); sections are outlined by this column when present
- `region`, `condition`, `timepoint` — Used for filter chips

Control display order of metadata values and section ordering via `metadata_value_order`:

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

If `adata.obsp` contains a neighbor graph (`spatial_connectivities`, `connectivities`, `neighbors`, or `neighbor_graph`), KaroSpace exposes graph overlay and neighbor-hover controls.

To enable contact-conditioned interaction markers, pass `interaction_markers_groupby=[...]` (typically the same categorical column used for coloring).

To enable pairwise cluster DE in `Insights → Genes → Compare`, pass `cluster_de_groupby=[...]`.

## Examples

See [`examples/`](examples/) for complete dataset-specific export scripts.

## Viewer Controls

### Grid View
- **Click a section** — Open detailed modal view
- **Color dropdown** — Switch between annotation columns
- **Gene input** — Fuzzy search with keyboard navigation, recent genes, saved panels, and marker suggestions
- **Size slider** — Adjust spot size
- **Filter chips** — Filter sections by metadata
- **Legend items** — Toggle categories; spotlight one across grid and UMAP
- **Insights button** — Toggle the insights panel (`Summary`, `Compare`, `Genes`, `Neighborhood`, `Regions`)
- **Screenshot / Theme** — Download snapshot or toggle dark/light mode

### Modal View
- **Scroll** — Zoom in/out
- **Drag** — Pan; `Space` + drag to pan while Select or Annotate is active
- **Split button** — Open the A/B comparison panel; drag the divider line directly on the canvas
- **Magic Wand** — Draw lasso selection; opens `Genes in selection` ranked gene panel
- **Annotate** — Draw persistent polygons; export as JSON with cell indices
- **Pick type** — Click a cell to spotlight its category
- **Focused view** — Open the current selection as a filtered subview
- **Hide tools** — Collapse the toolbar for an unobstructed canvas
- **Escape or click outside** — Close modal

### UMAP View
- **UMAP button** — Toggle panel; dock to any corner
- **Magic Wand** — Lasso selection synced to spatial views
- **Scroll / drag** — Zoom and pan the embedding

## Deployment and Sharing

- **Single-file HTML** — Default `gene_storage="embedded"` writes one standalone file; works offline
- **Sidecar mode** — `gene_storage="sidecar"` writes HTML + `<name>.genes.json` + `<name>.genes/` shards; requires HTTP(S) serving
- **`.karospace` packages** — Self-contained ZIP bundles combining the viewer HTML and all sidecar assets; open via the hosted loader or a local `loader.html`. All processing happens in the browser — data never leaves the user's machine
- **Static hosting** — Both embedded and sidecar modes work on GitHub Pages, S3, or any lab intranet
- **Sharing** — Send the HTML directly (embedded), or share HTML + manifest + shard directory via a hosted location (sidecar)

### Package an existing sidecar into `.karospace`

```bash
# Short form (auto-detects sidecar paths from the HTML)
karospace package-sidecar BALO.html --output BALO.karospace

# Explicit form
karospace package-sidecar BALO.html \
  --output BALO.karospace \
  --gene-aux-path BALO.genes.json \
  --gene-shard-dir BALO.genes \
  --loader-output BALO.loader.html
```

This wraps existing sidecar assets into a `.karospace` archive without recomputing analytics or rewriting the viewer HTML.

### Integrate polygon annotations back into AnnData

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

## License

MIT License

## Author

Christoffer Mattsson Langseth — Karolinska Institutet
