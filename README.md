# KaroSpace

<p align="center">
  <img src="assets/logo.png" alt="KaroSpace Logo" width="300">
</p>

**KaroSpace** is an interactive HTML viewer for exploring spatial transcriptomics data. It generates standalone HTML files from h5ad files that can be shared and viewed in any web browser without requiring a server or Python installation.

Originally developed at Karolinska Institutet for visualizing Xenium spatial transcriptomics data across multiple tissue sections.

## Features

- **Multi-section grid view** - Display dozens or hundreds of tissue sections in a responsive grid layout
- **Interactive zoom and pan** - Click any section to open a detailed view with mouse wheel zoom and drag-to-pan
- **UMAP view with Magic Wand selection** - Toggle UMAP panel, draw to select cells, highlights sync across views
- **Category toggling** - Click legend items to show/hide specific cell types or clusters; hidden cells appear as grey
- **Gene expression visualization** - Pre-load genes of interest and switch between them with a viridis colormap
- **Multiple color columns** - Switch between different annotation columns (e.g., cell types, clusters, conditions)
- **Metadata filtering** - Filter sections by metadata like experimental condition, timepoint, or region
- **Cell tooltips** - Hover over cells to see their type or expression value
- **Course-based borders** - Section panels are outlined with colors indicating their experimental course/condition
- **Dark/light mode** - Toggle between themes with preference saved to browser localStorage
- **Adjustable spot size** - Control cell/spot size in both grid view and detailed modal view
- **Standalone HTML** - Generated files are self-contained with embedded data and JavaScript

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

### Python API

```python
from karospace import load_spatial_data, export_to_html

# Load your h5ad file
dataset = load_spatial_data(
    "your_data.h5ad",
    groupby="sample_id",  # Column identifying each section
)

# Export to HTML
export_to_html(
    dataset,
    output_path="viewer.html",
    color="cell_type",           # Initial color column
    title="KaroSpace",
    min_panel_size=150,          # Min panel width (grid auto-adjusts)
    spot_size=2.0,               # Cell/spot size
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
)
```

### Command Line

```bash
karospace your_data.h5ad -o viewer.html --color leiden --cols 6
```

#### CLI Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --output` | Output HTML file path | `karospace.html` |
| `-c, --color` | Initial color column | `leiden` |
| `-g, --groupby` | Column to group sections by | `sample_id` |
| `--min-panel-size` | Minimum panel width in pixels (grid auto-adjusts) | `150` |
| `--spot-size` | Cell/spot size | `2.0` |
| `--downsample` | Max cells per section | None |
| `--theme` | Color theme (`light` or `dark`) | `light` |
| `--title` | Page title | `KaroSpace` |

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

## Example

See [example.py](example.py) for a complete working example.

```python
from karospace import load_spatial_data, export_to_html

dataset = load_spatial_data(
    "your_data.h5ad",
    groupby="sample_id",
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} cells")

export_to_html(
    dataset,
    output_path="viewer.html",
    color="anno_L2",
    title="KaroSpace",
    min_panel_size=120,
    spot_size=1.5,
    downsample=30000,
    additional_colors=['anno_L3', 'anno_L2', 'anno_L1', 'leiden'],
    genes=["Cd4", "Cd8a", "Gfap", "Mbp"],
)
```

## Viewer Controls

### Grid View
- **Click a section** - Open detailed modal view
- **Color dropdown** - Switch between annotation columns
- **Gene input** - Type a gene name to view expression (must be pre-loaded)
- **Size slider** - Adjust spot size
- **Filter chips** - Click to filter sections by metadata
- **Legend items** - Click to toggle categories on/off
- **Theme button** - Toggle dark/light mode

### Modal View (Detailed Section)
- **Mouse wheel** - Zoom in/out
- **Click and drag** - Pan around
- **Zoom buttons** - +/- zoom controls
- **Reset button** - Return to default zoom/pan
- **Size slider** - Adjust spot size for this view
- **Hover over cells** - See cell type or expression value
- **Escape or click outside** - Close modal

### UMAP View (if available)
If your h5ad file contains UMAP coordinates (`adata.obsm['X_umap']`), a UMAP toggle button appears:
- **UMAP button** - Toggle the UMAP panel on/off
- **Magic Wand** - Activate lasso selection mode
- **Draw selection** - Click and drag to draw a selection area
- **Clear** - Clear the current selection
- **Mouse wheel** - Zoom the UMAP view
- **Click and drag** (without Magic Wand) - Pan the UMAP view

Selected cells are highlighted with a yellow/gold outline in both UMAP and spatial views.

## Performance Tips

- Use `downsample` parameter for datasets with >50,000 cells per section
- Limit `genes` list to only essential genes (each adds to file size)
- Consider splitting very large datasets into multiple viewers

## License

MIT License

## Author

Christoffer Mattsson Langseth - Karolinska Institutet
