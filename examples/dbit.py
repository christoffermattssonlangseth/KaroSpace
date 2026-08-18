"""
Example usage of KaroSpace.

This script demonstrates how to load Xenium spatial transcriptomics data
and export it to an interactive HTML viewer.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from karospace import load_spatial_data, export_to_html

# Path to your h5ad file
# Update this path to point to your EAE/MANA data
H5AD_PATH = '/Volumes/processing2/dbit-nature/data/adata_multi_pp-mana.h5ad'#'/Volumes/processing2/RRmap/data/EAE_proseg_clustered_louvain_leiden_all_sections_annotated_rotated_scVI_mana_embedding_clustered.h5ad'

# Load the dataset
# - section_key: column in adata.obs that identifies each section
# - spatial_key: key in adata.obsm containing coordinates (default: 'spatial')
dataset = load_spatial_data(
    H5AD_PATH,
    section_key="sample_id",  # adjust to match your data
    # Choose which obs columns appear as filter chips in the viewer
    section_metadata=['condition'],
    metadata_value_order={
    },
    # metadata_max_columns=4,  # optional: limit number of metadata columns used
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available annotation columns: {dataset.obs_columns[:10]}...")  # first 10

# Choose gene source for expression:
# - True: use highly variable genes (if present, capped to 20)
# - False: use the explicit genes list below
USE_HVGS = True
OUTLINE_BY = "condition"

# Export to HTML with full features
# For your 107 sections with course/region metadata:
export_to_html(
    dataset,
    output_path="DBIT.html",
    main_cell_annotation='leiden_0p5',  # Initial annotation (categorical)
    title="KaroSpace",
    min_panel_size=120,  # minimum panel width in pixels, grid auto-adjusts
    spot_size="auto",  # adaptive default based on section density
    downsample=100000,  # limit cells per section to keep file manageable
    outline_by=OUTLINE_BY,  # metadata column for panel outline colors

    # Include additional color options for the dropdown
    cell_annotations=[
       'leiden_1p0','leiden_1p5'
    ],

    # Pre-load specific genes for expression visualization
    # These will be available in the gene input field
    features=[
        # Example marker genes - replace with your genes of interest
        "Arg1",
        #"C3",
        "Cd74",
        "Cldn11",
        "Col1a2",
        "Ctss",
        #"Ermn",
        "Foxp3",
        "Gfap",
        "Gpnmb",
        "Grn",
        "H2-Aa",
        "H2-Ab1",
        "H2-Eb1",
        #"Igf2",
       # "Klk6",
        #"Mag",
        "Mbp",
        "Meg3",
        "Mki67",
        "Ptgds",
        "Serpina3n",
        #"Serping1",
       # "Slc47a1",
       # "Snap25",
       #"Vtn"
    ],
    use_hvgs=USE_HVGS,
    hvg_limit=200,

    # Compute marker genes for these categorical annotation columns
    # (appears in the Color panel under "Marker genes")
    marker_gene_annotations=[
           'leiden_0p5','leiden_1p0','leiden_1p5'
    ],
    marker_genes_top_n=50,
    # Force permutation z-scores (auto mode disables permutations for very large datasets).
    neighbor_stats_permutations=25,
    neighbor_stats_seed=42,
    # Contact-conditioned interaction markers (source near target vs source not near target).
    interaction_marker_annotations=None,
    interaction_markers_top_targets=6,
    interaction_markers_top_genes=15,
    interaction_markers_min_cells=30,
    interaction_markers_min_neighbors=1,
)

# The viewer now supports:
# 1. Filter by course (peak_I, peak_II, peak_III) or other metadata
# 2. Switch between different annotation columns
# 3. View gene expression for pre-loaded genes
# 4. Click to expand sections with zoom/pan
# 5. Toggle categories on/off in the legend

print("\nDone! Open DBIT.html in a browser.")
print("Use the filter chips to show only specific courses (e.g., peak_III)")
print("Use the Annotation selector to switch between different annotations")
print("Type a gene name to view expression (must be in the genes list)")
