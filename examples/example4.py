"""
Example usage of KaroSpace.

This script demonstrates how to load Xenium spatial transcriptomics data
and export it to an interactive HTML viewer.
"""

from karospace import load_spatial_data, export_to_html

# Path to your h5ad file
# Update this path to point to your EAE/MANA data
H5AD_PATH = '/Volumes/processing2/RRmap/data/cell_paper_mana.h5ad'#'/Volumes/processing2/RRmap/data/EAE_proseg_clustered_louvain_leiden_all_sections_annotated_rotated_scVI_mana_embedding_clustered.h5ad'

# Load the dataset
# - groupby: column in adata.obs that identifies each section
# - spatial_key: key in adata.obsm containing coordinates (default: 'spatial')
dataset = load_spatial_data(
    H5AD_PATH,
    groupby="sample_id",  # adjust to match your data
    # Choose which obs columns appear as filter chips in the viewer
    metadata_columns=["time_type", "region",'SCORE','SEX'],
    # metadata_max_columns=4,  # optional: limit number of metadata columns used
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available color columns: {dataset.obs_columns[:10]}...")  # first 10

# Choose gene source for expression:
# - True: use highly variable genes (if present, capped to 20)
# - False: use the explicit genes list below
USE_HVGS = False

# Export to HTML with full features
# For your 107 sections with course/region metadata:
export_to_html(
    dataset,
    output_path="CELL-PAPER.html",
    color='Annotation 3 (medium with DA)',  # Initial color (categorical)
    title="KaroSpace",
    min_panel_size=120,  # minimum panel width in pixels, grid auto-adjusts
    spot_size=1.5,  # smaller spots for dense data
    downsample=100000,  # limit cells per section to keep file manageable
    theme="light",  # or "dark"
    outline_by="course",  # metadata column for panel outline colors

    # Include additional color options for the dropdown
    additional_colors=[
       'NON_DA', 'DA', 'Annotation 1 (broad)',
       'Annotation 2 (medium)',
       'Annotation 4 (high)', 'Annotation 3 (medium with DA)_colors',
       'compartment', 'Annotation 3.5','leiden_mana_0.3', 'leiden_mana_0.5',
       'leiden_mana_0.8', 'leiden_mana_1.0', 'leiden_mana_1.5',
       'leiden_mana_2'
    ],

    # Pre-load specific genes for expression visualization
    # These will be available in the gene input field
    genes=[
        # Example marker genes - replace with your genes of interest
         # A1
    "H2-D1", "B2m", "C4b", "Gfap", "Serpina3n","Cd74"

    ],
    use_hvgs=USE_HVGS,
    hvg_limit=500,

    # Compute marker genes for these categorical color columns
    # (appears in the Color panel under "Marker genes")
    marker_genes_groupby=[
       'Annotation 3 (medium with DA)','leiden_mana_1.5',
       'leiden_mana_2','compartment'
    ],
    marker_genes_top_n=50,
)

# The viewer now supports:
# 1. Filter by course (peak_I, peak_II, peak_III) or other metadata
# 2. Switch between different color columns
# 3. View gene expression for pre-loaded genes
# 4. Click to expand sections with zoom/pan
# 5. Toggle categories on/off in the legend

print("\nDone! Open eae_mana_viewer.html in a browser.")
print("Use the filter chips to show only specific courses (e.g., peak_III)")
print("Use the Color dropdown to switch between different annotations")
print("Type a gene name to view expression (must be in the genes list)")
