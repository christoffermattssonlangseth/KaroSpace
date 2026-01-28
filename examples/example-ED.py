"""
Example usage of KaroSpace.

This script demonstrates how to load Xenium spatial transcriptomics data
and export it to an interactive HTML viewer.
"""

from karospace import load_spatial_data, export_to_html

# Path to your h5ad file
# Update this path to point to your EAE/MANA data
H5AD_PATH = '/Volumes/processing2/erectile_dysfunction/data/adata/ED_5k_filtered_clustered_cytetype_cellcharter_metadata.h5ad'#'/Volumes/processing2/RRmap/data/EAE_proseg_clustered_louvain_leiden_all_sections_annotated_rotated_scVI_mana_embedding_clustered.h5ad'

# Load the dataset
# - groupby: column in adata.obs that identifies each section
# - spatial_key: key in adata.obsm containing coordinates (default: 'spatial')
dataset = load_spatial_data(
    H5AD_PATH,
    groupby="sample_name_updated",  # adjust to match your data
    # Choose which obs columns appear as filter chips in the viewer
    metadata_columns=[ 'ED', 'age',
       'score EHS ', 'status', 'Erection per op', 'Diabetes'],
    # metadata_max_columns=4,  # optional: limit number of metadata columns used
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available color columns: {dataset.obs_columns[:10]}...")  # first 10

# Choose gene source for expression:
# - True: use highly variable genes (if present, capped to 20)
# - False: use the explicit genes list below
USE_HVGS = True

# Export to HTML with full features
# For your 107 sections with course/region metadata:
export_to_html(
    dataset,
    output_path="erectile-dys-Göritz-lab.html",
    color='cytetype_annotation_leiden_3.5',  # Initial color (categorical)
    title="KaroSpace",
    min_panel_size=120,  # minimum panel width in pixels, grid auto-adjusts
    spot_size=1.5,  # smaller spots for dense data
    downsample=100000,  # limit cells per section to keep file manageable
    theme="light",  # or "dark"
    outline_by="course",  # metadata column for panel outline colors

    # Include additional color options for the dropdown
    additional_colors=[
       'leiden_2', 'leiden_2.5',
       'leiden_3', 'leiden_3.5', 'cytetype_annotation_leiden_3.5',
        'cell_class', 'cell_subclass','cluster_cellcharter_5',
       'cluster_cellcharter_10', 'cluster_cellcharter_15',
       'cluster_cellcharter_20', 'cluster_cellcharter_25',
       'cluster_cellcharter_30',
    ],

    # Pre-load specific genes for expression visualization
    # These will be available in the gene input field
    genes=[
        # Example marker genes - replace with your genes of interest
         # A1
    

    ],
    use_hvgs=USE_HVGS,
    hvg_limit=20,

    # Compute marker genes for these categorical color columns
    # (appears in the Color panel under "Marker genes")
    marker_genes_groupby=[
       'cytetype_annotation_leiden_3.5',
        'cell_class', 'cell_subclass',
       'cluster_cellcharter_30'
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
