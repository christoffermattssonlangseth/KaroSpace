"""
Example usage of KaroSpace.

This script demonstrates how to load Xenium spatial transcriptomics data
and export it to an interactive HTML viewer.
"""

from karospace import load_spatial_data, export_to_html

# Path to your h5ad file
# Update this path to point to your EAE/MANA data
H5AD_PATH = '/Volumes/processing2/RRmap/data/EAE_MANA_annotated_gmm_clust_with_scores.h5ad'#'/Volumes/processing2/RRmap/data/EAE_proseg_clustered_louvain_leiden_all_sections_annotated_rotated_scVI_mana_embedding_clustered.h5ad'

# Load the dataset
# - groupby: column in adata.obs that identifies each section
# - spatial_key: key in adata.obsm containing coordinates (default: 'spatial')
dataset = load_spatial_data(
    H5AD_PATH,
    groupby="sample_id",  # adjust to match your data
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
    output_path="eae_mana_viewer2.html",
    color="anno_L2",  # Initial color (categorical)
    title="KaroSpace",
    min_panel_size=120,  # minimum panel width in pixels, grid auto-adjusts
    spot_size=1.5,  # smaller spots for dense data
    downsample=100000,  # limit cells per section to keep file manageable
    theme="light",  # or "dark"

    # Include additional color options for the dropdown
    additional_colors=[
       'anno_L3', 'anno_L2', 'anno_L1','leiden_3.5','leiden_mana_1.0','leiden_mana_0.8','leiden_mana_0.3','leiden_mana_0.5','gmm_25_gauss'

        # Add any other obs columns you want to switch between
    ],

    # Pre-load specific genes for expression visualization
    # These will be available in the gene input field
    genes=[
        # Example marker genes - replace with your genes of interest
        "Arg1",
        "C3",
        "Cd74",
        "Cldn11",
        "Col1a2",
        "Ctss",
        "Ermn",
        "Foxp3",
        "Gfap",
        "Gpnmb",
        "Grn",
        "H2-Aa",
        "H2-Ab1",
        "H2-Eb1",
        "Igf2",
        "Klk6",
        "Mag",
        "Mbp",
        "Meg3",
        "Mki67",
        "Ptgds",
        "Serpina3n",
        "Serping1",
        "Slc47a1",
        "Snap25",
        "Vtn"
    ],
    use_hvgs=USE_HVGS,
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
