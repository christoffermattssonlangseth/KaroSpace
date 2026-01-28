"""
Example usage of KaroSpace.

This script demonstrates how to load Xenium spatial transcriptomics data
and export it to an interactive HTML viewer.
"""

from karospace import load_spatial_data, export_to_html

# Path to your h5ad file
# Update this path to point to your EAE/MANA data
H5AD_PATH = '/Volumes/processing2/ST_BRICHOS/data/ST_BRICHOS_region_subcluster.h5ad'#'/Volumes/processing2/RRmap/data/EAE_proseg_clustered_louvain_leiden_all_sections_annotated_rotated_scVI_mana_embedding_clustered.h5ad'

# Load the dataset
# - groupby: column in adata.obs that identifies each section
# - spatial_key: key in adata.obsm containing coordinates (default: 'spatial')
dataset = load_spatial_data(
    H5AD_PATH,
    groupby="sample_id",  # adjust to match your data
    # Choose which obs columns appear as filter chips in the viewer
    metadata_columns=["treatment"],
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
    output_path="ST-BRICHOS.html",
    color="re_annotation_regions",  # Initial color (categorical)
    title="KaroSpace",
    min_panel_size=120,  # minimum panel width in pixels, grid auto-adjusts
    spot_size=1.5,  # smaller spots for dense data
    downsample=100000,  # limit cells per section to keep file manageable
    theme="light",  # or "dark"
    outline_by="course",  # metadata column for panel outline colors

    # Include additional color options for the dropdown
    additional_colors=[
       'leiden_0.5',
       'leiden_0.75', 'leiden_1', 'leiden_1.5', 'leiden_2', 'leiden_2.5',
       'region_annotation',
    ],

    # Pre-load specific genes for expression visualization
    # These will be available in the gene input field
    genes=[
        # Example marker genes - replace with your genes of interest
         # A1
    "H2-D1", "B2m", "C4b", "Gfap", "Serpina3n",

    # DAM / ARM
    "Apoe", "Axl", "Cd63", "Cd63-ps", "Cd9", "Ctsb", "Ctsd", "Ctsl", "Ctsz",
    "H2-K1", "Hexa", "Lgals3bp", "Lyz2", "Npc2", "Trem2", "Tyrobp",

    # PIGs
    "Arpc1b", "C1qa", "C1qb", "C1qc", "C4a", "Cts3", "Ctsa", "Ctss", "Ctsh",
    "Clu", "Csf1r", "Cx3cr1", "Cyba", "Fcer1g", "Fcgr3", "Fcrls", "Grn",
    "Gusb", "Gns", "Gpx4", "Gpx4-ps", "Hexb", "Igfbp5", "Itgb5", "Itm2b",
    "Laptm5", "Lgmn", "Ly86", "Man2b1", "Mpeg1", "Olfml3", "Plek", "Prdx6",
    "Rpl18a", "S100a6", "Vsir"
    ],
    use_hvgs=USE_HVGS,
    hvg_limit=500,

    # Compute marker genes for these categorical color columns
    # (appears in the Color panel under "Marker genes")
    marker_genes_groupby=[
        'leiden_0.5',
       'leiden_0.75', 'leiden_1', 'leiden_1.5', 'leiden_2', 'leiden_2.5',
       'region_annotation','re_annotation_regions'
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
