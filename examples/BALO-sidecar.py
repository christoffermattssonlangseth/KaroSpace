"""
Example usage of KaroSpace with sidecar-based gene loading.

This script demonstrates how to load Xenium spatial transcriptomics data
and export it to an interactive HTML viewer plus an auxiliary gene JSON file.
"""

import os
import sys
from pathlib import Path

# Prefer the local repo checkout over any older site-packages install.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Keep numba JIT enabled for normal performance.
# If your local environment has scanpy/numba import issues, set NUMBA_DISABLE_JIT=1 manually.
# Suppress Intel/OpenMP info messages during long compute phases.
os.environ.setdefault("KMP_WARNINGS", "0")

from karospace import load_spatial_data, export_to_html

# Path to your h5ad file
H5AD_PATH = os.environ.get("BALO_H5AD_PATH","/tmp/baloMS_companion.ready.h5ad")# "/Users/chrislangseth/Downloads/baloMS_indep_clust_balo_MANA_balo_annot.h5ad")

if H5AD_PATH.startswith("/path/to/"):
    raise SystemExit("Set BALO_H5AD_PATH to your .h5ad file before running examples/BALO-sidecar.py.")

# Load the dataset
dataset = load_spatial_data(
    H5AD_PATH,
    groupby="sample_id",
    metadata_columns=['condition'],
    metadata_value_order={
        "stage": [],
    },
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available color columns: {dataset.obs_columns[:10]}...")

# Choose gene source for expression:
# - True: use highly variable genes (if present, capped to hvg_limit)
# - False: use the explicit genes list below for embedded startup genes only
USE_HVGS = False
OUTLINE_BY = "condition"
ENABLE_ANALYTICS = True

export_to_html(
    dataset,
    output_path="BALO.html",
    color="leiden_2",
    title="KaroSpace",
    min_panel_size=120,
    spot_size="auto",
    downsample=100000,
    theme="light",
    outline_by=OUTLINE_BY,
    additional_colors=[
        'leiden_0.5',
    ],
    genes=[
     #   "Arg1",
     #   "Cd74",
     #   "Cldn11",
     #   "Col1a2",
     #   "Ctss",
     ##   "Foxp3",
     #   "Gfap",
     #   "Gpnmb",
     #   "Grn",
     #   "H2-Aa",
     #   "H2-Ab1",
     #   "H2-Eb1",
     #   "Mbp",
     #   "Meg3",
     #   "Mki67",
     #   "Ptgds",
     #   "Serpina3n",
    ],
    use_hvgs=USE_HVGS,
    hvg_limit=200,
    gene_storage="sidecar",
    gene_aux_path="BALO.genes.json",
    marker_genes_groupby=['leiden_0.5'] if ENABLE_ANALYTICS else None,
    marker_genes_top_n=50,
    neighbor_stats_permutations=25 if ENABLE_ANALYTICS else 0,
    cluster_de_groupby=["leiden_2"],
    cluster_de_top_n=20,
    cluster_de_method="t-test",
    cluster_de_layer="normalized",
    cluster_de_min_cells=20,
    neighbor_stats_seed=42,
    interaction_markers_groupby=None,
    interaction_markers_top_targets=6,
    interaction_markers_top_genes=15,
    interaction_markers_min_cells=30,
    interaction_markers_min_neighbors=1,
)

print("\nDone! Open BALO.html through a local web server.")
print("This export also writes BALO.genes.json for lazy downstream gene loading.")
print("Example: python -m http.server 8765")
print("Then open: http://127.0.0.1:8765/BALO.html")
