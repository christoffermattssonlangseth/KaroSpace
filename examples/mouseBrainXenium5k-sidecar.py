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
H5AD_PATH = os.environ.get(
    "MOUSEBRAIN_H5AD_PATH",
    "/tmp/mouseBrain5k_cellcharter.companion.ready.h5ad",
)

if H5AD_PATH.startswith("/path/to/"):
    raise SystemExit(
        "Set MOUSEBRAIN_H5AD_PATH to your .h5ad file before running "
        "examples/mouseBrainXenium5k-sidecar.py."
    )

PRIMARY_CLUSTER = "CellCharter_10"
ANALYTICS_COLUMNS = [PRIMARY_CLUSTER, "CellCharter_5", "leiden_0.5"]
USE_HVGS = False
ENABLE_ANALYTICS = True

# Load the dataset
dataset = load_spatial_data(
    H5AD_PATH,
    groupby="sample_id",
    metadata_columns=[""],
    metadata_value_order={
       
    },
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available color columns: {dataset.obs_columns[:10]}...")

export_to_html(
    dataset,
    output_path="mouseBrainXenium5k.html",
    color=PRIMARY_CLUSTER,
    title="KaroSpace",
    min_panel_size=120,
    spot_size="auto",
    downsample=100000,
    theme="light",
    outline_by=None,
    additional_colors=[
        "CellCharter_5",
        "CellCharter_15",
        "CellCharter_20",
        "leiden_0.1",
        "leiden_0.5",
        "leiden_1.0",
    ],
    genes=[
        "Arg1",
        "Cd74",
        "Cldn11",
        "Col1a2",
        "Ctss",
        "Foxp3",
        "Gfap",
        "Gpnmb",
        "Grn",
        "H2-Aa",
        "H2-Ab1",
        "H2-Eb1",
        "Mbp",
        "Meg3",
        "Mki67",
        "Ptgds",
        "Serpina3n",
    ],
    use_hvgs=USE_HVGS,
    hvg_limit=200,
    gene_correlation_top_n=10,
    cluster_means_n_genes=500,
    gene_storage="sidecar",
    gene_aux_path="mouseBrainXenium5k.genes.json",
    marker_genes_groupby=ANALYTICS_COLUMNS if ENABLE_ANALYTICS else None,
    marker_genes_top_n=30,
    neighbor_stats_groupby=ANALYTICS_COLUMNS if ENABLE_ANALYTICS else None,
    neighbor_stats_permutations=0,
    neighbor_stats_seed=42,
    cluster_de_groupby=ANALYTICS_COLUMNS if ENABLE_ANALYTICS else None,
    cluster_de_top_n=20,
    cluster_de_method="t-test",
    cluster_de_layer="normalized",
    cluster_de_min_cells=20,
    interaction_markers_groupby=None,
    interaction_markers_top_targets=6,
    interaction_markers_top_genes=15,
    interaction_markers_min_cells=30,
    interaction_markers_min_neighbors=1,
)

print("\nDone! Open mouseBrainXenium5k.html through a local web server.")
print("This export also writes mouseBrainXenium5k.genes.json for lazy downstream gene loading.")
print("Example: python -m http.server 8765")
print("Then open: http://127.0.0.1:8765/mouseBrainXenium5k.html")
