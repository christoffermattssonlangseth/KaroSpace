"""
Example usage of KaroSpace with sidecar-based feature loading.

This script demonstrates how to load Xenium spatial transcriptomics data
and export it to an interactive HTML viewer plus an auxiliary feature JSON file.
"""

import os
import sys
from pathlib import Path

# Prefer the local repo checkout for this example.
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
ENABLE_ANALYTICS = True

# Load the dataset
dataset = load_spatial_data(
    H5AD_PATH,
    section_key="sample_id",
    section_metadata=[""],
    metadata_value_order={
       
    },
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available annotation columns: {dataset.obs_columns[:10]}...")

export_to_html(
    dataset,
    output_path="mouseBrainXenium5k.html",
    main_cell_annotation=PRIMARY_CLUSTER,
    title="KaroSpace",
    min_panel_size=120,
    spot_size="auto",
    downsample=100000,
    outline_by=None,
    cell_annotations=[
        "CellCharter_5",
        "CellCharter_15",
        "CellCharter_20",
        "leiden_0.1",
        "leiden_0.5",
        "leiden_1.0",
    ],
    features=[
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
    feature_correlation_top_n=10,
    feature_storage="sidecar",
    feature_manifest_path="mouseBrainXenium5k.features.json",
    neighbor_stats_annotations=ANALYTICS_COLUMNS if ENABLE_ANALYTICS else None,
    neighbor_stats_permutations=0,
    neighbor_stats_seed=42,
    statistics_additional_annotations=ANALYTICS_COLUMNS if ENABLE_ANALYTICS else None,
    pseudobulk_embed_top_n_per_comparison=20,
    pseudobulk_counts_layer="normalized",
    pseudobulk_min_cells_per_pseudobulk=20,
    interaction_markers_top_targets=6,
    interaction_markers_top_features=15,
    interaction_markers_min_cells=30,
    interaction_markers_min_neighbors=1,
)

print("\nDone! Open mouseBrainXenium5k.html through a local web server.")
print("This export also writes mouseBrainXenium5k.features.json for lazy downstream feature loading.")
print("Example: python -m http.server 8765")
print("Then open: http://127.0.0.1:8765/mouseBrainXenium5k.html")
