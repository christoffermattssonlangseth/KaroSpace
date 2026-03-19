"""
Example usage of KaroSpace with sidecar-based gene loading.

This script demonstrates how to load the Xenium mouse pup dataset
and export it to an interactive HTML viewer plus an auxiliary gene JSON file.
"""

import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("KMP_WARNINGS", "0")

from karospace import load_spatial_data, export_to_html

H5AD_PATH = os.environ.get(
    "XENIUM_PUP_H5AD_PATH",
    "/tmp/xenium_mouse_embryo.companion.ready.h5ad",
)

if H5AD_PATH.startswith("/path/to/"):
    raise SystemExit(
        "Set XENIUM_PUP_H5AD_PATH to your .h5ad file before running "
        "examples/xenium-pup-sidecar.py."
    )

PRIMARY_CLUSTER = "leiden_0.5"
ANALYTICS_COLUMNS = [PRIMARY_CLUSTER, "leiden_0.1", "leiden_1", "leiden_1.5", "leiden_2"]
additional_colors=ANALYTICS_COLUMNS[1:],
marker_genes_groupby=ANALYTICS_COLUMNS,
cluster_de_groupby=ANALYTICS_COLUMNS,
neighbor_stats_groupby=ANALYTICS_COLUMNS,
OUTPUT_PATH = "xenium-mouse-pup-sidecar.html"
GENE_AUX_PATH = "xenium-mouse-pup-sidecar.genes.json"

dataset = load_spatial_data(
    H5AD_PATH,
    groupby="sample_id",
    metadata_columns=[],
    metadata_value_order={
        "condition": [],
    },
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available color columns: {dataset.obs_columns[:10]}...")

export_to_html(
      dataset,
      output_path=OUTPUT_PATH,
      color=PRIMARY_CLUSTER,
      title="KaroSpace",
      min_panel_size=120,
      spot_size="auto",
      downsample=10_000_000,
      theme="light",
      outline_by=None,
      additional_colors=ANALYTICS_COLUMNS[1:],
      genes=[],
      use_hvgs=False,
      hvg_limit=50,
      gene_storage="sidecar",
      gene_aux_path=GENE_AUX_PATH,
      gene_sidecar_shard_size=8,
      marker_genes_groupby=ANALYTICS_COLUMNS,
      marker_genes_top_n=50,
      cluster_de_groupby=ANALYTICS_COLUMNS,
      cluster_de_top_n=20,
      cluster_de_method="t-test",
      cluster_de_layer="normalized",
      cluster_de_min_cells=20,
      neighbor_stats_groupby=ANALYTICS_COLUMNS,
      neighbor_stats_permutations=0,
      neighbor_stats_seed=42,
      interaction_markers_groupby=None,
  )

print(f"\nDone! Open {OUTPUT_PATH} through a local web server.")
print(f"This export also writes {GENE_AUX_PATH} for lazy downstream gene loading.")
print("Example: python -m http.server 8765")
print(f"Then open: http://127.0.0.1:8765/{OUTPUT_PATH}")
