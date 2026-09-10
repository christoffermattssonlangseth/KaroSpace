"""
Example usage of KaroSpace with sidecar-based feature loading.

This script demonstrates how to load the Xenium mouse pup dataset
and export it to an interactive HTML viewer plus an auxiliary feature JSON file.
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
cell_annotations=ANALYTICS_COLUMNS[1:],
pseudobulk_additional_annotations=ANALYTICS_COLUMNS,
neighbor_stats_annotations=ANALYTICS_COLUMNS,
OUTPUT_PATH = "xenium-mouse-pup-sidecar.html"
FEATURE_MANIFEST_PATH = "xenium-mouse-pup-sidecar.features.json"

dataset = load_spatial_data(
    H5AD_PATH,
    section_key="sample_id",
    section_metadata=[],
    metadata_value_order={
        "condition": [],
    },
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available annotation columns: {dataset.obs_columns[:10]}...")

export_to_html(
      dataset,
      output_path=OUTPUT_PATH,
      main_cell_annotation=PRIMARY_CLUSTER,
      title="KaroSpace",
      min_panel_size=120,
      spot_size="auto",
      downsample=10_000_000,
      outline_by=None,
      cell_annotations=ANALYTICS_COLUMNS[1:],
      features=[],
      feature_storage="sidecar",
      feature_manifest_path=FEATURE_MANIFEST_PATH,
      feature_sidecar_shard_size=8,
      pseudobulk_additional_annotations=ANALYTICS_COLUMNS,
      pseudobulk_embed_top_n_per_comparison=20,
      pseudobulk_counts_layer="normalized",
      pseudobulk_min_cells_per_pseudobulk=20,
      neighbor_stats_annotations=ANALYTICS_COLUMNS,
      neighbor_stats_permutations=0,
      neighbor_stats_seed=42,
  )

print(f"\nDone! Open {OUTPUT_PATH} through a local web server.")
print(f"This export also writes {FEATURE_MANIFEST_PATH} for lazy downstream feature loading.")
print("Example: python -m http.server 8765")
print(f"Then open: http://127.0.0.1:8765/{OUTPUT_PATH}")
