"""
Example usage of KaroSpace with both sidecar and .karospace export targets.

This mirrors the Xenium mouse pup setup from examples/xenium-pup-sidecar.py, but writes:
1. an unpacked sidecar viewer bundle
2. a packaged .karospace bundle with matching settings
"""

import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("KMP_WARNINGS", "0")

from karospace import export_to_html, load_spatial_data

H5AD_PATH = os.environ.get(
    "XENIUM_PUP_H5AD_PATH",
    "/tmp/xenium_mouse_embryo.companion.ready.h5ad",
)

if H5AD_PATH.startswith("/path/to/"):
    raise SystemExit(
        "Set XENIUM_PUP_H5AD_PATH to your .h5ad file before running "
        "examples/xenium-pup-sidecar-package.py."
    )

PRIMARY_CLUSTER = "leiden_0.5"
ANALYTICS_COLUMNS = [PRIMARY_CLUSTER, "leiden_0.1", "leiden_1", "leiden_1.5", "leiden_2"]
SIDECAR_OUTPUT = "xenium-mouse-pup-sidecar.html"
PACKAGE_OUTPUT = "xenium-mouse-pup.karospace"
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

common_kwargs = dict(
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

export_to_html(
    dataset,
    output_path=SIDECAR_OUTPUT,
    **common_kwargs,
)

export_to_html(
    dataset,
    output_path=PACKAGE_OUTPUT,
    **common_kwargs,
)

print(f"\nDone! Wrote unpacked sidecar viewer: {SIDECAR_OUTPUT}")
print(f"  - gene manifest: {GENE_AUX_PATH}")
print(f"  - shard directory: {Path(GENE_AUX_PATH).with_suffix('')}")
print(f"Wrote packaged viewer: {PACKAGE_OUTPUT}")
print(f"  - local opener: {Path(PACKAGE_OUTPUT).with_suffix('.loader.html')}")
print("Share either route:")
print(f"  - local web server flow: {SIDECAR_OUTPUT} + {GENE_AUX_PATH} + shard directory")
print(f"  - no-install local package flow: {PACKAGE_OUTPUT} + {Path(PACKAGE_OUTPUT).with_suffix('.loader.html')}")
