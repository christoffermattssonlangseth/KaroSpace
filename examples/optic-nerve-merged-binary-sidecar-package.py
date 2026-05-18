"""
Example usage of KaroSpace with binary sidecar and .karospace export targets.

This script is configured for the optic_nerve_merged companion-ready h5ad and writes:
1. an unpacked binary sidecar viewer bundle
2. a packaged .karospace bundle with matching settings
"""

import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("KMP_WARNINGS", "0")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/karospace-mpl-cache")

from karospace import export_to_html, load_spatial_data

H5AD_PATH = os.environ.get(
    "OPTIC_NERVE_MERGED_H5AD_PATH",
    "/Users/chrislangseth/Downloads/optic_nerve_merged.scanpy.companion.ready.h5ad",
)

PRIMARY_COLOR = "leiden"
ADDITIONAL_COLORS = [
    "leiden",
    "leiden_0_2",
    "leiden_0_4",
    "leiden_0_6",
    "leiden_0_8",
    "leiden_1_0",
    "leiden_1_5",
    "leiden_2_0",
    "leiden_2_5",
    "leiden_3_0",
    "leiden_3_5",
    "leiden_4_0",
    "CellCharter_6",
    "CellCharter_8",
    "CellCharter_10",
    "CellCharter_12",
    "CellCharter_15",
    "CellCharter_20",
    "CellCharter_25",
    "CellCharter_30",
    
]
SIDECAR_OUTPUT = "optic-nerve-merged-binary-sidecar.html"
PACKAGE_OUTPUT = "optic-nerve-merged-binary.karospace"
GENE_AUX_PATH = "optic-nerve-merged-binary.genes.json"

if not Path(H5AD_PATH).exists():
    raise SystemExit(
        "optic_nerve_merged h5ad not found. Set OPTIC_NERVE_MERGED_H5AD_PATH before running "
        "examples/optic-nerve-merged-binary-sidecar-package.py."
    )

dataset = load_spatial_data(
    H5AD_PATH,
    groupby="sample_id",
    spatial_key="spatial",
    metadata_columns=[
       
    ],
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available color columns: {dataset.obs_columns[:10]}...")

common_kwargs = dict(
    color=PRIMARY_COLOR,
    title="optic_nerve_merged",
    min_panel_size=120,
    spot_size="auto",
    downsample=10_000_000,
    theme="light",
    outline_by=None,
    additional_colors=ADDITIONAL_COLORS,
    genes=[],
    use_hvgs=False,
    hvg_limit=50,
    gene_storage="sidecar",
    gene_sidecar_format="binary-v1",
    gene_encoding="auto",
    gene_value_encoding="uint8",
    gene_aux_path=GENE_AUX_PATH,
    gene_sidecar_shard_size=128,
    marker_genes_groupby=[
        "leiden",
    "leiden_0_2",
    "leiden_0_4",
    "leiden_0_6",
    "leiden_0_8",
    "leiden_1_0",
    "leiden_1_5",
    "leiden_2_0",
    "leiden_2_5",
    "leiden_3_0",
    "leiden_3_5",
    "leiden_4_0",
    "CellCharter_6",
    "CellCharter_8",
    "CellCharter_10",
    "CellCharter_12",
    "CellCharter_15",
    "CellCharter_20",
    "CellCharter_25",
    "CellCharter_30",
    ],
    marker_genes_top_n=30,
    neighbor_stats_groupby=[
       "leiden",
    "leiden_0_2",
    "leiden_0_4",
    "leiden_0_6",
    "leiden_0_8",
    "leiden_1_0",
    "leiden_1_5",
    "leiden_2_0",
    "leiden_2_5",
    "leiden_3_0",
    "leiden_3_5",
    "leiden_4_0",
    "CellCharter_6",
    "CellCharter_8",
    "CellCharter_10",
    "CellCharter_12",
    "CellCharter_15",
    "CellCharter_20",
    "CellCharter_25",
    "CellCharter_30",
    ],
    neighbor_stats_permutations=0,
    neighbor_stats_seed=42,
    cluster_de_groupby=[
        "leiden",
    "leiden_0_2",
    "leiden_0_4",
    "leiden_0_6",
    "leiden_0_8",
    "leiden_1_0",
    "leiden_1_5",
    "leiden_2_0",
    "leiden_2_5",
    "leiden_3_0",
    "leiden_3_5",
    "leiden_4_0",
    "CellCharter_6",
    "CellCharter_8",
    "CellCharter_10",
    "CellCharter_12",
    "CellCharter_15",
    "CellCharter_20",
    "CellCharter_25",
    "CellCharter_30",
    ],
    cluster_de_top_n=20,
    cluster_de_method="t-test",
    cluster_de_layer=None,
    cluster_de_min_cells=20,
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

print(f"\nDone! Wrote unpacked binary sidecar viewer: {SIDECAR_OUTPUT}")
print(f"  - gene manifest: {GENE_AUX_PATH}")
print(f"  - shard directory: {Path(GENE_AUX_PATH).with_suffix('')}")
print(f"Wrote packaged binary viewer: {PACKAGE_OUTPUT}")
print(f"  - local opener: {Path(PACKAGE_OUTPUT).with_suffix('.loader.html')}")
print("Share either route:")
print(f"  - local web server flow: {SIDECAR_OUTPUT} + {GENE_AUX_PATH} + shard directory")
print(
    "  - no-install local package flow: "
    f"{PACKAGE_OUTPUT} + {Path(PACKAGE_OUTPUT).with_suffix('.loader.html')}"
)
