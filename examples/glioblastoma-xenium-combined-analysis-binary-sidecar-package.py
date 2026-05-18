"""
Example usage of KaroSpace with binary sidecar and .karospace export targets.

This script is configured for the Glioblastoma Xenium Combined Analysis companion-ready h5ad and writes:
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
    "GLIOBLASTOMA_XENIUM_COMBINED_ANALYSIS_H5AD_PATH",
    "/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/data/processed/glioblastoma-xenium/glioblastoma-xenium-combined-analysis.companion.ready_with_polygons_sample_id_remapped.h5ad",
)

PRIMARY_COLOR = "kmeans_cluster"
ADDITIONAL_COLORS = [
    "sample_id_original",
    "section_id",
    "run_name",
]
SIDECAR_OUTPUT = "glioblastoma-xenium-combined-analysis-binary-sidecar.html"
PACKAGE_OUTPUT = "glioblastoma-xenium-combined-analysis-binary.karospace"
GENE_AUX_PATH = "glioblastoma-xenium-combined-analysis-binary.genes.json"

if not Path(H5AD_PATH).exists():
    raise SystemExit(
        "Glioblastoma Xenium Combined Analysis h5ad not found. Set GLIOBLASTOMA_XENIUM_COMBINED_ANALYSIS_H5AD_PATH before running "
        "examples/glioblastoma-xenium-combined-analysis-binary-sidecar-package.py."
    )

dataset = load_spatial_data(
    H5AD_PATH,
    groupby="sample_id",
    spatial_key="spatial",
    metadata_columns=[
        "sample_id",
        "sample_id_original",
        "section_id",
        "run_name",
        "panel_name",
        "preservation_method",
    ],
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available color columns: {dataset.obs_columns[:10]}...")

common_kwargs = dict(
    color=PRIMARY_COLOR,
    title="Glioblastoma Xenium Combined Analysis",
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
        "kmeans_cluster",
    ],
    marker_genes_top_n=30,
    neighbor_stats_groupby=[
        "kmeans_cluster",
    ],
    neighbor_stats_permutations=0,
    neighbor_stats_seed=42,
    cluster_de_groupby=[
        "kmeans_cluster",
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
