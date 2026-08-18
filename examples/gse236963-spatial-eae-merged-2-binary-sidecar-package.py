"""
Example usage of KaroSpace with binary sidecar and .karospace export targets.

This script is configured for the GSE236963 Spatial EAE Merged Data 2 companion-ready h5ad and writes:
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
    "GSE236963_SPATIAL_EAE_MERGED_2_H5AD_PATH",
    "/Users/chrislangseth/Downloads/GSE236963_spatialEaeMergedDataGeo_2.companion.ready.h5ad",
)

PRIMARY_ANNOTATION = "group"
ADDITIONAL_ANNOTATIONS = [
    "sample",
]
SIDECAR_OUTPUT = "gse236963-spatial-eae-merged-2-binary-sidecar.html"
PACKAGE_OUTPUT = "gse236963-spatial-eae-merged-2-binary.karospace"
FEATURE_MANIFEST_PATH = "gse236963-spatial-eae-merged-2-binary.features.json"

if not Path(H5AD_PATH).exists():
    raise SystemExit(
        "GSE236963 Spatial EAE Merged Data 2 h5ad not found. Set GSE236963_SPATIAL_EAE_MERGED_2_H5AD_PATH before running "
        "examples/gse236963-spatial-eae-merged-2-binary-sidecar-package.py."
    )

dataset = load_spatial_data(
    H5AD_PATH,
    section_key="sample",
    spatial_key="spatial",
    section_metadata=[
        "sample",
        "group",
    ],
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available annotation columns: {dataset.obs_columns[:10]}...")

common_kwargs = dict(
    main_cell_annotation=PRIMARY_ANNOTATION,
    title="GSE236963 Spatial EAE Merged Data 2",
    min_panel_size=120,
    spot_size="auto",
    downsample=10_000_000,
    outline_by=None,
    cell_annotations=ADDITIONAL_ANNOTATIONS,
    features=[],
    use_hvgs=False,
    hvg_limit=50,
    feature_storage="sidecar",
    feature_encoding="auto",
    feature_value_encoding="uint16",
    feature_manifest_path=FEATURE_MANIFEST_PATH,
    feature_sidecar_shard_size=16,
    marker_gene_annotations=[
        "group",
        "sample",
    ],
    marker_genes_top_n=30,
    neighbor_stats_annotations=[
        "group",
        "sample",
    ],
    neighbor_stats_permutations=0,
    neighbor_stats_seed=42,
    pseudobulk_de_annotations=[
        "group",
        "sample",
    ],
    pseudobulk_de_top_n=20,
    pseudobulk_de_method="t-test",
    pseudobulk_de_layer=None,
    pseudobulk_de_min_cells=20,
    interaction_marker_annotations=None,
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
print(f"  - feature manifest: {FEATURE_MANIFEST_PATH}")
print(f"  - shard directory: {Path(FEATURE_MANIFEST_PATH).with_suffix('')}")
print(f"Wrote packaged binary viewer: {PACKAGE_OUTPUT}")
print(f"  - local opener: {Path(PACKAGE_OUTPUT).with_suffix('.loader.html')}")
print("Share either route:")
print(f"  - local web server flow: {SIDECAR_OUTPUT} + {FEATURE_MANIFEST_PATH} + shard directory")
print(
    "  - no-install local package flow: "
    f"{PACKAGE_OUTPUT} + {Path(PACKAGE_OUTPUT).with_suffix('.loader.html')}"
)
