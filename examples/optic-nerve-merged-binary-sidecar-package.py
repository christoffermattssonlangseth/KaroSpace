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
    "/Users/chrislangseth/Downloads/optic_nerve_merged.scanpy.companion.ready_with_polygons.h5ad",
)

PRIMARY_ANNOTATION = "leiden"
ADDITIONAL_ANNOTATIONS = [
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
FEATURE_MANIFEST_PATH = "optic-nerve-merged-binary.features.json"

if not Path(H5AD_PATH).exists():
    raise SystemExit(
        "optic_nerve_merged h5ad not found. Set OPTIC_NERVE_MERGED_H5AD_PATH before running "
        "examples/optic-nerve-merged-binary-sidecar-package.py."
    )

dataset = load_spatial_data(
    H5AD_PATH,
    section_key="karospace_polygon_labels",
    spatial_key="spatial",
    section_metadata=[
       'animal','condition','timepoint'
    ],
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available annotation columns: {dataset.obs_columns[:10]}...")

common_kwargs = dict(
    main_cell_annotation=PRIMARY_ANNOTATION,
    title="optic_nerve_merged",
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
    feature_value_encoding="uint8",
    feature_manifest_path=FEATURE_MANIFEST_PATH,
    feature_sidecar_shard_size=128,
    marker_gene_annotations=[
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
    neighbor_stats_annotations=[
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
    pseudobulk_de_annotations=[
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
