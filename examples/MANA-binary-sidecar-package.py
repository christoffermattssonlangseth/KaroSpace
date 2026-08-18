"""
Example usage of KaroSpace with binary sidecar and .karospace export targets.

This script mirrors the RRMap / MANA setup from examples/MANA.py, but writes:
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

from karospace import export_to_html, load_spatial_data

H5AD_PATH = os.environ.get(
    "MANA_H5AD_PATH",
    '/tmp/rrmap.companion.ready.h5ad',  # "/Volumes/processing2/RRmap/data/RRmap_metadata_fixed_update.h5ad",
)

USE_HVGS = False
OUTLINE_BY = "stage"
PRIMARY_ANNOTATION = "anno_L2"
SIDECAR_OUTPUT = "RRMap-binary-sidecar.html"
PACKAGE_OUTPUT = "RRMap-binary.karospace"
FEATURE_MANIFEST_PATH = "RRMap-binary.features.json"

dataset = load_spatial_data(
    H5AD_PATH,
    section_key="sample_id",
    section_metadata=["stage", "condition", "day_of_sacrifice", "score_sacrifice", "region", "sex", "model"],
    metadata_value_order={
        "stage": [
            "MOG CFA",
            "PLP CFA",
            "NONSYMPTOM",
            "OS1",
            "ONSET1",
            "ONSET2",
            "MONOPHASIC",
            "PEAK1",
            "REMISSION1",
            "PEAK2",
            "REMISSION2",
            "PEAK3",
            "LONG",
        ],
    },
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available annotation columns: {dataset.obs_columns[:10]}...")

common_kwargs = dict(
    main_cell_annotation=PRIMARY_ANNOTATION,
    title="KaroSpace",
    min_panel_size=120,
    spot_size="auto",
    downsample=100000,
    outline_by=OUTLINE_BY,
    cell_annotations=[
        "anno_L3",
        "anno_L2",
        "anno_L1",
        "leiden_3.5",
        "compartment_mana",
        "stage",
        "condition",
        "day_of_sacrifice",
        "score_sacrifice",
        "region",
        "sex",
        "model",
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
    use_hvgs=USE_HVGS,
    hvg_limit=50,
    feature_storage="sidecar",
    feature_encoding="auto",
    feature_value_encoding="uint8",
    feature_manifest_path=FEATURE_MANIFEST_PATH,
    marker_gene_annotations=[
        "anno_L3",
        "anno_L2",
        "anno_L1",
        "leiden_3.5",
        "compartment_mana",
        "stage",
        "condition",
    ],
    marker_genes_top_n=50,
    neighbor_stats_permutations=25,
    neighbor_stats_seed=42,
    interaction_marker_annotations=None,
    interaction_markers_top_targets=6,
    interaction_markers_top_genes=15,
    interaction_markers_min_cells=30,
    interaction_markers_min_neighbors=1,
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
print(f"  - no-install local package flow: {PACKAGE_OUTPUT} + {Path(PACKAGE_OUTPUT).with_suffix('.loader.html')}")
