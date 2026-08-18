"""
Example usage of KaroSpace with binary sidecar and .karospace export targets.

This script is configured for the Talbot Xenium Tumor CellCharter companion-ready h5ad and writes:
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
    "TALBOT_XENIUM_TUMOR_CELLCHARTER_H5AD_PATH",
    "/Users/chrislangseth/Downloads/talbot_xenium_tumor_annotated_updated_cellcharter.companion.ready.h5ad",
)

PRIMARY_ANNOTATION = "cytetype_annotation_leiden_4"
ADDITIONAL_ANNOTATIONS = [
      "CellCharter_20",
      "CellCharter_15",
      "condition",
      "genotype",
      "cytetype_cellState_leiden_4",
      "leiden_4"

]
SIDECAR_OUTPUT = "talbot-xenium-tumor-cellcharter-binary-sidecar.html"
PACKAGE_OUTPUT = "talbot-xenium-tumor-cellcharter-binary.karospace"
FEATURE_MANIFEST_PATH = "talbot-xenium-tumor-cellcharter-binary.features.json"

if not Path(H5AD_PATH).exists():
    raise SystemExit(
        "Talbot Xenium Tumor CellCharter h5ad not found. Set TALBOT_XENIUM_TUMOR_CELLCHARTER_H5AD_PATH before running "
        "examples/talbot-xenium-tumor-cellcharter-binary-sidecar-package.py."
    )

dataset = load_spatial_data(
    H5AD_PATH,
    section_key="sample_id",
    spatial_key="spatial",
    section_metadata=[
        "sample_id",
        "condition",
        "genotype",
    ],
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available annotation columns: {dataset.obs_columns[:10]}...")

common_kwargs = dict(
    main_cell_annotation=PRIMARY_ANNOTATION,
    title="Talbot Xenium Tumor CellCharter",
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
        "CellCharter_20",
      "CellCharter_15",
      "condition",
      "genotype",
      "cytetype_annotation_leiden_4",
      "cytetype_cellState_leiden_4",
      "leiden_4"
    ],
    marker_genes_top_n=30,
    neighbor_stats_annotations=[
        "CellCharter_20",
      "CellCharter_15",
      "condition",
      "genotype",
      "cytetype_annotation_leiden_4",
      "cytetype_cellState_leiden_4",
      "leiden_4"
    ],
    neighbor_stats_permutations=0,
    neighbor_stats_seed=42,
    pseudobulk_de_annotations=[
      "CellCharter_20",
      "CellCharter_15",
      "condition",
      "genotype",
      "cytetype_annotation_leiden_4",
      "cytetype_cellState_leiden_4",
      "leiden_4"
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
