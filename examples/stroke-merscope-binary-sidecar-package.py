"""
Example usage of KaroSpace with binary sidecar and .karospace export targets.

This script is configured for the stroke_all_clustered MERSCOPE companion-ready
h5ad (Stroke_merscop_Fan_CG) and writes:
1. an unpacked binary sidecar viewer bundle
2. a packaged .karospace bundle with matching settings

What the h5ad exposes for coloring:
  - leiden     (16 clusters; composition_cell_type, has stored leiden_colors)
  - all 500 panel genes (streamed to the binary feature sidecar)
  - the experimental metadata we added (per-cell obs columns), used both as
    filters and as categorical/continuous color tracks:
        animal_id, line, model, model_norm, treatment,
        timepoint, timepoint_days, replicate, injured, condition

Sections are split by `sample` (33 sections). The companion analytics were
prepared on leiden + condition + treatment + timepoint + model + line + injured.
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
    "STROKE_MERSCOPE_H5AD_PATH",
    "/Volumes/T7/Stroke_merscop_Fan_CG/h5ad/stroke_all_clustered.companion.ready.h5ad",
)

# The experimental metadata we added to obs.
METADATA_COLUMNS = [
    "animal_id",
    "line",
    "model",
    "model_norm",
    "treatment",
    "timepoint",
    "timepoint_days",
    "replicate",
    "injured",
    "condition",
]

# Clusters drive the primary color; the metadata columns are also selectable.
PRIMARY_ANNOTATION = "leiden"
ADDITIONAL_ANNOTATIONS = ["leiden", *METADATA_COLUMNS]

# Categorical annotations used for analytics (per uns/karospace_companion).
ANNOTATION_GROUPBYS = [
    "leiden",
    "condition",
    "treatment",
    "timepoint",
    "model",
    "line",
    "injured",
]

SIDECAR_OUTPUT = "stroke-merscope-binary-sidecar.html"
PACKAGE_OUTPUT = "stroke-merscope-binary.karospace"
FEATURE_MANIFEST_PATH = "stroke-merscope-binary.features.json"

if not Path(H5AD_PATH).exists():
    raise SystemExit(
        "stroke_all_clustered companion-ready h5ad not found. Set STROKE_MERSCOPE_H5AD_PATH "
        "before running examples/stroke-merscope-binary-sidecar-package.py."
    )

dataset = load_spatial_data(
    H5AD_PATH,
    section_key="sample",
    spatial_key="spatial",
    section_metadata=METADATA_COLUMNS,
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available annotation columns: {dataset.obs_columns[:10]}...")

common_kwargs = dict(
    main_cell_annotation=PRIMARY_ANNOTATION,
    title="Stroke MERSCOPE (Fan/CG)",
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
    marker_gene_annotations=ANNOTATION_GROUPBYS,
    marker_genes_top_n=30,
    neighbor_stats_annotations=ANNOTATION_GROUPBYS,
    neighbor_stats_permutations=0,
    neighbor_stats_seed=42,
    pseudobulk_de_annotations=ANNOTATION_GROUPBYS,
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
