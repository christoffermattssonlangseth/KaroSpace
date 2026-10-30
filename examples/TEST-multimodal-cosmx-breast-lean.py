"""
TEST RUN (2026-09-16) - MULTIMODAL case, LEAN analytics.

Same CosMx multiomic breast dataset and multimodal viewer (RNA + protein +
protein_channel) as TEST-multimodal-cosmx-breast.py, but with the statistics
recompute minimized so it fits in memory on this machine:
  - Wilcoxon runs on the MAIN annotation only (no extra annotations)
  - statistics restricted to the RNA modality
  - neighbor-stat permutations off

This exists because the full-config run OOMs during the Wilcoxon recompute of
the extra annotation across modalities. The proper fix is to reuse precomputed
analytics from the .companion.ready file instead of recomputing (tracked
separately). Use this to exercise the multimodal viewer UX now.

Writes into examples/_testruns/ so no existing example outputs are overwritten.
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
    "COSMX_BREAST_COMPANION_H5AD_PATH",
    "../KaroSpaceDataWrangling/notebooks/data/cosmx_multiomic_breast/cosmx_multiomic_breast.companion.ready.h5ad",
)

OUT_DIR = REPO_ROOT / "examples" / "_testruns"
OUT_DIR.mkdir(parents=True, exist_ok=True)
BASE_NAME = "TEST-multimodal-cosmx-breast-lean"
SIDECAR_OUTPUT = str(OUT_DIR / f"{BASE_NAME}-sidecar.html")
PACKAGE_OUTPUT = str(OUT_DIR / f"{BASE_NAME}.karospace")
SIDECAR_MANIFEST = str(OUT_DIR / f"{BASE_NAME}.features.json")
PACKAGE_MANIFEST = f"{BASE_NAME}.features.json"

if not Path(H5AD_PATH).exists():
    raise SystemExit(
        f"CosMx multiomic companion h5ad not found at {H5AD_PATH}. "
        "Set COSMX_BREAST_COMPANION_H5AD_PATH."
    )

print(f"Loading multimodal CosMx dataset: {H5AD_PATH}")
dataset = load_spatial_data(
    H5AD_PATH,
    section_key="sample",
    spatial_key="spatial",
    section_metadata=["sample", "tissue", "Run_Tissue_name", "Panel"],
)

print(f"Loaded {dataset.n_sections} section(s) with {dataset.n_cells:,} total cells")
print(f"Detected modalities: {list(dataset.modalities.keys())}")

common_kwargs = dict(
    main_cell_annotation="leiden_rna",
    title="TEST - CosMx Multiomic Breast (multimodal, lean stats)",
    min_panel_size=150,
    spot_size="auto",
    downsample=None,
    # leiden_protein stays viewable, but is NOT added to statistics (that extra
    # Wilcoxon pass across modalities is what OOMs).
    cell_annotations=["leiden_protein"],
    features=[],
    feature_storage="sidecar",
    feature_encoding="auto",
    feature_value_encoding="uint8",
    feature_sidecar_shard_size=128,
    # Lean statistics: main annotation only, RNA modality only, no permutations.
    statistics_additional_annotations=[],
    statistics_modalities=["rna"],
    neighbor_stats_annotations=["leiden_rna"],
    neighbor_stats_permutations=0,
    # Viewer stays fully multimodal.
    modalities=["rna", "protein", "protein_channel"],
)

print("\nExporting sidecar bundle...")
export_to_html(dataset, output_path=SIDECAR_OUTPUT, feature_manifest_path=SIDECAR_MANIFEST, **common_kwargs)
print("Exporting .karospace package...")
export_to_html(dataset, output_path=PACKAGE_OUTPUT, feature_manifest_path=PACKAGE_MANIFEST, **common_kwargs)

print("\nDone.")
print(f"  sidecar  : {SIDECAR_OUTPUT}")
print(f"  manifest : {SIDECAR_MANIFEST}")
print(f"  package  : {PACKAGE_OUTPUT}")
print(f"  loader   : {Path(PACKAGE_OUTPUT).with_suffix('.loader.html')}")
