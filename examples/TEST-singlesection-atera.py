"""
TEST RUN (2026-09-15) - SINGLE-SECTION case.

Atera human breast cancer companion-ready h5ad (a single-section dataset).
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
    "ATERA_HUMAN_BREAST_CANCER_H5AD_PATH",
    "/Users/chrislangseth/Downloads/atera_human-breast_cancer.companion.ready.h5ad",
)

OUT_DIR = REPO_ROOT / "examples" / "_testruns"
OUT_DIR.mkdir(parents=True, exist_ok=True)
BASE_NAME = "TEST-singlesection-atera"
SIDECAR_OUTPUT = str(OUT_DIR / f"{BASE_NAME}-sidecar.html")
PACKAGE_OUTPUT = str(OUT_DIR / f"{BASE_NAME}.karospace")
SIDECAR_MANIFEST = str(OUT_DIR / f"{BASE_NAME}.features.json")
PACKAGE_MANIFEST = f"{BASE_NAME}.features.json"

if not Path(H5AD_PATH).exists():
    raise SystemExit(
        f"Atera human breast cancer h5ad not found at {H5AD_PATH}. "
        "Set ATERA_HUMAN_BREAST_CANCER_H5AD_PATH."
    )

print(f"Loading Atera dataset: {H5AD_PATH}")
dataset = load_spatial_data(
    H5AD_PATH,
    section_key="section_id",
    spatial_key="spatial",
    section_metadata=["section_id"],
)

print(f"Loaded {dataset.n_sections} section(s) with {dataset.n_cells:,} total cells")

ANNOTATIONS = ["CellCharter_5", "CellCharter_10", "CellCharter_15", "CellCharter_20"]

common_kwargs = dict(
    main_cell_annotation="leiden",
    title="TEST - Atera Human Breast Cancer (single section)",
    min_panel_size=120,
    spot_size="auto",
    downsample=None,
    cell_annotations=ANNOTATIONS,
    features=[],
    feature_storage="sidecar",
    feature_encoding="auto",
    feature_value_encoding="uint8",
    feature_sidecar_shard_size=128,
    statistics_additional_annotations=["leiden"] + ANNOTATIONS,
    neighbor_stats_annotations=["leiden"] + ANNOTATIONS,
    neighbor_stats_permutations=0,
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
