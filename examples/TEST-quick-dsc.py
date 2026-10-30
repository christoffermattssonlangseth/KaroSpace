"""
TEST RUN (2026-09-15) - fast sanity check of the current export pipeline.

Uses the small D_SC companion-ready h5ad (~137 MB) to validate that the
post-merge API (statistics_* terminology, wilcoxon defaults) exports cleanly.
Writes everything into examples/_testruns/ so no existing example outputs are
overwritten.
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
    "D_SC_H5AD_PATH",
    "/Users/chrislangseth/Downloads/D_SC.scanpy.companion.ready.h5ad",
)

OUT_DIR = REPO_ROOT / "examples" / "_testruns"
OUT_DIR.mkdir(parents=True, exist_ok=True)
BASE_NAME = "TEST-quick-dsc"
SIDECAR_OUTPUT = str(OUT_DIR / f"{BASE_NAME}-sidecar.html")
PACKAGE_OUTPUT = str(OUT_DIR / f"{BASE_NAME}.karospace")
# Sidecar export takes a full path; the .karospace export requires a bare filename
# (the manifest is packed inside the bundle).
SIDECAR_MANIFEST = str(OUT_DIR / f"{BASE_NAME}.features.json")
PACKAGE_MANIFEST = f"{BASE_NAME}.features.json"

if not Path(H5AD_PATH).exists():
    raise SystemExit(f"D_SC h5ad not found at {H5AD_PATH}. Set D_SC_H5AD_PATH.")

dataset = load_spatial_data(
    H5AD_PATH,
    section_key="sample_id",
    spatial_key="spatial",
)

print(f"Loaded {dataset.n_sections} section(s) with {dataset.n_cells:,} total cells")

ANNOTATIONS = ["leiden_0_4", "leiden_1_0", "CellCharter_10"]

common_kwargs = dict(
    main_cell_annotation="leiden",
    title="TEST - D_SC (quick sanity)",
    min_panel_size=120,
    spot_size="auto",
    downsample=None,
    cell_annotations=ANNOTATIONS,
    features=[],
    feature_storage="sidecar",
    feature_encoding="auto",
    feature_value_encoding="uint16",
    feature_sidecar_shard_size=64,
    statistics_additional_annotations=ANNOTATIONS,
    neighbor_stats_annotations=ANNOTATIONS,
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
