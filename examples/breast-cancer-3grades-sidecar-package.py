"""
Example usage of KaroSpace with binary sidecar and .karospace export targets.

This script is configured for the breast cancer (3 grades) Illumina NovaSeq X
companion-ready h5ad and writes:
1. an unpacked binary sidecar viewer bundle
2. a packaged .karospace bundle with matching settings

What the h5ad exposes for coloring:
  - leiden               (21 transcriptional clusters)
  - cellcharter_domains  (26 spatial domains; the canonical CellCharter track)
  - all panel genes (streamed to the binary feature sidecar)
  - sample_id, which doubles as the per-section tumor/grade label:
        DCIS_IDC_Grade1, IDC_Grade2, IDC_Grade3

Sections are split by `sample_id` (3 tumors). The companion analytics were
prepared on leiden + cellcharter_domains, grouped by sample_id, so the export
reuses those precomputed tables (no recompute on the full ~33 GB matrix).

The file also carries extra CellCharter resolutions in obs
(cellcharter_k10/k12/k15/k18/k22/k26) if you want to swap the domain track.
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
    "BREAST_CANCER_H5AD_PATH",
    "/Volumes/processing/illumina_novaseqx_breast/breast_cancer_3grades.companion.ready.h5ad",
)

# Sections split on sample_id; that column also names the tumor / grade per
# section, so it is the metadata we surface as a filter chip + color track.
GROUPBY = "sample_id"
METADATA_COLUMNS = ["sample_id"]
METADATA_VALUE_ORDER = {
    "sample_id": ["DCIS_IDC_Grade1", "IDC_Grade2", "IDC_Grade3"],
}
METADATA_LABELS = {"sample_id": "Tumor / grade"}

# The two cluster/domain tracks the insights views run on.
PRIMARY_ANNOTATION = "leiden"
ADDITIONAL_ANNOTATIONS = ["leiden", "cellcharter_domains", "sample_id"]

# Categorical annotations used for analytics. These are exactly the columns the
# companion was prepared on, so the export reuses the precomputed marker genes /
# pseudobulk DE (instant, computed on the `normalized` layer) instead of recomputing
# on the raw-count X. `sample_id` is intentionally NOT here: it is the per-section
# tumor/grade label (metadata + color), not an insights DE section_key.
ANNOTATION_GROUPBYS = ["leiden", "cellcharter_domains"]

SIDECAR_OUTPUT = "breast-cancer-3grades-binary-sidecar.html"
PACKAGE_OUTPUT = "breast-cancer-3grades-binary.karospace"
FEATURE_MANIFEST_PATH = "breast-cancer-3grades-binary.features.json"

if not Path(H5AD_PATH).exists():
    raise SystemExit(
        "breast_cancer_3grades companion-ready h5ad not found. Set BREAST_CANCER_H5AD_PATH "
        "before running examples/breast-cancer-3grades-sidecar-package.py."
    )

dataset = load_spatial_data(
    H5AD_PATH,
    section_key=GROUPBY,
    spatial_key="spatial",
    section_metadata=METADATA_COLUMNS,
    metadata_value_order=METADATA_VALUE_ORDER,
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available annotation columns: {dataset.obs_columns[:10]}...")


def extract_he_images(dataset, out_dir="breast-cancer-3grades-he"):
    """Write each section's H&E hires image (uns['spatial'][sample]['images']
    ['he_hires']) to a PNG and return a {section_id: path} map for section_images.

    The viewer auto-fits each image over its section's data bounds at load; use
    the in-viewer Align panel to fine-tune, then Save Session to persist it.
    """
    from PIL import Image

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    spatial = dataset.adata.uns.get("spatial", {})
    images = {}
    for section_id in spatial:
        node = spatial[section_id]
        arr = node.get("images", {}).get("he_hires")
        if arr is None:
            print(f"  No he_hires image for section '{section_id}' — skipping.")
            continue
        png_path = out / f"he_{section_id}.png"
        Image.fromarray(arr).save(png_path)
        images[str(section_id)] = str(png_path)
        print(f"  Extracted H&E for '{section_id}': {arr.shape[1]}x{arr.shape[0]}px -> {png_path}")
    return images


SECTION_IMAGES = extract_he_images(dataset)

common_kwargs = dict(
    main_cell_annotation=PRIMARY_ANNOTATION,
    title="Breast Cancer (3 grades)",
    min_panel_size=120,
    spot_size="auto",
    downsample=10_000_000,
    outline_by=None,
    cell_annotations=ADDITIONAL_ANNOTATIONS,
    metadata_labels=METADATA_LABELS,
    section_images=SECTION_IMAGES,
    section_images_max_px=4096,
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
