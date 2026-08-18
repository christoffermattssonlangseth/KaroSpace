"""
KaroSpace export for dataset_2 (companion-ready) — brain/astrocyte spatial data
with all leiden resolutions, all GMM-CellCharter niche clusterings, and the cell
annotations. No embedded images.

Colour options include every clustering (leiden 0.2/0.5/1, gmm_CC 15–35) and the
annotation columns (celltype_major, astro_subclass, cytetype annotations). Because
the file is companion-ready, per-cluster analytics for every clustering are cheap
precomputed lookups.

Sections are split by sample_repl (sample × replicate). Writes a binary-sidecar
viewer bundle and a matching packaged .karospace.

Usage:
    python examples/dataset_2_karospace.py
    DATASET2_H5AD=/path/to/file.h5ad python examples/dataset_2_karospace.py
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
    "DATASET2_H5AD",
    "/Users/chrislangseth/Downloads/dataset_2_karospace.companion.ready.h5ad",
)

LEIDEN = ["leiden_0.2", "leiden_0.5", "leiden_1"]
GMM = ["gmm_CC_15", "gmm_CC_20", "gmm_CC_25", "gmm_CC_30", "gmm_CC_35"]
ANNOTATIONS = [
    "celltype_major",
    "astro_subclass",
    "cytetype_annotation_leiden_0.2",
    "cytetype_cellState_leiden_0.2",
    "cytetype_leiden_0.2_astrosub",
]

PRIMARY_ANNOTATION = "celltype_major"
ADDITIONAL_ANNOTATIONS = ANNOTATIONS + LEIDEN + GMM + ["total_counts", "n_genes_by_counts"]
# Run per-cluster analytics on every clustering + the main annotations.
CLUSTER_COLUMNS = ["celltype_major", "astro_subclass", "cytetype_annotation_leiden_0.2"] + LEIDEN + GMM

SIDECAR_OUTPUT = "dataset_2_karospace.html"
PACKAGE_OUTPUT = "dataset_2_karospace.karospace"
FEATURE_MANIFEST_PATH = "dataset_2_karospace.features.json"


def main() -> None:
    if not Path(H5AD_PATH).exists():
        raise SystemExit(
            f"H5AD not found: {H5AD_PATH}\nSet DATASET2_H5AD or edit the path above."
        )

    print("Loading spatial data (10.6 GB — this can take a while)...")
    dataset = load_spatial_data(
        H5AD_PATH,
        section_key="sample_repl",
        spatial_key="spatial",
        section_metadata=["sample_id", "sample_repl", "replicate", "batch"],
    )
    print(f"  Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")

    common_kwargs = dict(
        main_cell_annotation=PRIMARY_ANNOTATION,
        title="Dataset 2 — leiden / GMM-CellCharter / annotations",
        min_panel_size=120,
        spot_size="auto",
        downsample=10_000_000,
        outline_by=None,
        cell_annotations=ADDITIONAL_ANNOTATIONS,
        features=[],
        use_hvgs=False,
        feature_storage="sidecar",
        feature_encoding="auto",
        feature_value_encoding="uint8",
        feature_manifest_path=FEATURE_MANIFEST_PATH,
        feature_sidecar_shard_size=128,
        marker_gene_annotations=CLUSTER_COLUMNS,
        marker_genes_top_n=30,
        neighbor_stats_annotations=CLUSTER_COLUMNS,
        neighbor_stats_permutations=0,
        neighbor_stats_seed=42,
        pseudobulk_de_annotations=CLUSTER_COLUMNS,
        pseudobulk_de_top_n=20,
        pseudobulk_de_method="t-test",
        pseudobulk_de_layer=None,
        pseudobulk_de_min_cells=20,
        interaction_marker_annotations=None,
    )

    print("Exporting binary-sidecar viewer...")
    export_to_html(dataset, output_path=SIDECAR_OUTPUT, **common_kwargs)

    print("Packaging .karospace archive...")
    export_to_html(dataset, output_path=PACKAGE_OUTPUT, **common_kwargs)

    print(f"\nDone!\n  sidecar: {SIDECAR_OUTPUT} (+ {FEATURE_MANIFEST_PATH} + shard dir)")
    print(f"  package: {PACKAGE_OUTPUT} (+ {Path(PACKAGE_OUTPUT).with_suffix('.loader.html')})")
    print("Coloured by celltype_major; switch between all leiden / GMM-CellCharter / annotations in the colour selector.")


if __name__ == "__main__":
    main()
