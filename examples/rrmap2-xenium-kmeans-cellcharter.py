"""
KaroSpace export for the RRMAP2 Xenium "all samples" kmeans-separated object,
annotated + filtered + processed with CellCharter (companion-ready).

Colours by ALL leiden resolutions and ALL CellCharter resolutions; because the
file is companion-ready, the per-cluster analytics (marker genes / DE / neighbor
enrichment) for every clustering are cheap precomputed lookups.

~1.42M cells across 54 samples — writes a binary-sidecar viewer bundle and a
matching packaged .karospace.

Usage:
    python examples/rrmap2-xenium-kmeans-cellcharter.py
    RRMAP2_KMEANS_CELLCHARTER_H5AD=/path/to/file.h5ad python examples/rrmap2-xenium-kmeans-cellcharter.py
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
    "RRMAP2_KMEANS_CELLCHARTER_H5AD",
    "/Volumes/moldiassd/RRMAP2_xenium_adata/kmeans_separated/"
    "RRMAP2_xenium_all_samples.cellcharter.companion.ready.with_metadata.rerun.h5ad",
)

# Sample-level clinical/experimental metadata — exposed BOTH as section metadata
# and as colour options.
METADATA_COLS = ["stage", "condition", "region", "sex", "model"]

# All clusterings = every CellCharter resolution + every leiden resolution.
CELLCHARTER = [
    "CellCharter_5", "CellCharter_10", "CellCharter_15", "CellCharter_20",
    "CellCharter_25", "CellCharter_30", "CellCharter_35", "CellCharter_40",
    "CellCharter_45", "CellCharter_50",
]
LEIDEN = [
    "leiden_0.5", "leiden_1", "leiden_1.5", "leiden_2",
    "leiden_2.5", "leiden_3", "leiden_3.5", "leiden_4.0",
]
ALL_CLUSTERINGS = LEIDEN + CELLCHARTER

PRIMARY_ANNOTATION = "leiden_2.5"
ADDITIONAL_ANNOTATIONS = ALL_CLUSTERINGS + METADATA_COLS

SIDECAR_OUTPUT = "rrmap2-xenium-kmeans-cellcharter.html"
PACKAGE_OUTPUT = "rrmap2-xenium-kmeans-cellcharter.karospace"
FEATURE_MANIFEST_PATH = "rrmap2-xenium-kmeans-cellcharter.features.json"


def main() -> None:
    if not Path(H5AD_PATH).exists():
        raise SystemExit(
            f"H5AD not found: {H5AD_PATH}\nSet RRMAP2_KMEANS_CELLCHARTER_H5AD or edit the path above."
        )

    print("Loading spatial data (15 GB / ~1.4M cells — this can take a while)...")
    dataset = load_spatial_data(
        H5AD_PATH,
        section_key="kmeans_split_id",
        spatial_key="spatial",
        section_metadata=METADATA_COLS,
    )
    print(f"  Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")

    common_kwargs = dict(
        main_cell_annotation=PRIMARY_ANNOTATION,
        title="RRMAP2 Xenium All Samples — kmeans / CellCharter",
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
        # All clusterings get analytics (cheap precomputed companion lookups).
        marker_gene_annotations=ALL_CLUSTERINGS,
        marker_genes_top_n=30,
        neighbor_stats_annotations=ALL_CLUSTERINGS,
        neighbor_stats_permutations=0,
        neighbor_stats_seed=42,
        pseudobulk_de_annotations=ALL_CLUSTERINGS,
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
    print("Coloured by leiden_2.5; switch between all leiden / CellCharter resolutions in the colour selector.")


if __name__ == "__main__":
    main()
