"""
Export a KaroSpace viewer for the stroke (dMCAO) dataset stroke_all_clustered —
3 samples, ~261k cells, split into one section per sample.

No embedded morphology images (no uns/spatial); it's a straight spatial cell map
coloured by the leiden clustering, with cell morphology / QC metrics available as
alternative continuous colourings (area, density, elongation, counts, ...).

Usage:
    python examples/stroke_all_clustered.py
    STROKE_H5AD=/path/to/file.h5ad python examples/stroke_all_clustered.py
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
    "STROKE_H5AD",
    "/Users/chrislangseth/Downloads/stroke_all_clustered.companion.ready.h5ad",
)

OUTPUT_DIR = REPO_ROOT
SIDECAR_OUTPUT = str(OUTPUT_DIR / "stroke_all_clustered.html")
PACKAGE_OUTPUT = str(OUTPUT_DIR / "stroke_all_clustered.karospace")
FEATURE_MANIFEST_PATH = str(OUTPUT_DIR / "stroke_all_clustered.features.json")

PRIMARY_ANNOTATION = "leiden"
ADDITIONAL_ANNOTATIONS = [
    "leiden",
    "total_counts",
    "n_genes_by_counts",
    "density",
    "area",
    "elongation",
    "avg_confidence",
]
CLUSTER_COLUMNS = ["leiden"]


def main() -> None:
    if not Path(H5AD_PATH).exists():
        raise SystemExit(
            f"H5AD not found: {H5AD_PATH}\nSet STROKE_H5AD or edit the path above."
        )

    print("Loading spatial data (large file — this can take a while)...")
    dataset = load_spatial_data(
        H5AD_PATH,
        section_key="sample",
        spatial_key="spatial",
        section_metadata=["sample"],
    )
    print(f"  {dataset.n_sections} section(s), {dataset.n_cells:,} cells")
    print(f"  section IDs: {[s.section_id for s in dataset.sections]}")

    common_kwargs = dict(
        main_cell_annotation=PRIMARY_ANNOTATION,
        title="Stroke (dMCAO) — leiden clusters",
        min_panel_size=140,
        spot_size="auto",
        outline_by=None,
        cell_annotations=ADDITIONAL_ANNOTATIONS,
        features=[],
        use_hvgs=False,
        feature_storage="sidecar",
        feature_encoding="auto",
        feature_value_encoding="uint8",
        feature_sidecar_shard_size=128,
        marker_gene_annotations=CLUSTER_COLUMNS,
        marker_genes_top_n=30,
        neighbor_stats_annotations=CLUSTER_COLUMNS,
        neighbor_stats_permutations=0,
        pseudobulk_de_annotations=CLUSTER_COLUMNS,
        pseudobulk_de_top_n=20,
        pseudobulk_de_method="t-test",
        pseudobulk_de_layer=None,
        interaction_marker_annotations=None,
    )

    # The sidecar writes its feature manifest next to the HTML (full path is fine);
    # the .karospace packager requires a bare filename (it lives inside the archive).
    print("Exporting sidecar HTML...")
    export_to_html(
        dataset, output_path=SIDECAR_OUTPUT, feature_manifest_path=FEATURE_MANIFEST_PATH, **common_kwargs
    )

    print("Packaging .karospace archive...")
    export_to_html(
        dataset,
        output_path=PACKAGE_OUTPUT,
        feature_manifest_path=Path(FEATURE_MANIFEST_PATH).name,
        **common_kwargs,
    )

    print(f"Done!\n  {SIDECAR_OUTPUT}\n  {PACKAGE_OUTPUT}")
    print("Coloured by leiden; switch to morphology/QC metrics via the colour selector.")


if __name__ == "__main__":
    main()
