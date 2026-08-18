"""
Export a KaroSpace viewer for the detailed-annotation Xenium dataset
(detailedAnnot_incl_noAnnot_260602) — a 26-sample peripheral-nerve / ED study
(212k cells), split into one section per sample.

No embedded morphology images (no uns/spatial); it's a straight spatial cell map
shown by the detailed cell annotations:
  - cell_annotation_detailed (39 fine cell types) — default annotation
  - cell_class (5) / cell_subclass (16) — coarser hierarchies
  - cytetype_annotation_leiden0.2, cluster_cellcharter_10 (spatial niches)
Clinical variables (condition, status, diabetes, patient, slide) are exposed as
section metadata for filtering / the samples view.

NOTE: the source .h5ad is ~6.5 GB / 212k cells / 26 sections, so this export is
heavier and slower than the other examples.

Usage:
    python examples/detailedAnnot_incl_noAnnot.py
    XENIUM_DETAILED_ANNOT_H5AD=/path/to/file.h5ad python examples/detailedAnnot_incl_noAnnot.py
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
    "XENIUM_DETAILED_ANNOT_H5AD",
    "/Users/chrislangseth/Downloads/detailedAnnot_incl_noAnnot_260602.companion.ready.h5ad",
)

OUTPUT_DIR = REPO_ROOT
SIDECAR_OUTPUT = str(OUTPUT_DIR / "detailedAnnot_incl_noAnnot_260602.html")
PACKAGE_OUTPUT = str(OUTPUT_DIR / "detailedAnnot_incl_noAnnot_260602.karospace")
FEATURE_MANIFEST_PATH = str(OUTPUT_DIR / "detailedAnnot_incl_noAnnot_260602.features.json")

PRIMARY_ANNOTATION = "cell_annotation_detailed"
ADDITIONAL_ANNOTATIONS = [
    "cell_annotation_detailed",
    "cell_class",
    "cell_subclass",
    "cytetype_annotation_leiden0.2",
    "cluster_cellcharter_10",
    "leiden0.2",
    "n_counts",
    "n_genes",
]
# Annotation columns to run the cluster analytics (markers / DE / neighbor stats) on.
CLUSTER_COLUMNS = ["cell_annotation_detailed", "cell_class"]


def main() -> None:
    if not Path(H5AD_PATH).exists():
        raise SystemExit(
            f"H5AD not found: {H5AD_PATH}\nSet XENIUM_DETAILED_ANNOT_H5AD or edit the path above."
        )

    print("Loading spatial data (large file — this can take a while)...")
    dataset = load_spatial_data(
        H5AD_PATH,
        section_key="sample_id",
        spatial_key="spatial",
        section_metadata=[
            "sample_id",
            "sample_name",
            "condition",
            "condition_subtype",
            "status",
            "patient_id",
            "Diabetes",
            "ED",
            "slide",
        ],
    )
    print(f"  {dataset.n_sections} section(s), {dataset.n_cells:,} cells")
    print(f"  section IDs: {[s.section_id for s in dataset.sections][:6]} ...")

    common_kwargs = dict(
        main_cell_annotation=PRIMARY_ANNOTATION,
        title="Xenium — Detailed Annotation (26 samples)",
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
    print(
        "Coloured by cell_annotation_detailed; switch to cell_class / cell_subclass / "
        "cluster_cellcharter_10 (niches) via the colour selector. Clinical variables are "
        "available as section metadata."
    )


if __name__ == "__main__":
    main()
