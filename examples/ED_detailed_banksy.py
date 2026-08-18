"""
Export a KaroSpace viewer for the ED peripheral-nerve dataset with Banksy
clustering (deitaled_wBanky-lam08res03_pID) — 26 samples, ~212k cells, one
section per sample (by sample_NewNAME).

Per collaborator request:
  - sections grouped by sample_NewNAME (friendly names, e.g. CT1-045148)
  - patient identifiers exposed as metadata (pID = new IDs, patient_id = old)
  - condition / status removed
  - colour options limited to: cell_annotation_detailed, cell_class, banksy
    (cell_subclass, cellcharter, leiden removed)

No embedded morphology images (no uns/spatial).

Usage:
    python examples/ED_detailed_banksy.py
    ED_BANKSY_H5AD=/path/to/file.h5ad python examples/ED_detailed_banksy.py
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
    "ED_BANKSY_H5AD",
    "/Users/chrislangseth/Downloads/deitaled_wBanky-lam08res03_pID.companion.ready.h5ad",
)

OUTPUT_DIR = REPO_ROOT
SIDECAR_OUTPUT = str(OUTPUT_DIR / "ED_detailed_banksy.html")
PACKAGE_OUTPUT = str(OUTPUT_DIR / "ED_detailed_banksy.karospace")
FEATURE_MANIFEST_PATH = str(OUTPUT_DIR / "ED_detailed_banksy.features.json")

PRIMARY_ANNOTATION = "cell_annotation_detailed"
# Only the requested colour options: detailed annotation, coarse class, Banksy.
ADDITIONAL_ANNOTATIONS = [
    "cell_annotation_detailed",
    "cell_class",
    "banksy",
]
CLUSTER_COLUMNS = ["cell_annotation_detailed", "cell_class", "banksy"]


def main() -> None:
    if not Path(H5AD_PATH).exists():
        raise SystemExit(
            f"H5AD not found: {H5AD_PATH}\nSet ED_BANKSY_H5AD or edit the path above."
        )

    print("Loading spatial data (large file — this can take a while)...")
    dataset = load_spatial_data(
        H5AD_PATH,
        section_key="sample_NewNAME",
        spatial_key="spatial",
        # condition / status intentionally excluded; expose patient identifiers.
        section_metadata=["sample_NewNAME", "pID", "patient_id", "slide"],
    )
    print(f"  {dataset.n_sections} section(s), {dataset.n_cells:,} cells")
    print(f"  section IDs: {[s.section_id for s in dataset.sections][:6]} ...")

    common_kwargs = dict(
        main_cell_annotation=PRIMARY_ANNOTATION,
        title="ED Peripheral Nerve — Detailed Annotation + Banksy",
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
        "Coloured by cell_annotation_detailed; switch to cell_class or banksy via the "
        "colour selector. Sections grouped by sample_NewNAME."
    )


if __name__ == "__main__":
    main()
