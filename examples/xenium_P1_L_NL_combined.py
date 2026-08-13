"""
Export a KaroSpace viewer for the annotated P1 skin Xenium dataset, Lesional (L)
and Non-lesional (NL) combined into one object (two sections, by sample_id).

This dataset has NO embedded morphology images (no uns/spatial), so there is no
DAPI/H&E overlay — it's a straight spatial cell map shown by the Seurat
annotations:
  - predicted.celltype (19 skin cell types) — default annotation
  - niches (5 spatial niches)
  - seurat_clusters / Xenium_snn_res.0.2
  - an IFN-alpha REACTOME signature score (continuous)

Usage:
    python examples/xenium_P1_L_NL_combined.py
    XENIUM_P1_H5AD=/path/to/file.h5ad python examples/xenium_P1_L_NL_combined.py
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
    "XENIUM_P1_H5AD",
    "/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/"
    "data/processed/xenium-P1-annotated/P1_L_NL_combined.companion.ready.h5ad",
)

OUTPUT_DIR = REPO_ROOT
SIDECAR_OUTPUT = str(OUTPUT_DIR / "xenium_P1_L_NL_combined.html")
PACKAGE_OUTPUT = str(OUTPUT_DIR / "xenium_P1_L_NL_combined.karospace")
GENE_AUX_PATH = str(OUTPUT_DIR / "xenium_P1_L_NL_combined.genes.json")

PRIMARY_COLOR = "predicted.celltype"
ADDITIONAL_COLORS = [
    "predicted.celltype",
    "niches",
    "seurat_clusters",
    "segmentation_method",
    "nCount_Xenium",
    "nFeature_Xenium",
    "signature_1IFN_alpha_REACTOME",
]
# Annotation columns to run the cluster analytics (markers / DE / neighbor stats) on.
CLUSTER_COLUMNS = ["predicted.celltype", "niches"]


def main() -> None:
    if not Path(H5AD_PATH).exists():
        raise SystemExit(
            f"H5AD not found: {H5AD_PATH}\nSet XENIUM_P1_H5AD or edit the path above."
        )

    print("Loading spatial data...")
    dataset = load_spatial_data(
        H5AD_PATH,
        groupby="sample_id",
        spatial_key="spatial",
        metadata_section=["sample_id", "condition", "sample_batch"],
    )
    print(f"  {dataset.n_sections} section(s), {dataset.n_cells:,} cells")
    print(f"  section IDs: {[s.section_id for s in dataset.sections]}")

    common_kwargs = dict(
        main_cells_annotation=PRIMARY_COLOR,
        title="Xenium P1 Skin — Lesional + Non-lesional (annotated)",
        min_panel_size=140,
        spot_size="auto",
        outline_by=None,
        cells_annotations=ADDITIONAL_COLORS,
        genes=[],
        use_hvgs=False,
        gene_storage="sidecar",
        gene_encoding="auto",
        gene_value_encoding="uint8",
        gene_sidecar_shard_size=128,
        marker_genes_groupby=CLUSTER_COLUMNS,
        marker_genes_top_n=30,
        neighbor_stats_groupby=CLUSTER_COLUMNS,
        neighbor_stats_permutations=0,
        cluster_de_groupby=CLUSTER_COLUMNS,
        cluster_de_top_n=20,
        cluster_de_method="t-test",
        cluster_de_layer=None,
        interaction_markers_groupby=None,
    )

    # The sidecar writes its gene manifest next to the HTML (full path is fine);
    # the .karospace packager requires a bare filename (it lives inside the archive).
    print("Exporting sidecar HTML...")
    export_to_html(
        dataset, output_path=SIDECAR_OUTPUT, gene_aux_path=GENE_AUX_PATH, **common_kwargs
    )

    print("Packaging .karospace archive...")
    export_to_html(
        dataset,
        output_path=PACKAGE_OUTPUT,
        gene_aux_path=Path(GENE_AUX_PATH).name,
        **common_kwargs,
    )

    print(f"Done!\n  {SIDECAR_OUTPUT}\n  {PACKAGE_OUTPUT}")
    print(
        "Coloured by predicted.celltype; switch to 'niches', 'seurat_clusters', or the "
        "IFN-alpha signature via the colour selector."
    )


if __name__ == "__main__":
    main()
