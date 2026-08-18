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
    "BALO_H5AD_PATH",
    "/Users/yuk.kit.lor/bioinfo/balo/baloMS_indep_clust_balo_MANA_balo_noMT_annot_neutro_mana.h5ad",
)

# Optional: map section IDs to H&E image paths.
# Leave empty if you have no staining images yet.
SECTION_IMAGES = {
    "Balo-1": "/Users/yuk.kit.lor/bioinfo/XeniumAlign/Balo1_hne.tif",
    "Balo-2": "/Users/yuk.kit.lor/bioinfo/XeniumAlign/Balo2_hne.tif",
}

PRIMARY_ANNOTATION = "leiden_2_names"
ADDITIONAL_ANNOTATIONS = [
    "n_genes_by_counts",
    "n_counts",
    "n_genes",
    "leiden_0.5_names",
    "leiden_1_names",
    "leiden_2_names",
    "gmm_mana_5",
    "gmm_mana_8",
    "gmm_mana_10",
    "gmm_mana_12",
    "gmm_mana_15",
    "gmm_mana_20",
]

CLUSTER_COLUMNS = [
    "leiden_0.5_names",
    "leiden_1_names",
    "leiden_2_names",
    "gmm_mana_5",
    "gmm_mana_8",
    "gmm_mana_10",
    "gmm_mana_12",
    "gmm_mana_15",
    "gmm_mana_20",
]

SIDECAR_OUTPUT = "balo_may2026.html"
PACKAGE_OUTPUT = "balo_may2026.karospace"
FEATURE_MANIFEST_PATH  = "balo_may2026.features.json"

if not Path(H5AD_PATH).exists():
    raise SystemExit(f"H5AD not found: {H5AD_PATH}\nSet BALO_H5AD_PATH or update the path above.")

dataset = load_spatial_data(
    H5AD_PATH,
    section_key="sample_id",
    spatial_key="spatial",
    section_metadata=["sample_id"],
)

print(f"Loaded {dataset.n_sections} sections, {dataset.n_cells:,} cells")
print(f"Section IDs: {[s.section_id for s in dataset.sections]}")

common_kwargs = dict(
    main_cell_annotation=PRIMARY_ANNOTATION,
    title="Balo Lesion May 2026",
    min_panel_size=120,
    spot_size="auto",
    outline_by=None,
    cell_annotations=ADDITIONAL_ANNOTATIONS,
    features=[],
    use_hvgs=False,
    hvg_limit=50,
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
    section_images=SECTION_IMAGES or None,
    section_images_max_px=4096,
)

print("Exporting sidecar HTML...")
export_to_html(dataset, output_path=SIDECAR_OUTPUT, **common_kwargs)

print("Packaging .karospace archive...")
export_to_html(dataset, output_path=PACKAGE_OUTPUT, **common_kwargs)

print(f"Done! Outputs: {SIDECAR_OUTPUT}, {PACKAGE_OUTPUT}")
