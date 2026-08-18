import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("KMP_WARNINGS", "0")

from karospace import export_to_html, load_spatial_data

H5AD_PATH = os.environ.get(
    "BALO_H5AD_PATH",
    "/tmp/baloMS_companion.ready.h5ad",
)

PRIMARY_CLUSTER = "leiden_2_names_sub"
USE_HVGS = False
ENABLE_ANALYTICS = True
OUTPUT_PATH = "BALO.karospace"
FEATURE_MANIFEST_PATH = "BALO.features.json"

dataset = load_spatial_data(
    H5AD_PATH,
    section_key="sample_id",
    section_metadata=[
        "sample_id",
        "run",
    ],
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available annotation columns: {dataset.obs_columns[:10]}...")

export_to_html(
    dataset,
    output_path=OUTPUT_PATH,
    main_cell_annotation=PRIMARY_CLUSTER,
    title="KaroSpace",
    min_panel_size=120,
    spot_size="auto",
    downsample=100000,
    outline_by="sample_id",
    cell_annotations=[
        "leiden_2",
        "leiden_0.5",
        "run",
    ],
    features=[],
    use_hvgs=USE_HVGS,
    hvg_limit=200,
    feature_storage="sidecar",
    feature_manifest_path=FEATURE_MANIFEST_PATH,
    marker_gene_annotations=[PRIMARY_CLUSTER] if ENABLE_ANALYTICS else None,
    marker_genes_top_n=30,
    neighbor_stats_annotations=[PRIMARY_CLUSTER] if ENABLE_ANALYTICS else None,
    neighbor_stats_permutations=0,
    neighbor_stats_seed=42,
    pseudobulk_de_annotations=[PRIMARY_CLUSTER] if ENABLE_ANALYTICS else None,
    pseudobulk_de_top_n=20,
    pseudobulk_de_method="t-test",
    pseudobulk_de_layer="normalized",
    pseudobulk_de_min_cells=20,
    interaction_marker_annotations=None,
)

print(f"\nDone! Share {OUTPUT_PATH} together with BALO.loader.html.")
print("Collaborator workflow:")
print("  1. Double-click BALO.loader.html")
print(f"  2. Choose or drop {OUTPUT_PATH}")
