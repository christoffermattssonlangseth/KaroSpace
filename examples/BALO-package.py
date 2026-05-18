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
GENE_AUX_PATH = "BALO.genes.json"

dataset = load_spatial_data(
    H5AD_PATH,
    groupby="sample_id",
    metadata_columns=[
        "sample_id",
        "run",
    ],
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Available color columns: {dataset.obs_columns[:10]}...")

export_to_html(
    dataset,
    output_path=OUTPUT_PATH,
    color=PRIMARY_CLUSTER,
    title="KaroSpace",
    min_panel_size=120,
    spot_size="auto",
    downsample=100000,
    theme="light",
    outline_by="sample_id",
    additional_colors=[
        "leiden_2",
        "leiden_0.5",
        "run",
    ],
    genes=[],
    use_hvgs=USE_HVGS,
    hvg_limit=200,
    gene_storage="sidecar",
    gene_aux_path=GENE_AUX_PATH,
    marker_genes_groupby=[PRIMARY_CLUSTER] if ENABLE_ANALYTICS else None,
    marker_genes_top_n=30,
    neighbor_stats_groupby=[PRIMARY_CLUSTER] if ENABLE_ANALYTICS else None,
    neighbor_stats_permutations=0,
    neighbor_stats_seed=42,
    cluster_de_groupby=[PRIMARY_CLUSTER] if ENABLE_ANALYTICS else None,
    cluster_de_top_n=20,
    cluster_de_method="t-test",
    cluster_de_layer="normalized",
    cluster_de_min_cells=20,
    interaction_markers_groupby=None,
)

print(f"\nDone! Share {OUTPUT_PATH} together with BALO.loader.html.")
print("Collaborator workflow:")
print("  1. Double-click BALO.loader.html")
print(f"  2. Choose or drop {OUTPUT_PATH}")
