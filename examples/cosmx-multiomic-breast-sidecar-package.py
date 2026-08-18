"""
Example usage of KaroSpace with CosMx multiomic breast cancer data.
Uses binary sidecar storage and creates a .karospace package.
"""

import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Reduce noise and set up caching
os.environ.setdefault("KMP_WARNINGS", "0")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/karospace-mpl-cache")

from karospace import export_to_html, load_spatial_data

# Dataset path provided by user
H5AD_PATH = os.environ.get(
    "COSMX_BREAST_H5AD_PATH",
    "../KaroSpaceDataWrangling/notebooks/data/cosmx_multiomic_breast/cosmx_multiomic_breast.h5ad",
)

PRIMARY_ANNOTATION = "leiden_rna"
ADDITIONAL_ANNOTATIONS = [
    "leiden_protein",
]

# Output naming
BASE_NAME = "cosmx-multiomic-breast"
SIDECAR_OUTPUT = f"{BASE_NAME}-sidecar.html"
PACKAGE_OUTPUT = f"{BASE_NAME}.karospace"
FEATURE_MANIFEST_PATH = f"{BASE_NAME}.features.json"

if not Path(H5AD_PATH).exists():
    raise SystemExit(
        f"CosMx h5ad not found at {H5AD_PATH}. "
        "Please check the path or set COSMX_BREAST_H5AD_PATH."
    )

print(f"Loading CosMx dataset: {H5AD_PATH}")
dataset = load_spatial_data(
    H5AD_PATH,
    section_key="sample", # 'sample' identifies the sections in this dataset
    spatial_key="spatial",
    section_metadata=[
        "sample",
        "tissue",
        "Run_Tissue_name",
        "Panel",
    ],
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")
print(f"Detected modalities: {list(dataset.modalities.keys())}")

# Common export settings favoring large-scale multiomic performance
common_kwargs = dict(
    main_cell_annotation=PRIMARY_ANNOTATION,
    title="CosMx Multiomic Breast Cancer",
    min_panel_size=150,
    spot_size="auto",
    downsample=None, # Keep all cells for high-res exploration
    cell_annotations=ADDITIONAL_ANNOTATIONS,
    
    # Feature discovery
    features=[],        # Can specify initial genes if desired
    use_hvgs=True,
    hvg_limit=50,
    
    # Storage and encoding
    feature_storage="sidecar",
    feature_encoding="auto",
    feature_value_encoding="uint8", # Quantization for smaller sidecar shards
    feature_manifest_path=FEATURE_MANIFEST_PATH,
    feature_sidecar_shard_size=128,
    
    # Analytics
    marker_gene_annotations=[
        "leiden_rna",
        "leiden_protein",
    ],
    neighbor_stats_annotations=[
        "leiden_rna",
        "leiden_protein",
    ],
    pseudobulk_de_annotations=[
        "leiden_rna",
        "leiden_protein",
    ],
    
    # Use detected modalities for protein exploration
    # KaroSpace detects 'protein' from adata.obsm['protein'] and adata.uns['protein_var']
    modalities=["rna", "protein", "protein_channel"],
)

print("\nStarting export to HTML with sidecar...")
export_to_html(
    dataset,
    output_path=SIDECAR_OUTPUT,
    **common_kwargs,
)

print(f"Starting packaging to {PACKAGE_OUTPUT}...")
export_to_html(
    dataset,
    output_path=PACKAGE_OUTPUT,
    **common_kwargs,
)

print(f"\nSuccess! Wrote:")
print(f"  1. Sidecar Viewer: {SIDECAR_OUTPUT}")
print(f"  2. Gene Manifest: {FEATURE_MANIFEST_PATH}")
print(f"  3. Shard Directory: {Path(FEATURE_MANIFEST_PATH).with_suffix('')}/")
print(f"  4. KaroSpace Package: {PACKAGE_OUTPUT}")
print(f"  5. Local Loader: {Path(PACKAGE_OUTPUT).with_suffix('.loader.html')}")
