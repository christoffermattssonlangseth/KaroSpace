"""
Export a KaroSpace viewer for the mouse-brain Xenium dataset that carries BOTH a
DAPI morphology image and a post-Xenium H&E image.

The two images were resampled onto the same DAPI hires grid (they are pixel
aligned), so in the viewer they become two overlay LAYERS on the same section:
open a sample, then use the dropdown in the "H&E Overlay" controls to switch
between DAPI and H&E. They share one alignment — align once, then fine-tune and
"Save session" to persist it.

Usage:
    python examples/xenium_dapi_he_mouse_brain.py
    XENIUM_DAPI_HE_H5AD=/path/to/file.h5ad python examples/xenium_dapi_he_mouse_brain.py
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

import h5py
import numpy as np
from PIL import Image

from karospace import export_to_html, load_spatial_data

H5AD_PATH = os.environ.get(
    "XENIUM_DAPI_HE_H5AD",
    "/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/"
    "data/processed/xenium-dapi-he/mouse_brain_prime_ff.companion.ready.h5ad",
)

OUTPUT_DIR = REPO_ROOT
IMAGE_DIR = OUTPUT_DIR / "xenium_dapi_he_images"
SIDECAR_OUTPUT = str(OUTPUT_DIR / "xenium_dapi_he_mouse_brain.html")
PACKAGE_OUTPUT = str(OUTPUT_DIR / "xenium_dapi_he_mouse_brain.karospace")
FEATURE_MANIFEST_PATH = str(OUTPUT_DIR / "xenium_dapi_he_mouse_brain.features.json")

PRIMARY_ANNOTATION = "leiden"
ADDITIONAL_ANNOTATIONS = [
    "leiden",
    "total_counts",
    "transcript_counts",
    "detected_genes",
    "cell_area",
    "nucleus_area",
]
CLUSTER_COLUMNS = ["leiden"]


def _to_uint8_gray(arr: np.ndarray) -> np.ndarray:
    """Contrast-stretch a (possibly uint16) grayscale array to uint8 (2–99.5%)."""
    a = np.asarray(arr).astype(np.float32)
    lo, hi = np.percentile(a, (2.0, 99.5))
    if hi <= lo:
        hi = lo + 1.0
    a = np.clip((a - lo) / (hi - lo), 0.0, 1.0)
    return (a * 255.0 + 0.5).astype(np.uint8)


def extract_images(h5ad_path: str, out_dir: Path) -> dict:
    """Pull dapi_hires / he_hires out of uns/spatial and write them as PNGs.

    Returns {"DAPI": path, "H&E": path} for the layers that were found.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    layers: list[tuple[str, str]] = []
    with h5py.File(h5ad_path, "r") as f:
        spatial = f["uns/spatial"]
        sample_key = list(spatial.keys())[0]
        images = spatial[sample_key]["images"]

        if "dapi_hires" in images:
            dapi = _to_uint8_gray(images["dapi_hires"][()])
            dapi_path = out_dir / "dapi_hires.png"
            Image.fromarray(dapi, mode="L").save(dapi_path)
            layers.append(("DAPI", str(dapi_path)))
            print(f"  wrote {dapi_path.name} {dapi.shape}")

        if "he_hires" in images:
            he = np.asarray(images["he_hires"][()])
            if he.dtype != np.uint8:
                he = _to_uint8_gray(he) if he.ndim == 2 else he.astype(np.uint8)
            he_path = out_dir / "he_hires.png"
            Image.fromarray(he, mode="RGB" if he.ndim == 3 else "L").save(he_path)
            layers.append(("H&E", str(he_path)))
            print(f"  wrote {he_path.name} {he.shape}")

    return dict(layers)


def main() -> None:
    if not Path(H5AD_PATH).exists():
        raise SystemExit(
            f"H5AD not found: {H5AD_PATH}\nSet XENIUM_DAPI_HE_H5AD or edit the path above."
        )

    print("Extracting DAPI + H&E images...")
    image_layers = extract_images(H5AD_PATH, IMAGE_DIR)
    if not image_layers:
        print("  Warning: no dapi_hires / he_hires found in uns/spatial — exporting without overlays.")

    print("Loading spatial data...")
    dataset = load_spatial_data(
        H5AD_PATH,
        section_key="sample_id",
        spatial_key="spatial",
        section_metadata=["sample_id"],
    )
    print(f"  {dataset.n_sections} section(s), {dataset.n_cells:,} cells")
    print(f"  section IDs: {[s.section_id for s in dataset.sections]}")

    # Attach BOTH images (DAPI + H&E) to every section as named overlay layers.
    # They share the DAPI hires grid, so the same pair applies to each section.
    section_images = (
        {s.section_id: image_layers for s in dataset.sections} if image_layers else None
    )

    common_kwargs = dict(
        main_cell_annotation=PRIMARY_ANNOTATION,
        title="Mouse Brain Xenium — DAPI + H&E",
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
        section_images=section_images,
        section_images_max_px=4096,
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
        "Open a section, then use the overlay dropdown (H&E Overlay controls) to "
        "switch DAPI <-> H&E. Align once, then Save session to persist it."
    )


if __name__ == "__main__":
    main()
