"""
Export a KaroSpace viewer for the ANNOTATED two-sample breast-cancer Xenium
dataset (breast_cancer_rep1 + breast_cancer_rep2). Same data as
xenium_dapi_he_breast_cancer_2samples.py, but with a human-readable cell-type
annotation column (`leiden_label`) used as the default coloring.

Each sample carries BOTH a DAPI morphology image and a post-Xenium H&E image
(pixel aligned on the DAPI hires grid). Because there are two sections and EACH
has its own image pair, the overlays are attached PER SECTION. In the viewer,
open a section and use the dropdown in the "H&E Overlay" controls (or the global
"Overlay" switch in the top bar) to switch DAPI <-> H&E; they share one
alignment, or use ✨ Auto-align.

Usage:
    python examples/xenium_dapi_he_breast_cancer_2samples_annotated.py
    XENIUM_DAPI_HE_BREAST2_ANNOT_H5AD=/path/to/file.h5ad python examples/xenium_dapi_he_breast_cancer_2samples_annotated.py
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
    "XENIUM_DAPI_HE_BREAST2_ANNOT_H5AD",
    "/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/"
    "data/processed/xenium-dapi-he/breast_cancer_2samples.companion.annotated.ready.h5ad",
)

OUTPUT_DIR = REPO_ROOT
IMAGE_DIR = OUTPUT_DIR / "xenium_dapi_he_breast2_annotated_images"
SIDECAR_OUTPUT = str(OUTPUT_DIR / "xenium_dapi_he_breast_cancer_2samples_annotated.html")
PACKAGE_OUTPUT = str(OUTPUT_DIR / "xenium_dapi_he_breast_cancer_2samples_annotated.karospace")
GENE_AUX_PATH = str(OUTPUT_DIR / "xenium_dapi_he_breast_cancer_2samples_annotated.genes.json")

# The annotated label column is the most useful default; keep raw leiden available too.
PRIMARY_COLOR = "leiden_label"
ADDITIONAL_COLORS = [
    "leiden_label",
    "leiden",
    "total_counts",
    "transcript_counts",
    "detected_genes",
    "cell_area",
    "nucleus_area",
]
# Run the cluster analytics (markers / DE / neighbor stats) on both the annotation
# and the raw leiden so either can be explored in the Insights panel.
CLUSTER_COLUMNS = ["leiden_label", "leiden"]


def _to_uint8_gray(arr: np.ndarray) -> np.ndarray:
    """Contrast-stretch a (possibly uint16) grayscale array to uint8 (2–99.5%)."""
    a = np.asarray(arr).astype(np.float32)
    lo, hi = np.percentile(a, (2.0, 99.5))
    if hi <= lo:
        hi = lo + 1.0
    a = np.clip((a - lo) / (hi - lo), 0.0, 1.0)
    return (a * 255.0 + 0.5).astype(np.uint8)


def extract_images_per_sample(h5ad_path: str, out_dir: Path) -> dict:
    """Pull dapi_hires / he_hires for EVERY sample in uns/spatial.

    Returns {sample_key: {"DAPI": path, "H&E": path}} for the layers found.
    Filenames are sample-prefixed so the two samples don't collide.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    per_sample: dict[str, dict] = {}
    with h5py.File(h5ad_path, "r") as f:
        spatial = f["uns/spatial"]
        for sample_key in spatial.keys():
            images = spatial[sample_key].get("images")
            if images is None:
                continue
            layers: list[tuple[str, str]] = []

            if "dapi_hires" in images:
                dapi = _to_uint8_gray(images["dapi_hires"][()])
                dapi_path = out_dir / f"{sample_key}_dapi_hires.png"
                Image.fromarray(dapi, mode="L").save(dapi_path)
                layers.append(("DAPI", str(dapi_path)))
                print(f"  [{sample_key}] wrote {dapi_path.name} {dapi.shape}")

            if "he_hires" in images:
                he = np.asarray(images["he_hires"][()])
                if he.dtype != np.uint8:
                    he = _to_uint8_gray(he) if he.ndim == 2 else he.astype(np.uint8)
                he_path = out_dir / f"{sample_key}_he_hires.png"
                Image.fromarray(he, mode="RGB" if he.ndim == 3 else "L").save(he_path)
                layers.append(("H&E", str(he_path)))
                print(f"  [{sample_key}] wrote {he_path.name} {he.shape}")

            if layers:
                per_sample[str(sample_key)] = dict(layers)
    return per_sample


def main() -> None:
    if not Path(H5AD_PATH).exists():
        raise SystemExit(
            f"H5AD not found: {H5AD_PATH}\nSet XENIUM_DAPI_HE_BREAST2_ANNOT_H5AD or edit the path above."
        )

    print("Extracting DAPI + H&E images for each sample...")
    images_by_sample = extract_images_per_sample(H5AD_PATH, IMAGE_DIR)
    if not images_by_sample:
        print("  Warning: no dapi_hires / he_hires found in uns/spatial — exporting without overlays.")

    print("Loading spatial data...")
    dataset = load_spatial_data(
        H5AD_PATH,
        groupby="sample_id",
        spatial_key="spatial",
        metadata_section=["sample_id"],
    )
    print(f"  {dataset.n_sections} section(s), {dataset.n_cells:,} cells")
    print(f"  section IDs: {[s.section_id for s in dataset.sections]}")

    # Attach each section's OWN image pair (rep1 images -> rep1, rep2 -> rep2).
    # The section_id from groupby='sample_id' equals the uns/spatial sample key.
    section_images = {}
    for s in dataset.sections:
        layers = images_by_sample.get(str(s.section_id))
        if layers:
            section_images[s.section_id] = layers
        else:
            print(f"  Warning: no images for section '{s.section_id}' — it will export without an overlay.")
    if not section_images:
        section_images = None

    common_kwargs = dict(
        main_cells_annotation=PRIMARY_COLOR,
        title="Breast Cancer Xenium (2 samples, annotated) — DAPI + H&E",
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
        section_images=section_images,
        section_images_max_px=4096,
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
        "Open a section, then use the overlay dropdown (H&E Overlay controls) or the "
        "global 'Overlay' switch to toggle DAPI <-> H&E. The default annotation is "
        "leiden_label."
    )


if __name__ == "__main__":
    main()
