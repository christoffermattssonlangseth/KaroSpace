"""
KaroSpace export for the full stroke (dMCAO) dataset stroke_all (companion-ready)
— 49 samples, ~3.66M cells, one section per sample.

Each sample carries a grayscale `hires` morphology image (uns/spatial), attached
PER SECTION as an overlay (open a sample → use the H&E Overlay controls or
✨ Auto-align). Coloured by the coarse `cluster` annotation, with cell
morphology / QC metrics available as continuous colourings.

~3.66M cells / 500 genes — writes a binary-sidecar viewer bundle and a matching
packaged .karospace.

Usage:
    python examples/stroke_all.py
    STROKE_ALL_H5AD=/path/to/file.h5ad python examples/stroke_all.py
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
    "STROKE_ALL_H5AD",
    "/Users/chrislangseth/Downloads/stroke_all.companion.ready.h5ad",
)

OUTPUT_DIR = REPO_ROOT
IMAGE_DIR = OUTPUT_DIR / "stroke_all_images"
SIDECAR_OUTPUT = "stroke_all.html"
PACKAGE_OUTPUT = "stroke_all.karospace"
# Must be a BARE filename: it's written next to the sidecar HTML and packaged by
# name inside the .karospace archive (a full path is rejected for the package).
GENE_AUX_PATH = "stroke_all.genes.json"

PRIMARY_COLOR = "cluster"
ADDITIONAL_COLORS = [
    "cluster",
    "density",
    "area",
    "elongation",
    "n_transcripts",
    "avg_confidence",
    "max_cluster_frac",
    "lifespan",
]
CLUSTER_COLUMNS = ["cluster"]


def extract_images_per_sample(h5ad_path: str, out_dir: Path) -> dict:
    """Pull the grayscale 'hires' morphology image for every sample in uns/spatial.

    Returns {sample_key: {"Morphology": path}} for samples that have an image.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    per_sample: dict[str, dict] = {}
    with h5py.File(h5ad_path, "r") as f:
        spatial = f["uns/spatial"]
        for sample_key in spatial.keys():
            node = spatial[sample_key]
            images = node.get("images") if isinstance(node, h5py.Group) else None
            if images is None or "hires" not in images:
                continue
            arr = np.asarray(images["hires"][()])
            if arr.dtype != np.uint8:
                a = arr.astype(np.float32)
                lo, hi = np.percentile(a, (2.0, 99.5))
                hi = hi if hi > lo else lo + 1.0
                arr = (np.clip((a - lo) / (hi - lo), 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
            mode = "RGB" if arr.ndim == 3 else "L"
            path = out_dir / f"{sample_key}_hires.png"
            Image.fromarray(arr, mode=mode).save(path)
            per_sample[str(sample_key)] = {"Morphology": str(path)}
    return per_sample


def main() -> None:
    if not Path(H5AD_PATH).exists():
        raise SystemExit(
            f"H5AD not found: {H5AD_PATH}\nSet STROKE_ALL_H5AD or edit the path above."
        )

    print("Extracting per-sample morphology images...")
    images_by_sample = extract_images_per_sample(H5AD_PATH, IMAGE_DIR)
    print(f"  extracted {len(images_by_sample)} images")

    print("Loading spatial data (16 GB / ~3.66M cells — this can take a while)...")
    dataset = load_spatial_data(
        H5AD_PATH,
        groupby="sample",
        spatial_key="spatial",
        metadata_columns=["sample"],
    )
    print(f"  {dataset.n_sections} section(s), {dataset.n_cells:,} cells")

    # Attach each section's own morphology image (matched by sample key).
    section_images = {}
    for s in dataset.sections:
        layers = images_by_sample.get(str(s.section_id))
        if layers:
            section_images[s.section_id] = layers
    if not section_images:
        section_images = None

    common_kwargs = dict(
        color=PRIMARY_COLOR,
        title="Stroke (dMCAO) — all samples",
        min_panel_size=120,
        spot_size="auto",
        downsample=10_000_000,
        theme="light",
        outline_by=None,
        additional_colors=ADDITIONAL_COLORS,
        genes=[],
        use_hvgs=False,
        gene_storage="sidecar",
        gene_sidecar_format="binary-v1",
        gene_encoding="auto",
        gene_value_encoding="uint8",
        gene_aux_path=GENE_AUX_PATH,
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

    print("Exporting binary-sidecar viewer...")
    export_to_html(dataset, output_path=SIDECAR_OUTPUT, **common_kwargs)

    print("Packaging .karospace archive...")
    export_to_html(dataset, output_path=PACKAGE_OUTPUT, **common_kwargs)

    print(f"\nDone!\n  sidecar: {SIDECAR_OUTPUT}\n  package: {PACKAGE_OUTPUT}")
    print("Coloured by cluster; switch to morphology/QC metrics via the colour selector. Per-section morphology images attach as overlays.")


if __name__ == "__main__":
    main()
