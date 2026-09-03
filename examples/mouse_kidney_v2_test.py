"""
KaroSpace viewer for the mouse_kidney Illumina NovaSeq X spatial dataset
(companion-ready) — a single section, 276,925 cells × 56,748 genes.

Carries an RGB H&E `he_hires` image (uns/spatial/mouse_kidney) attached as an
overlay (open the sample → H&E Overlay controls or ✨ Auto-align). Gene
expression is shown from the `lognorm` layer; analytics / DE are enabled for the
two columns the companion precomputed stats for: `leiden` and
`cellcharter_domains`. The remaining cellcharter_k* columns are exposed as plain
categorical overlays (no precomputed analytics).

Analytics (marker genes, cluster-vs-cluster DE, neighbor composition stats) are
read straight from uns/karospace_companion — requesting the two analytics
columns reuses those precomputed results instead of recomputing them.

Usage:
    python examples/mouse_kidney.py
    MOUSE_KIDNEY_H5AD=/path/to/file.h5ad python examples/mouse_kidney.py
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
    "MOUSE_KIDNEY_H5AD",
    "/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/"
    "data/processed/illumina-novaseqx-spatial/mouse_kidney.companion.ready.h5ad",
)

OUTPUT_DIR = REPO_ROOT
IMAGE_DIR = OUTPUT_DIR / "mouse_kidney_images"
SIDECAR_OUTPUT = str(OUTPUT_DIR / "mouse_kidney_v2.html")
PACKAGE_OUTPUT = str(OUTPUT_DIR / "mouse_kidney_v2.karospace")
# Sidecar accepts a full path; the .karospace packager needs a bare filename
# (it lives inside the archive), so pass Path(...).name there.
GENE_AUX_PATH = str(OUTPUT_DIR / "mouse_kidney_v2.genes.json")

GROUPBY = "sample_id"
PRIMARY_COLOR = "leiden"

# Only leiden + cellcharter_domains have precomputed companion analytics; these
# get DE / marker / neighbor-stats panels.
ANALYTICS_COLUMNS = ["leiden", "cellcharter_domains"]

# Extra cellcharter resolutions: plain categorical overlays, NOT analytics.
PLAIN_CATEGORICAL = [
    "cellcharter_k8",
    "cellcharter_k10",
    "cellcharter_k12",
    "cellcharter_k15",
    "cellcharter_k18",
    "cellcharter_k22",
]

# Continuous QC / morphology metrics available as colourings.
CONTINUOUS_METRICS = [
    "total_counts",
    "n_genes",
    "area_um2",
    "nuc_area_um2",
    "nuc_cyto_ratio",
]

ADDITIONAL_COLORS = ANALYTICS_COLUMNS + PLAIN_CATEGORICAL + CONTINUOUS_METRICS


def _read_sample_ids(f: "h5py.File") -> np.ndarray:
    """Return obs/sample_id as a string array, handling categorical encoding."""
    node = f["obs/sample_id"]
    if isinstance(node, h5py.Group):  # categorical: codes + categories
        cats = np.asarray(node["categories"][()])
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
        codes = np.asarray(node["codes"][()])
        out = np.where(codes >= 0, cats[codes.clip(min=0)], "")
        return out.astype(str)
    arr = np.asarray(node[()])
    return np.array([c.decode() if isinstance(c, bytes) else str(c) for c in arr])


def extract_he_image(h5ad_path: str, out_dir: Path) -> dict:
    """Pull each sample's RGB `he_hires` H&E image, cropped to its cells.

    The embedded `he_hires` bitmap carries gray margins around the tissue, but
    the viewer auto-fits the WHOLE image to the cell bounding box — so the tissue
    ends up uniformly shrunk. We crop the image to the exact cell bounding box
    (converted from microns to hires pixels via the squidpy scalefactors) so the
    auto-fit lands 1:1. Conversion: hires_px = micron * he_hires_scalef /
    pixel_size_um (obsm/spatial is in microns; he_hires = full-res * scalef).

    Returns {sample_key: {"H&E": path}} for samples carrying an image.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    per_sample: dict[str, dict] = {}
    with h5py.File(h5ad_path, "r") as f:
        if "spatial" not in f.get("uns", {}):
            return per_sample
        spatial = f["uns/spatial"]
        coords = np.asarray(f["obsm/spatial"][()])[:, :2]
        sample_ids = _read_sample_ids(f)
        for sample_key in spatial.keys():
            node = spatial[sample_key]
            images = node.get("images") if isinstance(node, h5py.Group) else None
            if images is None or "he_hires" not in images:
                continue
            arr = np.asarray(images["he_hires"][()])
            if arr.dtype != np.uint8:
                a = arr.astype(np.float32)
                lo, hi = np.percentile(a, (2.0, 99.5))
                hi = hi if hi > lo else lo + 1.0
                arr = (np.clip((a - lo) / (hi - lo), 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
            mode = "RGB" if arr.ndim == 3 else "L"
            img = Image.fromarray(arr, mode=mode)

            # Crop to this sample's cell bounding box, in hires-image pixels.
            sf = node.get("scalefactors") if isinstance(node, h5py.Group) else None
            mask = sample_ids == str(sample_key)
            if sf is not None and "he_hires_scalef" in sf and "pixel_size_um" in sf and mask.any():
                scalef = float(sf["he_hires_scalef"][()])
                pixel_um = float(sf["pixel_size_um"][()])
                m2px = scalef / pixel_um  # microns -> hires pixels
                cc = coords[mask]
                W, H = img.size
                left = max(0, int(np.floor(cc[:, 0].min() * m2px)))
                right = min(W, int(np.ceil(cc[:, 0].max() * m2px)))
                upper = max(0, int(np.floor(cc[:, 1].min() * m2px)))
                lower = min(H, int(np.ceil(cc[:, 1].max() * m2px)))
                if right > left and lower > upper:
                    img = img.crop((left, upper, right, lower))
                    print(
                        f"  Cropped '{sample_key}' H&E to cell bbox "
                        f"({right - left}×{lower - upper}px from {W}×{H}px, "
                        f"{1.0 / m2px:.3f} µm/px)"
                    )
                else:
                    print(f"  Warning: degenerate crop box for '{sample_key}' — using full image")
            else:
                print(f"  Warning: no scalefactors for '{sample_key}' — using full (uncropped) image")

            path = out_dir / f"{sample_key}_he_hires.png"
            img.save(path)
            per_sample[str(sample_key)] = {"H&E": str(path)}
    return per_sample


def main() -> None:
    if not Path(H5AD_PATH).exists():
        raise SystemExit(
            f"H5AD not found: {H5AD_PATH}\nSet MOUSE_KIDNEY_H5AD or edit the path above."
        )

    print("Extracting H&E image(s)...")
    images_by_sample = extract_he_image(H5AD_PATH, IMAGE_DIR)
    print(f"  extracted {len(images_by_sample)} image(s)")

    print("Loading spatial data (276k cells × 56.7k genes — this can take a while)...")
    dataset = load_spatial_data(
        H5AD_PATH,
        section_key=GROUPBY,          # was: groupby (renamed in the merge)
        spatial_key="spatial",
        section_metadata=[GROUPBY],   # was: metadata_columns
    )
    print(f"  {dataset.n_sections} section(s), {dataset.n_cells:,} cells")

    # The viewer's gene-expression path reads a layer literally named
    # "normalized"; this file stores normalized counts under "lognorm". Alias it
    # so both gene display and cluster DE use the lognorm values.
    if "normalized" not in dataset.adata.layers and "lognorm" in dataset.adata.layers:
        dataset.adata.layers["normalized"] = dataset.adata.layers["lognorm"]
        print("  Aliased layers['lognorm'] -> layers['normalized'] for gene display/DE")

    # Attach the H&E overlay to its matching section.
    section_images = {
        s.section_id: images_by_sample[str(s.section_id)]
        for s in dataset.sections
        if str(s.section_id) in images_by_sample
    } or None

    common_kwargs = dict(
        color=PRIMARY_COLOR,
        title="Mouse Kidney — NovaSeq X Spatial",
        min_panel_size=140,
        spot_size="auto",
        theme="light",
        outline_by=None,
        additional_colors=ADDITIONAL_COLORS,
        genes=[],
        use_hvgs=False,
        gene_storage="sidecar",
        gene_sidecar_format="binary-v1",
        gene_encoding="auto",
        gene_value_encoding="uint8",
        gene_sidecar_shard_size=128,
        # Analytics enabled ONLY for the two companion-precomputed columns;
        # requesting them reuses uns/karospace_companion instead of recomputing.
        marker_genes_groupby=ANALYTICS_COLUMNS,
        marker_genes_top_n=30,
        neighbor_stats_groupby=ANALYTICS_COLUMNS,
        neighbor_stats_permutations=0,
        cluster_de_groupby=ANALYTICS_COLUMNS,
        cluster_de_top_n=20,
        cluster_de_method="t-test",
        cluster_de_layer="normalized",
        interaction_markers_groupby=None,
        section_images=section_images,
        section_images_max_px=4096,
    )

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

    print(f"\nDone!\n  sidecar: {SIDECAR_OUTPUT}\n  package: {PACKAGE_OUTPUT}")
    print(
        "Coloured by leiden; DE/analytics enabled for leiden & cellcharter_domains. "
        "cellcharter_k* are plain categorical overlays. Open the sample for the H&E overlay."
    )


if __name__ == "__main__":
    main()
