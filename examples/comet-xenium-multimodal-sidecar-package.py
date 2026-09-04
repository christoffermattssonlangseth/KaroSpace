"""
Export a KaroSpace viewer for the COMET multimodal Xenium dataset.

This is a 4-core tissue microarray (hepatocellular carcinoma / non-tumour HCC
liver / tonsil / hepatocellular adenoma) profiled with BOTH:
  * RNA   — 5,001-gene Xenium panel (adata.X)
  * Protein — 16-channel COMET immunofluorescence (adata.obsm['protein'],
    channel names in adata.uns['protein_var'], e.g. DAPI, CK8, CD45, HepPar1)

Each core also carries a post-stain H&E image. The H&E PNGs are NOT embedded in
the h5ad — `uns/section_images` holds a relative PATH per section and the actual
files live in `he_images/` next to the h5ad. We resolve those paths and attach
each core's own H&E as an overlay, so in the viewer you can open a core and
toggle / auto-align the H&E under the cells.

Because this carries an extra (protein) modality, the exporter requires
feature_storage="sidecar" (a single embedded HTML can't carry a 2nd modality).

Usage:
    python examples/comet-xenium-multimodal-sidecar-package.py
    COMET_XENIUM_H5AD=/path/to/comet_xenium_multimodal.h5ad \
        python examples/comet-xenium-multimodal-sidecar-package.py
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
import scipy.sparse as sp
from scipy.spatial import cKDTree

from karospace import export_to_html, load_spatial_data

# k for the spatial kNN graph. Typical Xenium neighbourhoods use ~6 neighbours.
SPATIAL_GRAPH_K = 6

H5AD_PATH = os.environ.get(
    "COMET_XENIUM_H5AD",
    "/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/"
    "notebooks/data/comet_xenium_multimodal/comet_xenium_multimodal.h5ad",
)

OUTPUT_DIR = REPO_ROOT
BASE_NAME = "comet-xenium-multimodal"
SIDECAR_OUTPUT = str(OUTPUT_DIR / f"{BASE_NAME}-sidecar.html")
PACKAGE_OUTPUT = str(OUTPUT_DIR / f"{BASE_NAME}.karospace")
FEATURE_MANIFEST_PATH = str(OUTPUT_DIR / f"{BASE_NAME}.features.json")

# leiden_rna is the most familiar default; keep the protein-based clustering and
# the tissue label available too.
PRIMARY_ANNOTATION = "leiden_rna"
ADDITIONAL_ANNOTATIONS = [
    "leiden_rna",
    "leiden_protein",
    "tissue",
    "area_px",
]
# Run cluster analytics (markers / DE / neighbor stats) on both clusterings so
# either can be explored in the Insights panel.
CLUSTER_COLUMNS = ["leiden_rna", "leiden_protein"]


def resolve_section_images(h5ad_path: str) -> dict:
    """Map each section to its H&E overlay by resolving the paths in
    uns/section_images. Paths are stored relative to the data repo root; fall
    back to the he_images/ dir sitting next to the h5ad if that fails.

    Returns {section_id: {"H&E": absolute_png_path}} for images found on disk.
    """
    h5ad = Path(h5ad_path).resolve()
    h5ad_dir = h5ad.parent
    # .../<repo>/notebooks/data/<dataset>/file.h5ad -> repo root is parents[3]
    repo_root = h5ad.parents[3] if len(h5ad.parents) > 3 else h5ad_dir

    per_section: dict[str, dict] = {}
    with h5py.File(h5ad_path, "r") as f:
        if "section_images" not in f["uns"]:
            return per_section
        for section_key in f["uns/section_images"].keys():
            raw = f["uns/section_images"][section_key][()]
            rel = raw.decode() if isinstance(raw, bytes) else str(raw)
            candidates = [
                repo_root / rel,
                h5ad_dir / rel,
                h5ad_dir / "he_images" / Path(rel).name,
            ]
            png = next((c for c in candidates if c.exists()), None)
            if png is None:
                print(f"  Warning: H&E for '{section_key}' not found (tried: {rel})")
                continue
            per_section[str(section_key)] = {"H&E": str(png.resolve())}
            print(f"  [{section_key}] H&E -> {png.name}")
    return per_section


def build_spatial_graph(adata, section_col: str, spatial_key: str = "spatial",
                        k: int = SPATIAL_GRAPH_K) -> None:
    """Build a per-section symmetric kNN spatial graph into
    adata.obsp['spatial_connectivities'].

    The exporter's spatially-variable-genes step (Moran's I) only *reads* a graph
    from obsp — unlike neighbor stats, it never builds one from coordinates. This
    plain h5ad ships no graph, so we build one here. Edges are added within each
    section only, so Moran's I never draws spurious links between the 4 cores.
    """
    coords = np.asarray(adata.obsm[spatial_key])[:, :2].astype(float)
    n = coords.shape[0]
    sections = np.asarray(adata.obs[section_col].astype(str))

    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    for sec in np.unique(sections):
        idx = np.flatnonzero(sections == sec)
        if idx.size < 2:
            continue
        kk = int(min(k, idx.size - 1))
        tree = cKDTree(coords[idx])
        # k+1 because the first neighbour is the point itself.
        _, nbr = tree.query(coords[idx], k=kk + 1)
        nbr = np.atleast_2d(nbr)[:, 1:]  # drop self column
        src = np.repeat(idx, nbr.shape[1])
        dst = idx[nbr.ravel()]
        rows.append(src)
        cols.append(dst)

    if not rows:
        print("  Warning: could not build a spatial graph (too few cells).")
        return

    r = np.concatenate(rows)
    c = np.concatenate(cols)
    data = np.ones(r.shape[0], dtype=np.float32)
    W = sp.coo_matrix((data, (r, c)), shape=(n, n)).tocsr()
    W = W.maximum(W.T)  # symmetrize (mutual + one-directional kNN edges)
    # Self-neighbours were already dropped, so the diagonal is zero by construction.
    adata.obsp["spatial_connectivities"] = W
    print(f"  Built spatial graph: {n:,} cells, {W.nnz:,} edges "
          f"(k={k}, per-section, symmetric).")


def clean_protein_channel_names(dataset) -> None:
    """Strip the fluorophore suffix from protein channel names so they're
    addressable by antibody target alone.

    The COMET channels ship as 'CD45 - TRITC', 'CD66 - FITC', etc. In the viewer
    that means colouring requires typing the full 'CD45 - TRITC'; bare 'CD45'
    matches nothing. The suffix is the detection channel, not the target, so we
    drop it — but only if every stripped name stays unique (no collisions).
    """
    mod = dataset.modalities.get("protein")
    if mod is None:
        return
    orig = list(mod.var.index)
    cleaned = [str(n).split(" - ")[0].strip() for n in orig]
    if len(set(cleaned)) != len(cleaned):
        print("  Skipping protein channel rename — stripping suffixes would collide.")
        return
    mod.var.index = cleaned
    renamed = [f"{o} -> {c}" for o, c in zip(orig, cleaned) if o != c]
    print(f"  Cleaned {len(renamed)} protein channel name(s), e.g. {renamed[:3]}")


def main() -> None:
    if not Path(H5AD_PATH).exists():
        raise SystemExit(
            f"H5AD not found: {H5AD_PATH}\nSet COMET_XENIUM_H5AD or edit the path above."
        )

    print("Resolving per-core H&E overlays...")
    images_by_section = resolve_section_images(H5AD_PATH)
    if not images_by_section:
        print("  Warning: no H&E overlays resolved — exporting without images.")

    print("Loading spatial data...")
    dataset = load_spatial_data(
        H5AD_PATH,
        section_key="sample_id",
        spatial_key="spatial",
        section_metadata=["sample_id", "tissue"],
    )
    print(f"  {dataset.n_sections} section(s), {dataset.n_cells:,} cells")
    print(f"  section IDs: {[s.section_id for s in dataset.sections]}")
    print(f"  detected modalities: {list(dataset.modalities.keys())}")

    # Make protein channels addressable by antibody target (CD45, not CD45 - TRITC).
    print("Cleaning protein channel names...")
    clean_protein_channel_names(dataset)

    # The exporter's Moran's I step needs a spatial graph in obsp (it won't build
    # one). Add a per-section kNN graph so spatially variable genes populate.
    if "spatial_connectivities" not in dataset.adata.obsp:
        print("Building spatial neighbour graph for Moran's I...")
        build_spatial_graph(dataset.adata, section_col="sample_id", spatial_key="spatial")

    # Attach each core's OWN H&E (section_id from section_key='sample_id' matches
    # the uns/section_images keys).
    section_images = {}
    for s in dataset.sections:
        layers = images_by_section.get(str(s.section_id))
        if layers:
            section_images[s.section_id] = layers
        else:
            print(f"  Warning: no H&E for section '{s.section_id}' — no overlay.")
    section_images = section_images or None

    common_kwargs = dict(
        main_cell_annotation=PRIMARY_ANNOTATION,
        title="COMET Multimodal Xenium (RNA + 16-plex Protein) — H&E",
        min_panel_size=100,
        spot_size="auto",
        outline_by=None,
        cell_annotations=ADDITIONAL_ANNOTATIONS,
        features=[],
        use_hvgs=False,
        # Extra (protein) modality requires sidecar storage.
        feature_storage="sidecar",
        feature_encoding="auto",
        feature_value_encoding="uint8",
        feature_sidecar_shard_size=128,
        # Cluster analytics (computed here — this is a plain, non-companion h5ad).
        marker_gene_annotations=CLUSTER_COLUMNS,
        marker_genes_top_n=30,
        neighbor_stats_annotations=CLUSTER_COLUMNS,
        neighbor_stats_permutations=0,
        pseudobulk_de_annotations=CLUSTER_COLUMNS,
        pseudobulk_de_top_n=20,
        pseudobulk_de_method="t-test",
        pseudobulk_de_layer=None,
        interaction_marker_annotations=None,
        # RNA panel + 16-channel COMET protein.
        modalities=["rna", "protein"],
        # H&E overlays.
        section_images=section_images,
        section_images_max_px=4096,
    )

    # Sidecar manifest is written next to the HTML (full path OK); the .karospace
    # packager needs a bare filename (the manifest lives inside the archive).
    print("\nExporting sidecar HTML...")
    export_to_html(
        dataset,
        output_path=SIDECAR_OUTPUT,
        feature_manifest_path=FEATURE_MANIFEST_PATH,
        **common_kwargs,
    )

    print(f"Packaging {PACKAGE_OUTPUT}...")
    export_to_html(
        dataset,
        output_path=PACKAGE_OUTPUT,
        feature_manifest_path=Path(FEATURE_MANIFEST_PATH).name,
        **common_kwargs,
    )

    print("\nDone! Wrote multimodal viewer:")
    print(f"  1. Sidecar Viewer:  {SIDECAR_OUTPUT}")
    print(f"  2. Feature Manifest:{FEATURE_MANIFEST_PATH}")
    print(f"  3. Shard Directory: {Path(FEATURE_MANIFEST_PATH).with_suffix('')}/")
    print(f"  4. KaroSpace Package:{PACKAGE_OUTPUT}")
    print(f"  5. Local Loader:    {Path(PACKAGE_OUTPUT).with_suffix('.loader.html')}")
    print(
        "\nOpen a core, switch the modality (RNA <-> Protein) to color by any of the "
        "16 antibodies, and use the H&E overlay controls to toggle / auto-align the stain."
    )


if __name__ == "__main__":
    main()
