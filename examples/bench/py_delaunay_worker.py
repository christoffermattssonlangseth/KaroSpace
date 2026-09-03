"""Python scipy Delaunay baseline — mirrors KaroSpaceCompanion's graph.rs.

Reads an h5ad, builds a Delaunay spatial graph with 99th-percentile long-link
trimming, writes obsp["spatial_connectivities"] (1.0 weights) and
obsp["spatial_distances"] (edge lengths), and saves a new h5ad.

End-to-end: read + build + write, matching what the Rust CLI does.

Usage: python py_delaunay_worker.py INPUT.h5ad OUTPUT.h5ad [PERCENTILE]
"""

import sys
import numpy as np
import scipy.sparse as sp
from scipy.spatial import Delaunay
import anndata as ad


def build_delaunay_graph(coords: np.ndarray, percentile: float):
    """Return (connectivities, distances) CSR matrices mirroring graph.rs."""
    n = coords.shape[0]
    tri = Delaunay(coords)
    simp = tri.simplices  # (M, 3)

    # Every triangle contributes its 3 undirected edges; sort endpoints so i<j.
    e = np.vstack([simp[:, [0, 1]], simp[:, [1, 2]], simp[:, [0, 2]]])
    e = np.sort(e, axis=1)
    e = np.unique(e, axis=0)

    d = np.linalg.norm(coords[e[:, 0]] - coords[e[:, 1]], axis=1)

    # Trim long links: threshold = <percentile> of positive edge lengths.
    # numpy's default linear interpolation matches graph.rs::percentile_value.
    pos = d[d > 0.0]
    if pos.size:
        thr = np.percentile(pos, percentile)
        keep = d <= thr
        e, d = e[keep], d[keep]

    # Symmetric CSR (both directions), like csr_from_neighbors.
    rows = np.concatenate([e[:, 0], e[:, 1]])
    cols = np.concatenate([e[:, 1], e[:, 0]])
    conn = sp.csr_matrix(
        (np.ones(rows.shape[0], dtype=np.float32), (rows, cols)), shape=(n, n)
    )
    dist = sp.csr_matrix(
        (np.concatenate([d, d]).astype(np.float32), (rows, cols)), shape=(n, n)
    )
    return conn, dist


def main():
    inp, out = sys.argv[1], sys.argv[2]
    percentile = float(sys.argv[3]) if len(sys.argv) > 3 else 99.0

    adata = ad.read_h5ad(inp)
    coords = np.asarray(adata.obsm["spatial"])[:, :2].astype(np.float64)

    conn, dist = build_delaunay_graph(coords, percentile)
    adata.obsp["spatial_connectivities"] = conn
    adata.obsp["spatial_distances"] = dist

    adata.write_h5ad(out)


if __name__ == "__main__":
    main()
