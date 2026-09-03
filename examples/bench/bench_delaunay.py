"""End-to-end Delaunay graph benchmark: KaroSpaceCompanion (Rust) vs scipy.

Both paths do the same end-to-end work: read an h5ad, build a Delaunay spatial
graph with 99th-percentile long-link trimming, write a new h5ad. Each is run as
a fresh subprocess so the wall clock reflects the real "time to enriched file"
(the Python side pays interpreter + import startup, which is part of that cost).

After timing, both output graphs are compared by nnz for parity — if the edge
counts diverge a lot, the two aren't doing the same work and the times aren't
comparable.
"""

import os
import subprocess
import sys
import time
from statistics import median

import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
WORKER = os.path.join(HERE, "py_delaunay_worker.py")
RUST_BIN = os.path.expanduser(
    "~/work/karolinska_institutet/projects/KaroSpaceCompanion/target/release/karospace-companion"
)
PERCENTILE = 99.0
REPEATS = 3
TMPDIR = os.path.join(HERE, "_bench_tmp")

DATASETS = [
    ("VisiumHD_96k", os.path.expanduser("~/Downloads/Human_VisiumHD_compressed_v1.h5ad")),
    ("MERFISH_401k", os.path.expanduser("~/Downloads/GSE284005_merfish_all.h5ad")),
    (
        "XeniumPup_1.36M",
        os.path.expanduser(
            "~/work/karolinska_institutet/projects/xenium_mouse_embryo/derived_scanpy/xenium_mouse_pup_gene_only.h5ad"
        ),
    ),
]


def time_cmd(cmd, repeats):
    """Run cmd `repeats` times (after one warmup), return list of wall times."""
    subprocess.run(cmd, check=True, capture_output=True)  # warmup
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        subprocess.run(cmd, check=True, capture_output=True)
        times.append(time.perf_counter() - t0)
    return times


def graph_nnz(path):
    a = ad.read_h5ad(path, backed="r")
    conn = a.obsp["spatial_connectivities"]
    return int(conn.nnz)


def main():
    os.makedirs(TMPDIR, exist_ok=True)
    print(f"{'dataset':>16} | {'rust (s)':>9} | {'python (s)':>10} | "
          f"{'speedup':>7} | {'rust edges':>11} | {'py edges':>10} | {'Δ%':>6}")
    print("-" * 92)

    for name, path in DATASETS:
        if not os.path.exists(path):
            print(f"{name:>16} | MISSING: {path}")
            continue

        rust_out = os.path.join(TMPDIR, f"{name}.rust.h5ad")
        py_out = os.path.join(TMPDIR, f"{name}.py.h5ad")

        rust_cmd = [
            RUST_BIN, "prepare", path, "--output", rust_out,
            "--delaunay", "--remove-long-links-percentile", str(PERCENTILE),
            "--skip-normalized-layer", "--skip-aggregation",
        ]
        py_cmd = [sys.executable, WORKER, path, py_out, str(PERCENTILE)]

        rust_t = time_cmd(rust_cmd, REPEATS)
        py_t = time_cmd(py_cmd, REPEATS)

        rn, pn = graph_nnz(rust_out), graph_nnz(py_out)
        delta = 100.0 * (pn - rn) / rn if rn else float("nan")
        rmed, pmed = median(rust_t), median(py_t)
        speedup = pmed / rmed if rmed else float("nan")

        print(f"{name:>16} | {rmed:9.2f} | {pmed:10.2f} | {speedup:6.2f}x | "
              f"{rn:>11,} | {pn:>10,} | {delta:+5.1f}%")

    print("\nnnz = 2 x undirected edges (symmetric). Δ% is python vs rust edge count.")


if __name__ == "__main__":
    main()
