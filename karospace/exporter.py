"""
Export spatial data to standalone HTML viewer.

Creates self-contained HTML files with embedded data and JavaScript
for interactive visualization of spatial transcriptomics data.
"""

import base64
import json
import os
import re
import shutil
import struct
import tempfile
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Mapping, Optional, List, Tuple, Union
import numpy as np
import pandas as pd
try:
    from tqdm.auto import tqdm
except Exception:  # pragma: no cover - optional dependency fallback
    tqdm = None

from .data_loader import SpatialDataset

GENE_SIDECAR_SHARD_SIZE = 256
KAROSPACE_PACKAGE_MANIFEST = "karospace-package.json"
KAROSPACE_PACKAGE_LOADER_FILENAME = "karospace-package-loader.html"
GENE_SIDECAR_FORMAT_JSON_V2 = "json-v2"
GENE_SIDECAR_FORMAT_BINARY_V1 = "binary-v1"
GENE_SIDECAR_BINARY_MAGIC = b"KSB1"
GENE_SIDECAR_BINARY_VERSION = 1


def _load_logo_base64() -> Optional[str]:
    """Load logo from assets as base64 string."""
    logo_path = Path(__file__).parent.parent / "assets" / "logo.png"
    if logo_path.exists():
        with open(logo_path, "rb") as f:
            return base64.b64encode(f.read()).decode("utf-8")
    return None


def _chunked(values: List[str], size: int) -> List[List[str]]:
    if size < 1:
        raise ValueError("chunk size must be >= 1")
    return [values[i:i + size] for i in range(0, len(values), size)]


def _isoformat_utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _get_karospace_version() -> str:
    try:
        from . import __version__

        return str(__version__)
    except Exception:
        return "unknown"


def _guess_package_media_type(path: Union[str, Path]) -> str:
    suffix = Path(path).suffix.lower()
    if suffix == ".html":
        return "text/html"
    if suffix == ".json":
        return "application/json"
    if suffix == ".bin":
        return "application/octet-stream"
    return "application/octet-stream"


def _validate_gene_sidecar_format(gene_sidecar_format: str) -> str:
    normalized = str(gene_sidecar_format or GENE_SIDECAR_FORMAT_JSON_V2).strip().lower()
    if normalized not in {GENE_SIDECAR_FORMAT_JSON_V2, GENE_SIDECAR_FORMAT_BINARY_V1}:
        raise ValueError(
            "gene_sidecar_format must be one of: "
            f"'{GENE_SIDECAR_FORMAT_JSON_V2}', '{GENE_SIDECAR_FORMAT_BINARY_V1}'"
        )
    return normalized


def _u32_bytes(values: List[int]) -> bytes:
    if not values:
        return b""
    arr = np.asarray(values, dtype="<u4")
    return np.ascontiguousarray(arr).tobytes(order="C")


def _decode_sidecar_packed_bytes(value, *, dtype: str) -> bytes:
    if value is None:
        return b""
    if isinstance(value, str):
        return base64.b64decode(value.encode("ascii"))
    if isinstance(value, (bytes, bytearray, memoryview)):
        return bytes(value)
    arr = np.asarray(value, dtype=dtype)
    return np.ascontiguousarray(arr).tobytes(order="C")


def _binary_payload_kind_from_sections(section_entries: Mapping[str, dict]) -> int:
    kinds = set()
    for section_entry in section_entries.values():
        if not section_entry:
            continue
        if "sparse" in section_entry:
            kinds.add("sparse")
        else:
            kinds.add("dense")
    if kinds == {"dense"}:
        return 1
    if kinds == {"sparse"}:
        return 2
    return 0


def _build_binary_gene_payload(
    *,
    gene_entry: Mapping[str, object],
    section_order: List[str],
    section_cell_counts: Mapping[str, int],
) -> Tuple[int, bytes]:
    section_entries = gene_entry.get("sections") or {}
    payload_kind = _binary_payload_kind_from_sections(section_entries)
    chunks: List[bytes] = [struct.pack("<I", len(section_order))]

    for section_index, section_id in enumerate(section_order):
        section_entry = section_entries.get(section_id) or {}
        cell_count = int(section_cell_counts.get(section_id, 0))
        section_encoding = 0
        nnz = 0
        nan_values: List[int] = []
        section_data = b""

        if "sparse" in section_entry:
            sparse = section_entry.get("sparse") or {}
            index_bytes = _decode_sidecar_packed_bytes(sparse.get("ib64"), dtype="<u4")
            if "vq8b64" in sparse:
                value_bytes = _decode_sidecar_packed_bytes(sparse.get("vq8b64"), dtype="u1")
                section_encoding = 2
            elif "vq16b64" in sparse:
                value_bytes = _decode_sidecar_packed_bytes(sparse.get("vq16b64"), dtype="<u2")
                section_encoding = 4
            else:
                raise ValueError("binary-v1 requires quantized sparse sidecar values")
            nnz = len(index_bytes) // 4
            nan_values = [int(v) for v in (sparse.get("nan") or [])]
            section_data = index_bytes + value_bytes
        elif section_entry:
            if "dq8b64" in section_entry:
                value_bytes = _decode_sidecar_packed_bytes(section_entry.get("dq8b64"), dtype="u1")
                section_encoding = 1
            elif "dq16b64" in section_entry:
                value_bytes = _decode_sidecar_packed_bytes(section_entry.get("dq16b64"), dtype="<u2")
                section_encoding = 3
            else:
                raise ValueError("binary-v1 requires quantized dense sidecar values")
            nan_values = [int(v) for v in (section_entry.get("nan") or [])]
            section_data = value_bytes

        chunks.append(
            struct.pack(
                "<HBBIII",
                int(section_index),
                int(section_encoding),
                0,
                int(cell_count),
                int(nnz),
                int(len(nan_values)),
            )
        )
        chunks.append(section_data)
        chunks.append(_u32_bytes(nan_values))

    return payload_kind, b"".join(chunks)


def _write_binary_gene_shard(
    *,
    shard_path: Path,
    shard_genes: List[str],
    shard_data: Mapping[str, object],
    section_order: List[str],
    section_cell_counts: Mapping[str, int],
) -> None:
    genes_payload = shard_data.get("genes") or {}
    payload_blobs: List[Tuple[bytes, int, bytes]] = []
    index_size = 0

    for gene in shard_genes:
        gene_name_bytes = str(gene).encode("utf-8")
        payload_kind, payload_bytes = _build_binary_gene_payload(
            gene_entry=genes_payload.get(gene) or {},
            section_order=section_order,
            section_cell_counts=section_cell_counts,
        )
        payload_blobs.append((gene_name_bytes, payload_kind, payload_bytes))
        index_size += 2 + len(gene_name_bytes) + 1 + 1 + 8 + 8

    header = struct.pack(
        "<4sHHII",
        GENE_SIDECAR_BINARY_MAGIC,
        GENE_SIDECAR_BINARY_VERSION,
        0,
        len(payload_blobs),
        0,
    )
    payload_offset = len(header) + index_size
    out = bytearray()
    out.extend(header)

    current_offset = payload_offset
    for gene_name_bytes, payload_kind, payload_bytes in payload_blobs:
        out.extend(struct.pack("<H", len(gene_name_bytes)))
        out.extend(gene_name_bytes)
        out.extend(struct.pack("<BBQQ", payload_kind, 0, current_offset, len(payload_bytes)))
        current_offset += len(payload_bytes)

    for _, _, payload_bytes in payload_blobs:
        out.extend(payload_bytes)

    shard_path.parent.mkdir(parents=True, exist_ok=True)
    with open(shard_path, "wb") as handle:
        handle.write(out)


def _build_karospace_package_manifest(
    *,
    source_root: Path,
    entry_html: str,
    gene_manifest_path: str,
    gene_shard_dir: str,
    title: str,
    n_sections: int,
    total_cells: int,
) -> dict:
    files = {}
    for file_path in sorted(p for p in source_root.rglob("*") if p.is_file()):
        rel_path = file_path.relative_to(source_root).as_posix()
        files[rel_path] = {
            "media_type": _guess_package_media_type(rel_path),
            "size_bytes": int(file_path.stat().st_size),
        }

    return {
        "format": "karospace-package-v1",
        "package_version": 1,
        "entry_html": entry_html,
        "created_at": _isoformat_utc_now(),
        "producer": {
            "name": "karospace",
            "version": _get_karospace_version(),
        },
        "title": title,
        "n_sections": int(n_sections),
        "total_cells": int(total_cells),
        "viewer": {
            "mode": "sidecar-package",
            "gene_storage": "sidecar",
            "gene_manifest_path": gene_manifest_path,
            "gene_shard_dir": gene_shard_dir,
        },
        "files": files,
    }


def _write_karospace_package(
    *,
    package_path: Path,
    source_root: Path,
    entry_html: str,
    gene_manifest_path: str,
    gene_shard_dir: str,
    title: str,
    n_sections: int,
    total_cells: int,
) -> None:
    if not (source_root / entry_html).exists():
        raise FileNotFoundError(f"package entry HTML missing: {entry_html}")
    if not (source_root / gene_manifest_path).exists():
        raise FileNotFoundError(f"package gene manifest missing: {gene_manifest_path}")
    if not (source_root / gene_shard_dir).exists():
        raise FileNotFoundError(f"package gene shard directory missing: {gene_shard_dir}")

    package_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = _build_karospace_package_manifest(
        source_root=source_root,
        entry_html=entry_html,
        gene_manifest_path=gene_manifest_path,
        gene_shard_dir=gene_shard_dir,
        title=title,
        n_sections=n_sections,
        total_cells=total_cells,
    )

    def _package_sort_key(path: Path) -> tuple[int, str]:
        rel_path = path.relative_to(source_root).as_posix()
        if rel_path == entry_html:
            return (0, rel_path)
        if rel_path == gene_manifest_path:
            return (1, rel_path)
        if rel_path.startswith(f"{gene_shard_dir.rstrip('/')}/"):
            return (3, rel_path)
        return (2, rel_path)

    with zipfile.ZipFile(package_path, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(
            KAROSPACE_PACKAGE_MANIFEST,
            json.dumps(manifest, separators=(",", ":")),
            compress_type=zipfile.ZIP_DEFLATED,
        )
        for file_path in sorted((p for p in source_root.rglob("*") if p.is_file()), key=_package_sort_key):
            rel_path = file_path.relative_to(source_root).as_posix()
            compress_type = zipfile.ZIP_STORED if file_path.suffix.lower() == ".bin" else zipfile.ZIP_DEFLATED
            zf.write(file_path, arcname=rel_path, compress_type=compress_type)


def _find_karospace_package_loader_template() -> Optional[Path]:
    candidates = [
        Path(__file__).resolve().with_name("package_loader.html"),
        Path(__file__).resolve().parent.parent / KAROSPACE_PACKAGE_LOADER_FILENAME,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _write_karospace_package_loader(
    *,
    loader_path: Path,
) -> Optional[Path]:
    template_path = _find_karospace_package_loader_template()
    if template_path is None:
        return None
    loader_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(template_path, loader_path)
    return loader_path


def _extract_embedded_viewer_data(html_text: str) -> dict:
    match = re.search(
        r'<script id="karospace-data" type="application/json">(.*?)</script>',
        html_text,
        re.DOTALL,
    )
    if not match:
        raise ValueError("embedded karospace-data script not found in HTML")
    return json.loads(match.group(1).replace("<\\/", "</"))


def _extract_html_title(html_text: str) -> str:
    match = re.search(r"<title>(.*?)</title>", html_text, re.DOTALL | re.IGNORECASE)
    if not match:
        return "KaroSpace"
    return str(match.group(1)).strip() or "KaroSpace"


def package_sidecar_viewer(
    html_path: Union[str, Path],
    *,
    output_path: Optional[Union[str, Path]] = None,
    gene_manifest_path: Optional[Union[str, Path]] = None,
    gene_shard_dir: Optional[Union[str, Path]] = None,
    loader_output_path: Optional[Union[str, Path]] = None,
) -> str:
    """
    Package an existing sidecar viewer bundle into a `.karospace` archive.

    Parameters
    ----------
    html_path : str or Path
        Path to an existing KaroSpace HTML viewer generated with `gene_storage="sidecar"`.
    output_path : str or Path, optional
        Output `.karospace` path. Defaults to `<html stem>.karospace`.
    gene_manifest_path : str or Path, optional
        Actual filesystem path to the sidecar gene manifest JSON if it differs from the
        path referenced in the HTML.
    gene_shard_dir : str or Path, optional
        Actual filesystem path to the sidecar shard directory if it differs from the
        manifest stem directory.
    loader_output_path : str or Path, optional
        Optional output path for the companion `.loader.html` file. Defaults to
        `<output stem>.loader.html`.

    Returns
    -------
    str
        Absolute path to the written `.karospace` package.
    """
    source_html_path = Path(html_path).expanduser().resolve()
    if not source_html_path.exists():
        raise FileNotFoundError(f"sidecar HTML not found: {source_html_path}")

    print(f"  - reading existing viewer HTML: {source_html_path}")
    html_text = source_html_path.read_text(encoding="utf-8")
    data = _extract_embedded_viewer_data(html_text)
    package_title = _extract_html_title(html_text)

    gene_aux_url = str(data.get("gene_aux_url") or "").strip()
    if not gene_aux_url:
        raise ValueError("HTML viewer does not reference a sidecar gene manifest")

    html_parent = source_html_path.parent
    package_gene_manifest_rel = Path(gene_aux_url).as_posix()
    actual_gene_manifest_path = (
        Path(gene_manifest_path).expanduser().resolve()
        if gene_manifest_path is not None
        else (html_parent / package_gene_manifest_rel).resolve()
    )
    if not actual_gene_manifest_path.exists():
        raise FileNotFoundError(f"sidecar gene manifest not found: {actual_gene_manifest_path}")

    package_gene_shard_rel = Path(package_gene_manifest_rel).with_suffix("").as_posix()
    actual_gene_shard_dir = (
        Path(gene_shard_dir).expanduser().resolve()
        if gene_shard_dir is not None
        else (html_parent / package_gene_shard_rel).resolve()
    )
    if not actual_gene_shard_dir.exists():
        raise FileNotFoundError(f"sidecar shard directory not found: {actual_gene_shard_dir}")

    resolved_output_path = (
        Path(output_path).expanduser().resolve()
        if output_path is not None
        else source_html_path.with_suffix(".karospace")
    )
    if resolved_output_path.suffix.lower() != ".karospace":
        raise ValueError("output_path must end with .karospace")

    resolved_loader_output_path = (
        Path(loader_output_path).expanduser().resolve()
        if loader_output_path is not None
        else resolved_output_path.with_suffix(".loader.html")
    )

    n_sections = int(data.get("n_sections") or 0)
    total_cells = int(data.get("total_cells") or 0)

    print(f"  - resolved gene manifest: {actual_gene_manifest_path}")
    print(f"  - resolved shard directory: {actual_gene_shard_dir}")
    print(f"  - package title: {package_title}")
    print(f"  - package target: {resolved_output_path}")
    print(
        f"  - packaging {n_sections} sections and {total_cells:,} cells "
        f"from existing sidecar assets"
    )

    with tempfile.TemporaryDirectory(prefix="karospace-package-") as tmpdir:
        source_root = Path(tmpdir)
        print("  - staging package assets...")
        entry_html = "index.html"
        staged_entry_html = source_root / entry_html
        staged_entry_html.write_text(html_text, encoding="utf-8")

        staged_gene_manifest = source_root / package_gene_manifest_rel
        staged_gene_manifest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(actual_gene_manifest_path, staged_gene_manifest)

        staged_gene_shard_dir = source_root / package_gene_shard_rel
        staged_gene_shard_dir.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(actual_gene_shard_dir, staged_gene_shard_dir, dirs_exist_ok=True)

        print("  - writing .karospace archive...")
        _write_karospace_package(
            package_path=resolved_output_path,
            source_root=source_root,
            entry_html=entry_html,
            gene_manifest_path=package_gene_manifest_rel,
            gene_shard_dir=package_gene_shard_rel,
            title=package_title,
            n_sections=n_sections,
            total_cells=total_cells,
        )

    print("  - writing local package loader...")
    written_loader_path = _write_karospace_package_loader(
        loader_path=resolved_loader_output_path,
    )

    print(f"Exported KaroSpace package to: {resolved_output_path}")
    print(f"  - entry HTML: index.html (from {source_html_path.name})")
    print(f"  - packaged gene manifest: {package_gene_manifest_rel}")
    print(f"  - packaged shard directory: {package_gene_shard_rel}")
    if written_loader_path is not None:
        print(f"  - local package loader: {written_loader_path}")
    else:
        print(f"  - local package loader: unavailable (template not found)")
    print(f"  - {n_sections} sections")
    print(f"  - {total_cells:,} cells")

    return str(resolved_output_path)


def _select_top_variable_genes(adata, n: int) -> List[str]:
    """Return up to n gene names ranked by expression variance across cells."""
    X = adata.X
    if hasattr(X, "toarray"):
        # For sparse matrices compute variance without densifying the whole thing.
        # Var(X) = E[X^2] - E[X]^2
        mean_sq = np.asarray(X.power(2).mean(axis=0)).ravel()
        sq_mean = np.asarray(X.mean(axis=0)).ravel() ** 2
        variances = mean_sq - sq_mean
    else:
        variances = np.var(np.asarray(X, dtype=float), axis=0)
    top_idx = np.argsort(variances)[::-1][:n]
    return [adata.var_names[i] for i in top_idx]


def _compute_cluster_gene_means(adata, genes: List[str], obs_col: str) -> Optional[dict]:
    """Compute per-cluster and overall mean expression for each gene."""
    if obs_col not in adata.obs.columns:
        return None
    col = adata.obs[obs_col]
    if not hasattr(col, "cat"):
        try:
            col = col.astype("category")
        except Exception:
            return None
    categories = [str(c) for c in col.cat.categories]
    if not categories:
        return None

    X = adata[:, genes].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=float)

    means = {}
    col_vals = col.to_numpy()
    for cat in categories:
        mask = col_vals == cat
        if mask.sum() == 0:
            means[cat] = [0.0] * len(genes)
        else:
            means[cat] = np.round(X[mask].mean(axis=0), 4).tolist()

    background = np.round(X.mean(axis=0), 4).tolist()
    return {"categories": categories, "means": means, "background": background}


def _compute_morans_i(adata, genes: List[str], n_genes: int = 200) -> list:
    """Compute Moran's I spatial autocorrelation for the top variable genes.

    Returns a list of {gene, I} dicts sorted descending by I value, or [] if
    no spatial weight matrix is found in adata.obsp.

    W is kept sparse throughout. All genes are processed in a single sparse
    matmul (W @ Z), so cost is O(G * nnz) rather than O(G * N^2).
    """
    import scipy.sparse as sp

    # Find spatial weight matrix — keep sparse
    W = None
    for key in ("spatial_connectivities", "connectivities", "neighbors", "neighbor_graph"):
        if key in adata.obsp:
            W = adata.obsp[key]
            break
    if W is None:
        return []

    if not sp.issparse(W):
        W = sp.csr_matrix(W)
    else:
        W = W.tocsr().astype(float)

    N = W.shape[0]

    # Row-normalize in place using sparse diagonal scaling
    row_sums = np.asarray(W.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0
    W = sp.diags(1.0 / row_sums) @ W
    S0 = float(W.sum())
    if S0 == 0:
        return []

    # Select top variable genes and densify only those (N × G, manageable)
    selected = _select_top_variable_genes(adata, n_genes)
    selected = [g for g in selected if g in set(genes)]
    if not selected:
        return []

    X = adata[:, selected].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=float)          # (N, G)

    # Center each gene
    Z = X - X.mean(axis=0)                  # (N, G)

    # One sparse matmul for all genes: W @ Z  →  (N, G)
    WZ = W.dot(Z)                            # sparse × dense = dense (N, G)

    num = N * (Z * WZ).sum(axis=0)          # (G,)
    denom = S0 * (Z * Z).sum(axis=0)        # (G,)

    valid = denom > 0
    I_vals = np.where(valid, num / np.where(valid, denom, 1.0), 0.0)
    I_vals = np.clip(I_vals, -1.0, 1.0)

    results = [
        {"gene": gene, "I": round(float(I_vals[j]), 4)}
        for j, gene in enumerate(selected)
        if valid[j]
    ]
    results.sort(key=lambda d: d["I"], reverse=True)
    return results


def _compute_gene_correlations(adata, genes: List[str], top_n: int = 10) -> dict:
    """For each gene, return the top_n most positively correlated genes."""
    if not genes or top_n < 1 or len(genes) < 2:
        return {gene: [] for gene in genes} if genes else {}
    X = adata[:, genes].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=float)
    corr = np.corrcoef(X.T)  # shape (n_genes, n_genes)
    result = {}
    for i, gene in enumerate(genes):
        scores = corr[i].copy()
        scores[i] = -2.0  # exclude self
        top_idx = np.argsort(scores)[::-1][:top_n]
        result[gene] = [
            {"gene": genes[j], "r": round(float(corr[i, j]), 3)}
            for j in top_idx
            if j != i and corr[i, j] > 0
        ]
    return result


def _estimate_auto_spot_size(dataset: SpatialDataset, min_panel_size: int) -> float:
    """Estimate a reasonable default spot radius from section density."""
    panel_px = max(48.0, float(min_panel_size) - 16.0)  # Account for canvas padding.
    spacing_px_samples: List[float] = []

    for section in dataset.sections:
        coords = np.asarray(section.coordinates)
        if coords.ndim != 2 or coords.shape[0] < 2 or coords.shape[1] < 2:
            continue

        x = coords[:, 0]
        y = coords[:, 1]
        valid = np.isfinite(x) & np.isfinite(y)
        n = int(valid.sum())
        if n < 2:
            continue

        x = x[valid]
        y = y[valid]
        dx = float(np.max(x) - np.min(x))
        dy = float(np.max(y) - np.min(y))
        extent = max(dx, dy, 1e-12)

        # Use observed bbox area, but guard against near-line geometries.
        area = max(dx * dy, extent * extent * 0.05, 1e-12)
        spacing_data = float(np.sqrt(area / n))
        spacing_px = spacing_data * (panel_px / extent)
        if np.isfinite(spacing_px) and spacing_px > 0:
            spacing_px_samples.append(spacing_px)

    if not spacing_px_samples:
        return 2.0

    # Target a radius that keeps neighboring cells separable without looking sparse.
    est = 0.38 * float(np.median(spacing_px_samples))
    return float(np.clip(est, 0.35, 4.0))


def _resolve_spot_size(
    dataset: SpatialDataset,
    spot_size: Union[float, str, None],
    min_panel_size: int,
) -> Tuple[float, bool]:
    """
    Resolve spot size input to a numeric value.

    Returns
    -------
    Tuple[float, bool]
        (resolved_size, used_auto_mode)
    """
    if spot_size is None:
        return _estimate_auto_spot_size(dataset, min_panel_size), True

    if isinstance(spot_size, str):
        token = spot_size.strip().lower()
        if token in {"auto", "adaptive", "density"}:
            return _estimate_auto_spot_size(dataset, min_panel_size), True
        try:
            value = float(token)
        except ValueError as exc:
            raise ValueError("spot_size must be a positive number or 'auto'.") from exc
    else:
        value = float(spot_size)

    if not np.isfinite(value) or value <= 0:
        raise ValueError("spot_size must be a positive finite number or 'auto'.")
    return float(value), False


def _normalize_section_rotation(value: Union[int, float]) -> float:
    angle = float(value)
    if not np.isfinite(angle):
        raise ValueError("section rotation angles must be finite numbers")
    normalized = angle % 360.0
    if np.isclose(normalized, 360.0):
        normalized = 0.0
    return float(normalized)


def _resolve_section_rotations(
    dataset: SpatialDataset,
    section_rotations: Optional[Mapping[str, Union[int, float]]],
) -> Dict[str, float]:
    if section_rotations is None:
        return {}
    if not isinstance(section_rotations, Mapping):
        raise TypeError("section_rotations must be a mapping of section_id -> angle")

    valid_ids = {section.section_id for section in dataset.sections}
    resolved: Dict[str, float] = {}
    unknown_ids = sorted(str(section_id) for section_id in section_rotations if str(section_id) not in valid_ids)
    if unknown_ids:
        raise ValueError(
            "section_rotations contains unknown section_id values: " + ", ".join(unknown_ids)
        )

    for raw_section_id, raw_angle in section_rotations.items():
        section_id = str(raw_section_id)
        try:
            resolved[section_id] = _normalize_section_rotation(raw_angle)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"section_rotations[{section_id!r}] must be a finite numeric angle"
            ) from exc

    return resolved


def _resolve_marker_genes_groupby(
    dataset: SpatialDataset,
    color: str,
    marker_genes_groupby: Optional[List[str]],
) -> Optional[List[str]]:
    if not marker_genes_groupby:
        return marker_genes_groupby
    resolved = [str(value) for value in marker_genes_groupby if str(value).strip()]
    if color not in dataset.adata.obs.columns:
        return resolved or None
    if color in resolved:
        return resolved or None
    if pd.api.types.is_numeric_dtype(dataset.adata.obs[color]):
        return resolved or None
    return [str(color), *resolved]


# Default color palettes
DEFAULT_CATEGORICAL_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    "#393b79", "#5254a3", "#6b6ecf", "#9c9ede", "#637939",
    "#8ca252", "#b5cf6b", "#cedb9c", "#8c6d31", "#bd9e39",
    "#e7ba52", "#e7cb94", "#843c39", "#ad494a", "#d6616b",
    "#e7969c", "#7b4173", "#a55194", "#ce6dbd", "#de9ed6",
]

HTML_TEMPLATE = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    {favicon_link}
    <style>
        :root {{
            --background: {background};
            --text-color: {text_color};
            --header-bg: {header_bg};
            --panel-bg: {panel_bg};
            --border-color: {border_color};
            --input-bg: {input_bg};
            --muted-color: {muted_color};
            --hover-bg: {hover_bg};
            --graph-color: {graph_color};
            --accent: #870052;
            --accent-strong: #4F0433;
            --accent-warm: #FF876F;
            --accent-soft: #FEEEEB;
            --accent-cool: #EDF4F4;
        }}
        :root.dark {{
            --background: #1a1a1a;
            --text-color: #e0e0e0;
            --header-bg: #2a2a2a;
            --panel-bg: #2a2a2a;
            --border-color: #404040;
            --input-bg: #333333;
            --muted-color: #888888;
            --hover-bg: #3a3a3a;
            --graph-color: rgba(255, 255, 255, 0.12);
        }}
        :root.light {{
            --background: #f5f5f5;
            --text-color: #1a1a1a;
            --header-bg: #ffffff;
            --panel-bg: #ffffff;
            --border-color: #e0e0e0;
            --input-bg: #ffffff;
            --muted-color: #666666;
            --hover-bg: #f0f0f0;
            --graph-color: rgba(0, 0, 0, 0.12);
        }}
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            background:
                radial-gradient(800px 500px at 10% 0%, rgba(255, 135, 111, 0.08), rgba(0, 0, 0, 0)),
                radial-gradient(900px 600px at 100% 20%, rgba(135, 0, 82, 0.08), rgba(0, 0, 0, 0)),
                var(--background);
            color: var(--text-color);
            min-height: 100vh;
            display: flex;
            flex-direction: column;
            transition: background 0.3s, color 0.3s;
        }}
        .header {{
            padding: 8px 16px;
            background:
                linear-gradient(90deg, rgba(255, 135, 111, 0.12), rgba(135, 0, 82, 0.08)),
                var(--header-bg);
            border-bottom: 1px solid var(--border-color);
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
            gap: 8px;
            transition: background 0.3s, border-color 0.3s;
        }}
        .header-title {{
            display: flex;
            align-items: center;
            gap: 8px;
            position: relative;
        }}
        .header h1 {{ font-size: 16px; font-weight: 600; }}
        .info-trigger {{
            padding: 4px 8px;
            border: 1px solid var(--border-color);
            border-radius: 999px;
            background: var(--panel-bg);
            color: var(--text-color);
            font-size: 10px;
            letter-spacing: 0.03em;
            text-transform: uppercase;
            cursor: pointer;
            transition: background 0.2s, border-color 0.2s, color 0.2s;
        }}
        .info-trigger:hover {{
            background: var(--hover-bg);
            border-color: var(--accent-strong);
        }}
        .controls {{ display: flex; align-items: center; gap: 8px; flex-wrap: wrap; }}
        .control-group {{ display: flex; align-items: center; gap: 4px; }}
        #expression-scale-section {{
            flex-direction: column;
            align-items: flex-start;
        }}
        .control-group label {{ font-size: 11px; color: var(--muted-color); }}
        select, input[type="text"] {{
            padding: 5px 8px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--text-color);
            font-size: 12px;
            transition: background 0.3s, border-color 0.3s, color 0.3s;
        }}
        select {{ min-width: 120px; }}
        input[type="text"] {{ width: 140px; }}
        select:focus, input:focus {{ outline: none; border-color: var(--accent-strong); box-shadow: 0 0 0 2px rgba(135, 0, 82, 0.15); }}
        .gene-control-group {{
            position: relative;
            align-items: flex-start;
        }}
        .gene-input-shell {{
            position: relative;
            min-width: 180px;
        }}
        .gene-input-shell input[type="text"] {{
            width: 180px;
        }}
        .gene-discovery-panel {{
            position: absolute;
            top: calc(100% + 6px);
            left: 0;
            width: min(360px, calc(100vw - 32px));
            max-height: 420px;
            overflow: auto;
            padding: 10px;
            border: 1px solid var(--border-color);
            border-radius: 10px;
            background:
                linear-gradient(180deg, rgba(255, 135, 111, 0.08), rgba(135, 0, 82, 0.03)),
                var(--panel-bg);
            box-shadow: 0 14px 36px rgba(0, 0, 0, 0.16);
            display: none;
            z-index: 40;
        }}
        .gene-discovery-panel.active {{
            display: block;
        }}
        .gene-discovery-header {{
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 8px;
            margin-bottom: 8px;
        }}
        .gene-discovery-title {{
            font-size: 11px;
            font-weight: 600;
            letter-spacing: 0.04em;
            text-transform: uppercase;
            color: var(--muted-color);
        }}
        .gene-discovery-actions {{
            display: flex;
            align-items: center;
            gap: 6px;
        }}
        .gene-discovery-section + .gene-discovery-section {{
            margin-top: 10px;
            padding-top: 10px;
            border-top: 1px solid rgba(127, 127, 127, 0.18);
        }}
        .gene-discovery-label {{
            font-size: 11px;
            font-weight: 600;
            margin-bottom: 6px;
            color: var(--text-color);
        }}
        .gene-discovery-empty {{
            font-size: 11px;
            color: var(--muted-color);
            line-height: 1.4;
        }}
        .gene-token-grid {{
            display: flex;
            flex-wrap: wrap;
            gap: 6px;
        }}
        .gene-token-btn {{
            display: inline-flex;
            align-items: center;
            gap: 6px;
            padding: 5px 8px;
            border: 1px solid var(--border-color);
            border-radius: 999px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            font-size: 11px;
            line-height: 1.2;
            transition: background 0.2s, border-color 0.2s, transform 0.2s;
        }}
        .gene-token-btn:hover {{
            background: var(--hover-bg);
            border-color: var(--accent-strong);
            transform: translateY(-1px);
        }}
        .gene-token-btn.active {{
            border-color: var(--accent-strong);
            background: rgba(135, 0, 82, 0.10);
        }}
        .gene-token-btn.search-active {{
            box-shadow: 0 0 0 2px rgba(135, 0, 82, 0.18);
        }}
        .gene-token-btn.unloaded {{
            border-style: dashed;
        }}
        .gene-token-btn.disabled {{
            cursor: default;
            opacity: 0.78;
        }}
        .gene-token-btn.disabled:hover {{
            background: var(--input-bg);
            border-color: var(--border-color);
            transform: none;
        }}
        .gene-token-meta {{
            font-size: 10px;
            color: var(--muted-color);
        }}
        .gene-suggestion-group + .gene-suggestion-group {{
            margin-top: 8px;
        }}
        .gene-suggestion-group-title {{
            font-size: 10px;
            color: var(--muted-color);
            margin-bottom: 5px;
            text-transform: uppercase;
            letter-spacing: 0.03em;
        }}
        .gene-panel-card {{
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 8px;
            background: rgba(255, 255, 255, 0.04);
        }}
        .gene-panel-card + .gene-panel-card {{
            margin-top: 8px;
        }}
        .gene-panel-header {{
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 8px;
            margin-bottom: 6px;
        }}
        .gene-panel-name {{
            font-size: 12px;
            font-weight: 600;
        }}
        .gene-panel-actions {{
            display: flex;
            align-items: center;
            gap: 6px;
        }}
        .gene-panel-btn {{
            padding: 3px 7px;
            border: 1px solid var(--border-color);
            border-radius: 999px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            font-size: 10px;
            transition: background 0.2s, border-color 0.2s;
        }}
        .gene-panel-btn:hover {{
            background: var(--hover-bg);
        }}
        .gene-panel-genes {{
            display: flex;
            flex-wrap: wrap;
            gap: 6px;
        }}
        .stats {{ font-size: 11px; color: var(--muted-color); }}
        .size-control {{
            display: inline-flex;
            align-items: center;
            gap: 6px;
        }}
        .size-step {{
            width: 18px;
            height: 18px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            font-size: 12px;
            line-height: 1;
            display: inline-flex;
            align-items: center;
            justify-content: center;
            transition: background 0.2s, border-color 0.2s;
        }}
        .size-step:hover {{ background: var(--hover-bg); }}

        /* Theme toggle button */
        .theme-toggle {{
            background: var(--input-bg);
            border: 1px solid var(--border-color);
            border-radius: 4px;
            padding: 5px 10px;
            cursor: pointer;
            font-size: 14px;
            transition: background 0.3s, border-color 0.3s;
        }}
        .theme-toggle:hover {{ background: var(--hover-bg); }}
        .export-btn {{
            background: var(--input-bg);
            border: 1px solid var(--border-color);
            border-radius: 4px;
            padding: 5px 10px;
            cursor: pointer;
            font-size: 12px;
            transition: background 0.3s, border-color 0.3s;
        }}
        .export-btn:hover {{ background: var(--hover-bg); }}

        /* Filter bar */
        .filter-bar {{
            padding: 6px 16px;
            background: var(--header-bg);
            border-bottom: 1px solid var(--border-color);
            display: flex;
            align-items: center;
            gap: 12px;
            flex-wrap: wrap;
            transition: background 0.3s, border-color 0.3s;
        }}
        .filter-bar:empty {{ display: none; }}
        .filter-group {{ display: flex; align-items: center; gap: 4px; }}
        .filter-group label {{ font-size: 10px; color: var(--muted-color); text-transform: capitalize; }}
        .filter-chips {{ display: flex; gap: 3px; flex-wrap: wrap; }}
        .filter-chip {{
            padding: 2px 6px;
            font-size: 10px;
            border: 1px solid var(--border-color);
            border-radius: 10px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            transition: all 0.15s;
        }}
        .filter-chip:hover {{ background: var(--hover-bg); }}
        .filter-chip.active {{ background: var(--accent-strong); color: white; border-color: var(--accent-strong); }}
        .filter-chip.inactive {{ opacity: 0.4; }}
        .filter-reset-btn {{
            padding: 2px 8px;
            font-size: 10px;
            border: 1px solid var(--border-color);
            border-radius: 10px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            transition: background 0.15s, opacity 0.15s;
        }}
        .filter-reset-btn:hover {{ background: var(--hover-bg); }}
        .filter-reset-btn:disabled {{
            opacity: 0.45;
            cursor: default;
        }}

        .main-container {{ display: flex; flex: 1; min-height: 0; }}
        .content-column {{
            flex: 1;
            min-height: 0;
            display: flex;
            flex-direction: column;
            position: relative;
            --umap-reserve-left: 0px;
            --umap-reserve-right: 0px;
        }}
        .grid-container {{
            flex: 1;
            min-height: 0;
            padding-top: 8px;
            padding-bottom: 8px;
            padding-left: calc(8px + var(--umap-reserve-left));
            padding-right: calc(8px + var(--umap-reserve-right));
            overflow: auto;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(clamp({min_panel_size}px, 18vw, {max_panel_size}px), 1fr));
            gap: 8px;
            align-content: start;
            transition: padding 0.2s ease;
        }}
        .grid-container.single-section-layout {{
            grid-template-columns: minmax(320px, min(100%, 900px));
            justify-content: center;
        }}
        .grid-container.single-section-layout .section-panel {{
            width: 100%;
            max-width: 900px;
            margin: 0 auto;
        }}
        .grid-container.single-visible-filter-layout {{
            grid-template-columns: minmax(320px, min(50vw, 900px));
            justify-content: center;
        }}
        .grid-container.single-visible-filter-layout .section-panel {{
            width: 100%;
            max-width: 900px;
            margin: 0 auto;
        }}
        .section-panel {{
            background: var(--panel-bg);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            overflow: hidden;
            cursor: pointer;
            transition: box-shadow 0.2s, transform 0.2s, background 0.3s, border-color 0.3s;
        }}
        .section-panel:hover {{
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
            transform: translateY(-2px);
        }}
        .section-panel.filtered-out {{ display: none; }}
        .section-header {{
            padding: 4px 8px;
            background: var(--header-bg);
            border-bottom: 1px solid var(--border-color);
            font-size: 10px;
            font-weight: 500;
            display: flex;
            justify-content: space-between;
            align-items: center;
            transition: background 0.3s, border-color 0.3s;
        }}
        .section-header-main {{
            flex: 1;
            min-width: 0;
        }}
        .section-header-actions {{
            display: flex;
            align-items: center;
            gap: 4px;
            flex-shrink: 0;
        }}
        .section-rotation-label {{
            font-size: 8px;
            color: var(--muted-color);
            min-width: 54px;
            text-align: right;
        }}
        .section-rotate-group {{
            display: inline-flex;
            align-items: center;
            gap: 2px;
        }}
        .section-rotate-btn {{
            min-width: 28px;
            padding: 2px 4px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            font-size: 9px;
            line-height: 1.2;
            transition: background 0.2s, border-color 0.2s;
        }}
        .section-rotate-btn:hover {{
            background: var(--hover-bg);
        }}
        .section-header .expand-icon {{ font-size: 10px; opacity: 0.5; }}
        .section-meta {{ font-size: 8px; color: var(--muted-color); margin-top: 1px; }}
        .section-canvas {{ display: block; width: 100%; aspect-ratio: 1; }}

        .legend-container {{
            width: 200px;
            padding: 12px;
            background: var(--panel-bg);
            border-left: 1px solid var(--border-color);
            overflow-y: auto;
            font-size: 12px;
            transition: background 0.3s, border-color 0.3s, width 0.3s, padding 0.3s;
        }}
        .legend-container.collapsed {{
            width: 0;
            padding: 0;
            overflow: hidden;
            border-left: none;
        }}
        .color-panel {{
            width: 240px;
            padding: 12px;
            background: var(--panel-bg);
            border-left: 1px solid var(--border-color);
            overflow-y: auto;
            font-size: 12px;
            transition: background 0.3s, border-color 0.3s, width 0.3s, padding 0.3s;
        }}
        .color-panel.collapsed {{
            width: 0;
            padding: 0;
            overflow: hidden;
            border-left: none;
        }}
        .color-panel-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 8px;
            padding-bottom: 6px;
            border-bottom: 1px solid var(--border-color);
        }}
        .color-panel-title {{ font-size: 13px; font-weight: 600; }}
        .color-panel-section {{
            margin-bottom: 10px;
            display: flex;
            flex-direction: column;
            gap: 6px;
        }}
        .color-panel-section label {{ font-size: 10px; color: var(--muted-color); }}
        .color-search {{
            padding: 6px 8px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--text-color);
            font-size: 12px;
        }}
        .color-tabs {{
            display: flex;
            gap: 6px;
            flex-wrap: wrap;
        }}
        .color-tab {{
            flex: 1 1 96px;
            min-width: 0;
            padding: 4px 6px;
            font-size: 10px;
            line-height: 1.2;
            white-space: normal;
            overflow-wrap: anywhere;
            text-align: center;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            transition: background 0.2s, border-color 0.2s, color 0.2s;
        }}
        .color-tab.active {{
            background: var(--accent-strong);
            color: white;
            border-color: var(--accent-strong);
        }}
        .color-tab-content {{
            display: none;
            flex-direction: column;
            gap: 8px;
        }}
        .color-tab-content.active {{
            display: flex;
        }}
        .scale-controls {{
            display: flex;
            align-items: center;
            gap: 6px;
            flex-wrap: wrap;
        }}
        .scale-controls input[type="number"] {{
            width: 80px;
        }}
        .scale-sep {{
            font-size: 11px;
            color: var(--muted-color);
        }}
        .scale-hint {{
            font-size: 10px;
            color: var(--muted-color);
        }}
        .info-content {{
            display: grid;
            gap: 10px;
        }}
        .info-block {{
            padding: 10px;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            background: var(--panel-bg);
        }}
        .info-title {{
            font-size: 11px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.06em;
            color: var(--muted-color);
            margin-bottom: 6px;
        }}
        .info-text {{
            font-size: 12px;
            line-height: 1.5;
        }}
        .info-list {{
            display: grid;
            gap: 4px;
            font-size: 12px;
        }}
        .info-link {{
            color: inherit;
            text-decoration: none;
            border-bottom: 1px dotted var(--border-color);
        }}
        .info-link:hover {{
            color: var(--accent-strong);
            border-bottom-color: var(--accent-strong);
        }}
        .info-shortcuts-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 11px;
        }}
        .info-shortcuts-table th,
        .info-shortcuts-table td {{
            padding: 4px 0;
            border-top: 1px solid var(--border-color);
            text-align: left;
            vertical-align: top;
        }}
        .info-shortcuts-table tr:first-child th,
        .info-shortcuts-table tr:first-child td {{
            border-top: none;
            padding-top: 0;
        }}
        .info-shortcuts-key {{
            width: 76px;
            white-space: nowrap;
            color: var(--muted-color);
            font-variant-numeric: tabular-nums;
        }}
        .info-shortcuts-table kbd {{
            display: inline-block;
            padding: 1px 5px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--text-color);
            font: inherit;
        }}
        .info-popover {{
            position: absolute;
            top: calc(100% + 8px);
            left: 0;
            z-index: 10;
            width: 280px;
            max-width: calc(100vw - 32px);
            padding: 10px;
            border: 1px solid var(--border-color);
            border-radius: 10px;
            background: var(--panel-bg);
            box-shadow: 0 10px 30px rgba(0, 0, 0, 0.12);
            display: none;
        }}
        .info-popover.active {{
            display: block;
        }}
        .color-list {{
            display: flex;
            flex-direction: column;
            gap: 4px;
            max-height: 220px;
            overflow-y: auto;
            padding-right: 2px;
        }}
        .color-item {{
            padding: 5px 8px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            font-size: 11px;
            transition: background 0.2s, border-color 0.2s;
        }}
        .color-item:hover {{ background: var(--hover-bg); }}
        .color-item.active {{
            background: var(--accent-strong);
            color: white;
            border-color: var(--accent-strong);
        }}
        .color-aggregation {{
            display: flex;
            flex-direction: column;
            gap: 8px;
            min-width: 0;
            font-size: 11px;
        }}
        .color-aggregation.collapsed {{
            display: none;
        }}
        .marker-genes {{
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 8px;
            background: var(--input-bg);
            max-height: 420px;
            overflow: auto;
            display: flex;
            flex-direction: column;
            gap: 8px;
            resize: vertical;
        }}
        .marker-search {{
            padding: 6px 8px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--text-color);
            font-size: 11px;
        }}
        .marker-group-title {{ font-size: 11px; font-weight: 600; }}
        .marker-genes .gene-token-grid {{ gap: 4px; }}
        .marker-genes .gene-token-btn {{
            gap: 4px;
            padding: 3px 6px;
            font-size: 10px;
        }}
        .marker-genes .gene-token-meta {{ font-size: 9px; }}
        .marker-genes-list {{
            font-size: 10px;
            color: var(--muted-color);
            line-height: 1.4;
            word-break: break-word;
        }}
        .agg-group {{
            padding: 6px;
            border-radius: 6px;
            background: rgba(135, 0, 82, 0.06);
            border: 1px solid var(--border-color);
            min-width: 0;
            overflow: hidden;
        }}
        .agg-group-title {{ font-weight: 600; margin-bottom: 4px; }}
        .agg-group-meta {{ font-size: 10px; color: var(--muted-color); margin-bottom: 4px; }}
        .agg-row {{
            display: flex;
            align-items: center;
            gap: 6px;
            margin: 2px 0;
        }}
        .agg-dot {{
            width: 8px;
            height: 8px;
            border-radius: 50%;
            flex-shrink: 0;
        }}
        .agg-label {{ flex: 1; }}
        .agg-value {{ font-variant-numeric: tabular-nums; }}
        .comparison-stack {{
            display: flex;
            flex-direction: column;
            gap: 6px;
        }}
        .comparison-card {{
            padding: 6px;
            border-radius: 6px;
            background: var(--input-bg);
            border: 1px solid var(--border-color);
            min-width: 0;
        }}
        .comparison-card-title {{
            display: flex;
            align-items: flex-start;
            gap: 6px;
            font-weight: 600;
            margin-bottom: 4px;
            min-width: 0;
        }}
        .gene-search-btn {{
            flex: 0 0 auto;
            min-width: 0;
            padding: 2px 7px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--muted-color);
            cursor: pointer;
            font-size: 9px;
            font-weight: 600;
            line-height: 1.2;
            white-space: nowrap;
        }}
        .gene-search-btn:hover {{
            background: var(--hover-bg);
            color: var(--text-color);
        }}
        .comparison-card-title span:last-child {{
            min-width: 0;
            overflow-wrap: anywhere;
            word-break: break-word;
        }}
        .comparison-metric-grid {{
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 6px 8px;
        }}
        .comparison-metric {{
            min-width: 0;
        }}
        .comparison-metric-label {{
            display: block;
            font-size: 9px;
            color: var(--muted-color);
            margin-bottom: 2px;
        }}
        .comparison-metric-value {{
            display: block;
            font-variant-numeric: tabular-nums;
            overflow-wrap: anywhere;
            word-break: break-word;
        }}
        .trend-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 10px;
        }}
        .trend-table th, .trend-table td {{
            text-align: left;
            padding: 4px 6px;
            border-bottom: 1px solid var(--border-color);
            font-variant-numeric: tabular-nums;
            vertical-align: top;
        }}
        .trend-table th {{ color: var(--muted-color); font-weight: 600; }}
        .interaction-target {{
            display: flex;
            align-items: flex-start;
            gap: 6px;
            min-width: 0;
        }}
        .interaction-target span:last-child {{
            min-width: 0;
            overflow-wrap: anywhere;
            word-break: break-word;
        }}
        .interaction-markers {{
            font-size: 10px;
            color: var(--muted-color);
            line-height: 1.35;
        }}

        /* Dotplot */
        .dotplot-controls {{
            display: flex;
            flex-direction: column;
            gap: 8px;
        }}
        .dotplot-grid {{
            border: 1px solid var(--border-color);
            border-radius: 8px;
            background: var(--input-bg);
            overflow: auto;
            max-height: 420px;
        }}
        .dotplot-row {{
            display: grid;
            grid-template-columns: 180px repeat(var(--dotplot-cols, 1), 24px);
            align-items: center;
            gap: 6px;
            padding: 6px 8px;
            border-bottom: 1px solid var(--border-color);
        }}
        .dotplot-row:last-child {{ border-bottom: none; }}
        .dotplot-row.dotplot-header {{
            position: sticky;
            top: 0;
            background: var(--panel-bg);
            z-index: 1;
        }}
        .dotplot-label {{
            font-size: 10px;
            color: var(--text-color);
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }}
        .dotplot-gene {{
            font-size: 9px;
            color: var(--muted-color);
            text-align: center;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }}
        .dotplot-dot {{
            width: 20px;
            height: 20px;
            display: flex;
            align-items: center;
            justify-content: center;
        }}
        .cluster-de-controls {{
            display: flex;
            flex-direction: column;
            gap: 8px;
        }}
        .cluster-de-select-row {{
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 8px;
        }}
        .cluster-de-result {{
            border: 1px solid var(--border-color);
            border-radius: 8px;
            background: var(--input-bg);
            overflow: auto;
            max-height: 420px;
        }}
        .cluster-de-table {{
            width: 100%;
            border-collapse: collapse;
            table-layout: fixed;
            font-size: 10px;
        }}
        .cluster-de-table th, .cluster-de-table td {{
            text-align: left;
            padding: 6px 8px;
            border-bottom: 1px solid var(--border-color);
            font-variant-numeric: tabular-nums;
            vertical-align: top;
            overflow-wrap: anywhere;
            word-break: break-word;
        }}
        .cluster-de-table th {{
            position: sticky;
            top: 0;
            z-index: 1;
            background: var(--panel-bg);
            color: var(--muted-color);
            font-weight: 600;
        }}
        .cluster-de-table tr:last-child td {{ border-bottom: none; }}
        .cluster-de-gene-btn {{
            padding: 0;
            border: none;
            background: none;
            color: var(--accent-strong);
            cursor: pointer;
            font: inherit;
            text-align: left;
        }}
        .cluster-de-gene-btn:hover {{ text-decoration: underline; }}
        .cluster-de-meta {{
            padding: 8px;
            font-size: 10px;
            color: var(--muted-color);
            border-bottom: 1px solid var(--border-color);
        }}
        .legend-title {{
            font-size: 13px;
            font-weight: 600;
            margin-bottom: 8px;
            padding-bottom: 6px;
            border-bottom: 1px solid var(--border-color);
        }}
        .legend-actions {{ display: flex; gap: 6px; margin-bottom: 8px; }}
        .legend-btn {{
            flex: 1;
            padding: 3px 6px;
            font-size: 10px;
            border: 1px solid var(--border-color);
            border-radius: 3px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            transition: background 0.3s, border-color 0.3s, color 0.3s;
        }}
        .legend-btn:hover {{ background: var(--hover-bg); }}
        .legend-btn.active {{
            background: var(--accent-strong);
            color: #fff;
            border-color: var(--accent-strong);
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 6px;
            padding: 3px 6px;
            margin: 1px 0;
            font-size: 11px;
            cursor: pointer;
            border-radius: 3px;
            transition: background 0.15s;
        }}
        .legend-item:hover {{ background: var(--hover-bg); }}
        .legend-item.hidden {{ opacity: 0.3; text-decoration: line-through; }}
        .legend-item.spotlight {{
            background: color-mix(in srgb, var(--accent-strong) 10%, transparent);
            box-shadow: inset 0 0 0 1px color-mix(in srgb, var(--accent-strong) 35%, transparent);
        }}
        .legend-item.dimmed {{ opacity: 0.45; }}
        .legend-color {{
            width: 12px;
            height: 12px;
            border-radius: 50%;
            flex-shrink: 0;
            border: 2px solid transparent;
        }}
        .legend-item.hidden .legend-color {{ border-color: var(--muted-color); background: transparent !important; }}
        .legend-item.selected {{
            border-color: var(--accent-strong);
            box-shadow: 0 0 0 2px rgba(135, 0, 82, 0.2);
        }}
        .legend-item.selected .legend-color {{
            border-color: var(--accent-strong);
        }}
        .split-legend-list {{
            display: flex;
            flex-direction: column;
            gap: 8px;
        }}
        .split-legend-item {{
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 7px 8px;
            background: color-mix(in srgb, var(--input-bg) 88%, transparent);
        }}
        .split-legend-row {{
            display: flex;
            flex-direction: column;
            gap: 4px;
        }}
        .split-legend-side {{
            display: inline-flex;
            align-items: center;
            width: fit-content;
            font-size: 10px;
            font-weight: 700;
            letter-spacing: 0.04em;
            text-transform: uppercase;
            border-radius: 999px;
            padding: 1px 6px;
            border: 1px solid transparent;
        }}
        .split-legend-side-a {{
            color: #1d4ed8;
            background: color-mix(in srgb, #1d4ed8 12%, transparent);
            border-color: color-mix(in srgb, #1d4ed8 28%, transparent);
        }}
        .split-legend-side-b {{
            color: #be185d;
            background: color-mix(in srgb, #be185d 12%, transparent);
            border-color: color-mix(in srgb, #be185d 28%, transparent);
        }}
        .split-legend-label {{
            font-size: 11px;
            color: var(--text-color);
            line-height: 1.35;
            word-break: break-word;
        }}
        .split-legend-key {{
            margin-top: 6px;
            display: flex;
            align-items: center;
            gap: 6px;
            font-size: 10px;
            color: var(--muted-color);
            min-height: 12px;
        }}
        .split-legend-swatch {{
            display: inline-block;
            flex-shrink: 0;
            border: 1px solid color-mix(in srgb, var(--border-color) 75%, transparent);
            background: #9ca3af;
        }}
        .split-legend-swatch-dot {{
            width: 12px;
            height: 12px;
            border-radius: 50%;
        }}
        .split-legend-swatch-bar {{
            width: 40px;
            height: 10px;
            border-radius: 999px;
        }}
        .split-legend-categories {{
            margin-top: 6px;
            padding-top: 6px;
            border-top: 1px dashed color-mix(in srgb, var(--border-color) 72%, transparent);
            display: flex;
            flex-direction: column;
            gap: 4px;
            max-height: 190px;
            overflow-y: auto;
        }}
        .split-legend-cat {{
            display: flex;
            align-items: center;
            gap: 6px;
            font-size: 10px;
            min-width: 0;
        }}
        .split-legend-cat-label {{
            min-width: 0;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            color: var(--text-color);
        }}
        .colorbar {{ width: 16px; height: 150px; margin: 8px auto; border-radius: 2px; }}
        .colorbar-labels {{
            display: flex;
            flex-direction: column;
            justify-content: space-between;
            height: 150px;
            font-size: 10px;
            color: var(--muted-color);
            margin-left: 6px;
        }}
        .colorbar-container {{ display: flex; align-items: stretch; justify-content: center; }}

        /* Modal styles */
        .modal-overlay {{
            display: none;
            position: fixed;
            top: 0; left: 0; right: 0; bottom: 0;
            background: rgba(0,0,0,0.7);
            z-index: 1000;
            align-items: center;
            justify-content: center;
        }}
        .modal-overlay.active {{ display: flex; }}
        .modal-content {{
            background: var(--panel-bg);
            border-radius: 12px;
            width: 90vw;
            height: 90vh;
            max-width: 1400px;
            display: flex;
            flex-direction: column;
            overflow: hidden;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            transition: background 0.3s;
        }}
        .modal-header {{
            padding: 12px 16px;
            background: var(--header-bg);
            border-bottom: 1px solid var(--border-color);
            display: flex;
            justify-content: space-between;
            align-items: center;
            transition: background 0.3s, border-color 0.3s;
        }}
        .modal-header-actions {{
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        .modal-header h2 {{ font-size: 15px; font-weight: 600; }}
        .modal-header .modal-meta {{ font-size: 11px; color: var(--muted-color); margin-left: 10px; }}
        .modal-controls-toggle {{
            padding: 4px 8px;
            border: 1px solid var(--border-color);
            border-radius: 999px;
            background: var(--panel-bg);
            color: var(--text-color);
            font-size: 10px;
            letter-spacing: 0.03em;
            text-transform: uppercase;
            cursor: pointer;
            transition: background 0.2s, border-color 0.2s, color 0.2s;
        }}
        .modal-controls-toggle:hover {{
            background: var(--hover-bg);
            border-color: var(--accent-strong);
        }}
        .modal-close {{
            background: none;
            border: none;
            font-size: 22px;
            cursor: pointer;
            color: var(--text-color);
            padding: 2px 6px;
            border-radius: 4px;
        }}
        .modal-close:hover {{ background: var(--hover-bg); }}
        .modal-body {{ flex: 1; display: flex; overflow: hidden; }}
        .modal-canvas-container {{ flex: 1; position: relative; overflow: hidden; background: var(--panel-bg); }}
        .modal-canvas {{ position: absolute; top: 0; left: 0; width: 100%; height: 100%; transform-origin: 0 0; will-change: transform; }}
        .modal-controls {{
            position: absolute;
            bottom: 12px;
            left: 50%;
            transform: translateX(-50%);
            display: flex;
            flex-wrap: wrap;
            justify-content: flex-start;
            align-items: stretch;
            gap: 3px;
            width: max-content;
            max-width: min(960px, calc(100% - 72px));
            padding: 3px 4px;
            border: 1px solid color-mix(in srgb, var(--border-color) 82%, transparent);
            border-radius: 12px;
            background: color-mix(in srgb, var(--header-bg) 90%, transparent);
            backdrop-filter: blur(12px);
            box-shadow: 0 10px 28px rgba(0,0,0,0.2);
            transition: background 0.3s, border-color 0.3s, box-shadow 0.3s;
        }}
        .modal-controls.dragging {{
            box-shadow: 0 14px 34px rgba(0,0,0,0.28);
        }}
        .modal-control-group {{
            display: flex;
            align-items: center;
            gap: 3px;
            padding: 2px 4px;
            border: 1px solid color-mix(in srgb, var(--border-color) 76%, transparent);
            border-radius: 8px;
            background: color-mix(in srgb, var(--panel-bg) 92%, transparent);
            min-width: 0;
            cursor: grab;
        }}
        .modal-control-group-title {{
            flex-shrink: 0;
            font-size: 9px;
            font-weight: 700;
            letter-spacing: 0.08em;
            text-transform: uppercase;
            color: var(--muted-color);
            cursor: inherit;
            user-select: none;
        }}
        .modal-controls.dragging .modal-control-group,
        .modal-controls.dragging .modal-control-group-title {{ cursor: grabbing; }}
        .modal-control-group-body {{
            display: flex;
            flex-wrap: wrap;
            align-items: center;
            gap: 3px;
            min-width: 0;
        }}
        .modal-controls button {{
            min-height: 22px;
            padding: 1px 7px;
            border: 1px solid var(--border-color);
            border-radius: 999px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            font-size: 11px;
            font-weight: 600;
            line-height: 1.4;
            transition: background 0.3s, border-color 0.3s, color 0.3s, box-shadow 0.3s, opacity 0.3s;
        }}
        .modal-controls button:hover {{ background: var(--hover-bg); }}
        .modal-controls button.modal-btn-primary {{
            background: color-mix(in srgb, var(--accent-strong) 92%, var(--panel-bg));
            border-color: color-mix(in srgb, var(--accent-strong) 94%, white 6%);
            color: white;
            box-shadow: 0 0 0 1px color-mix(in srgb, var(--accent-strong) 24%, transparent);
        }}
        .modal-controls button.modal-btn-primary:hover {{
            background: color-mix(in srgb, var(--accent-strong) 100%, black 0%);
        }}
        .modal-controls button.modal-btn-muted {{
            color: var(--muted-color);
            border-color: color-mix(in srgb, var(--border-color) 62%, transparent);
            background: color-mix(in srgb, var(--input-bg) 72%, transparent);
            opacity: 0.86;
        }}
        .modal-controls button.modal-btn-muted:hover {{
            color: var(--text-color);
            opacity: 1;
        }}
        .modal-controls button.modal-btn-danger {{
            color: color-mix(in srgb, #ff7b7b 72%, var(--text-color));
            border-color: color-mix(in srgb, #ff7b7b 35%, var(--border-color));
            background: color-mix(in srgb, #ff7b7b 10%, var(--input-bg));
        }}
        .modal-controls button.modal-btn-danger:hover {{
            background: color-mix(in srgb, #ff7b7b 16%, var(--hover-bg));
        }}
        .modal-controls select {{
            min-height: 22px;
            padding: 1px 7px;
            border: 1px solid var(--border-color);
            border-radius: 999px;
            background: var(--input-bg);
            color: var(--text-color);
            font-size: 11px;
            line-height: 1.4;
        }}
        .modal-size-block {{
            display: flex;
            align-items: center;
            gap: 4px;
        }}
        .modal-size-value {{
            min-width: 26px;
            font-size: 10px;
            font-variant-numeric: tabular-nums;
            color: var(--text-color);
        }}
        .modal-controls.hidden {{ display: none; }}
        @media (max-width: 960px) {{
            .modal-controls {{
                left: 12px;
                right: 12px;
                width: auto;
                transform: none;
            }}
            .modal-control-group {{
                width: 100%;
                justify-content: space-between;
            }}
            .modal-control-group-body {{
                justify-content: flex-end;
                flex: 1;
            }}
        }}
        .modal-legend {{
            width: 180px;
            padding: 12px;
            border-left: 1px solid var(--border-color);
            overflow-y: auto;
            transition: border-color 0.3s, width 0.2s ease;
        }}
        .modal-legend.split-expanded {{ width: 280px; }}
        .zoom-info {{ font-size: 10px; color: var(--muted-color); margin-left: 6px; }}
        .modal-blend-panel {{
            position: absolute;
            top: 8px;
            right: 8px;
            width: min(420px, calc(100% - 16px));
            padding: 8px;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            background: color-mix(in srgb, var(--header-bg) 92%, transparent);
            backdrop-filter: blur(6px);
            z-index: 3;
            display: none;
            pointer-events: auto;
        }}
        .modal-blend-panel.visible {{ display: block; }}
        .modal-blend-title {{
            font-size: 10px;
            text-transform: uppercase;
            letter-spacing: 0.04em;
            color: var(--muted-color);
            margin-bottom: 6px;
        }}
        .modal-blend-row {{
            display: grid;
            grid-template-columns: 18px 112px 1fr;
            gap: 6px;
            align-items: center;
            margin-bottom: 6px;
        }}
        .modal-blend-row:last-child {{ margin-bottom: 0; }}
        .modal-blend-row select,
        .modal-blend-row input[type="text"],
        .modal-blend-row input[type="range"] {{
            min-width: 0;
            width: 100%;
            font-size: 11px;
        }}
        .modal-blend-side {{
            font-size: 11px;
            font-weight: 600;
            color: var(--text-color);
        }}
        .modal-blend-fields {{
            display: flex;
            gap: 4px;
            min-width: 0;
        }}
        .modal-blend-fields > * {{
            flex: 1 1 0;
            min-width: 0;
        }}
        .modal-blend-scale-row {{
            margin-top: -2px;
        }}
        .modal-blend-scale-label {{
            font-size: 10px;
            color: var(--muted-color);
        }}
        .modal-blend-scale-controls {{
            display: flex;
            align-items: center;
            gap: 4px;
            min-width: 0;
        }}
        .modal-blend-scale-controls input[type="number"] {{
            width: 64px;
            min-width: 0;
            padding: 4px 6px;
            font-size: 10px;
        }}
        .modal-blend-scale-controls .legend-btn {{
            flex: 0 0 auto;
            font-size: 10px;
            padding: 3px 6px;
            min-width: 84px;
        }}
        .modal-blend-meta {{
            margin-top: 4px;
            font-size: 10px;
            color: var(--muted-color);
            line-height: 1.35;
            white-space: pre-line;
            max-height: 180px;
            overflow: auto;
        }}
        .modal-blend-loading {{
            display: none;
            align-items: center;
            gap: 7px;
            margin: 6px 0 4px;
            padding: 6px 8px;
            border-radius: 8px;
            border: 1px solid var(--border-color);
            background: rgba(255, 255, 255, 0.05);
            color: var(--text-color);
            font-size: 10px;
            line-height: 1.3;
        }}
        .modal-blend-loading.visible {{
            display: flex;
        }}
        .modal-blend-loading::before {{
            content: '';
            width: 7px;
            height: 7px;
            border-radius: 999px;
            background: var(--accent-strong);
            animation: modalBlendPulse 1s ease-in-out infinite;
            flex: 0 0 auto;
        }}
        @keyframes modalBlendPulse {{
            0%, 100% {{ transform: scale(0.85); opacity: 0.45; }}
            50% {{ transform: scale(1.15); opacity: 1; }}
        }}
        .modal-blend-summary {{
            color: var(--text-color);
            font-weight: 500;
        }}
        .modal-blend-markers-toggle {{
            margin-top: 4px;
            padding: 0;
            border: none;
            background: none;
            color: var(--accent-strong);
            font-size: 10px;
            font-weight: 600;
            cursor: pointer;
        }}
        .modal-blend-markers-toggle:hover {{
            text-decoration: underline;
        }}
        .modal-blend-markers {{
            margin-top: 4px;
        }}
        .modal-blend-markers.collapsed {{
            display: none;
        }}
        .modal-blend-marker-header {{
            margin-top: 4px;
            font-weight: 600;
            color: var(--text-color);
        }}
        .modal-blend-marker-row {{
            margin-top: 2px;
        }}
        .modal-blend-marker-name {{
            color: var(--text-color);
            font-weight: 600;
        }}
        .modal-selection-info {{
            position: absolute;
            top: 8px;
            left: 8px;
            font-size: 11px;
            color: var(--muted-color);
            padding: 4px 8px;
            border: 1px solid var(--border-color);
            border-radius: 999px;
            background: color-mix(in srgb, var(--header-bg) 88%, transparent);
            backdrop-filter: blur(6px);
            pointer-events: none;
            z-index: 2;
        }}
        .selection-summary {{
            position: absolute;
            left: 8px;
            max-width: min(360px, calc(100% - 16px));
            max-height: min(70%, 420px);
            overflow-y: auto;
            padding: 8px 10px;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            background: color-mix(in srgb, var(--header-bg) 88%, transparent);
            backdrop-filter: blur(6px);
            font-size: 11px;
            color: var(--text-color);
            pointer-events: auto;
            z-index: 2;
        }}
        .selection-summary.expanded {{
            max-height: min(82%, 560px);
        }}
        .selection-summary-title {{
            margin: 2px 0 6px;
            font-size: 10px;
            text-transform: uppercase;
            letter-spacing: 0.04em;
            color: var(--muted-color);
        }}
        .selection-summary-row {{
            display: flex;
            justify-content: space-between;
            gap: 10px;
            line-height: 1.3;
        }}
        .selection-summary-row + .selection-summary-row {{ margin-top: 2px; }}
        .selection-summary-label {{
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            max-width: 210px;
        }}
        .selection-summary-count {{
            color: var(--muted-color);
            font-variant-numeric: tabular-nums;
        }}
        .selection-summary-meta {{
            margin-top: 6px;
            font-size: 10px;
            color: var(--muted-color);
            line-height: 1.35;
        }}
        .selection-summary-toggle {{
            margin-top: 6px;
            padding: 0;
            border: none;
            background: none;
            color: var(--accent-strong);
            font-size: 10px;
            text-decoration: underline;
            cursor: pointer;
            pointer-events: auto;
        }}
        .selection-summary-toggle:hover {{
            color: var(--accent);
        }}
        .selection-summary-expr {{
            margin-top: 8px;
            border-top: 1px solid var(--border-color);
            padding-top: 6px;
        }}
        .selection-summary-expr-row {{
            display: flex;
            align-items: center;
            gap: 5px;
            margin-bottom: 3px;
        }}
        .selection-summary-expr-gene {{
            width: 80px;
            flex-shrink: 0;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            font-size: 10px;
        }}
        .selection-summary-expr-bars {{
            flex: 1;
            display: flex;
            flex-direction: column;
            gap: 1px;
        }}
        .selection-summary-expr-bar {{
            height: 4px;
            border-radius: 2px;
            min-width: 1px;
        }}
        .selection-summary-expr-bar.sel {{ background: var(--accent-strong); }}
        .selection-summary-expr-bar.rest {{ background: var(--muted-color); opacity: 0.5; }}
        .selection-summary-expr-bar.region-b {{ background: #4cc9f0; opacity: 0.85; }}
        .selection-summary-expr-fc {{
            width: 34px;
            flex-shrink: 0;
            text-align: right;
            font-size: 10px;
            font-variant-numeric: tabular-nums;
            color: var(--muted-color);
        }}
        .selection-summary-compare-btn {{
            margin-top: 6px;
            width: 100%;
            padding: 4px 8px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--panel-bg);
            color: var(--text-color);
            cursor: pointer;
            font-size: 11px;
        }}
        .selection-summary-compare-btn:hover {{ background: var(--hover-bg); }}
        .selection-summary-actions {{
            display: flex;
            flex-direction: column;
            gap: 6px;
            margin-top: 6px;
        }}
        .selection-summary-actions .selection-summary-compare-btn {{
            margin-top: 0;
        }}
        .selection-summary-compare-header {{
            display: flex;
            gap: 4px;
            margin-bottom: 4px;
            font-size: 11px;
            font-weight: 600;
        }}
        .selection-summary-compare-label {{
            flex: 1;
            text-align: center;
            padding: 2px 4px;
            border-radius: 3px;
        }}
        .selection-summary-compare-label.region-a {{ background: var(--accent-strong); color: #fff; }}
        .selection-summary-compare-label.region-b {{ background: #4cc9f0; color: #000; }}
        .selection-summary-compare-row {{
            display: grid;
            grid-template-columns: 1fr auto auto auto;
            gap: 4px;
            align-items: center;
            font-size: 11px;
            padding: 1px 0;
        }}
        .selection-summary-compare-type {{ color: var(--text-color); font-weight: 500; }}
        .selection-summary-compare-a {{ color: var(--accent-strong); font-variant-numeric: tabular-nums; }}
        .selection-summary-compare-b {{ color: #4cc9f0; font-variant-numeric: tabular-nums; }}
        .selection-summary-compare-sep {{ color: var(--muted-color); }}
        .annot-comp-row {{ margin-bottom: 8px; }}
        .annot-comp-header {{ display: flex; align-items: center; gap: 6px; font-size: 11px; font-weight: 600; margin-bottom: 3px; }}
        .annot-comp-dot {{ width: 10px; height: 10px; border-radius: 50%; display: inline-block; flex-shrink: 0; }}
        .annot-comp-count {{ color: var(--muted-color); font-weight: normal; margin-left: auto; }}
        .annot-comp-bar-track {{ display: flex; height: 8px; border-radius: 4px; overflow: hidden; background: var(--border-color); }}
        .annot-comp-segment {{ height: 100%; }}
        .annot-comp-legend {{ display: flex; flex-wrap: wrap; gap: 4px; margin-top: 4px; font-size: 10px; }}
        .annot-comp-legend-item {{ display: flex; align-items: center; gap: 3px; color: var(--text-color); }}
        .annot-expr-row {{ display: flex; align-items: center; gap: 6px; margin: 2px 0; }}
        .annot-expr-bars {{ flex: 1; display: flex; flex-direction: column; gap: 1px; }}
        .annot-expr-bar {{ height: 4px; border-radius: 2px; min-width: 1px; }}
        .annot-de-controls {{
            display: flex;
            flex-direction: column;
            gap: 8px;
            margin-top: 10px;
        }}
        .modal-annotation-panel {{
            position: absolute;
            right: 8px;
            bottom: 82px;
            width: min(360px, calc(100% - 16px));
            max-height: min(56%, 380px);
            overflow-y: auto;
            padding: 8px;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            background: color-mix(in srgb, var(--header-bg) 88%, transparent);
            backdrop-filter: blur(6px);
            z-index: 2;
            pointer-events: auto;
        }}
        .modal-annotation-panel.dragging {{
            user-select: none;
        }}
        .modal-annotation-title {{
            font-size: 10px;
            text-transform: uppercase;
            letter-spacing: 0.04em;
            color: var(--muted-color);
            margin: 0 0 6px;
            cursor: grab;
            user-select: none;
            touch-action: none;
        }}
        .modal-annotation-panel.dragging .modal-annotation-title {{
            cursor: grabbing;
        }}
        .modal-annotation-actions {{
            display: flex;
            flex-wrap: wrap;
            gap: 6px;
            margin-bottom: 6px;
        }}
        .modal-annotation-actions button {{
            padding: 4px 8px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            font-size: 11px;
            transition: background 0.3s, border-color 0.3s, color 0.3s;
        }}
        .modal-annotation-actions button:hover {{ background: var(--hover-bg); }}
        .modal-annotation-list {{
            display: flex;
            flex-direction: column;
            gap: 6px;
        }}
        .modal-annotation-empty {{
            font-size: 10px;
            color: var(--muted-color);
            line-height: 1.35;
        }}
        .modal-annotation-row {{
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 6px;
            background: color-mix(in srgb, var(--panel-bg) 85%, transparent);
        }}
        .modal-annotation-row-main {{
            display: grid;
            grid-template-columns: 12px 1fr auto;
            gap: 6px;
            align-items: center;
        }}
        .modal-annotation-color {{
            width: 10px;
            height: 10px;
            border-radius: 50%;
            border: 1px solid rgba(255, 255, 255, 0.65);
        }}
        .modal-annotation-label {{
            min-width: 0;
            width: 100%;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--text-color);
            font-size: 11px;
            padding: 3px 6px;
        }}
        .modal-annotation-count {{
            font-size: 10px;
            color: var(--muted-color);
            white-space: nowrap;
        }}
        .modal-annotation-row-actions {{
            margin-top: 6px;
            display: flex;
            gap: 6px;
        }}
        .modal-annotation-row-actions button {{
            flex: 1;
            padding: 4px 6px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--text-color);
            font-size: 10px;
            cursor: pointer;
            transition: background 0.3s, border-color 0.3s, color 0.3s;
        }}
        .modal-annotation-row-actions button:hover {{ background: var(--hover-bg); }}

        .no-results {{
            grid-column: 1 / -1;
            text-align: center;
            padding: 40px;
            color: var(--muted-color);
            font-size: 14px;
        }}

        /* Tooltip styles */
        .cell-tooltip {{
            position: fixed;
            background: var(--header-bg);
            border: 1px solid var(--border-color);
            border-radius: 4px;
            padding: 6px 10px;
            font-size: 11px;
            pointer-events: none;
            z-index: 2000;
            display: none;
            box-shadow: 0 2px 8px rgba(0,0,0,0.15);
            max-width: 250px;
            transition: background 0.3s, border-color 0.3s;
        }}
        .cell-tooltip.visible {{ display: block; }}
        .cell-tooltip-color {{
            display: inline-block;
            width: 10px;
            height: 10px;
            border-radius: 50%;
            margin-right: 6px;
            vertical-align: middle;
        }}
        .cell-tooltip-label {{ font-weight: 500; }}
        .cell-tooltip-value {{ color: var(--muted-color); margin-left: 4px; }}
        .help-tooltip {{
            position: fixed;
            left: 0;
            top: 0;
            max-width: 280px;
            padding: 8px 10px;
            border: 1px solid var(--border-color);
            border-radius: 6px;
            background: color-mix(in srgb, var(--header-bg) 92%, transparent);
            color: var(--text-color);
            font-size: 11px;
            line-height: 1.35;
            box-shadow: 0 6px 20px rgba(0, 0, 0, 0.2);
            pointer-events: none;
            z-index: 2100;
            opacity: 0;
            transform: translateY(2px);
            transition: opacity 0.12s ease, transform 0.12s ease, background 0.3s, border-color 0.3s, color 0.3s;
        }}
        .help-tooltip.visible {{
            opacity: 1;
            transform: translateY(0);
        }}
        :root.dark button {{
            color: #ffffff;
        }}

        /* UMAP panel styles */
        .umap-panel {{
            width: 320px;
            max-width: min(560px, 90vw);
            aspect-ratio: 1 / 1;
            height: auto;
            margin: 0;
            position: absolute;
            z-index: 20;
            background: var(--panel-bg);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            box-shadow: 0 8px 20px rgba(0, 0, 0, 0.12);
            overflow: hidden;
            display: none;
            flex-direction: column;
            transition: background 0.3s, border-color 0.3s;
        }}
        .umap-panel.dock-top-right {{ top: 8px; right: 8px; left: auto; bottom: auto; }}
        .umap-panel.dock-top-left {{ top: 8px; left: 8px; right: auto; bottom: auto; }}
        .umap-panel.dock-bottom-right {{ bottom: 8px; right: 8px; left: auto; top: auto; }}
        .umap-panel.dock-bottom-left {{ bottom: 8px; left: 8px; right: auto; top: auto; }}
        .umap-panel.visible {{ display: flex; }}
        .umap-header {{
            padding: 8px 12px;
            background: var(--header-bg);
            border-bottom: 1px solid var(--border-color);
            display: flex;
            justify-content: space-between;
            align-items: center;
            transition: background 0.3s, border-color 0.3s;
        }}
        .umap-header h3 {{ font-size: 13px; font-weight: 600; }}
        .umap-header-actions {{
            display: flex;
            align-items: center;
            gap: 4px;
        }}
        .umap-canvas-container {{
            flex: 1;
            position: relative;
            min-height: 0;
        }}
        .umap-canvas {{
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
        }}
        .umap-controls {{
            position: absolute;
            left: 8px;
            right: 8px;
            bottom: 8px;
            padding: 6px;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            background: color-mix(in srgb, var(--header-bg) 88%, transparent);
            backdrop-filter: blur(6px);
            display: flex;
            flex-wrap: wrap;
            gap: 6px;
            align-items: center;
            transition: border-color 0.3s, background 0.3s;
        }}
        .umap-btn {{
            padding: 5px 10px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            font-size: 11px;
            transition: background 0.3s, border-color 0.3s, color 0.3s;
        }}
        .umap-btn:hover {{ background: var(--hover-bg); }}
        .umap-btn.active {{
            background: var(--accent-strong);
            color: white;
            border-color: var(--accent-strong);
        }}
        .umap-selection-info {{
            position: absolute;
            top: 8px;
            left: 8px;
            font-size: 11px;
            color: var(--muted-color);
            padding: 4px 8px;
            border: 1px solid var(--border-color);
            border-radius: 999px;
            background: color-mix(in srgb, var(--header-bg) 88%, transparent);
            backdrop-filter: blur(6px);
        }}
        .umap-selection-summary {{
            top: 38px;
        }}
        .modal-selection-summary {{
            top: 38px;
        }}
        .umap-toggle, .legend-toggle, .graph-toggle {{
            background: var(--input-bg);
            border: 1px solid var(--border-color);
            border-radius: 4px;
            padding: 5px 10px;
            cursor: pointer;
            font-size: 12px;
            transition: background 0.3s, border-color 0.3s;
        }}
        .umap-toggle:hover, .legend-toggle:hover, .graph-toggle:hover {{ background: var(--hover-bg); }}
        .umap-toggle.active, .legend-toggle.active, .graph-toggle.active {{
            background: var(--accent-strong);
            color: white;
            border-color: var(--accent-strong);
        }}
        .graph-toggle:disabled {{
            opacity: 0.5;
            cursor: not-allowed;
        }}
        .color-toggle {{
            background: var(--input-bg);
            border: 1px solid var(--border-color);
            border-radius: 4px;
            padding: 5px 10px;
            cursor: pointer;
            font-size: 12px;
            transition: background 0.3s, border-color 0.3s;
        }}
        .color-toggle:hover {{ background: var(--hover-bg); }}
        .color-toggle.active {{
            background: var(--accent-strong);
            color: white;
            border-color: var(--accent-strong);
        }}

        /* Selection highlight */
        .selection-highlight {{
            stroke: #ffd700;
            stroke-width: 2px;
        }}

        /* Footer logo */
        .footer-logo {{
            position: fixed;
            bottom: 10px;
            right: 10px;
            opacity: 0.7;
            transition: opacity 0.2s, transform 0.2s;
            z-index: 100;
            font-size: 11px;
            letter-spacing: 0.18em;
            text-transform: uppercase;
            padding: 6px 10px;
            border-radius: 999px;
            background: rgba(135, 0, 82, 0.12);
            color: var(--text-color);
            border: 1px solid var(--border-color);
            display: inline-flex;
            align-items: center;
            gap: 10px;
        }}
        .footer-logo:hover {{
            opacity: 1;
            transform: translateY(-1px);
        }}
        .footer-link {{
            letter-spacing: 0.12em;
            text-decoration: none;
            color: var(--text-color);
            opacity: 0.75;
            border-left: 1px solid var(--border-color);
            padding-left: 10px;
        }}
        .footer-link:hover {{ opacity: 1; }}
        .loading-overlay {{
            position: fixed;
            inset: 0;
            background:
                radial-gradient(140px 140px at 50% 28%, rgba(135, 0, 82, 0.35), rgba(0, 0, 0, 0)),
                radial-gradient(520px 360px at 50% 60%, rgba(255, 135, 111, 0.18), rgba(0, 0, 0, 0)),
                linear-gradient(180deg, #0b0508 0%, #14070f 60%, #1a0a14 100%);
            color: #f5dbe7;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            gap: 10px;
            z-index: 1000;
            text-transform: uppercase;
            letter-spacing: 0.18em;
        }}
        .loading-cloud {{
            position: relative;
            width: 140px;
            height: 90px;
            filter: drop-shadow(0 0 12px rgba(135, 0, 82, 0.6));
        }}
        .loading-dot {{
            position: absolute;
            width: 8px;
            height: 8px;
            border-radius: 50%;
            background: radial-gradient(circle at 30% 30%, #ffe3ef 0%, #ffb39f 45%, #870052 100%);
            opacity: 0.9;
            box-shadow: 0 0 8px rgba(135, 0, 82, 0.9);
            animation: drift 2.6s ease-in-out infinite;
        }}
        .loading-dot:nth-child(1) {{ left: 10px; top: 18px; animation-delay: 0s; }}
        .loading-dot:nth-child(2) {{ left: 34px; top: 46px; animation-delay: 0.2s; }}
        .loading-dot:nth-child(3) {{ left: 62px; top: 20px; animation-delay: 0.4s; }}
        .loading-dot:nth-child(4) {{ left: 88px; top: 52px; animation-delay: 0.1s; }}
        .loading-dot:nth-child(5) {{ left: 114px; top: 28px; animation-delay: 0.3s; }}
        .loading-dot:nth-child(6) {{ left: 22px; top: 72px; animation-delay: 0.5s; }}
        .loading-dot:nth-child(7) {{ left: 72px; top: 72px; animation-delay: 0.6s; }}
        .loading-dot:nth-child(8) {{ left: 48px; top: 6px; animation-delay: 0.7s; }}
        .loading-dot:nth-child(9) {{ left: 96px; top: 8px; animation-delay: 0.8s; }}
        .loading-dot:nth-child(10) {{ left: 6px; top: 54px; animation-delay: 0.9s; }}
        .loading-dot:nth-child(11) {{ left: 126px; top: 62px; animation-delay: 1.0s; }}
        .loading-dot:nth-child(12) {{ left: 58px; top: 40px; animation-delay: 1.1s; }}
        @keyframes drift {{
            0% {{ transform: translate(0, 0) scale(1); opacity: 0.6; }}
            50% {{ transform: translate(0, -10px) scale(1.25); opacity: 1; }}
            100% {{ transform: translate(0, 0) scale(1); opacity: 0.6; }}
        }}
        .loading-text {{
            font-size: 11px;
            color: rgba(245, 219, 231, 0.7);
            font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
        }}
    </style>
</head>
<body>
    <div class="loading-overlay" id="loading-overlay">
        <div class="loading-cloud" aria-hidden="true">
            <span class="loading-dot"></span>
            <span class="loading-dot"></span>
            <span class="loading-dot"></span>
            <span class="loading-dot"></span>
            <span class="loading-dot"></span>
            <span class="loading-dot"></span>
            <span class="loading-dot"></span>
            <span class="loading-dot"></span>
            <span class="loading-dot"></span>
            <span class="loading-dot"></span>
            <span class="loading-dot"></span>
            <span class="loading-dot"></span>
        </div>
        <div class="loading-text">Loading data...</div>
    </div>
    <div class="header">
        <div class="header-title">
            <h1>{title}</h1>
            <button class="info-trigger" id="info-trigger" type="button" title="Viewer info" data-help="Open viewer notes with dataset context, control tips, and usage guidance.">Info</button>
            <div class="info-popover" id="info-popover" aria-hidden="true">
                <div class="info-content">
                    {viewer_info_html}
                    <div class="info-block">
                        <div class="info-title">Keyboard Shortcuts</div>
                        <table class="info-shortcuts-table">
                            <tr><td class="info-shortcuts-key"><kbd>Esc</kbd></td><td>Close the modal, or close gene discovery if the modal is not open.</td></tr>
                            <tr><td class="info-shortcuts-key"><kbd>/</kbd></td><td>Focus the gene search input and open gene discovery.</td></tr>
                            <tr><td class="info-shortcuts-key"><kbd>←</kbd> <kbd>→</kbd></td><td>Move to the previous or next section while the modal is open.</td></tr>
                            <tr><td class="info-shortcuts-key"><kbd>T</kbd></td><td>Toggle theme.</td></tr>
                            <tr><td class="info-shortcuts-key"><kbd>U</kbd></td><td>Toggle the UMAP panel.</td></tr>
                            <tr><td class="info-shortcuts-key"><kbd>L</kbd></td><td>Toggle the legend panel.</td></tr>
                            <tr><td class="info-shortcuts-key"><kbd>F</kbd></td><td>Fit the current modal section to view.</td></tr>
                            <tr><td class="info-shortcuts-key"><kbd>+</kbd> <kbd>=</kbd></td><td>Zoom in inside the modal.</td></tr>
                            <tr><td class="info-shortcuts-key"><kbd>-</kbd></td><td>Zoom out inside the modal.</td></tr>
                        </table>
                    </div>
                </div>
            </div>
        </div>
        <div class="controls">
            <div class="control-group">
                <label>Color:</label>
                <select id="color-select"></select>
            </div>
            <div class="control-group gene-control-group">
                <label>Gene:</label>
                <div class="gene-input-shell" id="gene-input-shell">
                    <input type="text" id="gene-input" placeholder="e.g. Cd4, Gfap..." list="gene-list" autocomplete="off" spellcheck="false" aria-expanded="false" aria-controls="gene-discovery-panel">
                    <datalist id="gene-list"></datalist>
                    <div class="gene-discovery-panel" id="gene-discovery-panel" aria-hidden="true">
                        <div class="gene-discovery-header">
                            <div class="gene-discovery-title">Gene discovery</div>
                            <div class="gene-discovery-actions">
                                <button class="gene-panel-btn" id="gene-panel-new" type="button">New panel</button>
                            </div>
                        </div>
                        <div id="gene-discovery-content"></div>
                    </div>
                </div>
            </div>
            <div class="control-group" id="expression-scale-section" style="display: none;">
                <label>Scale:</label>
                <div class="scale-controls">
                    <input type="number" id="expr-vmin" step="0.001" placeholder="min">
                    <span class="scale-sep">to</span>
                    <input type="number" id="expr-vmax" step="0.001" placeholder="max">
                    <button class="legend-btn" id="expr-auto" type="button">Auto (1-99%)</button>
                </div>
                <div class="scale-hint" id="expr-scale-hint">Auto scale: 1-99 percentile.</div>
            </div>
            <div class="control-group">
                <label>Size:</label>
                <div class="size-control">
                    <button class="size-step" id="spot-size-dec" type="button">−</button>
                    <input type="range" id="spot-size" min="0.01" max="5" step="0.01" value="{spot_size}" style="width:80px">
                    <button class="size-step" id="spot-size-inc" type="button">+</button>
                </div>
            </div>
            <button class="umap-toggle" id="umap-toggle" title="Toggle UMAP view" data-help="Show or hide the UMAP panel for a global embedding view with shared selection tools." style="display: none;">
                UMAP
            </button>
            <button class="legend-toggle active" id="legend-toggle" title="Toggle legend panel" data-help="Show or hide the legend panel with color keys, category toggles, and spotlight controls.">
                Legend
            </button>
            <button class="color-toggle" id="color-toggle" title="Toggle insights panel" data-help="Insights opens analysis tools for the current view: summary, comparison, gene panels, neighborhood patterns, and regions.">
                Insights
            </button>
            <button class="export-btn" id="screenshot-btn" title="Download screenshot" data-help="Download a PNG screenshot of the current main viewer layout.">
                Screenshot
            </button>
            <button class="theme-toggle" id="theme-toggle" title="Toggle dark/light mode" data-help="Switch between light and dark themes.">
                <span id="theme-icon">{theme_icon}</span>
            </button>
        </div>
        <div class="stats"><span id="stats-text"></span></div>
    </div>

    <div class="filter-bar" id="filter-bar"></div>

    <div class="main-container">
        <div class="content-column" id="content-column">
            <div class="grid-container" id="grid"></div>
            <div class="umap-panel dock-top-right" id="umap-panel">
                <div class="umap-header">
                    <h3>UMAP</h3>
                    <div class="umap-header-actions">
                        <button class="umap-btn" id="umap-dock-btn" title="Cycle panel corner">TR</button>
                        <button class="umap-btn" id="umap-panel-smaller" title="Smaller panel">−</button>
                        <button class="umap-btn" id="umap-panel-larger" title="Larger panel">+</button>
                    </div>
                </div>
                <div class="umap-canvas-container" id="umap-canvas-container">
                    <canvas class="umap-canvas" id="umap-canvas"></canvas>
                    <div class="umap-controls">
                        <button class="umap-btn" id="magic-wand-btn" title="Draw to select cells">Magic Wand</button>
                        <button class="umap-btn" id="clear-selection-btn" title="Clear selection">Clear</button>
                        <span style="margin-left: 6px; font-size: 11px; color: var(--muted-color);">Size:</span>
                        <div class="size-control">
                            <button class="size-step" id="umap-spot-size-dec" type="button">−</button>
                            <input type="range" id="umap-spot-size" min="0.01" max="6" step="0.01" value="2" style="width: 60px;">
                            <button class="size-step" id="umap-spot-size-inc" type="button">+</button>
                        </div>
                        <span id="umap-spot-size-label" style="font-size: 11px; min-width: 20px;">2</span>
                    </div>
                    <div class="umap-selection-info" id="umap-selection-info">No cells selected</div>
                    <div class="selection-summary umap-selection-summary" id="umap-selection-summary">
                        <div class="selection-summary-meta">Draw a region with Magic Wand to inspect selected cells.</div>
                    </div>
                </div>
            </div>
        </div>
        <div class="color-panel collapsed" id="color-panel"></div>
        <div class="legend-container" id="legend"></div>
    </div>

    <div class="modal-overlay" id="modal">
        <div class="modal-content">
            <div class="modal-header">
                <div style="display: flex; align-items: center;">
                    <h2 id="modal-title">Section</h2>
                    <span class="modal-meta" id="modal-meta"></span>
                    <span class="zoom-info" id="zoom-info">100%</span>
                    <span class="zoom-info" id="modal-rotation-info">Rot 0 deg</span>
                </div>
                <div class="modal-header-actions">
                    <button class="modal-controls-toggle" id="modal-controls-toggle" type="button">Hide tools</button>
                    <button class="modal-close" id="modal-close">&times;</button>
                </div>
            </div>
            <div class="modal-body">
                <div class="modal-canvas-container" id="modal-canvas-container">
                    <canvas class="modal-canvas" id="modal-canvas"></canvas>
                    <div class="modal-selection-info" id="modal-selection-info">No cells selected</div>
                    <div class="selection-summary modal-selection-summary" id="modal-selection-summary">
                        <div class="selection-summary-meta">Draw a region with Magic Wand to inspect selected cells.</div>
                    </div>
                    <div class="modal-blend-panel" id="modal-blend-panel">
                        <div class="modal-blend-title">Variable Slider</div>
                        <div class="modal-blend-loading" id="modal-blend-loading" role="status" aria-live="polite"></div>
                        <div class="modal-blend-row">
                            <span class="modal-blend-side">A</span>
                            <select id="modal-blend-a-kind">
                                <option value="cell">Cell type</option>
                                <option value="gene">Gene</option>
                            </select>
                            <div class="modal-blend-fields">
                                <select id="modal-blend-a-color"></select>
                                <select id="modal-blend-a-category"></select>
                                <input type="text" id="modal-blend-a-gene" list="gene-list" placeholder="Gene symbol">
                            </div>
                        </div>
                        <div class="modal-blend-row modal-blend-scale-row" id="modal-blend-a-scale-row" style="display: none;">
                            <span></span>
                            <span class="modal-blend-scale-label">Scale</span>
                            <div class="modal-blend-scale-controls">
                                <input type="number" id="modal-blend-a-vmin" step="0.001" placeholder="min">
                                <span class="scale-sep">to</span>
                                <input type="number" id="modal-blend-a-vmax" step="0.001" placeholder="max">
                                <button class="legend-btn" id="modal-blend-a-auto" type="button">Auto (1-99%)</button>
                            </div>
                        </div>
                        <div class="modal-blend-row">
                            <span class="modal-blend-side">B</span>
                            <select id="modal-blend-b-kind">
                                <option value="cell">Cell type</option>
                                <option value="gene">Gene</option>
                            </select>
                            <div class="modal-blend-fields">
                                <select id="modal-blend-b-color"></select>
                                <select id="modal-blend-b-category"></select>
                                <input type="text" id="modal-blend-b-gene" list="gene-list" placeholder="Gene symbol">
                            </div>
                        </div>
                        <div class="modal-blend-row modal-blend-scale-row" id="modal-blend-b-scale-row" style="display: none;">
                            <span></span>
                            <span class="modal-blend-scale-label">Scale</span>
                            <div class="modal-blend-scale-controls">
                                <input type="number" id="modal-blend-b-vmin" step="0.001" placeholder="min">
                                <span class="scale-sep">to</span>
                                <input type="number" id="modal-blend-b-vmax" step="0.001" placeholder="max">
                                <button class="legend-btn" id="modal-blend-b-auto" type="button">Auto (1-99%)</button>
                            </div>
                        </div>
                        <div class="modal-blend-row">
                            <span class="modal-blend-side">Split</span>
                            <span id="modal-blend-mix-label" style="font-size: 10px; color: var(--muted-color);">A left 50% • B right 50%</span>
                            <input type="range" id="modal-blend-mix" min="0" max="100" step="1" value="50">
                        </div>
                        <div class="modal-blend-meta" id="modal-blend-meta"></div>
                    </div>
                    <div class="modal-annotation-panel" id="modal-annotation-panel">
                        <div class="modal-annotation-title">Polygon Annotations</div>
                        <div class="modal-annotation-actions">
                            <button id="modal-annotations-export" type="button" title="Export all polygon annotations as JSON">Export JSON</button>
                            <button id="modal-annotations-clear-section" type="button" title="Delete annotations in this section">Clear section</button>
                            <button id="modal-annotations-clear-all" type="button" title="Delete all annotations">Clear all</button>
                        </div>
                        <div class="modal-annotation-list" id="modal-annotation-list">
                            <div class="modal-annotation-empty">Enable Annotate and draw a region to create polygon annotations.</div>
                        </div>
                    </div>
                    <div class="modal-controls">
                        <div class="modal-control-group" data-modal-group="view">
                            <div class="modal-control-group-title">View</div>
                            <div class="modal-control-group-body">
                                <button id="zoom-in" type="button" title="Zoom in">+</button>
                                <button id="zoom-out" type="button" title="Zoom out">−</button>
                                <button id="zoom-reset" type="button" title="Reset zoom and pan">Fit</button>
                                <button id="modal-rotate-left" type="button" title="Rotate section −45°">↺</button>
                                <button id="modal-rotate-right" type="button" title="Rotate section +45°">↻</button>
                                <button id="modal-rotate-reset" type="button" title="Reset section rotation">0°</button>
                            </div>
                        </div>
                        <div class="modal-control-group" data-modal-group="tools">
                            <div class="modal-control-group-title">Tools</div>
                            <div class="modal-control-group-body">
                                <button class="graph-toggle" id="modal-blend-toggle" type="button" title="Split the view between two selected variables">Split</button>
                                <button class="graph-toggle" id="modal-magic-wand-btn" type="button" title="Draw to select cells in this section">Select</button>
                                <button class="graph-toggle" id="modal-annotate-btn" type="button" title="Draw persistent polygon annotations in this section">Annotate</button>
                                <button class="graph-toggle" id="modal-type-toggle" type="button" title="Select a category by clicking a cell">Pick type</button>
                                <button class="graph-toggle" id="modal-type-clear" type="button" title="Clear selected type" hidden>Clear type</button>
                            </div>
                        </div>
                        <div class="modal-control-group" data-modal-group="graph" hidden>
                            <div class="modal-control-group-title">Graph</div>
                            <div class="modal-control-group-body">
                                <button class="graph-toggle" id="modal-graph-toggle" type="button" title="Toggle neighborhood graph">Graph</button>
                                <button class="graph-toggle" id="modal-neighbor-hover-toggle" type="button" title="Toggle neighbor rings on hover">Neighbors</button>
                                <select id="modal-neighbor-hop-select" title="Neighbor hop display" hidden>
                                    <option value="1">1-hop</option>
                                    <option value="2">2-hop</option>
                                    <option value="3">3-hop</option>
                                    <option value="all" selected>All hops</option>
                                </select>
                            </div>
                        </div>
                        <div class="modal-control-group" data-modal-group="actions">
                            <div class="modal-control-group-title">Actions</div>
                            <div class="modal-control-group-body">
                                <button class="graph-toggle" id="modal-screenshot-btn" type="button" title="Download screenshot of this sample view">Screenshot</button>
                                <button class="graph-toggle" id="modal-clear-selection-btn" type="button" title="Clear selected cells" hidden>Clear sel.</button>
                                <button class="graph-toggle" id="modal-exit-subview-btn" type="button" title="Return to the full section view" hidden>Back view</button>
                                <div class="modal-size-block">
                                    <div class="size-control">
                                        <button class="size-step" id="modal-spot-size-dec" type="button">−</button>
                                        <input type="range" id="modal-spot-size" min="0.01" max="5" step="0.01" value="{spot_size}" style="width: 60px;">
                                        <button class="size-step" id="modal-spot-size-inc" type="button">+</button>
                                    </div>
                                    <span class="modal-size-value" id="modal-spot-size-label">{spot_size}</span>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                <div class="modal-legend" id="modal-legend"></div>
            </div>
        </div>
    </div>

    <div class="cell-tooltip" id="cell-tooltip"></div>
    <div class="help-tooltip" id="help-tooltip" aria-hidden="true"></div>

    <script>
    (function() {{
        const LOADING_WARNING_ID = 'karospace-loading-warning';
        const LOADING_WARNING_AUTO_DISMISS_MS = 12000;
        let loadingWarningTimer = null;

        function hide() {{
            const loader = document.getElementById('loading-overlay');
            if (loader) loader.style.display = 'none';
        }}

        function removeLoadingWarning() {{
            if (loadingWarningTimer) {{
                clearTimeout(loadingWarningTimer);
                loadingWarningTimer = null;
            }}
            const warning = document.getElementById(LOADING_WARNING_ID);
            if (warning) warning.remove();
        }}

        function showError(msg, id = null) {{
            hide();
            if (id) {{
                const existing = document.getElementById(id);
                if (existing) existing.remove();
            }}
            const el = document.createElement('div');
            if (id) el.id = id;
            el.style.position = 'fixed';
            el.style.left = '16px';
            el.style.right = '16px';
            el.style.bottom = '16px';
            el.style.padding = '10px 12px';
            el.style.background = 'rgba(20, 7, 15, 0.92)';
            el.style.border = '1px solid rgba(255,255,255,0.15)';
            el.style.borderRadius = '10px';
            el.style.color = '#f5dbe7';
            el.style.fontSize = '12px';
            el.style.zIndex = '2000';
            el.textContent = msg;
            document.body.appendChild(el);
        }}
        function showLoadingWarning(msg) {{
            showError(msg, LOADING_WARNING_ID);
            loadingWarningTimer = setTimeout(removeLoadingWarning, LOADING_WARNING_AUTO_DISMISS_MS);
        }}
        function formatAsyncError(reason) {{
            if (!reason) return 'Unknown async error';
            if (typeof reason === 'string') return reason;
            if (reason && typeof reason.message === 'string' && reason.message) return reason.message;
            try {{
                return String(reason);
            }} catch (_) {{
                return 'Unknown async error';
            }}
        }}
        window.__karospaceShowLoadingWarning = showLoadingWarning;
        window.__karospaceRemoveLoadingWarning = removeLoadingWarning;
        window.__karospaceShowStartupError = showError;
        window.addEventListener('error', (e) => {{
            showError(`KaroSpace failed to start: ${{e.message || 'Unknown error'}} (open DevTools console).`);
        }});
        window.addEventListener('unhandledrejection', (e) => {{
            const detail = formatAsyncError(e.reason);
            console.error('Unhandled promise rejection:', e.reason);
            showError(`KaroSpace async error: ${{detail}} (open DevTools console).`);
        }});
        // Fallback: never keep the loader up forever.
        setTimeout(() => {{
            const loader = document.getElementById('loading-overlay');
            if (loader && loader.style.display !== 'none') {{
                showLoadingWarning('Still loading… If this persists, open DevTools console for errors.');
            }} else {{
                removeLoadingWarning();
            }}
        }}, 4000);
    }})();
    </script>

    <script id="karospace-data" type="application/json">{data_json}</script>
    <script>
    const DATA = JSON.parse(document.getElementById('karospace-data').textContent);
    DATA.genes_meta = DATA.genes_meta || {{}};
    DATA.gene_encodings = DATA.gene_encodings || {{}};
    DATA.gene_value_encodings = DATA.gene_value_encodings || {{}};
    const PALETTE = {palette_json};
    const METADATA_LABELS = {metadata_labels_json};
    const OUTLINE_BY = {outline_by_json};
    const VIEWER_INFO_HTML = {viewer_info_html_json};
    const AVAILABLE_GENE_SET = new Set(DATA.available_genes || []);
    const GENE_NAME_BY_LOWER = new Map(
        (DATA.available_genes || []).map(gene => [String(gene).toLowerCase(), gene])
    );

    const USER_AGENT = navigator.userAgent || '';
    const IS_SAFARI = /Safari/i.test(USER_AGENT) &&
        !/Chrome|Chromium|Edg|OPR|CriOS|FxiOS|Android/i.test(USER_AGENT);
    const SAFARI_DPR_CAP = 1.0;

    function getRenderDpr() {{
        const dpr = window.devicePixelRatio || 1;
        if (IS_SAFARI) return Math.max(1, Math.min(dpr, SAFARI_DPR_CAP));
        return dpr;
    }}

    const ROTATION_STEP_DEG = 45;
    const ROTATION_EPSILON = 1e-6;

    function normalizeRotationDeg(value) {{
        const numeric = Number(value);
        if (!Number.isFinite(numeric)) return 0;
        let normalized = numeric % 360;
        if (normalized < 0) normalized += 360;
        if (Math.abs(normalized - 360) < ROTATION_EPSILON) normalized = 0;
        return normalized;
    }}

    function getSignedRotationDeg(value) {{
        let signed = normalizeRotationDeg(value);
        if (signed > 180) signed -= 360;
        if (Math.abs(signed) < ROTATION_EPSILON) signed = 0;
        return signed;
    }}

    function formatRotationAngle(value) {{
        const signed = getSignedRotationDeg(value);
        if (Math.abs(signed - Math.round(signed)) < ROTATION_EPSILON) {{
            return String(Math.round(signed));
        }}
        return signed.toFixed(2).replace(/\.?0+$/, '');
    }}

    function formatRotationLabel(value) {{
        return `Rot ${{formatRotationAngle(value)}} deg`;
    }}

    function getSectionRotationDeg(section) {{
        if (!section) return 0;
        const normalized = normalizeRotationDeg(section.rotation_deg ?? 0);
        section.rotation_deg = normalized;
        return normalized;
    }}

    function getSectionRotationMetrics(section) {{
        return getSectionRotationMetricsWithBounds(section, null);
    }}

    function getSectionRotationMetricsWithBounds(section, boundsOverride = null) {{
        const bounds = boundsOverride || ensureSectionBounds(section) || {{ xmin: 0, xmax: 0, ymin: 0, ymax: 0 }};
        const dataCenterX = (bounds.xmin + bounds.xmax) / 2;
        const dataCenterY = (bounds.ymin + bounds.ymax) / 2;
        const halfWidth = Math.max((bounds.xmax - bounds.xmin) / 2, ROTATION_EPSILON);
        const halfHeight = Math.max((bounds.ymax - bounds.ymin) / 2, ROTATION_EPSILON);
        const angleDeg = getSectionRotationDeg(section);
        const angleRad = angleDeg * Math.PI / 180;
        const cos = Math.cos(angleRad);
        const sin = Math.sin(angleRad);
        const rotatedHalfWidth = Math.max(Math.abs(cos) * halfWidth + Math.abs(sin) * halfHeight, ROTATION_EPSILON);
        const rotatedHalfHeight = Math.max(Math.abs(sin) * halfWidth + Math.abs(cos) * halfHeight, ROTATION_EPSILON);
        return {{
            angleDeg,
            angleRad,
            cos,
            sin,
            dataCenterX,
            dataCenterY,
            rotatedWidth: rotatedHalfWidth * 2,
            rotatedHeight: rotatedHalfHeight * 2,
        }};
    }}

    function createSectionViewTransform(section, options = {{}}) {{
        if (!section) return null;
        ensureSectionXY(section);
        const width = Math.max(1, Number(options.width) || 0);
        const height = Math.max(1, Number(options.height) || 0);
        const padding = Math.max(0, Number(options.padding) || 0);
        const zoom = Math.max(ROTATION_EPSILON, Number(options.zoom) || 1);
        const panX = Number(options.panX) || 0;
        const panY = Number(options.panY) || 0;
        const boundsOverride = options.boundsOverride || null;
        const metrics = getSectionRotationMetricsWithBounds(section, boundsOverride);
        const fitWidth = Math.max(ROTATION_EPSILON, width - 2 * padding);
        const fitHeight = Math.max(ROTATION_EPSILON, height - 2 * padding);
        const baseScale = Math.min(
            fitWidth / Math.max(metrics.rotatedWidth, ROTATION_EPSILON),
            fitHeight / Math.max(metrics.rotatedHeight, ROTATION_EPSILON),
        );
        const scale = Math.max(ROTATION_EPSILON, baseScale * zoom);
        const centerX = width / 2 + panX;
        const centerY = height / 2 + panY;

        const transform = {{
            width,
            height,
            padding,
            baseScale,
            scale,
            centerX,
            centerY,
            angleDeg: metrics.angleDeg,
            angleRad: metrics.angleRad,
            cos: metrics.cos,
            sin: metrics.sin,
            dataCenterX: metrics.dataCenterX,
            dataCenterY: metrics.dataCenterY,
            isModal: !!options.isModal,
        }};

        transform.dataToScreen = (x, y) => {{
            const dx = x - transform.dataCenterX;
            const dy = y - transform.dataCenterY;
            const rotatedX = dx * transform.cos - dy * transform.sin;
            const rotatedY = dx * transform.sin + dy * transform.cos;
            return {{
                x: transform.centerX + rotatedX * transform.scale,
                y: transform.centerY - rotatedY * transform.scale,
            }};
        }};
        transform.screenToData = (screenX, screenY) => {{
            const rotatedX = (screenX - transform.centerX) / transform.scale;
            const rotatedY = -(screenY - transform.centerY) / transform.scale;
            return {{
                x: transform.dataCenterX + rotatedX * transform.cos + rotatedY * transform.sin,
                y: transform.dataCenterY - rotatedX * transform.sin + rotatedY * transform.cos,
            }};
        }};
        transform.isPointVisible = (screenX, screenY, radius = 0) => {{
            return !(
                screenX < -radius ||
                screenX > transform.width + radius ||
                screenY < -radius ||
                screenY > transform.height + radius
            );
        }};
        return transform;
    }}

    // Outline color overrides (used for course by default)
    const OUTLINE_COLOR_OVERRIDES = {{
        'peak_I': 'rgba(228, 26, 28, 0.5)',
        'peak_II': 'rgba(55, 126, 184, 0.5)',
        'peak_III': 'rgba(77, 175, 74, 0.5)',
        'naive': 'rgba(152, 78, 163, 0.5)',
        'remission': 'rgba(255, 127, 0, 0.5)',
        'chronic': 'rgba(255, 255, 51, 0.5)',
        'acute': 'rgba(166, 86, 40, 0.5)',
        'control': 'rgba(153, 153, 153, 0.5)',
    }};

    function getOutlineColor(value) {{
        if (!value || !OUTLINE_BY) return null;
        if (OUTLINE_BY === 'course') {{
            if (OUTLINE_COLOR_OVERRIDES[value]) return OUTLINE_COLOR_OVERRIDES[value];
            const lowerValue = value.toLowerCase();
            for (const [key, color] of Object.entries(OUTLINE_COLOR_OVERRIDES)) {{
                if (key.toLowerCase() === lowerValue) return color;
            }}
        }}
        // Generate a consistent color for unknown values
        let hash = 0;
        for (let i = 0; i < value.length; i++) {{
            hash = value.charCodeAt(i) + ((hash << 5) - hash);
        }}
        const hue = Math.abs(hash) % 360;
        return `hsla(${{hue}}, 65%, 50%, 0.5)`;
    }}

    // State
    let currentColor = DATA.initial_color;
    let currentGene = null;
    const geneScaleOverrides = {{}};
    const geneScaleAuto = {{}};
    const GENE_SCALE_PMIN = 1;
    const GENE_SCALE_PMAX = 99;
    const GENE_SCALE_MAX_SAMPLES = 200000;
    const GENE_DISCOVERY_MAX_RESULTS = 10;
    const GENE_DISCOVERY_RECENT_LIMIT = 8;
    const GENE_DISCOVERY_SUGGESTION_GROUP_LIMIT = 6;
    const GENE_DISCOVERY_SUGGESTION_GENES_PER_GROUP = 4;
    let hiddenCategories = new Set();
    let linkedSpotlightEnabled = false;
    let spotlightPinnedCategory = null;
    let spotlightHoverCategory = null;
    let spotSize = {spot_size};
    let activeFilters = {{}};  // e.g. {{ course: new Set(['peak_I', 'peak_III']) }}
    let currentTheme = '{initial_theme}';
    let showGraph = false;
    let hoverNeighbors = null;
    let neighborHoverEnabled = false;
    let neighborHopMode = 'all';
    const MAX_HOVER_HOPS = 3;
    const HOVER_COLORS = [
        'rgba(255, 165, 0, 0.9)',
        'rgba(0, 200, 255, 0.9)',
        'rgba(255, 105, 180, 0.9)',
    ];
    const expandedAggGroups = new Set();
    const expandedNeighborGroups = new Set();
    const comparisonCountsCache = new Map();
    let celltypeTrendTarget = null;
    let interactionSourceCategory = null;
    let dotplotRenderToken = 0;
    let clusterDeGroupby = null;
    let clusterDeSourceCategory = null;
    let clusterDeReferenceCategory = null;
    let geneAuxManifest = null;
    let geneAuxManifestPromise = null;
    const geneAuxShardCache = new Map();
    const geneAuxShardPromises = new Map();
    const geneAuxBinaryShardIndexCache = new Map();
    const geneAuxBinaryShardIndexPromises = new Map();
    const geneAuxBinaryShardBufferCache = new Map();
    const geneAuxBinaryShardBufferPromises = new Map();
    const modalBlendGeneLoads = new Set();
    let geneAuxLoadingMessage = '';
    let geneDiscoveryOpen = false;
    let geneDiscoveryResults = [];
    let geneDiscoveryActiveIndex = -1;
    let recentGenes = [];
    let savedGenePanels = [];

    // Modal state
    let modalSection = null;
    let modalZoom = 1;
    let modalPanX = 0, modalPanY = 0;
    let modalSpotSize = {spot_size};
    let isDragging = false;
    let isSplitDragging = false;
    let dragStartX = 0, dragStartY = 0;
    let lastPanX = 0, lastPanY = 0;
    let modalPointerMoved = false;
    let modalTypeSelectEnabled = false;
    let modalSelectedCategory = null;
    let modalMagicWandActive = false;
    let isDrawingModalLasso = false;
    let modalLassoPath = [];  // Array of {{x, y}} points in modal canvas space
    let modalAnnotationModeActive = false;
    let isDrawingModalAnnotation = false;
    let modalAnnotationPath = [];  // Array of {{x, y}} points in modal canvas space
    let modalAnnotations = [];
    let modalNextAnnotationId = 1;
    let modalBlendEnabled = false;
    let modalBlendMix = 0.5;  // Split position from left: 0 = all B, 1 = all A
    let modalBlendMarkersCollapsed = true;
    let modalControlsCollapsed = false;
    let modalControlsCustomPosition = null;
    let isDraggingModalControls = false;
    let modalControlsDragOffsetX = 0;
    let modalControlsDragOffsetY = 0;
    let modalAnnotationPanelCustomPosition = null;
    let isDraggingModalAnnotationPanel = false;
    let modalAnnotationPanelDragOffsetX = 0;
    let modalAnnotationPanelDragOffsetY = 0;
    let modalRenderedView = null;
    let modalInteractionCommitTimer = null;
    let modalBlendRenderFrame = null;
    let modalBlendRenderCache = null;
    let modalSubview = null;
    const BLEND_ALL_CATEGORIES = '__ALL__';
    let modalBlendSpec = {{
        a: {{ kind: 'cell', color: null, category: null, gene: '' }},
        b: {{ kind: 'cell', color: null, category: null, gene: '' }},
    }};
    const MODAL_ANNOTATION_COLORS = [
        '#ff7f50', '#2a9d8f', '#ffd166', '#ef476f', '#06d6a0',
        '#118ab2', '#f4a261', '#4cc9f0', '#e63946', '#43aa8b'
    ];

    // UMAP state
    let umapVisible = false;
    let umapZoom = 1;
    let umapPanX = 0, umapPanY = 0;
    let umapSpotSize = 2;
    let umapPanelDock = 'top-right';
    let umapPanelSize = 320;
    const UMAP_PANEL_DOCKS = ['top-right', 'bottom-right', 'bottom-left', 'top-left'];
    const UMAP_PANEL_SIZE_STEP = 24;
    const UMAP_PANEL_MIN_SIZE = 220;
    const UMAP_PANEL_MAX_SIZE = 560;
    const UMAP_GRID_RESERVE_MARGIN = 16;
    const UMAP_PANEL_STORAGE_KEY = 'spatial-viewer-umap-panel';
    let isUmapDragging = false;
    let umapDragStartX = 0, umapDragStartY = 0;
    let umapLastPanX = 0, umapLastPanY = 0;

    // Selection state
    let magicWandActive = false;
    let isDrawingLasso = false;
    let lassoPath = [];  // Array of {{x, y}} points
    let selectedCells = new Set();  // Set of "sectionId:cellIdx" strings
    let selectedCellsB = new Set();  // Second region for comparison
    let lassoModeB = false;  // Next lasso draw fills region B
    let selectionSummaryColor = DATA.initial_color;
    let selectionSummaryExpanded = false;
    let annotationDeSourceId = null;
    let annotationDeReferenceId = null;
    const ANNOTATION_DE_TOP_N = 12;
    let annotationDeTopN = ANNOTATION_DE_TOP_N;
    let annotationDeFullRunToken = 0;
    let annotationDeFullRun = null;
    const annotationDeFullCache = new Map();
    const sectionById = new Map((DATA.sections || []).map(section => [section.id, section]));

    function updateSectionRotationIndicators(sectionId = null) {{
        document.querySelectorAll('.section-panel').forEach((panel) => {{
            if (sectionId && panel.dataset.sectionId !== sectionId) return;
            const section = sectionById.get(panel.dataset.sectionId);
            const label = panel.querySelector('[data-section-rotation-label]');
            if (!section || !label) return;
            label.textContent = formatRotationLabel(getSectionRotationDeg(section));
        }});
        const modalRotationInfo = document.getElementById('modal-rotation-info');
        if (modalRotationInfo) {{
            modalRotationInfo.textContent = modalSection
                ? formatRotationLabel(getSectionRotationDeg(modalSection))
                : 'Rot 0 deg';
        }}
    }}

    function getSelectionIndicesForSection(sectionId, cells = selectedCells) {{
        if (!sectionId || !cells || cells.size === 0) return null;
        const indices = [];
        let otherSectionFound = false;
        cells.forEach((key) => {{
            const sep = key.lastIndexOf(':');
            if (sep <= 0) return;
            const keySectionId = key.slice(0, sep);
            const idx = Number(key.slice(sep + 1));
            if (!Number.isInteger(idx) || idx < 0) return;
            if (keySectionId !== sectionId) {{
                otherSectionFound = true;
                return;
            }}
            indices.push(idx);
        }});
        if (otherSectionFound || !indices.length) return null;
        indices.sort((a, b) => a - b);
        return indices;
    }}

    function expandSubviewBounds(bounds, referenceBounds = null) {{
        if (!bounds) return null;
        const ref = referenceBounds || bounds;
        const spanX = Math.max(bounds.xmax - bounds.xmin, 0);
        const spanY = Math.max(bounds.ymax - bounds.ymin, 0);
        const refSpanX = Math.max(ref.xmax - ref.xmin, ROTATION_EPSILON);
        const refSpanY = Math.max(ref.ymax - ref.ymin, ROTATION_EPSILON);
        const padX = Math.max(spanX * 0.18, refSpanX * 0.01, ROTATION_EPSILON * 1000);
        const padY = Math.max(spanY * 0.18, refSpanY * 0.01, ROTATION_EPSILON * 1000);
        return {{
            xmin: bounds.xmin - padX,
            xmax: bounds.xmax + padX,
            ymin: bounds.ymin - padY,
            ymax: bounds.ymax + padY,
        }};
    }}

    function getSectionBoundsFromIndices(section, indices) {{
        if (!section || !Array.isArray(indices) || !indices.length) return null;
        ensureSectionXY(section);
        let xmin = Infinity;
        let xmax = -Infinity;
        let ymin = Infinity;
        let ymax = -Infinity;
        for (let k = 0; k < indices.length; k++) {{
            const idx = indices[k];
            const x = Number(section.x[idx]);
            const y = Number(section.y[idx]);
            if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
            if (x < xmin) xmin = x;
            if (x > xmax) xmax = x;
            if (y < ymin) ymin = y;
            if (y > ymax) ymax = y;
        }}
        if (!Number.isFinite(xmin) || !Number.isFinite(xmax) || !Number.isFinite(ymin) || !Number.isFinite(ymax)) {{
            return null;
        }}
        return {{ xmin, xmax, ymin, ymax }};
    }}

    function buildModalSubviewStateFromSelection(section = modalSection) {{
        if (!section || selectedCellsB.size > 0 || lassoModeB) return null;
        const indices = getSelectionIndicesForSection(section.id);
        if (!indices || !indices.length) return null;
        const bounds = getSectionBoundsFromIndices(section, indices);
        if (!bounds) return null;
        return {{
            sectionId: section.id,
            indices,
            indexSet: new Set(indices),
            bounds: expandSubviewBounds(bounds, ensureSectionBounds(section)),
        }};
    }}

    function getModalSubviewRenderToken() {{
        if (!modalSubview || !modalSubview.bounds) return '';
        const bounds = modalSubview.bounds;
        return [
            modalSubview.sectionId,
            modalSubview.indices.length,
            bounds.xmin.toFixed(3),
            bounds.xmax.toFixed(3),
            bounds.ymin.toFixed(3),
            bounds.ymax.toFixed(3),
        ].join(':');
    }}

    function getModalActiveBounds() {{
        if (!modalSection || !modalSubview || modalSubview.sectionId !== modalSection.id) return null;
        return modalSubview.bounds || null;
    }}

    function getModalActiveCellIndexSet(sectionId = modalSection?.id) {{
        if (!sectionId || !modalSubview || modalSubview.sectionId !== sectionId) return null;
        return modalSubview.indexSet || null;
    }}

    function filterModalCandidateIndices(section, candidateIndices = null) {{
        const activeIndexSet = getModalActiveCellIndexSet(section?.id);
        if (!activeIndexSet) return candidateIndices;
        if (!Array.isArray(candidateIndices)) return modalSubview?.indices ? modalSubview.indices.slice() : Array.from(activeIndexSet);
        return candidateIndices.filter((idx) => activeIndexSet.has(idx));
    }}

    function updateModalHeader() {{
        const titleEl = document.getElementById('modal-title');
        const metaEl = document.getElementById('modal-meta');
        if (!titleEl || !metaEl) return;
        if (!modalSection) {{
            titleEl.textContent = 'Section';
            metaEl.textContent = '';
            return;
        }}
        titleEl.textContent = modalSubview ? `${{modalSection.id}} • Focused view` : modalSection.id;
        const metaParts = Object.entries(modalSection.metadata || {{}})
            .map(([k, v]) => `${{formatMetadataLabel(k)}}: ${{v}}`);
        if (modalSubview?.indices?.length) {{
            metaParts.unshift(`Focused selection: ${{modalSubview.indices.length.toLocaleString()}} cells`);
        }}
        metaEl.textContent = metaParts.join(' | ');
    }}

    function activateModalSubviewFromSelection() {{
        if (!modalSection) return false;
        const nextSubview = buildModalSubviewStateFromSelection(modalSection);
        if (!nextSubview) return false;
        selectedCellsB.clear();
        lassoModeB = false;
        modalSubview = nextSubview;
        modalZoom = 1;
        modalPanX = 0;
        modalPanY = 0;
        invalidateModalRenderedView();
        updateModalHeader();
        updateSelectionInfo();
        renderModalSection();
        return true;
    }}

    function exitModalSubview() {{
        if (!modalSubview) return false;
        modalSubview = null;
        modalZoom = 1;
        modalPanX = 0;
        modalPanY = 0;
        invalidateModalRenderedView();
        updateModalHeader();
        updateSelectionInfo();
        renderModalSection();
        return true;
    }}

    function zoomModalToFit() {{
        if (!modalSection) return;
        modalZoom = 1;
        modalPanX = 0;
        modalPanY = 0;
        renderModalSection();
    }}

    function toggleLegendPanel() {{
        const legend = document.getElementById('legend');
        const btn = document.getElementById('legend-toggle');
        if (!legend || !btn) return;
        legend.classList.toggle('collapsed');
        btn.classList.toggle('active');
        requestAnimationFrame(renderAllSections);
    }}

    function toggleInsightsPanel() {{
        const colorToggle = document.getElementById('color-toggle');
        const colorPanel = document.getElementById('color-panel');
        if (!colorToggle || !colorPanel) return;
        colorPanel.classList.toggle('collapsed');
        colorToggle.classList.toggle('active');
        if (!colorPanel.classList.contains('collapsed')) {{
            if (document.getElementById('color-tab-compare')?.classList.contains('active')) {{
                renderClusterDE();
            }} else if (document.getElementById('color-tab-neighbors')?.classList.contains('active')) {{
                renderNeighborStats();
                renderInteractionBrowser();
            }} else if (document.getElementById('color-tab-genes')?.classList.contains('active')) {{
                if (document.getElementById('genes-tab-dotplot')?.classList.contains('active')) renderDotplot();
                if (document.getElementById('genes-tab-markers')?.classList.contains('active')) renderMarkerGenes();
            }} else if (document.getElementById('color-tab-regions')?.classList.contains('active')) {{
                renderAnnotationComparison();
            }} else {{
                renderColorAggregation();
                renderCellTypeTrend();
            }}
        }}
        requestAnimationFrame(renderAllSections);
    }}

    function navigateModalSection(step) {{
        if (!modalSection || !Array.isArray(DATA.sections) || !DATA.sections.length) return false;
        const currentIndex = DATA.sections.findIndex((section) => section.id === modalSection.id);
        if (currentIndex < 0) return false;
        const nextIndex = currentIndex + Number(step || 0);
        if (nextIndex < 0 || nextIndex >= DATA.sections.length) return false;
        openModal(DATA.sections[nextIndex].id);
        return true;
    }}

    function isEditableKeyboardTarget(target) {{
        if (!target || !(target instanceof Element)) return false;
        if (target.closest('input, textarea, select')) return true;
        if (target.closest('[contenteditable="true"]')) return true;
        return !!target.isContentEditable;
    }}

    function initKeyboardShortcuts() {{
        document.addEventListener('keydown', (event) => {{
            if (event.defaultPrevented) return;
            if (event.metaKey || event.ctrlKey || event.altKey) return;
            if (isEditableKeyboardTarget(event.target)) return;

            const key = event.key;

            if (key === 'Escape') {{
                if (modalSection) {{
                    event.preventDefault();
                    closeModal();
                    return;
                }}
                if (geneDiscoveryOpen) {{
                    event.preventDefault();
                    setGeneDiscoveryOpen(false);
                    return;
                }}
                return;
            }}

            if (key === '/') {{
                const geneInput = document.getElementById('gene-input');
                if (!geneInput) return;
                event.preventDefault();
                geneInput.focus();
                geneInput.select?.();
                setGeneDiscoveryOpen(true);
                renderGeneDiscoveryPanel();
                return;
            }}

            if (key === 'ArrowLeft') {{
                if (!modalSection) return;
                event.preventDefault();
                navigateModalSection(-1);
                return;
            }}

            if (key === 'ArrowRight') {{
                if (!modalSection) return;
                event.preventDefault();
                navigateModalSection(1);
                return;
            }}

            if (key === 't' || key === 'T') {{
                event.preventDefault();
                toggleTheme();
                return;
            }}

            if (key === 'u' || key === 'U') {{
                if (!DATA.has_umap) return;
                event.preventDefault();
                toggleUMAP();
                return;
            }}

            if (key === 'l' || key === 'L') {{
                event.preventDefault();
                toggleLegendPanel();
                return;
            }}

            if (key === 'f' || key === 'F') {{
                if (!modalSection) return;
                event.preventDefault();
                zoomModalToFit();
                return;
            }}

            if (key === '+' || key === '=') {{
                if (!modalSection) return;
                event.preventDefault();
                modalZoom = Math.min(modalZoom * 1.5, 20);
                renderModalSection();
                return;
            }}

            if (key === '-') {{
                if (!modalSection) return;
                event.preventDefault();
                modalZoom = Math.max(modalZoom / 1.5, 0.1);
                renderModalSection();
            }}
        }});
    }}

    function setModalButtonPriority(button, priority = 'default') {{
        if (!button) return;
        button.classList.remove('modal-btn-primary', 'modal-btn-muted', 'modal-btn-danger');
        if (priority === 'primary') {{
            button.classList.add('modal-btn-primary');
        }} else if (priority === 'muted') {{
            button.classList.add('modal-btn-muted');
        }} else if (priority === 'danger') {{
            button.classList.add('modal-btn-danger');
        }}
    }}

    function setModalControlsCollapsed(collapsed) {{
        modalControlsCollapsed = !!collapsed;
        const controls = document.querySelector('#modal .modal-controls');
        if (controls) controls.classList.toggle('hidden', modalControlsCollapsed);
        const toggle = document.getElementById('modal-controls-toggle');
        if (toggle) toggle.textContent = modalControlsCollapsed ? 'Show tools' : 'Hide tools';
        if (modalSection) {{
            layoutModalControls();
            layoutModalAnnotationPanel();
        }}
    }}

    function getModalControlsBounds() {{
        const container = document.getElementById('modal-canvas-container');
        const controls = container?.querySelector('.modal-controls');
        if (!container || !controls) return null;
        const containerRect = container.getBoundingClientRect();
        if (containerRect.width <= 0 || containerRect.height <= 0) return null;
        const margin = 8;
        const controlsWidth = controls.offsetWidth || Math.min(920, Math.max(220, containerRect.width - margin * 2));
        const controlsHeight = controls.offsetHeight || 72;
        return {{
            minX: margin,
            maxX: Math.max(margin, containerRect.width - controlsWidth - margin),
            minY: margin,
            maxY: Math.max(margin, containerRect.height - controlsHeight - margin),
        }};
    }}

    function clampModalControlsPosition(x, y) {{
        const bounds = getModalControlsBounds();
        if (!bounds) return {{ x, y }};
        return {{
            x: Math.min(bounds.maxX, Math.max(bounds.minX, x)),
            y: Math.min(bounds.maxY, Math.max(bounds.minY, y)),
        }};
    }}

    function applyModalControlsPosition(position) {{
        const controls = document.querySelector('#modal .modal-controls');
        if (!controls) return;
        if (!position) {{
            controls.style.left = '50%';
            controls.style.right = '';
            controls.style.top = '';
            controls.style.bottom = '12px';
            controls.style.transform = 'translateX(-50%)';
            return;
        }}
        controls.style.left = `${{position.x}}px`;
        controls.style.top = `${{position.y}}px`;
        controls.style.right = 'auto';
        controls.style.bottom = 'auto';
        controls.style.transform = 'none';
    }}

    function layoutModalControls(forceReset = false) {{
        const controls = document.querySelector('#modal .modal-controls');
        if (!controls) return;
        if (forceReset) modalControlsCustomPosition = null;
        const bounds = getModalControlsBounds();
        if (!bounds) return;
        if (modalControlsCollapsed) return;
        if (modalControlsCustomPosition) {{
            modalControlsCustomPosition = clampModalControlsPosition(
                modalControlsCustomPosition.x,
                modalControlsCustomPosition.y
            );
            applyModalControlsPosition(modalControlsCustomPosition);
            return;
        }}
        applyModalControlsPosition(null);
    }}

    function initModalControlsDragging() {{
        const controls = document.querySelector('#modal .modal-controls');
        const container = document.getElementById('modal-canvas-container');
        if (!controls || !container) return;

        const stopDragging = () => {{
            if (!isDraggingModalControls) return;
            isDraggingModalControls = false;
            controls.classList.remove('dragging');
        }};

        controls.addEventListener('mousedown', (e) => {{
            if (e.button !== 0 || modalControlsCollapsed) return;
            const target = e.target;
            if (!target || !(target instanceof Element)) return;
            if (target.closest('button, select, input, option, .size-control')) return;
            const group = target.closest('.modal-control-group');
            if (!group) return;
            const containerRect = container.getBoundingClientRect();
            const controlsRect = controls.getBoundingClientRect();
            const next = clampModalControlsPosition(
                controlsRect.left - containerRect.left,
                controlsRect.top - containerRect.top
            );
            modalControlsCustomPosition = next;
            applyModalControlsPosition(next);
            modalControlsDragOffsetX = e.clientX - controlsRect.left;
            modalControlsDragOffsetY = e.clientY - controlsRect.top;
            isDraggingModalControls = true;
            controls.classList.add('dragging');
            e.preventDefault();
        }});

        document.addEventListener('mousemove', (e) => {{
            if (!isDraggingModalControls) return;
            const containerRect = container.getBoundingClientRect();
            const next = clampModalControlsPosition(
                e.clientX - containerRect.left - modalControlsDragOffsetX,
                e.clientY - containerRect.top - modalControlsDragOffsetY
            );
            modalControlsCustomPosition = next;
            applyModalControlsPosition(next);
            layoutModalAnnotationPanel();
        }});

        document.addEventListener('mouseup', stopDragging);
        window.addEventListener('blur', stopDragging);
    }}

    function updateModalToolbarState() {{
        const controls = document.querySelector('#modal .modal-controls');
        if (!controls) return;

        const fitBtn = document.getElementById('zoom-reset');
        const zoomInBtn = document.getElementById('zoom-in');
        const zoomOutBtn = document.getElementById('zoom-out');
        const rotateLeftBtn = document.getElementById('modal-rotate-left');
        const rotateRightBtn = document.getElementById('modal-rotate-right');
        const rotateResetBtn = document.getElementById('modal-rotate-reset');
        const splitBtn = document.getElementById('modal-blend-toggle');
        const selectBtn = document.getElementById('modal-magic-wand-btn');
        const annotateBtn = document.getElementById('modal-annotate-btn');
        const pickTypeBtn = document.getElementById('modal-type-toggle');
        const clearTypeBtn = document.getElementById('modal-type-clear');
        const screenshotBtn = document.getElementById('modal-screenshot-btn');
        const clearSelectionBtn = document.getElementById('modal-clear-selection-btn');
        if (clearSelectionBtn) clearSelectionBtn.hidden = selectedCells.size === 0;
        const exitSubviewBtn = document.getElementById('modal-exit-subview-btn');
        if (exitSubviewBtn) exitSubviewBtn.hidden = !modalSubview;

        const graphGroup = controls.querySelector('[data-modal-group="graph"]');
        if (graphGroup) graphGroup.hidden = !DATA.has_neighbors;

        const typeToggleBtn = document.getElementById('modal-type-toggle');
        const typeClearBtn = document.getElementById('modal-type-clear');
        const config = getColorConfig();
        const blendActive = !!getModalBlendRuntimes(modalSection);
        const typeSelectionDisabled = !modalSection || config.is_continuous || blendActive || modalAnnotationModeActive;
        if (typeSelectionDisabled) {{
            modalSelectedCategory = null;
            modalTypeSelectEnabled = false;
        }} else if (modalSelectedCategory && !config.categories?.includes(modalSelectedCategory)) {{
            modalSelectedCategory = null;
        }}
        if (typeToggleBtn) {{
            typeToggleBtn.disabled = typeSelectionDisabled;
            typeToggleBtn.classList.toggle('active', modalTypeSelectEnabled && !typeSelectionDisabled);
        }}
        if (typeClearBtn) {{
            typeClearBtn.disabled = typeSelectionDisabled;
            typeClearBtn.hidden = typeSelectionDisabled || !modalSelectedCategory;
        }}

        const modalGraphBtn = document.getElementById('modal-graph-toggle');
        const modalNeighborBtn = document.getElementById('modal-neighbor-hover-toggle');
        const modalHopSelect = document.getElementById('modal-neighbor-hop-select');
        if (modalGraphBtn) modalGraphBtn.hidden = !DATA.has_neighbors;
        if (modalNeighborBtn) modalNeighborBtn.hidden = !DATA.has_neighbors;
        if (modalHopSelect) {{
            const showHopSelect = DATA.has_neighbors && neighborHoverEnabled;
            modalHopSelect.hidden = !showHopSelect;
            modalHopSelect.disabled = !showHopSelect;
        }}

        [
            fitBtn,
            zoomInBtn,
            zoomOutBtn,
            rotateLeftBtn,
            rotateRightBtn,
            rotateResetBtn,
            splitBtn,
            selectBtn,
            annotateBtn,
            pickTypeBtn,
            clearTypeBtn,
            screenshotBtn,
            clearSelectionBtn,
            exitSubviewBtn,
            modalGraphBtn,
            modalNeighborBtn,
        ].forEach((btn) => setModalButtonPriority(btn));

        [rotateLeftBtn, rotateRightBtn, rotateResetBtn, screenshotBtn, modalGraphBtn, modalNeighborBtn].forEach((btn) => {{
            setModalButtonPriority(btn, 'muted');
        }});

        if (clearSelectionBtn && !clearSelectionBtn.hidden) {{
            setModalButtonPriority(clearSelectionBtn, 'danger');
        }}
        if (typeClearBtn && !typeClearBtn.hidden) {{
            setModalButtonPriority(typeClearBtn, 'muted');
        }}

        if (exitSubviewBtn && !exitSubviewBtn.hidden) {{
            setModalButtonPriority(exitSubviewBtn, 'primary');
        }} else if (modalAnnotationModeActive) {{
            setModalButtonPriority(annotateBtn, 'primary');
            [selectBtn, splitBtn, pickTypeBtn].forEach((btn) => setModalButtonPriority(btn, 'muted'));
            setModalButtonPriority(fitBtn, 'muted');
        }} else if (modalMagicWandActive) {{
            setModalButtonPriority(selectBtn, 'primary');
            [annotateBtn, splitBtn, pickTypeBtn].forEach((btn) => setModalButtonPriority(btn, 'muted'));
            setModalButtonPriority(fitBtn, 'muted');
        }} else if (modalTypeSelectEnabled && !typeSelectionDisabled) {{
            setModalButtonPriority(pickTypeBtn, 'primary');
            setModalButtonPriority(splitBtn, 'muted');
        }} else if (blendActive) {{
            setModalButtonPriority(splitBtn, 'primary');
            setModalButtonPriority(fitBtn, 'muted');
            setModalButtonPriority(pickTypeBtn, 'muted');
        }} else {{
            setModalButtonPriority(fitBtn, 'primary');
            setModalButtonPriority(splitBtn, 'muted');
        }}

        if (modalSection) {{
            layoutModalControls();
            layoutModalAnnotationPanel();
        }}
    }}

    function rotateSectionBy(sectionId, deltaDeg) {{
        const section = sectionById.get(sectionId);
        if (!section) return;
        section.rotation_deg = normalizeRotationDeg(getSectionRotationDeg(section) + deltaDeg);
        updateSectionRotationIndicators(sectionId);
        renderAllSections();
        if (modalSection && modalSection.id === sectionId) renderModalSection();
    }}

    function resetSectionRotation(sectionId) {{
        const section = sectionById.get(sectionId);
        if (!section) return;
        section.rotation_deg = 0;
        updateSectionRotationIndicators(sectionId);
        renderAllSections();
        if (modalSection && modalSection.id === sectionId) renderModalSection();
    }}

    function createViewerStorage() {{
        try {{
            const storage = window.localStorage;
            if (!storage) throw new Error('localStorage unavailable');
            const probeKey = '__karospace_storage_probe__';
            storage.setItem(probeKey, '1');
            storage.removeItem(probeKey);
            return storage;
        }} catch (error) {{
            console.warn('Persistent browser storage unavailable; using in-memory fallback.', error);
            const memoryStorage = new Map();
            return {{
                getItem(key) {{
                    const storageKey = String(key);
                    return memoryStorage.has(storageKey) ? memoryStorage.get(storageKey) : null;
                }},
                setItem(key, value) {{
                    memoryStorage.set(String(key), String(value));
                }},
                removeItem(key) {{
                    memoryStorage.delete(String(key));
                }},
            }};
        }}
    }}

    const VIEWER_STORAGE = createViewerStorage();

    // Theme toggle
    function toggleTheme() {{
        currentTheme = currentTheme === 'light' ? 'dark' : 'light';
        document.documentElement.classList.remove('light', 'dark');
        document.documentElement.classList.add(currentTheme);
        document.getElementById('theme-icon').textContent = currentTheme === 'dark' ? '☀️' : '🌙';
        VIEWER_STORAGE.setItem('spatial-viewer-theme', currentTheme);
        // Re-render canvases with new background
        renderAllSections();
        if (modalSection) renderModalSection();
        if (umapVisible) renderUMAP();
    }}

    function initTheme() {{
        // Check for saved preference or use initial theme
        const saved = VIEWER_STORAGE.getItem('spatial-viewer-theme');
        if (saved && (saved === 'light' || saved === 'dark')) {{
            currentTheme = saved;
        }}
        document.documentElement.classList.add(currentTheme);
        document.getElementById('theme-icon').textContent = currentTheme === 'dark' ? '☀️' : '🌙';
        document.getElementById('theme-toggle').addEventListener('click', toggleTheme);
    }}

    function getScreenshotTimestamp() {{
        return new Date().toISOString().replace(/[:.]/g, '-');
    }}

    let html2canvasPromise = null;
    function ensureHtml2CanvasLoaded() {{
        if (typeof html2canvas === 'function') return Promise.resolve();
        if (html2canvasPromise) return html2canvasPromise;
        html2canvasPromise = new Promise((resolve, reject) => {{
            const script = document.createElement('script');
            script.src = 'https://cdn.jsdelivr.net/npm/html2canvas@1.4.1/dist/html2canvas.min.js';
            script.async = true;
            script.onload = () => resolve();
            script.onerror = () => reject(new Error('Failed to load html2canvas'));
            document.head.appendChild(script);
        }});
        return html2canvasPromise;
    }}

    function downloadCanvasImage(canvas, filename) {{
        if (!canvas) return;
        const link = document.createElement('a');
        link.href = canvas.toDataURL('image/png');
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        link.remove();
    }}

    function downloadJsonFile(data, filename) {{
        const blob = new Blob([JSON.stringify(data, null, 2)], {{ type: 'application/json' }});
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        link.remove();
        URL.revokeObjectURL(url);
    }}

    function downloadTextFile(text, filename, mimeType = 'text/plain;charset=utf-8') {{
        const blob = new Blob([String(text ?? '')], {{ type: mimeType }});
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        link.remove();
        URL.revokeObjectURL(url);
    }}

    function replaceCanvasesWithImages(root) {{
        const originals = document.querySelectorAll('canvas');
        const clones = root.querySelectorAll('canvas');
        originals.forEach((canvas, idx) => {{
            const cloneCanvas = clones[idx];
            if (!cloneCanvas || !cloneCanvas.parentNode) return;
            const img = document.createElement('img');
            img.src = canvas.toDataURL('image/png');
            const rect = canvas.getBoundingClientRect();
            img.style.width = `${{rect.width}}px`;
            img.style.height = `${{rect.height}}px`;
            img.style.display = 'block';
            img.setAttribute('width', `${{canvas.width}}`);
            img.setAttribute('height', `${{canvas.height}}`);
            cloneCanvas.parentNode.replaceChild(img, cloneCanvas);
        }});
    }}

    function screenshotFullPage() {{
        const name = `spatial-viewer-${{getScreenshotTimestamp()}}.png`;
        ensureHtml2CanvasLoaded()
            .then(() => html2canvas(document.body, {{
                backgroundColor: null,
                scale: getRenderDpr(),
                useCORS: true
            }}))
            .then(canvas => {{
                downloadCanvasImage(canvas, name);
            }})
            .catch(() => {{
                alert('Screenshot failed (offline? blocked CDN?).');
            }});
    }}

    function sanitizeFilenamePart(value) {{
        if (!value) return 'section';
        const s = String(value).trim().replace(/[^a-zA-Z0-9._-]+/g, '_').replace(/^_+|_+$/g, '');
        return s || 'section';
    }}

    function screenshotModalView() {{
        if (!modalSection) {{
            alert('Open a section first.');
            return;
        }}
        const canvas = document.getElementById('modal-canvas');
        if (!canvas) {{
            alert('Sample canvas not available.');
            return;
        }}
        const sectionName = sanitizeFilenamePart(modalSection.id);
        const name = `spatial-viewer-${{sectionName}}-${{getScreenshotTimestamp()}}.png`;
        downloadCanvasImage(canvas, name);
    }}

    // Color utilities
    function magmaRgb(t) {{
        const colors = [
            [0.001, 0.000, 0.015], [0.092, 0.047, 0.256], [0.235, 0.073, 0.386],
            [0.388, 0.100, 0.451], [0.531, 0.136, 0.430], [0.651, 0.188, 0.392],
            [0.741, 0.259, 0.331], [0.813, 0.354, 0.255], [0.870, 0.477, 0.171],
            [0.918, 0.624, 0.110], [0.987, 0.855, 0.185]
        ];
        const idx = Math.min(Math.floor(t * 10), 9);
        const frac = (t * 10) - idx;
        const c1 = colors[idx], c2 = colors[idx + 1];
        const r = c1[0] + frac * (c2[0] - c1[0]);
        const g = c1[1] + frac * (c2[1] - c1[1]);
        const b = c1[2] + frac * (c2[2] - c1[2]);
        return [Math.round(r * 255), Math.round(g * 255), Math.round(b * 255)];
    }}

    function magma(t) {{
        const rgb = magmaRgb(t);
        return `rgb(${{rgb[0]}}, ${{rgb[1]}}, ${{rgb[2]}})`;
    }}

    function getCategoryColor(idx) {{ return PALETTE[idx % PALETTE.length]; }}

    const cssColorRgbCache = new Map();
    function cssColorToRgb(color) {{
        if (!color) return [128, 128, 128];
        if (cssColorRgbCache.has(color)) return cssColorRgbCache.get(color);
        let rgb = null;
        const value = String(color).trim();
        const hex = value.match(/^#([0-9a-fA-F]{{3}}|[0-9a-fA-F]{{6}})$/);
        if (hex) {{
            const body = hex[1];
            if (body.length === 3) {{
                rgb = [
                    parseInt(body[0] + body[0], 16),
                    parseInt(body[1] + body[1], 16),
                    parseInt(body[2] + body[2], 16),
                ];
            }} else {{
                rgb = [
                    parseInt(body.slice(0, 2), 16),
                    parseInt(body.slice(2, 4), 16),
                    parseInt(body.slice(4, 6), 16),
                ];
            }}
        }} else {{
            const m = value.match(/^rgb\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*\)$/i);
            if (m) {{
                rgb = [
                    Math.max(0, Math.min(255, Math.round(parseFloat(m[1])))),
                    Math.max(0, Math.min(255, Math.round(parseFloat(m[2])))),
                    Math.max(0, Math.min(255, Math.round(parseFloat(m[3])))),
                ];
            }}
        }}
        if (!rgb) rgb = [128, 128, 128];
        cssColorRgbCache.set(color, rgb);
        return rgb;
    }}

    function rgbToCss(rgb) {{
        if (!Array.isArray(rgb) || rgb.length < 3) return 'rgb(128, 128, 128)';
        return `rgb(${{rgb[0]}}, ${{rgb[1]}}, ${{rgb[2]}})`;
    }}

    function clamp01(v) {{
        if (!Number.isFinite(v)) return 0;
        return Math.max(0, Math.min(1, v));
    }}

    function formatMetadataLabel(key) {{
        return METADATA_LABELS[key] || key.replace(/_/g, ' ');
    }}

    function formatScaleNumber(value) {{
        if (!Number.isFinite(value)) return 'n/a';
        const abs = Math.abs(value);
        if (abs >= 1000 || (abs > 0 && abs < 0.001)) return value.toExponential(2);
        return value.toFixed(3);
    }}

    function getGeneScaleRange(gene) {{
        const base = DATA.genes_meta?.[gene] || {{}};
        const autoScale = geneScaleAuto[gene];
        const overrideScale = geneScaleOverrides[gene];
        let vmin = Number.isFinite(overrideScale?.vmin) ? overrideScale.vmin
            : (Number.isFinite(autoScale?.vmin) ? autoScale.vmin : base.vmin);
        let vmax = Number.isFinite(overrideScale?.vmax) ? overrideScale.vmax
            : (Number.isFinite(autoScale?.vmax) ? autoScale.vmax : base.vmax);
        if (!Number.isFinite(vmin)) vmin = 0;
        if (!Number.isFinite(vmax)) vmax = 1;
        if (!(vmax > vmin)) {{
            if (vmax === vmin) {{
                vmin -= 1e-6;
                vmax += 1e-6;
            }} else {{
                const lo = Math.min(vmin, vmax);
                const hi = Math.max(vmin, vmax);
                vmin = lo;
                vmax = hi;
            }}
        }}
        return {{ vmin, vmax }};
    }}

    // Get current color config
    function getColorConfig() {{
        if (currentGene && DATA.genes_meta[currentGene]) {{
            const scale = getGeneScaleRange(currentGene);
            return {{
                is_continuous: true,
                categories: null,
                vmin: scale.vmin,
                vmax: scale.vmax
            }};
        }}
        return DATA.colors_meta[currentColor] || {{ is_continuous: false, categories: [], vmin: 0, vmax: 1 }};
    }}

    function getLinkedSpotlightCategory(config = getColorConfig()) {{
        if (!linkedSpotlightEnabled) return null;
        if (!config || config.is_continuous) return null;
        const categories = config.categories || [];
        const candidate = spotlightHoverCategory || spotlightPinnedCategory;
        if (!candidate) return null;
        if (!categories.includes(candidate)) return null;
        if (hiddenCategories.has(candidate)) return null;
        return candidate;
    }}

    function updateLegendSpotlightClasses(targetId = 'legend') {{
        const legend = document.getElementById(targetId);
        if (!legend) return;
        const config = getColorConfig();
        const activeSpotlight = getLinkedSpotlightCategory(config);
        legend.querySelectorAll('.legend-item').forEach(item => {{
            const cat = item.dataset.category;
            const isSpotlight = !!activeSpotlight && cat === activeSpotlight;
            const isDimmed = !!activeSpotlight && cat !== activeSpotlight;
            item.classList.toggle('spotlight', isSpotlight);
            item.classList.toggle('dimmed', isDimmed);
        }});
        const toggleBtn = document.getElementById(`${{targetId}}-spotlight-toggle`);
        if (toggleBtn) toggleBtn.classList.toggle('active', linkedSpotlightEnabled);
    }}

    function updateAllLegendSpotlightClasses() {{
        updateLegendSpotlightClasses('legend');
        updateLegendSpotlightClasses('modal-legend');
    }}

    function rerenderForSpotlightChange() {{
        renderAllSections();
        if (modalSection) renderModalSection();
        if (umapVisible) renderUMAP();
    }}

    function rerenderForExpressionScaleChange() {{
        updateExpressionScaleUI();
        renderLegend('legend');
        renderLegend('modal-legend');
        renderAllSections();
        if (modalSection) renderModalSection();
        if (umapVisible) renderUMAP();
    }}

    function formatNeighborCount(value) {{
        if (!Number.isFinite(value)) return '0';
        if (Math.abs(value - Math.round(value)) < 1e-6) return Math.round(value).toLocaleString();
        return value.toFixed(2);
    }}

    function computeGenePercentiles(gene, pmin = GENE_SCALE_PMIN, pmax = GENE_SCALE_PMAX) {{
        const samples = [];
        let seenNonZero = 0;
        let totalCells = 0;
        let totalNonZero = 0;
        DATA.sections.forEach(section => {{
            const sparse = section.genes_sparse?.[gene];
            if (sparse && typeof sparse.vb64 === 'string') {{
                const sectionCells = section.n_cells ?? section.x?.length ?? 0;
                const vals = base64ToFloat32Array(sparse.vb64);
                totalCells += sectionCells;
                totalNonZero += vals.length;
                for (let i = 0; i < vals.length; i++) {{
                    const v = vals[i];
                    if (v === null || v === undefined || Number.isNaN(v) || v === 0) continue;
                    seenNonZero += 1;
                    if (samples.length < GENE_SCALE_MAX_SAMPLES) {{
                        samples.push(v);
                    }} else {{
                        const j = Math.floor(Math.random() * seenNonZero);
                        if (j < GENE_SCALE_MAX_SAMPLES) samples[j] = v;
                    }}
                }}
                return;
            }}
            if (sparse && Array.isArray(sparse.v)) {{
                const sectionCells = section.n_cells ?? section.x?.length ?? 0;
                totalCells += sectionCells;
                totalNonZero += Array.isArray(sparse.i) ? sparse.i.length : sparse.v.length;
                for (let i = 0; i < sparse.v.length; i++) {{
                    const v = sparse.v[i];
                    if (v === null || v === undefined || Number.isNaN(v) || v === 0) continue;
                    seenNonZero += 1;
                    if (samples.length < GENE_SCALE_MAX_SAMPLES) {{
                        samples.push(v);
                    }} else {{
                        const j = Math.floor(Math.random() * seenNonZero);
                        if (j < GENE_SCALE_MAX_SAMPLES) samples[j] = v;
                    }}
                }}
                return;
            }}

            const vals = section.genes?.[gene];
            if (!vals) return;
            totalCells += vals.length;
            for (let i = 0; i < vals.length; i++) {{
                const v = vals[i];
                if (v === null || v === undefined || Number.isNaN(v)) continue;
                if (v !== 0) {{
                    totalNonZero += 1;
                    seenNonZero += 1;
                    if (samples.length < GENE_SCALE_MAX_SAMPLES) {{
                        samples.push(v);
                    }} else {{
                        const j = Math.floor(Math.random() * seenNonZero);
                        if (j < GENE_SCALE_MAX_SAMPLES) samples[j] = v;
                    }}
                }}
            }}
        }});
        if (totalCells === 0) return null;
        if (totalNonZero === 0) return {{ vmin: 0, vmax: 0, pmin, pmax }};
        if (samples.length === 0) return null;
        samples.sort((a, b) => a - b);
        const nonZeroFrac = Math.max(0, Math.min(1, totalNonZero / totalCells));
        const zeroFrac = 1 - nonZeroFrac;
        const denom = Math.max(1e-12, 1 - zeroFrac);
        const qLo = Math.max(0, Math.min(1, pmin / 100));
        const qHi = Math.max(0, Math.min(1, pmax / 100));
        const nonZeroQuantile = (q) => {{
            const idx = Math.max(0, Math.floor(q * (samples.length - 1)));
            return samples[idx];
        }};
        let vmin = qLo <= zeroFrac ? 0 : nonZeroQuantile((qLo - zeroFrac) / denom);
        let vmax = qHi <= zeroFrac ? 0 : nonZeroQuantile((qHi - zeroFrac) / denom);
        if (!Number.isFinite(vmin) || !Number.isFinite(vmax)) return null;
        if (vmin === vmax) {{
            const minAll = samples[0];
            const maxAll = samples[samples.length - 1];
            if (minAll !== maxAll) {{
                vmin = minAll;
                vmax = maxAll;
            }}
        }}
        return {{ vmin, vmax, pmin, pmax }};
    }}

    const geneDenseCache = new Map(); // key: sectionId::gene -> Float32Array

    function base64ToBytes(b64) {{
        const bin = atob(b64);
        const bytes = new Uint8Array(bin.length);
        for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
        return bytes;
    }}

    function base64ToFloat32Array(b64) {{
        const bytes = base64ToBytes(b64);
        return new Float32Array(bytes.buffer, bytes.byteOffset, Math.floor(bytes.byteLength / 4));
    }}

    function base64ToUint32Array(b64) {{
        const bytes = base64ToBytes(b64);
        return new Uint32Array(bytes.buffer, bytes.byteOffset, Math.floor(bytes.byteLength / 4));
    }}

    function base64ToUint16Array(b64) {{
        const bytes = base64ToBytes(b64);
        return new Uint16Array(bytes.buffer, bytes.byteOffset, Math.floor(bytes.byteLength / 2));
    }}

    function base64ToUint8Array(b64) {{
        const bytes = base64ToBytes(b64);
        return new Uint8Array(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    }}

    function decodeQuantizedValues(qvals, qmin, qmax, maxCode) {{
        const resolvedMin = Number.isFinite(Number(qmin)) ? Number(qmin) : 0;
        const resolvedMax = Number.isFinite(Number(qmax)) ? Number(qmax) : resolvedMin;
        const scale = resolvedMax > resolvedMin ? (resolvedMax - resolvedMin) / Number(maxCode) : 0;
        const arr = new Float32Array(qvals.length);
        for (let i = 0; i < qvals.length; i++) {{
            arr[i] = scale > 0 ? (resolvedMin + qvals[i] * scale) : resolvedMin;
        }}
        return arr;
    }}

    function hydratePackedSections() {{
        // Keep initial load fast: don't eagerly base64-decode large arrays here.
        // Decode on-demand when a section is rendered (grid/modal/UMAP).
        if (!DATA || !Array.isArray(DATA.sections)) return;
        DATA.sections.forEach(section => {{
            if (!section.colors) section.colors = {{}};
            if (!section.colors_b64) section.colors_b64 = {{}};
            if (!section._colorCache) section._colorCache = {{}};
            if (!section._edgesCache) section._edgesCache = null;
            section.rotation_deg = normalizeRotationDeg(section.rotation_deg ?? 0);
        }});
    }}

    function ensureSectionXY(section) {{
        if (!section) return false;
        if ((section.x === null || section.x === undefined) && typeof section.xb64 === 'string') {{
            section.x = base64ToFloat32Array(section.xb64);
            delete section.xb64;
        }}
        if ((section.y === null || section.y === undefined) && typeof section.yb64 === 'string') {{
            section.y = base64ToFloat32Array(section.yb64);
            delete section.yb64;
        }}
        if (section.x === null || section.x === undefined) section.x = [];
       if (section.y === null || section.y === undefined) section.y = [];
        return true;
    }}

    function ensureSectionBounds(section) {{
        if (!section) return null;
        if (
            section.bounds &&
            Number.isFinite(section.bounds.xmin) &&
            Number.isFinite(section.bounds.xmax) &&
            Number.isFinite(section.bounds.ymin) &&
            Number.isFinite(section.bounds.ymax)
        ) {{
            return section.bounds;
        }}
        ensureSectionXY(section);
        let xmin = Infinity;
        let xmax = -Infinity;
        let ymin = Infinity;
        let ymax = -Infinity;
        const n = Math.min(section.x.length, section.y.length);
        for (let i = 0; i < n; i++) {{
            const x = Number(section.x[i]);
            const y = Number(section.y[i]);
            if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
            if (x < xmin) xmin = x;
            if (x > xmax) xmax = x;
            if (y < ymin) ymin = y;
            if (y > ymax) ymax = y;
        }}
        if (!Number.isFinite(xmin) || !Number.isFinite(xmax) || !Number.isFinite(ymin) || !Number.isFinite(ymax)) {{
            section.bounds = {{ xmin: 0, xmax: 0, ymin: 0, ymax: 0 }};
            return section.bounds;
        }}
        section.bounds = {{ xmin, xmax, ymin, ymax }};
        return section.bounds;
    }}

    function clampNumber(value, minValue, maxValue) {{
        return Math.min(maxValue, Math.max(minValue, value));
    }}

    function ensureSectionSpatialIndex(section) {{
        if (!section) return null;
        if (section._spatialIndex) return section._spatialIndex;
        ensureSectionXY(section);
        const bounds = ensureSectionBounds(section);
        if (!bounds) return null;

        const n = Math.min(section.x.length, section.y.length);
        const width = Math.max(bounds.xmax - bounds.xmin, ROTATION_EPSILON);
        const height = Math.max(bounds.ymax - bounds.ymin, ROTATION_EPSILON);
        const targetPointsPerBucket = 192;
        const targetBuckets = Math.max(16, Math.ceil(n / targetPointsPerBucket));
        const aspect = width / height;
        const cols = clampNumber(Math.round(Math.sqrt(targetBuckets * aspect)), 4, 256);
        const rows = clampNumber(Math.round(targetBuckets / Math.max(cols, 1)), 4, 256);
        const cellWidth = Math.max(width / cols, ROTATION_EPSILON);
        const cellHeight = Math.max(height / rows, ROTATION_EPSILON);
        const invCellWidth = 1 / cellWidth;
        const invCellHeight = 1 / cellHeight;
        const buckets = new Map();

        for (let i = 0; i < n; i++) {{
            const x = Number(section.x[i]);
            const y = Number(section.y[i]);
            if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
            const col = clampNumber(Math.floor((x - bounds.xmin) * invCellWidth), 0, cols - 1);
            const row = clampNumber(Math.floor((y - bounds.ymin) * invCellHeight), 0, rows - 1);
            const key = row * cols + col;
            let bucket = buckets.get(key);
            if (!bucket) {{
                bucket = [];
                buckets.set(key, bucket);
            }}
            bucket.push(i);
        }}

        section._spatialIndex = {{
            xmin: bounds.xmin,
            xmax: bounds.xmax,
            ymin: bounds.ymin,
            ymax: bounds.ymax,
            cols,
            rows,
            invCellWidth,
            invCellHeight,
            buckets,
        }};
        return section._spatialIndex;
    }}

    function getBoundsFromPoints(points) {{
        if (!Array.isArray(points) || points.length === 0) return null;
        let xmin = Infinity;
        let xmax = -Infinity;
        let ymin = Infinity;
        let ymax = -Infinity;
        for (let i = 0; i < points.length; i++) {{
            const point = points[i];
            const x = Number(point?.x);
            const y = Number(point?.y);
            if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
            if (x < xmin) xmin = x;
            if (x > xmax) xmax = x;
            if (y < ymin) ymin = y;
            if (y > ymax) ymax = y;
        }}
        if (!Number.isFinite(xmin) || !Number.isFinite(xmax) || !Number.isFinite(ymin) || !Number.isFinite(ymax)) {{
            return null;
        }}
        return {{ xmin, xmax, ymin, ymax }};
    }}

    function getDataBoundsFromScreenRect(transform, left, top, right, bottom) {{
        if (!transform) return null;
        const points = [
            transform.screenToData(left, top),
            transform.screenToData(right, top),
            transform.screenToData(left, bottom),
            transform.screenToData(right, bottom),
        ];
        return getBoundsFromPoints(points);
    }}

    function querySectionSpatialIndex(section, bbox) {{
        if (!section || !bbox) return null;
        const index = ensureSectionSpatialIndex(section);
        if (!index) return null;
        if (
            bbox.xmax < index.xmin || bbox.xmin > index.xmax ||
            bbox.ymax < index.ymin || bbox.ymin > index.ymax
        ) {{
            return [];
        }}

        const minCol = clampNumber(Math.floor((bbox.xmin - index.xmin) * index.invCellWidth), 0, index.cols - 1);
        const maxCol = clampNumber(Math.floor((bbox.xmax - index.xmin) * index.invCellWidth), 0, index.cols - 1);
        const minRow = clampNumber(Math.floor((bbox.ymin - index.ymin) * index.invCellHeight), 0, index.rows - 1);
        const maxRow = clampNumber(Math.floor((bbox.ymax - index.ymin) * index.invCellHeight), 0, index.rows - 1);
        const matches = [];
        for (let row = minRow; row <= maxRow; row++) {{
            for (let col = minCol; col <= maxCol; col++) {{
                const bucket = index.buckets.get(row * index.cols + col);
                if (!bucket || !bucket.length) continue;
                for (let i = 0; i < bucket.length; i++) {{
                    matches.push(bucket[i]);
                }}
            }}
        }}
        return matches;
    }}

    function getSectionVisibleScreenCandidates(section, transform, padding = 0) {{
        if (!section || !transform) return null;
        const bbox = getDataBoundsFromScreenRect(
            transform,
            -padding,
            -padding,
            transform.width + padding,
            transform.height + padding,
        );
        return querySectionSpatialIndex(section, bbox);
    }}

    function ensureSectionUMAP(section) {{
        if (!section) return false;
        if ((section.umap_x === null || section.umap_x === undefined) && typeof section.umap_xb64 === 'string') {{
            section.umap_x = base64ToFloat32Array(section.umap_xb64);
            delete section.umap_xb64;
        }}
        if ((section.umap_y === null || section.umap_y === undefined) && typeof section.umap_yb64 === 'string') {{
            section.umap_y = base64ToFloat32Array(section.umap_yb64);
            delete section.umap_yb64;
        }}
        return true;
    }}

    function ensureSectionObsIndices(section) {{
        if (!section) return false;
        if ((section.obs_idx === null || section.obs_idx === undefined) && typeof section.obs_idxb64 === 'string') {{
            section.obs_idx = base64ToUint32Array(section.obs_idxb64);
            delete section.obs_idxb64;
        }}
        if (section.obs_idx === null || section.obs_idx === undefined) section.obs_idx = [];
        return true;
    }}

    function getSectionEdgesPacked(section) {{
        if (Array.isArray(section.edges)) return section.edges;
        if (section._edgesCache) return section._edgesCache;
        if (typeof section.edges_b64 !== 'string') return null;
        const pairs = base64ToUint32Array(section.edges_b64);
        section._edgesCache = pairs;
        return pairs;
    }}

    function getSectionColorValues(section, color) {{
        const dense = section.colors?.[color];
        if (dense) return dense;
        const b64 = section.colors_b64?.[color];
        if (typeof b64 !== 'string') return null;
        section._colorCache = section._colorCache || {{}};
        if (section._colorCache[color]) return section._colorCache[color];
        const decoded = base64ToFloat32Array(b64);
        section._colorCache[color] = decoded;
        return decoded;
    }}

    function getSectionGeneValues(section, gene) {{
        const dense = section.genes?.[gene];
        if (dense) return dense;

        const sparse = section.genes_sparse?.[gene];
        if (!sparse) return null;

        const key = `${{section.id}}::${{gene}}`;
        const cached = geneDenseCache.get(key);
        if (cached) return cached;

        const n = section.n_cells ?? section.x?.length ?? 0;
        const arr = new Float32Array(n);
        const geneMeta = DATA.genes_meta?.[gene] || null;
        if (typeof sparse.ib64 === 'string' && typeof sparse.vq16b64 === 'string') {{
            const idxs = base64ToUint32Array(sparse.ib64);
            const vals = base64ToUint16Array(sparse.vq16b64);
            const m = Math.min(idxs.length, vals.length);
            const qmin = Number(geneMeta?.vmin ?? 0);
            const qmax = Number(geneMeta?.vmax ?? qmin);
            const scale = qmax > qmin ? (qmax - qmin) / 65535.0 : 0;
            for (let k = 0; k < m; k++) {{
                const idx = idxs[k];
                if (idx < n) arr[idx] = scale > 0 ? (qmin + vals[k] * scale) : qmin;
            }}
        }} else if (typeof sparse.ib64 === 'string' && typeof sparse.vq8b64 === 'string') {{
            const idxs = base64ToUint32Array(sparse.ib64);
            const vals = base64ToUint8Array(sparse.vq8b64);
            const m = Math.min(idxs.length, vals.length);
            const qmin = Number(geneMeta?.vmin ?? 0);
            const qmax = Number(geneMeta?.vmax ?? qmin);
            const scale = qmax > qmin ? (qmax - qmin) / 255.0 : 0;
            for (let k = 0; k < m; k++) {{
                const idx = idxs[k];
                if (idx < n) arr[idx] = scale > 0 ? (qmin + vals[k] * scale) : qmin;
            }}
        }} else if (typeof sparse.ib64 === 'string' && typeof sparse.vb64 === 'string') {{
            const idxs = base64ToUint32Array(sparse.ib64);
            const vals = base64ToFloat32Array(sparse.vb64);
            const m = Math.min(idxs.length, vals.length);
            for (let k = 0; k < m; k++) {{
                const idx = idxs[k];
                if (idx < n) arr[idx] = vals[k];
            }}
        }} else if (
            (Array.isArray(sparse.i) || ArrayBuffer.isView(sparse.i)) &&
            (Array.isArray(sparse.v) || ArrayBuffer.isView(sparse.v))
        ) {{
            const m = Math.min(sparse.i.length, sparse.v.length);
            for (let k = 0; k < m; k++) {{
                const idx = sparse.i[k];
                if (idx >= 0 && idx < n) arr[idx] = sparse.v[k];
            }}
        }} else {{
            return null;
        }}
        if (Array.isArray(sparse.nan) || ArrayBuffer.isView(sparse.nan)) {{
            for (let k = 0; k < sparse.nan.length; k++) {{
                const idx = sparse.nan[k];
                if (idx >= 0 && idx < n) arr[idx] = NaN;
            }}
        }}
        geneDenseCache.set(key, arr);
        return arr;
    }}

    function decodeDenseSidecarSection(sectionEntry) {{
        if (!sectionEntry) return null;
        let arr = null;
        if (typeof sectionEntry.dq16b64 === 'string') {{
            const qvals = base64ToUint16Array(sectionEntry.dq16b64);
            arr = decodeQuantizedValues(qvals, sectionEntry.qmin, sectionEntry.qmax, 65535);
        }} else if (typeof sectionEntry.dq8b64 === 'string') {{
            const qvals = base64ToUint8Array(sectionEntry.dq8b64);
            arr = decodeQuantizedValues(qvals, sectionEntry.qmin, sectionEntry.qmax, 255);
        }} else if (typeof sectionEntry.db64 === 'string') {{
            arr = base64ToFloat32Array(sectionEntry.db64);
        }}
        if (!arr) return null;
        if (Array.isArray(sectionEntry.nan)) {{
            for (let k = 0; k < sectionEntry.nan.length; k++) {{
                const idx = sectionEntry.nan[k];
                if (idx >= 0 && idx < arr.length) arr[idx] = NaN;
            }}
        }}
        return arr;
    }}

    function setGeneLoadingState(isLoading, message = '') {{
        geneAuxLoadingMessage = isLoading ? (message || 'Loading gene expression…') : '';
        const geneInput = document.getElementById('gene-input');
        if (geneInput) {{
            geneInput.disabled = !!isLoading;
            if (isLoading) {{
                geneInput.dataset.prevPlaceholder = geneInput.placeholder || '';
                geneInput.placeholder = geneAuxLoadingMessage;
            }} else if (geneInput.dataset.prevPlaceholder !== undefined) {{
                geneInput.placeholder = geneInput.dataset.prevPlaceholder || 'e.g. Cd4, Gfap...';
                delete geneInput.dataset.prevPlaceholder;
            }}
        }}
        ['modal-blend-a-gene', 'modal-blend-b-gene'].forEach((id) => {{
            const control = document.getElementById(id);
            if (control) control.disabled = !!isLoading;
        }});
        const modalBlendLoading = document.getElementById('modal-blend-loading');
        if (modalBlendLoading) {{
            if (isLoading) {{
                modalBlendLoading.textContent = geneAuxLoadingMessage;
                modalBlendLoading.classList.add('visible');
            }} else {{
                modalBlendLoading.textContent = '';
                modalBlendLoading.classList.remove('visible');
            }}
        }}
    }}

    function hydrateGeneFromAux(gene, auxData) {{
        const geneEntry = auxData?.genes?.[gene];
        const geneMeta = auxData?.genes_meta?.[gene] || geneAuxManifest?.genes_meta?.[gene];
        if (!geneEntry || !geneMeta) return false;
        DATA.genes_meta[gene] = geneMeta;
        const encoding = auxData?.gene_encodings?.[gene] || geneAuxManifest?.gene_encodings?.[gene];
        if (encoding) {{
            DATA.gene_encodings[gene] = encoding;
        }}
        const valueEncoding = auxData?.gene_value_encodings?.[gene] || geneAuxManifest?.gene_value_encodings?.[gene];
        if (valueEncoding) {{
            DATA.gene_value_encodings[gene] = valueEncoding;
        }}

        (DATA.sections || []).forEach((section) => {{
            const sectionEntry = geneEntry.sections?.[section.id];
            if (!sectionEntry) return;
            section.genes = section.genes || {{}};
            section.genes_sparse = section.genes_sparse || {{}};
            if (typeof sectionEntry.db64 === 'string' || typeof sectionEntry.dq16b64 === 'string' || typeof sectionEntry.dq8b64 === 'string') {{
                const denseValues = decodeDenseSidecarSection({{
                    ...sectionEntry,
                    qmin: geneMeta?.vmin,
                    qmax: geneMeta?.vmax,
                }});
                if (!denseValues) return;
                section.genes[gene] = denseValues;
                if (section.genes_sparse[gene]) delete section.genes_sparse[gene];
            }} else if (Array.isArray(sectionEntry.dense)) {{
                section.genes[gene] = sectionEntry.dense;
                if (section.genes_sparse[gene]) delete section.genes_sparse[gene];
            }} else if (sectionEntry.sparse) {{
                section.genes_sparse[gene] = sectionEntry.sparse;
                if (section.genes[gene]) delete section.genes[gene];
            }}
        }});
        return true;
    }}

    function getGeneAuxSidecarFormat(manifest = null) {{
        const resolved = manifest || geneAuxManifest;
        return resolved?.gene_sidecar_format || 'json-v2';
    }}

    function getPackageSession() {{
        return window.__karospacePackageSession || null;
    }}

    function parseBinaryShardIndexBuffer(buffer) {{
        if (!(buffer instanceof ArrayBuffer) || buffer.byteLength < 16) return null;
        const view = new DataView(buffer);
        const magic = String.fromCharCode(
            view.getUint8(0),
            view.getUint8(1),
            view.getUint8(2),
            view.getUint8(3)
        );
        if (magic !== 'KSB1') {{
            throw new Error('Unsupported binary gene shard format');
        }}
        const version = view.getUint16(4, true);
        if (version !== 1) {{
            throw new Error(`Unsupported binary gene shard version: ${{version}}`);
        }}
        const geneCount = view.getUint32(8, true);
        const decoder = new TextDecoder('utf-8');
        const entries = Object.create(null);
        let offset = 16;
        for (let index = 0; index < geneCount; index++) {{
            if (offset + 2 > buffer.byteLength) return null;
            const nameLength = view.getUint16(offset, true);
            offset += 2;
            if (offset + nameLength + 18 > buffer.byteLength) return null;
            const geneName = decoder.decode(new Uint8Array(buffer, offset, nameLength));
            offset += nameLength;
            const payloadKind = view.getUint8(offset);
            offset += 1;
            offset += 1;
            const payloadOffset = Number(view.getBigUint64(offset, true));
            offset += 8;
            const payloadLength = Number(view.getBigUint64(offset, true));
            offset += 8;
            entries[geneName] = {{ payloadKind, payloadOffset, payloadLength }};
        }}
        return {{ entries, indexBytes: offset }};
    }}

    async function loadBinaryShardWholeBuffer(shardUrl) {{
        if (geneAuxBinaryShardBufferCache.has(shardUrl)) return geneAuxBinaryShardBufferCache.get(shardUrl);
        if (geneAuxBinaryShardBufferPromises.has(shardUrl)) return geneAuxBinaryShardBufferPromises.get(shardUrl);
        const promise = fetch(shardUrl, {{ credentials: 'same-origin' }})
            .then((response) => {{
                if (!response.ok) {{
                    throw new Error(`HTTP ${{response.status}} while loading binary gene shard`);
                }}
                return response.arrayBuffer();
            }})
            .then((buffer) => {{
                geneAuxBinaryShardBufferCache.set(shardUrl, buffer);
                return buffer;
            }})
            .catch((error) => {{
                geneAuxBinaryShardBufferPromises.delete(shardUrl);
                throw error;
            }});
        geneAuxBinaryShardBufferPromises.set(shardUrl, promise);
        return promise;
    }}

    async function loadBinaryShardIndex(shardUrl) {{
        if (geneAuxBinaryShardIndexCache.has(shardUrl)) return geneAuxBinaryShardIndexCache.get(shardUrl);
        if (geneAuxBinaryShardIndexPromises.has(shardUrl)) return geneAuxBinaryShardIndexPromises.get(shardUrl);

        const promise = (async () => {{
            const session = getPackageSession();
            if (window.__karospacePackageMode && session && typeof session.readRange === 'function' && session.supportsRange?.(shardUrl)) {{
                const stat = session.stat?.(shardUrl) || null;
                const maxBytes = Math.max(16, Number(stat?.uncompressedSize || 262144));
                let prefixBytes = Math.min(maxBytes, 262144);
                while (prefixBytes <= maxBytes) {{
                    const buffer = await session.readRange(shardUrl, 0, prefixBytes);
                    const parsed = parseBinaryShardIndexBuffer(buffer);
                    if (parsed) {{
                        const info = {{ ...parsed, rangeCapable: true, buffer: null }};
                        geneAuxBinaryShardIndexCache.set(shardUrl, info);
                        return info;
                    }}
                    if (prefixBytes === maxBytes) break;
                    prefixBytes = Math.min(maxBytes, prefixBytes * 2);
                }}
                throw new Error('Binary gene shard index could not be parsed from package range reads');
            }}

            const buffer = await loadBinaryShardWholeBuffer(shardUrl);
            const parsed = parseBinaryShardIndexBuffer(buffer);
            if (!parsed) {{
                throw new Error('Binary gene shard index could not be parsed');
            }}
            const info = {{ ...parsed, rangeCapable: false, buffer }};
            geneAuxBinaryShardIndexCache.set(shardUrl, info);
            return info;
        }})().catch((error) => {{
            geneAuxBinaryShardIndexPromises.delete(shardUrl);
            throw error;
        }});

        geneAuxBinaryShardIndexPromises.set(shardUrl, promise);
        return promise;
    }}

    async function readBinaryGenePayload(shardUrl, gene, indexInfo = null) {{
        const resolvedIndex = indexInfo || await loadBinaryShardIndex(shardUrl);
        const entry = resolvedIndex?.entries?.[gene];
        if (!entry) return null;
        if (resolvedIndex.rangeCapable) {{
            const session = getPackageSession();
            if (session && typeof session.readRange === 'function') {{
                return session.readRange(shardUrl, entry.payloadOffset, entry.payloadOffset + entry.payloadLength);
            }}
        }}
        const buffer = resolvedIndex.buffer || await loadBinaryShardWholeBuffer(shardUrl);
        return buffer.slice(entry.payloadOffset, entry.payloadOffset + entry.payloadLength);
    }}

    function parseBinaryGenePayload(buffer, sectionOrder, geneMeta = null) {{
        if (!(buffer instanceof ArrayBuffer)) return null;
        if (!Array.isArray(sectionOrder) || !sectionOrder.length) {{
            throw new Error('Binary gene sidecar manifest is missing section_order');
        }}
        if (buffer.byteLength < 4) {{
            throw new Error('Binary gene payload is truncated');
        }}
        const view = new DataView(buffer);
        let offset = 0;
        const sectionCount = view.getUint32(offset, true);
        offset += 4;
        const sections = {{}};
        for (let i = 0; i < sectionCount; i++) {{
            if (offset + 16 > buffer.byteLength) {{
                throw new Error('Binary gene payload section header is truncated');
            }}
            const sectionIndex = view.getUint16(offset, true);
            offset += 2;
            const sectionEncoding = view.getUint8(offset);
            offset += 1;
            offset += 1;
            const cellCount = view.getUint32(offset, true);
            offset += 4;
            const nnz = view.getUint32(offset, true);
            offset += 4;
            const nanCount = view.getUint32(offset, true);
            offset += 4;
            const sectionId = sectionOrder[sectionIndex];
            if (!sectionId) {{
                throw new Error(`Binary gene payload references unknown section index ${{sectionIndex}}`);
            }}

            if (sectionEncoding === 0) {{
                sections[sectionId] = {{ sparse: {{ i: new Uint32Array(0), v: new Float32Array(0), nan: new Uint32Array(0) }} }};
                continue;
            }}

            if (sectionEncoding === 1 || sectionEncoding === 3) {{
                const valueBytes = sectionEncoding === 1 ? cellCount : cellCount * 2;
                if (offset + valueBytes + (nanCount * 4) > buffer.byteLength) {{
                    throw new Error('Binary dense gene payload is truncated');
                }}
                const valuesBuffer = buffer.slice(offset, offset + valueBytes);
                offset += valueBytes;
                const values = sectionEncoding === 1
                    ? decodeQuantizedValues(new Uint8Array(valuesBuffer), geneMeta?.vmin, geneMeta?.vmax, 255)
                    : decodeQuantizedValues(new Uint16Array(valuesBuffer), geneMeta?.vmin, geneMeta?.vmax, 65535);
                const nanBuffer = buffer.slice(offset, offset + (nanCount * 4));
                const nanIdxs = nanCount ? new Uint32Array(nanBuffer) : new Uint32Array(0);
                offset += nanCount * 4;
                for (let k = 0; k < nanIdxs.length; k++) {{
                    const idx = nanIdxs[k];
                    if (idx < values.length) values[idx] = NaN;
                }}
                sections[sectionId] = {{ dense: values }};
                continue;
            }}

            if (sectionEncoding === 2 || sectionEncoding === 4) {{
                const indexBytes = nnz * 4;
                const valueBytes = sectionEncoding === 2 ? nnz : nnz * 2;
                const totalBytes = indexBytes + valueBytes + (nanCount * 4);
                if (offset + totalBytes > buffer.byteLength) {{
                    throw new Error('Binary sparse gene payload is truncated');
                }}
                const idxBuffer = buffer.slice(offset, offset + indexBytes);
                const idxs = nnz ? new Uint32Array(idxBuffer) : new Uint32Array(0);
                offset += indexBytes;
                const valueBuffer = buffer.slice(offset, offset + valueBytes);
                offset += valueBytes;
                const values = sectionEncoding === 2
                    ? decodeQuantizedValues(new Uint8Array(valueBuffer), geneMeta?.vmin, geneMeta?.vmax, 255)
                    : decodeQuantizedValues(new Uint16Array(valueBuffer), geneMeta?.vmin, geneMeta?.vmax, 65535);
                const nanBuffer = buffer.slice(offset, offset + (nanCount * 4));
                const nanIdxs = nanCount ? new Uint32Array(nanBuffer) : new Uint32Array(0);
                offset += nanCount * 4;
                sections[sectionId] = {{ sparse: {{ i: idxs, v: values, nan: nanIdxs }} }};
                continue;
            }}

            throw new Error(`Unsupported binary section encoding: ${{sectionEncoding}}`);
        }}
        return {{ sections }};
    }}

    function hydrateGeneFromBinary(gene, geneEntry) {{
        const geneMeta = geneAuxManifest?.genes_meta?.[gene];
        if (!geneEntry || !geneMeta) return false;
        DATA.genes_meta[gene] = geneMeta;
        const encoding = geneAuxManifest?.gene_encodings?.[gene];
        if (encoding) {{
            DATA.gene_encodings[gene] = encoding;
        }}
        const valueEncoding = geneAuxManifest?.gene_value_encodings?.[gene];
        if (valueEncoding) {{
            DATA.gene_value_encodings[gene] = valueEncoding;
        }}

        (DATA.sections || []).forEach((section) => {{
            const sectionEntry = geneEntry?.sections?.[section.id];
            if (!sectionEntry) return;
            section.genes = section.genes || {{}};
            section.genes_sparse = section.genes_sparse || {{}};
            if (sectionEntry.dense) {{
                section.genes[gene] = sectionEntry.dense;
                if (section.genes_sparse[gene]) delete section.genes_sparse[gene];
            }} else if (sectionEntry.sparse) {{
                section.genes_sparse[gene] = sectionEntry.sparse;
                if (section.genes[gene]) delete section.genes[gene];
            }}
        }});
        return true;
    }}

    async function loadGeneAuxManifest() {{
        if (geneAuxManifest) return geneAuxManifest;
        if (geneAuxManifestPromise) return geneAuxManifestPromise;
        if (!DATA.gene_aux_url) return null;
        if (window.location.protocol === 'file:' && !window.__karospacePackageMode) {{
            window.__karospaceShowStartupError?.(
                'This viewer was exported with sidecar gene loading. Open it over HTTP(S) to load additional genes.'
            );
            return null;
        }}

        setGeneLoadingState(true, 'Loading gene manifest…');
        window.__karospaceShowLoadingWarning?.('Loading gene sidecar manifest…');
        geneAuxManifestPromise = fetch(DATA.gene_aux_url, {{ credentials: 'same-origin' }})
            .then((response) => {{
                if (!response.ok) {{
                    throw new Error(`HTTP ${{response.status}} while loading gene sidecar manifest`);
                }}
                return response.json();
            }})
            .then((payload) => {{
                const format = payload?.format;
                if (!payload || (format !== 'karospace-gene-sidecar-manifest-v2' && format !== 'karospace-gene-sidecar-manifest-v3')) {{
                    throw new Error('Unsupported gene sidecar manifest format');
                }}
                if (getGeneAuxSidecarFormat(payload) === 'binary-v1' && !Array.isArray(payload.section_order)) {{
                    throw new Error('Binary gene sidecar manifest is missing section_order');
                }}
                geneAuxManifest = payload;
                return payload;
            }})
            .catch((error) => {{
                console.error('Failed to load gene sidecar manifest:', error);
                window.__karospaceShowStartupError?.(
                    `Failed to load auxiliary gene data: ${{error?.message || 'Unknown error'}}`
                );
                geneAuxManifestPromise = null;
                return null;
            }})
            .finally(() => {{
                setGeneLoadingState(false);
                window.__karospaceRemoveLoadingWarning?.();
            }});
        return geneAuxManifestPromise;
    }}

    async function loadGeneAuxShard(shardUrl) {{
        if (!shardUrl) return null;
        if (geneAuxShardCache.has(shardUrl)) return geneAuxShardCache.get(shardUrl);
        if (geneAuxShardPromises.has(shardUrl)) return geneAuxShardPromises.get(shardUrl);

        setGeneLoadingState(true, 'Loading gene expression…');
        window.__karospaceShowLoadingWarning?.('Loading requested gene expression…');
        const promise = fetch(shardUrl, {{ credentials: 'same-origin' }})
            .then((response) => {{
                if (!response.ok) {{
                    throw new Error(`HTTP ${{response.status}} while loading gene shard`);
                }}
                return response.json();
            }})
            .then((payload) => {{
                if (!payload || payload.format !== 'karospace-gene-sidecar-shard-v2') {{
                    throw new Error('Unsupported gene sidecar shard format');
                }}
                geneAuxShardCache.set(shardUrl, payload);
                return payload;
            }})
            .catch((error) => {{
                console.error('Failed to load gene shard:', error);
                window.__karospaceShowStartupError?.(
                    `Failed to load requested gene data: ${{error?.message || 'Unknown error'}}`
                );
                geneAuxShardPromises.delete(shardUrl);
                return null;
            }})
            .finally(() => {{
                setGeneLoadingState(false);
                window.__karospaceRemoveLoadingWarning?.();
            }});
        geneAuxShardPromises.set(shardUrl, promise);
        return promise;
    }}

    async function ensureGeneAvailable(gene, options = {{}}) {{
        const token = String(gene || '').trim();
        const showErrors = options.showErrors !== false;
        if (!token) return false;
        if (DATA.genes_meta?.[token]) return true;
        if (!AVAILABLE_GENE_SET.has(token)) {{
            if (showErrors) {{
                alert(`Gene "${{token}}" was not found in this dataset.`);
            }}
            return false;
        }}
        const manifest = await loadGeneAuxManifest();
        if (!manifest) return false;
        const shardUrl = manifest?.gene_to_shard?.[token];
        if (!shardUrl) {{
            if (showErrors) {{
                alert(`Gene "${{token}}" is listed in the dataset but was not indexed in the auxiliary manifest.`);
            }}
            return false;
        }}
        let hydrated = false;
        if (getGeneAuxSidecarFormat(manifest) === 'binary-v1') {{
            try {{
                const payloadBuffer = await readBinaryGenePayload(shardUrl, token);
                if (!payloadBuffer) return false;
                const geneEntry = parseBinaryGenePayload(payloadBuffer, manifest.section_order || [], manifest?.genes_meta?.[token] || null);
                hydrated = hydrateGeneFromBinary(token, geneEntry);
            }} catch (error) {{
                console.error('Failed to load binary gene payload:', error);
                if (showErrors) {{
                    window.__karospaceShowStartupError?.(
                        `Failed to load requested gene data: ${{error?.message || 'Unknown error'}}`
                    );
                }}
                return false;
            }}
        }} else {{
            const shardData = await loadGeneAuxShard(shardUrl);
            if (!shardData) return false;
            hydrated = hydrateGeneFromAux(token, shardData);
        }}
        if (!hydrated && showErrors) {{
            alert(`Gene "${{token}}" is listed in the dataset but was not found in the auxiliary gene file.`);
        }}
        return hydrated;
    }}

    async function runAsyncUIAction(label, fn) {{
        try {{
            return await fn();
        }} catch (error) {{
            console.error(`${{label}} failed`, error);
            const detail = (error && error.message) ? error.message : String(error);
            window.__karospaceShowStartupError?.(`${{label}} failed: ${{detail}} (open DevTools console).`);
            return null;
        }}
    }}

    function requestModalBlendGene(gene) {{
        const token = String(gene || '').trim();
        if (!token || DATA.genes_meta?.[token] || modalBlendGeneLoads.has(token)) return;
        modalBlendGeneLoads.add(token);
        runAsyncUIAction(`Modal split gene load (${{token}})`, async () => {{
            const ok = await ensureGeneAvailable(token, {{ showErrors: false }});
            if (ok) {{
                ensureGeneAutoScale(token);
                renderLegend('modal-legend');
                if (modalSection) renderModalSection();
            }}
        }}).finally(() => {{
            modalBlendGeneLoads.delete(token);
        }});
    }}

    function setGeneDiscoveryOpen(isOpen) {{
        geneDiscoveryOpen = !!isOpen;
        const panel = document.getElementById('gene-discovery-panel');
        const input = document.getElementById('gene-input');
        if (panel) {{
            panel.classList.toggle('active', geneDiscoveryOpen);
            panel.setAttribute('aria-hidden', geneDiscoveryOpen ? 'false' : 'true');
        }}
        if (input) {{
            if (!geneDiscoveryOpen) {{
                input.value = currentGene || '';
            }}
            input.setAttribute('aria-expanded', geneDiscoveryOpen ? 'true' : 'false');
        }}
    }}

    function renderGeneTokenButton(gene, options = {{}}) {{
        const rawToken = String(gene || '').trim();
        const token = resolveCanonicalGeneName(rawToken);
        if (!token && options.allowUnknown !== true) return '';
        const label = token || rawToken;
        if (!label) return '';
        const canActivate = !!token;
        const classes = ['gene-token-btn'];
        if (options.isActive) classes.push('active');
        if (options.isSearchActive) classes.push('search-active');
        if (!canActivate) classes.push('disabled');
        else if (!DATA.genes_meta?.[token]) classes.push('unloaded');
        const showMeta = options.showMeta !== false;
        const metaLabel = options.metaLabel !== undefined
            ? options.metaLabel
            : (!canActivate ? 'name only' : (DATA.genes_meta?.[token] ? 'loaded' : 'sidecar'));
        const metaHtml = showMeta && metaLabel
            ? `<span class="gene-token-meta">${{escapeHtml(metaLabel)}}</span>`
            : '';
        return `
            <button
                type="button"
                class="${{classes.join(' ')}}"
                ${{canActivate ? `data-gene-activate="${{escapeHtml(token)}}"` : ''}}
                title="${{escapeHtml(options.title || label)}}"
                ${{canActivate ? '' : 'disabled'}}
            >
                <span>${{escapeHtml(label)}}</span>
                ${{metaHtml}}
            </button>
        `;
    }}

    function renderGeneGoogleSearchButton(gene, options = {{}}) {{
        const label = String(gene || '').trim();
        if (!label) return '';
        const query = options.query || `${{label}} gene`;
        return `
            <button
                type="button"
                class="gene-search-btn"
                data-gene-google-search="${{escapeHtml(query)}}"
                title="${{escapeHtml(options.title || `Search Google for ${{label}}`)}}"
            >
                Google
            </button>
        `;
    }}

    function bindGeneActivateButtons(container, rerenderFn = null) {{
        if (!container) return;
        container.querySelectorAll('[data-gene-activate]').forEach((btn) => {{
            btn.addEventListener('click', async () => {{
                const gene = btn.getAttribute('data-gene-activate') || '';
                if (!gene) return;
                const ok = await activateViewerGene(gene, {{ showErrors: true }});
                if (ok && typeof rerenderFn === 'function') rerenderFn();
            }});
        }});
    }}

    function bindGeneGoogleSearchButtons(container) {{
        if (!container) return;
        container.querySelectorAll('[data-gene-google-search]').forEach((btn) => {{
            btn.addEventListener('click', () => {{
                const query = btn.getAttribute('data-gene-google-search') || '';
                if (!query) return;
                const url = `https://www.google.com/search?q=${{encodeURIComponent(query)}}`;
                window.open(url, '_blank', 'noopener,noreferrer');
            }});
        }});
    }}

    function renderGeneDiscoveryPanel() {{
        const content = document.getElementById('gene-discovery-content');
        const input = document.getElementById('gene-input');
        if (!content || !input) return;

        const query = String(input.value || '').trim();
        geneDiscoveryResults = getGeneSearchResults(query);
        if (!geneDiscoveryResults.length) {{
            geneDiscoveryActiveIndex = -1;
        }} else if (geneDiscoveryActiveIndex < 0 || geneDiscoveryActiveIndex >= geneDiscoveryResults.length) {{
            geneDiscoveryActiveIndex = 0;
        }}

        const sections = [];
        if (query) {{
            const searchRows = geneDiscoveryResults.length
                ? `<div class="gene-token-grid">${{geneDiscoveryResults.map((gene, idx) => renderGeneTokenButton(gene, {{
                    isActive: gene === currentGene,
                    isSearchActive: idx === geneDiscoveryActiveIndex,
                    title: 'Load gene into the viewer',
                }})).join('')}}</div>`
                : `<div class="gene-discovery-empty">No genes matched "${{escapeHtml(query)}}". Try a shorter token or browse suggestions.</div>`;
            sections.push(`
                <div class="gene-discovery-section">
                    <div class="gene-discovery-label">Search results</div>
                    ${{searchRows}}
                </div>
            `);
        }}

        const recentRows = recentGenes.length
            ? `<div class="gene-token-grid">${{recentGenes.map((gene) => renderGeneTokenButton(gene, {{
                isActive: gene === currentGene,
                title: 'Recently viewed gene',
            }})).join('')}}</div>`
            : '<div class="gene-discovery-empty">Recent genes will appear here after you load them.</div>';
        sections.push(`
            <div class="gene-discovery-section">
                <div class="gene-discovery-label">Recent genes</div>
                ${{recentRows}}
            </div>
        `);

        if (currentGene && DATA.gene_correlations?.[currentGene]?.length) {{
            const corrData = DATA.gene_correlations[currentGene];
            const corrRows = `<div class="gene-token-grid">${{corrData.map(({{gene, r}}) =>
                renderGeneTokenButton(gene, {{
                    isActive: gene === currentGene,
                    metaLabel: `r=${{r.toFixed(2)}}`,
                    title: `Pearson r = ${{r.toFixed(2)}} with ${{escapeHtml(currentGene)}}`,
                }})
            ).join('')}}</div>`;
            sections.push(`
                <div class="gene-discovery-section">
                    <div class="gene-discovery-label">Correlated with ${{escapeHtml(currentGene)}}</div>
                    ${{corrRows}}
                </div>
            `);
        }}

        if (DATA.spatial_variable_genes?.length) {{
            const topSVG = DATA.spatial_variable_genes.slice(0, 12);
            const svgRows = `<div class="gene-token-grid">${{topSVG.map((item) =>
                renderGeneTokenButton(item.gene, {{
                    isActive: item.gene === currentGene,
                    metaLabel: `I=${{item.I.toFixed(2)}}`,
                    title: `Moran's I = ${{item.I.toFixed(4)}}`,
                }})
            ).join('')}}</div>`;
            sections.push(`
                <div class="gene-discovery-section">
                    <div class="gene-discovery-label">Spatially variable genes</div>
                    ${{svgRows}}
                </div>
            `);
        }}

        const suggestionInfo = getGeneSuggestionGroups();
        const suggestionRows = suggestionInfo.groups.length
            ? suggestionInfo.groups.map(group => `
                <div class="gene-suggestion-group">
                    <div class="gene-suggestion-group-title">${{escapeHtml(group.category)}}</div>
                    <div class="gene-token-grid">
                        ${{group.genes.map((gene) => renderGeneTokenButton(gene, {{
                            isActive: gene === currentGene,
                            title: `Suggested marker for ${{group.category}}`,
                        }})).join('')}}
                    </div>
                </div>
            `).join('')
            : `<div class="gene-discovery-empty">${{escapeHtml(suggestionInfo.subtitle || 'No suggestions available.')}}</div>`;
        const hiddenSuggestionNote = suggestionInfo.hiddenCount > 0
            ? `<div class="gene-discovery-empty" style="margin-top: 6px;">+${{suggestionInfo.hiddenCount}} more categories available in Insights → Markers.</div>`
            : '';
        sections.push(`
            <div class="gene-discovery-section">
                <div class="gene-discovery-label">${{escapeHtml(suggestionInfo.title)}}</div>
                ${{suggestionRows}}
                ${{hiddenSuggestionNote}}
            </div>
        `);

        const panelRows = savedGenePanels.length
            ? savedGenePanels.map((panel) => `
                <div class="gene-panel-card">
                    <div class="gene-panel-header">
                        <div class="gene-panel-name">${{escapeHtml(panel.name)}}</div>
                        <div class="gene-panel-actions">
                            <button type="button" class="gene-panel-btn" data-gene-panel-add="${{escapeHtml(panel.name)}}" title="Add the current gene to this panel">+ current</button>
                            <button type="button" class="gene-panel-btn" data-gene-panel-delete="${{escapeHtml(panel.name)}}" title="Delete this panel">Delete</button>
                        </div>
                    </div>
                    <div class="gene-panel-genes">
                        ${{panel.genes.length
                            ? panel.genes.map((gene) => renderGeneTokenButton(gene, {{
                                isActive: gene === currentGene,
                                title: `Load ${{gene}} from ${{panel.name}}`,
                            }})).join('')
                            : '<div class="gene-discovery-empty">Panel is empty. Use + current to add a gene.</div>'
                        }}
                    </div>
                </div>
            `).join('')
            : '<div class="gene-discovery-empty">No saved panels yet. Create one from the current active gene or an exact gene in the input.</div>';
        sections.push(`
            <div class="gene-discovery-section">
                <div class="gene-discovery-label">Saved panels</div>
                ${{panelRows}}
            </div>
        `);

        content.innerHTML = sections.join('');
    }}

    async function activateViewerGene(gene, options = {{}}) {{
        const rawToken = String(gene || '').trim();
        const token = resolveCanonicalGeneName(rawToken);
        const showErrors = options.showErrors !== false;
        const geneInput = document.getElementById('gene-input');

        if (!rawToken) {{
            currentGene = null;
            hiddenCategories.clear();
            if (geneInput) geneInput.value = '';
            updateExpressionScaleUI();
            renderLegend('legend');
            renderLegend('modal-legend');
            renderAllSections();
            if (modalSection) renderModalSection();
            if (umapVisible) renderUMAP();
            updateSelectionInfo();
            recentGenes = loadRecentGenes();
            savedGenePanels = loadSavedGenePanels();
            renderGeneDiscoveryPanel();
            return true;
        }}

        if (!token) {{
            if (showErrors) {{
                alert(`Gene "${{rawToken}}" was not found in this dataset.`);
            }}
            renderGeneDiscoveryPanel();
            return false;
        }}

        if (geneInput) geneInput.value = token;
        const ok = await runAsyncUIAction('Gene selection', async () => {{
            if (!(await ensureGeneAvailable(token, {{ showErrors }}))) {{
                return false;
            }}
            currentGene = token;
            geneDenseCache.clear();
            modalSelectedCategory = null;
            modalTypeSelectEnabled = false;
            hiddenCategories.clear();
            ensureGeneAutoScale(currentGene);
            updateExpressionScaleUI();
            renderLegend('legend');
            renderLegend('modal-legend');
            renderAllSections();
            if (modalSection) renderModalSection();
            if (umapVisible) renderUMAP();
            updateSelectionInfo();
            return true;
        }});
        if (ok) {{
            recordRecentGene(token);
            geneDiscoveryResults = [];
            geneDiscoveryActiveIndex = -1;
            renderGeneDiscoveryPanel();
            setGeneDiscoveryOpen(false);
            return true;
        }}
        renderGeneDiscoveryPanel();
        return false;
    }}

    function ensureGeneAutoScale(gene) {{
        if (!gene) return;
        if (!geneScaleAuto[gene]) {{
            const autoScale = computeGenePercentiles(gene);
            if (autoScale) geneScaleAuto[gene] = autoScale;
        }}
    }}

    function updateExpressionScaleUI() {{
        const section = document.getElementById('expression-scale-section');
        const vminInput = document.getElementById('expr-vmin');
        const vmaxInput = document.getElementById('expr-vmax');
        const hint = document.getElementById('expr-scale-hint');
        if (!section || !vminInput || !vmaxInput || !hint) return;
        if (!currentGene) {{
            section.style.display = 'none';
            return;
        }}
        section.style.display = 'block';
        const config = getColorConfig();
        vminInput.value = Number.isFinite(config.vmin) ? config.vmin.toFixed(3) : '';
        vmaxInput.value = Number.isFinite(config.vmax) ? config.vmax.toFixed(3) : '';
        if (geneScaleOverrides[currentGene]) {{
            hint.textContent = 'Custom scale (manual).';
        }} else if (geneScaleAuto[currentGene]) {{
            hint.textContent = `Auto scale: ${{GENE_SCALE_PMIN}}-${{GENE_SCALE_PMAX}} percentile.`;
        }} else {{
            hint.textContent = 'Auto scale unavailable; using data range.';
        }}
    }}

    function parseGeneList(text) {{
        if (!text) return [];
        const parts = String(text)
            .split(/[,\s]+/)
            .map(s => s.trim())
            .filter(Boolean);
        const seen = new Set();
        const genes = [];
        parts.forEach(g => {{
            if (seen.has(g)) return;
            if (AVAILABLE_GENE_SET.has(g)) {{
                seen.add(g);
                genes.push(g);
            }}
        }});
        return genes;
    }}

    function autocompleteDotplotGeneToken(inputEl) {{
        if (!inputEl) return false;
        const value = inputEl.value || '';
        const cursor = Number.isFinite(inputEl.selectionStart) ? inputEl.selectionStart : value.length;
        const left = value.slice(0, cursor);
        const right = value.slice(cursor);
        const tokenMatch = left.match(/^(.*?)([^,\s]*)$/);
        if (!tokenMatch) return false;
        const prefix = tokenMatch[1] || '';
        const token = tokenMatch[2] || '';
        if (!token) return false;

        const lower = token.toLowerCase();
        const match = (DATA.available_genes || []).find(g => g.toLowerCase().startsWith(lower));
        if (!match || match === token) return false;

        inputEl.value = `${{prefix}}${{match}}${{right}}`;
        const nextCursor = (prefix + match).length;
        inputEl.setSelectionRange(nextCursor, nextCursor);
        return true;
    }}

    function getCategoricalColorColumns() {{
        const cols = (DATA.available_colors || []).filter(col => {{
            const meta = DATA.colors_meta?.[col];
            return meta && !meta.is_continuous && Array.isArray(meta.categories);
        }});
        cols.sort((a, b) => a.localeCompare(b));
        return cols;
    }}

    function getCategoriesForColorColumn(colorCol) {{
        const meta = DATA.colors_meta?.[colorCol];
        if (!meta || meta.is_continuous || !Array.isArray(meta.categories)) return [];
        return meta.categories;
    }}

    function getModalBlendVariableLabel(spec) {{
        if (!spec) return 'Unknown';
        if (spec.kind === 'gene') return spec.gene ? `Gene: ${{spec.gene}}` : 'Gene';
        const colLabel = spec.color ? formatMetadataLabel(spec.color) : 'Cell type';
        if (spec.category === BLEND_ALL_CATEGORIES) return `${{colLabel}}: All`;
        return spec.category ? `${{colLabel}}: ${{spec.category}}` : colLabel;
    }}

    function getModalSplitLegendEntry(side) {{
        const spec = modalBlendSpec?.[side];
        if (!spec) return null;
        const sideLabel = side === 'a' ? 'A (left)' : 'B (right)';
        if (spec.kind === 'gene') {{
            const geneLabel = spec.gene ? `Gene: ${{spec.gene}}` : 'Gene';
            const scale = getGeneScaleRange(spec.gene);
            return {{
                side,
                sideLabel,
                label: geneLabel,
                detail: `Expression (${{formatScaleNumber(scale.vmin)}} to ${{formatScaleNumber(scale.vmax)}})`,
                swatchClass: 'split-legend-swatch-bar',
                swatchBackground: `linear-gradient(90deg, ${{magma(0)}}, ${{magma(0.5)}}, ${{magma(1)}})`,
                categories: null,
            }};
        }}

        const colorCol = spec.color;
        const categories = getCategoriesForColorColumn(colorCol);
        const colLabel = colorCol ? formatMetadataLabel(colorCol) : 'Cell type';
        const categoryEntries = categories.map((cat, idx) => ({{
            label: cat,
            color: getCategoryColor(idx),
        }}));
        if (spec.category && spec.category !== BLEND_ALL_CATEGORIES) {{
            const catIdx = categories.indexOf(spec.category);
            if (catIdx >= 0) {{
                return {{
                    side,
                    sideLabel,
                    label: `${{colLabel}}: ${{spec.category}}`,
                    detail: 'Selected category',
                    swatchClass: 'split-legend-swatch-dot',
                    swatchBackground: getCategoryColor(catIdx),
                    categories: [categoryEntries[catIdx]],
                }};
            }}
        }}

        const gradientStops = categories
            .slice(0, 5)
            .map((_, idx) => getCategoryColor(idx));
        return {{
            side,
            sideLabel,
            label: `${{colLabel}}: All`,
            detail: categories.length ? `${{categories.length}} categories` : 'All categories',
            swatchClass: 'split-legend-swatch-bar',
            swatchBackground: gradientStops.length > 1
                ? `linear-gradient(90deg, ${{gradientStops.join(', ')}})`
                : (gradientStops[0] || '#9ca3af'),
            categories: categoryEntries,
        }};
    }}

    function getModalBlendRuntime(section, spec) {{
        if (!section || !spec) return null;
        if (spec.kind === 'gene') {{
            const gene = (spec.gene || '').trim();
            if (!gene) return null;
            if (!DATA.genes_meta?.[gene]) {{
                requestModalBlendGene(gene);
                return null;
            }}
            const values = getSectionGeneValues(section, gene);
            if (!values) return null;
            const scale = getGeneScaleRange(gene);
            return {{ kind: 'gene', values, vmin: scale.vmin, vmax: scale.vmax }};
        }}

        const colorCol = spec.color;
        const categories = getCategoriesForColorColumn(colorCol);
        if (!colorCol || categories.length === 0) return null;
        const values = getSectionColorValues(section, colorCol);
        if (!values) return null;
        if (spec.category === BLEND_ALL_CATEGORIES) {{
            return {{
                kind: 'cell-all',
                values,
                paletteRgb: categories.map((_, idx) => cssColorToRgb(getCategoryColor(idx))),
            }};
        }}
        const catIdx = categories.indexOf(spec.category);
        if (catIdx < 0) return null;
        return {{
            kind: 'cell',
            values,
            catIdx,
            activeRgb: cssColorToRgb(getCategoryColor(catIdx)),
            inactiveRgb: [176, 176, 176],
        }};
    }}

    function getModalBlendRuntimes(section) {{
        if (!modalBlendEnabled || !section) return null;
        const a = getModalBlendRuntime(section, modalBlendSpec.a);
        const b = getModalBlendRuntime(section, modalBlendSpec.b);
        if (!a || !b) return null;
        return {{ a, b }};
    }}

    function getModalBlendCellRgb(runtime, idx) {{
        if (!runtime) return [128, 128, 128];
        const raw = runtime.values?.[idx];
        if (runtime.kind === 'cell-all') {{
            if (!Number.isFinite(raw)) return [160, 160, 160];
            const catIdx = Math.round(raw);
            return runtime.paletteRgb?.[catIdx] || [160, 160, 160];
        }}
        if (runtime.kind === 'cell') {{
            if (!Number.isFinite(raw)) return runtime.inactiveRgb;
            return Math.round(raw) === runtime.catIdx ? runtime.activeRgb : runtime.inactiveRgb;
        }}
        if (!Number.isFinite(raw)) return [140, 140, 140];
        const t = clamp01((raw - runtime.vmin) / (runtime.vmax - runtime.vmin));
        return magmaRgb(t);
    }}

    function updateDotplotAggregateValueOptions() {{
        const aggSelect = document.getElementById('dotplot-aggregate-by');
        const valueWrap = document.getElementById('dotplot-aggregate-value-wrap');
        const valueSelect = document.getElementById('dotplot-aggregate-value');
        if (!aggSelect || !valueWrap || !valueSelect) return;
        const key = aggSelect.value;
        if (!key) {{
            valueWrap.style.display = 'none';
            valueSelect.innerHTML = '';
            return;
        }}
        const values = (DATA.metadata_filters && DATA.metadata_filters[key]) ? DATA.metadata_filters[key] : [];
        valueWrap.style.display = 'block';
        const opts = ['<option value="__ALL__">All</option>']
            .concat(values.map(v => `<option value="${{v}}">${{v}}</option>`));
        valueSelect.innerHTML = opts.join('');
    }}

    async function renderDotplot() {{
        const status = document.getElementById('dotplot-status');
        const grid = document.getElementById('dotplot-grid');
        const groupbySelect = document.getElementById('dotplot-groupby');
        const genesInput = document.getElementById('dotplot-genes');
        const aggSelect = document.getElementById('dotplot-aggregate-by');
        const aggValueSelect = document.getElementById('dotplot-aggregate-value');
        if (!status || !grid || !groupbySelect || !genesInput || !aggSelect) return;

        const groupbyColor = groupbySelect.value;
        if (!groupbyColor) {{
            status.textContent = 'Pick a categorical color to group by.';
            grid.innerHTML = '';
            return;
        }}
        const meta = DATA.colors_meta?.[groupbyColor];
        if (!meta || meta.is_continuous || !Array.isArray(meta.categories)) {{
            status.textContent = 'Dotplot requires a categorical color column.';
            grid.innerHTML = '';
            return;
        }}

        const genes = parseGeneList(genesInput.value);
        if (genes.length === 0) {{
            status.textContent = 'Enter one or more genes (comma-separated).';
            grid.innerHTML = '';
            return;
        }}

        const aggKey = aggSelect.value || '';
        const aggValue = (aggKey && aggValueSelect) ? (aggValueSelect.value || '__ALL__') : '__ALL__';
        const aggLabel = aggKey ? `${{formatMetadataLabel(aggKey)}}=${{aggValue}}` : 'All sections';

        dotplotRenderToken += 1;
        const token = dotplotRenderToken;
        status.textContent = `Loading genes for dotplot (${{genes.length}})…`;
        grid.innerHTML = '';

        const missing = [];
        for (const gene of genes) {{
            if (token !== dotplotRenderToken) return;
            const ok = await ensureGeneAvailable(gene, {{ showErrors: false }});
            if (!ok) missing.push(gene);
        }}
        if (token !== dotplotRenderToken) return;
        if (missing.length) {{
            status.textContent = `Gene data unavailable: ${{missing.join(', ')}}`;
            return;
        }}

        status.textContent = `Computing dotplot (${{groupbyColor}}, genes=${{genes.length}}, ${{aggLabel}})…`;

        setTimeout(() => {{
            if (token !== dotplotRenderToken) return;

            const categories = meta.categories;
            const k = categories.length;

            // Eligible sections + pre-count totals per category once.
            const eligible = [];
            const totals = new Uint32Array(k);
            for (let s = 0; s < DATA.sections.length; s++) {{
                if (token !== dotplotRenderToken) return;
                const section = DATA.sections[s];
                if (!sectionPassesFilter(section)) continue;
                if (aggKey) {{
                    const val = section.metadata?.[aggKey] || 'unknown';
                    if (aggValue !== '__ALL__' && val !== aggValue) continue;
                }}
                const groupVals = getSectionColorValues(section, groupbyColor);
                if (!groupVals || !groupVals.length) continue;
                eligible.push({{ section, groupVals }});
                for (let i = 0; i < groupVals.length; i++) {{
                    const gv = groupVals[i];
                    if (gv === null || gv === undefined || Number.isNaN(gv)) continue;
                    const ci = Math.round(gv);
                    if (ci >= 0 && ci < k) totals[ci] += 1;
                }}
            }}

            if (eligible.length === 0) {{
                status.textContent = 'No sections match the current filters.';
                return;
            }}

            const sums = genes.map(() => new Float64Array(k));
            const nnz = genes.map(() => new Uint32Array(k));
            let usedDenseFallback = false;

            for (let g = 0; g < genes.length; g++) {{
                if (token !== dotplotRenderToken) return;
                const gene = genes[g];
                for (let e = 0; e < eligible.length; e++) {{
                    const {{ section, groupVals }} = eligible[e];
                    const sparse = section.genes_sparse?.[gene];
                    if (sparse) {{
                        if (typeof sparse.ib64 === 'string' && typeof sparse.vb64 === 'string') {{
                            const idxs = base64ToUint32Array(sparse.ib64);
                            const vals = base64ToFloat32Array(sparse.vb64);
                            const m = Math.min(idxs.length, vals.length);
                            for (let j = 0; j < m; j++) {{
                                const idx = idxs[j];
                                if (idx >= groupVals.length) continue;
                                const gv = groupVals[idx];
                                if (!Number.isFinite(gv)) continue;
                                const ci = Math.round(gv);
                                if (ci < 0 || ci >= k) continue;
                                const v = vals[j];
                                if (!Number.isFinite(v) || v === 0) continue;
                                sums[g][ci] += v;
                                nnz[g][ci] += 1;
                            }}
                            continue;
                        }}
                        if (Array.isArray(sparse.i) && Array.isArray(sparse.v)) {{
                            const m = Math.min(sparse.i.length, sparse.v.length);
                            for (let j = 0; j < m; j++) {{
                                const idx = sparse.i[j];
                                if (idx === null || idx === undefined) continue;
                                if (idx < 0 || idx >= groupVals.length) continue;
                                const gv = groupVals[idx];
                                if (!Number.isFinite(gv)) continue;
                                const ci = Math.round(gv);
                                if (ci < 0 || ci >= k) continue;
                                const v = sparse.v[j];
                                if (!Number.isFinite(v) || v === 0) continue;
                                sums[g][ci] += v;
                                nnz[g][ci] += 1;
                            }}
                            continue;
                        }}
                    }}

                    const dense = section.genes?.[gene];
                    if (dense && dense.length) {{
                        usedDenseFallback = true;
                        const n = Math.min(dense.length, groupVals.length);
                        for (let i = 0; i < n; i++) {{
                            const v = dense[i];
                            if (!Number.isFinite(v) || v === 0) continue;
                            const gv = groupVals[i];
                            if (!Number.isFinite(gv)) continue;
                            const ci = Math.round(gv);
                            if (ci < 0 || ci >= k) continue;
                            sums[g][ci] += v;
                            nnz[g][ci] += 1;
                        }}
                    }}
                }}
            }}

            if (token !== dotplotRenderToken) return;
            grid.style.setProperty('--dotplot-cols', String(genes.length));
            const header = `
                <div class="dotplot-row dotplot-header">
                    <div class="dotplot-label">Cell type</div>
                    ${{genes.map(g => `<div class="dotplot-gene" title="${{g}}">${{g}}</div>`).join('')}}
                </div>
            `;

            const rows = categories.map((cat, ci) => {{
                const total = totals[ci];
                const cells = genes.map((gene, gi) => {{
                    if (!total) return `<div class="dotplot-dot" title="No cells"></div>`;
                    const mean = sums[gi][ci] / total;
                    const frac = nnz[gi][ci] / total;
                    const vmax = (DATA.genes_meta?.[gene]?.vmax ?? 0) || 0;
                    const tRaw = vmax > 0 ? (mean / vmax) : 0;
                    const t = Math.max(0, Math.min(1, tRaw));
                    const color = magma(0.1 + 0.9 * t);
                    const r = Math.max(0.5, Math.min(8, 8 * Math.sqrt(frac)));
                    const title = `${{gene}} · mean=${{mean.toFixed(3)}} · %expr=${{(frac*100).toFixed(1)}} · n=${{total.toLocaleString()}}`;
                    return `
                        <div class="dotplot-dot" title="${{title}}">
                            <svg width="20" height="20" viewBox="0 0 20 20">
                                <circle cx="10" cy="10" r="${{r}}" fill="${{color}}" stroke="rgba(0,0,0,0.10)" stroke-width="1"></circle>
                            </svg>
                        </div>
                    `;
                }}).join('');
                return `
                    <div class="dotplot-row">
                        <div class="dotplot-label" title="${{cat}}">${{cat}}</div>
                        ${{cells}}
                    </div>
                `;
            }}).join('');

            grid.innerHTML = header + rows;
            const denseNote = usedDenseFallback ? ' (some genes were dense; may be slower)' : '';
            status.textContent = `Dotplot ready (${{eligible.length}} sections, ${{aggLabel}})${{denseNote}}.`;
        }}, 0);
    }}

    // Get values for a section
    function getSectionValues(section) {{
        if (currentGene) {{
            const geneVals = getSectionGeneValues(section, currentGene);
            if (geneVals) return geneVals;
        }}
        const colVals = getSectionColorValues(section, currentColor);
        return colVals || [];
    }}

    // Check if section passes filters
    function sectionPassesFilter(section) {{
        for (const [key, values] of Object.entries(activeFilters)) {{
            if (values.size === 0) continue;
            const sectionVal = section.metadata[key];
            if (!sectionVal || !values.has(sectionVal)) return false;
        }}
        return true;
    }}

    // Get current panel background color from CSS variable
    function getPanelBg() {{
        return getComputedStyle(document.documentElement).getPropertyValue('--panel-bg').trim();
    }}

    function getGraphColor() {{
        const color = getComputedStyle(document.documentElement).getPropertyValue('--graph-color').trim();
        return color || 'rgba(0, 0, 0, 0.12)';
    }}

    function getSectionAdjacency(section) {{
        if (section._adj) return section._adj;
        ensureSectionXY(section);
        const n = section.x.length;
        const adj = Array.from({{ length: n }}, () => []);
        const edges = getSectionEdgesPacked(section);
        if (Array.isArray(edges)) {{
            edges.forEach(edge => {{
                const i = edge[0];
                const j = edge[1];
                if (i >= 0 && j >= 0 && i < n && j < n) {{
                    adj[i].push(j);
                    adj[j].push(i);
                }}
            }});
        }} else if (edges instanceof Uint32Array) {{
            for (let e = 0; e + 1 < edges.length; e += 2) {{
                const i = edges[e];
                const j = edges[e + 1];
                if (i < n && j < n) {{
                    adj[i].push(j);
                    adj[j].push(i);
                }}
            }}
        }}
        section._adj = adj;
        return adj;
    }}

    function computeNeighborRings(section, startIdx, maxHops) {{
        const edges = getSectionEdgesPacked(section);
        if (!edges || edges.length === 0) return [];
        const adj = getSectionAdjacency(section);
        const visited = new Set([startIdx]);
        let frontier = [startIdx];
        const rings = [];

        for (let hop = 1; hop <= maxHops; hop++) {{
            const nextSet = new Set();
            frontier.forEach(node => {{
                adj[node].forEach(nb => {{
                    if (!visited.has(nb)) {{
                        visited.add(nb);
                        nextSet.add(nb);
                    }}
                }});
            }});
            if (nextSet.size === 0) break;
            const ring = Array.from(nextSet);
            rings.push(ring);
            frontier = ring;
        }}

        return rings;
    }}

    function updateHoverNeighbors(section, cellIdx) {{
        if (!neighborHoverEnabled) return false;
        const edges = section ? getSectionEdgesPacked(section) : null;
        if (!section || !edges || edges.length === 0) {{
            hoverNeighbors = null;
            return false;
        }}
        if (hoverNeighbors &&
            hoverNeighbors.sectionId === section.id &&
            hoverNeighbors.centerIdx === cellIdx) {{
            return false;
        }}
        const rings = computeNeighborRings(section, cellIdx, MAX_HOVER_HOPS);
        hoverNeighbors = {{
            sectionId: section.id,
            centerIdx: cellIdx,
            rings,
        }};
        return true;
    }}

    function drawNeighborHighlights(ctx, section, transform, adjustedSpotSize) {{
        if (!hoverNeighbors || hoverNeighbors.sectionId !== section.id) return;
        let rings = hoverNeighbors.rings || [];
        const centerIdx = hoverNeighbors.centerIdx;
        if (centerIdx === null || centerIdx === undefined) return;
        const activeIndexSet = getModalActiveCellIndexSet(section.id);
        if (activeIndexSet && !activeIndexSet.has(centerIdx)) return;

        const config = getColorConfig();
        const values = getSectionValues(section);

        if (neighborHopMode !== 'all') {{
            const hopIdx = Math.max(1, Math.min(MAX_HOVER_HOPS, parseInt(neighborHopMode, 10))) - 1;
            rings = rings[hopIdx] ? [rings[hopIdx]] : [];
        }}

        const centerPoint = transform.dataToScreen(section.x[centerIdx], section.y[centerIdx]);
        const xCenter = centerPoint.x;
        const yCenter = centerPoint.y;

        ctx.save();
        ctx.strokeStyle = 'rgba(255, 255, 255, 0.9)';
        ctx.lineWidth = Math.max(1.5, adjustedSpotSize * 0.35);
        ctx.beginPath();
        ctx.arc(xCenter, yCenter, adjustedSpotSize + 2, 0, Math.PI * 2);
        ctx.stroke();

        rings.forEach((ring, idx) => {{
            const color = HOVER_COLORS[idx] || HOVER_COLORS[HOVER_COLORS.length - 1];
            ctx.strokeStyle = color;
            ctx.lineWidth = Math.max(1, adjustedSpotSize * 0.25);
            ring.forEach(cellIdx => {{
                if (activeIndexSet && !activeIndexSet.has(cellIdx)) return;
                const val = values[cellIdx];
                if (val === null || val === undefined) return;
                if (!config.is_continuous) {{
                    const catIdx = Math.round(val);
                    const catName = config.categories[catIdx];
                    if (hiddenCategories.has(catName)) return;
                }}
                const point = transform.dataToScreen(section.x[cellIdx], section.y[cellIdx]);
                const x = point.x;
                const y = point.y;
                if (!transform.isPointVisible(x, y, adjustedSpotSize)) return;
                ctx.beginPath();
                ctx.arc(x, y, adjustedSpotSize + 1 + idx, 0, Math.PI * 2);
                ctx.stroke();
            }});
        }});

        rings.forEach((ring, idx) => {{
            const color = HOVER_COLORS[idx] || HOVER_COLORS[HOVER_COLORS.length - 1];
            ctx.strokeStyle = color;
            ctx.lineWidth = Math.max(1, adjustedSpotSize * 0.2);
            ctx.globalAlpha = 0.8;
            ctx.beginPath();
            ring.forEach(cellIdx => {{
                if (activeIndexSet && !activeIndexSet.has(cellIdx)) return;
                const val = values[cellIdx];
                if (val === null || val === undefined) return;
                if (!config.is_continuous) {{
                    const catIdx = Math.round(val);
                    const catName = config.categories[catIdx];
                    if (hiddenCategories.has(catName)) return;
                }}
                const point = transform.dataToScreen(section.x[cellIdx], section.y[cellIdx]);
                const x = point.x;
                const y = point.y;
                if (!transform.isPointVisible(x, y, adjustedSpotSize)) return;
                ctx.moveTo(xCenter, yCenter);
                ctx.lineTo(x, y);
            }});
            ctx.stroke();
            ctx.globalAlpha = 1;
        }});
        ctx.restore();
    }}

    // Check if a cell is selected
    function isCellSelected(sectionId, cellIdx) {{
        return selectedCells.has(`${{sectionId}}:${{cellIdx}}`);
    }}

    // Point-in-polygon test using ray casting algorithm
    function pointInPolygon(x, y, polygon) {{
        if (polygon.length < 3) return false;
        let inside = false;
        for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {{
            const xi = polygon[i].x, yi = polygon[i].y;
            const xj = polygon[j].x, yj = polygon[j].y;
            if (((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {{
                inside = !inside;
            }}
        }}
        return inside;
    }}

    function getSelectionSummaryColorColumn() {{
        const isValidCategorical = (col) => {{
            if (!col) return false;
            const meta = DATA.colors_meta?.[col];
            return !!(meta && !meta.is_continuous && Array.isArray(meta.categories));
        }};

        if (isValidCategorical(selectionSummaryColor)) return selectionSummaryColor;
        if (isValidCategorical(currentColor)) return currentColor;
        if (isValidCategorical(DATA.initial_color)) return DATA.initial_color;

        const fallbackCols = getCategoricalColorColumns();
        return fallbackCols.length ? fallbackCols[0] : null;
    }}

    function computeSelectionSummary(cells = selectedCells) {{
        const summary = {{
            total: cells.size,
            sections: [],
            typeColumn: null,
            types: [],
            missingTypeValues: 0,
            _cells: cells,
        }};
        if (cells.size === 0) return summary;

        const sectionCounts = new Map();
        const typeColumn = getSelectionSummaryColorColumn();
        summary.typeColumn = typeColumn;
        const typeCounts = new Map();
        const typeMeta = typeColumn ? DATA.colors_meta?.[typeColumn] : null;
        const typeCategories = Array.isArray(typeMeta?.categories) ? typeMeta.categories : null;

        cells.forEach((key) => {{
            const sep = key.lastIndexOf(':');
            if (sep <= 0) return;
            const sectionId = key.slice(0, sep);
            const idx = Number(key.slice(sep + 1));
            if (!Number.isInteger(idx) || idx < 0) return;

            sectionCounts.set(sectionId, (sectionCounts.get(sectionId) || 0) + 1);

            if (!typeColumn || !typeCategories) return;
            const section = sectionById.get(sectionId);
            if (!section) return;
            const values = getSectionColorValues(section, typeColumn);
            if (!values || idx >= values.length) {{
                summary.missingTypeValues += 1;
                return;
            }}
            const raw = values[idx];
            if (!Number.isFinite(raw)) {{
                summary.missingTypeValues += 1;
                return;
            }}
            const catIdx = Math.round(raw);
            const catName = typeCategories[catIdx] ?? 'unknown';
            typeCounts.set(catName, (typeCounts.get(catName) || 0) + 1);
        }});

        summary.sections = Array.from(sectionCounts.entries()).sort((a, b) => b[1] - a[1]);
        summary.types = Array.from(typeCounts.entries()).sort((a, b) => b[1] - a[1]);
        return summary;
    }}

    function computeSelectionMeanExpression(summary) {{
        // Prefer cluster-means lookup: covers all genes including sidecar ones.
        const clusterMeansData = DATA.cluster_gene_means;
        if (clusterMeansData && summary?.typeColumn) {{
            const colData = clusterMeansData.columns?.[summary.typeColumn];
            if (colData && summary.types.length) {{
                const genes = clusterMeansData.genes;
                const n = genes.length;
                const selMeans = new Float64Array(n);
                const totalSelected = summary.total - summary.missingTypeValues;
                if (totalSelected > 0) {{
                    for (const [cat, count] of summary.types) {{
                        const catMeans = colData.means[cat];
                        if (!catMeans) continue;
                        const weight = count / totalSelected;
                        for (let i = 0; i < n; i++) {{
                            selMeans[i] += catMeans[i] * weight;
                        }}
                    }}
                }}
                const bgMeans = colData.background;
                const result = genes.map((gene, i) => ({{
                    gene,
                    meanSel: selMeans[i],
                    meanRest: bgMeans[i],
                }}));
                result.sort((a, b) => b.meanSel - a.meanSel);
                return result;
            }}
        }}

        // Fallback: compute directly from loaded gene vectors (exact, embedded genes only).
        const loadedGenes = Object.keys(DATA.genes_meta || {{}});
        if (!loadedGenes.length || !summary.total) return [];
        // Rebuild cell set from summary sections for O(1) lookup
        const cellSetKeys = summary._cells || selectedCells;
        const result = [];
        for (const gene of loadedGenes) {{
            let selSum = 0, selCount = 0, restSum = 0, restCount = 0;
            for (const section of (DATA.sections || [])) {{
                const vals = getSectionGeneValues(section, gene);
                if (!vals) continue;
                for (let i = 0; i < vals.length; i++) {{
                    const v = vals[i] ?? 0;
                    if (cellSetKeys.has(`${{section.id}}:${{i}}`)) {{
                        selSum += v; selCount++;
                    }} else {{
                        restSum += v; restCount++;
                    }}
                }}
            }}
            const meanSel = selCount > 0 ? selSum / selCount : 0;
            const meanRest = restCount > 0 ? restSum / restCount : 0;
            result.push({{ gene, meanSel, meanRest }});
        }}
        result.sort((a, b) => b.meanSel - a.meanSel);
        return result;
    }}

    function getAnnotationCellSet(annotation) {{
        const cells = new Set();
        if (!annotation || !annotation.sectionId) return cells;
        (annotation.localCellIndices || []).forEach((idx) => {{
            if (Number.isInteger(idx) && idx >= 0) {{
                cells.add(`${{annotation.sectionId}}:${{idx}}`);
            }}
        }});
        return cells;
    }}

    function getAnnotationDisplayLabel(annotation) {{
        if (!annotation) return 'Unknown annotation';
        const label = String(annotation.label || `Annotation ${{annotation.id}}`).trim();
        const sectionLabel = annotation.sectionId ? ` · ${{annotation.sectionId}}` : '';
        return `${{label || `Annotation ${{annotation.id}}`}}${{sectionLabel}}`;
    }}

    function syncAnnotationDEState(annotations = modalAnnotations) {{
        const available = (annotations || []).filter((annotation) => Number.isFinite(Number(annotation?.id)));
        const availableIds = available.map((annotation) => Number(annotation.id));
        if (!availableIds.length) {{
            annotationDeSourceId = null;
            annotationDeReferenceId = null;
            return {{ annotations: available, source: null, reference: null }};
        }}

        if (!availableIds.includes(annotationDeSourceId)) {{
            annotationDeSourceId = availableIds[0];
        }}
        if (
            !availableIds.includes(annotationDeReferenceId) ||
            annotationDeReferenceId === annotationDeSourceId
        ) {{
            annotationDeReferenceId = availableIds.find((id) => id !== annotationDeSourceId) || null;
        }}

        const source = available.find((annotation) => Number(annotation.id) === annotationDeSourceId) || null;
        const reference = available.find((annotation) => Number(annotation.id) === annotationDeReferenceId) || null;
        return {{ annotations: available, source, reference }};
    }}

    function computeRegionAnnotationDE(annotationA, annotationB, options = {{}}) {{
        const topN = Math.max(1, Number(options.topN) || ANNOTATION_DE_TOP_N);
        if (!annotationA || !annotationB) {{
            return {{ available: false, reason: 'missing_annotations', results: [] }};
        }}
        if (Number(annotationA.id) === Number(annotationB.id)) {{
            return {{ available: false, reason: 'same_annotation', results: [] }};
        }}

        const sectionA = sectionById.get(annotationA.sectionId);
        const sectionB = sectionById.get(annotationB.sectionId);
        const indicesA = Array.from(new Set((annotationA.localCellIndices || []).filter((idx) => Number.isInteger(idx) && idx >= 0)));
        const indicesB = Array.from(new Set((annotationB.localCellIndices || []).filter((idx) => Number.isInteger(idx) && idx >= 0)));
        if (!sectionA || !sectionB || !indicesA.length || !indicesB.length) {{
            return {{
                available: false,
                reason: 'empty_region',
                nA: indicesA.length,
                nB: indicesB.length,
                results: [],
            }};
        }}

        const loadedGenes = Object.keys(DATA.genes_meta || {{}}).sort((a, b) => a.localeCompare(b));
        const totalGenes = Array.isArray(DATA.available_genes) ? DATA.available_genes.length : loadedGenes.length;
        if (!loadedGenes.length) {{
            return {{
                available: false,
                reason: 'no_loaded_genes',
                nA: indicesA.length,
                nB: indicesB.length,
                loadedGeneCount: 0,
                totalGeneCount: totalGenes,
                results: [],
            }};
        }}

        const sameSection = annotationA.sectionId === annotationB.sectionId;
        const results = [];
        const eps = 1e-6;

        loadedGenes.forEach((gene) => {{
            const valsA = getSectionGeneValues(sectionA, gene);
            const valsB = sameSection ? valsA : getSectionGeneValues(sectionB, gene);
            if (!valsA || !valsB) return;

            let sumA = 0;
            let sumSqA = 0;
            let nnzA = 0;
            for (let i = 0; i < indicesA.length; i++) {{
                const idx = indicesA[i];
                const value = Number(valsA[idx] ?? 0);
                if (!Number.isFinite(value)) continue;
                sumA += value;
                sumSqA += value * value;
                if (value > 0) nnzA += 1;
            }}

            let sumB = 0;
            let sumSqB = 0;
            let nnzB = 0;
            for (let i = 0; i < indicesB.length; i++) {{
                const idx = indicesB[i];
                const value = Number(valsB[idx] ?? 0);
                if (!Number.isFinite(value)) continue;
                sumB += value;
                sumSqB += value * value;
                if (value > 0) nnzB += 1;
            }}

            const nA = indicesA.length;
            const nB = indicesB.length;
            if (!nA || !nB) return;

            const meanA = sumA / nA;
            const meanB = sumB / nB;
            const varA = nA > 1 ? Math.max(0, (sumSqA - (sumA * sumA) / nA) / (nA - 1)) : 0;
            const varB = nB > 1 ? Math.max(0, (sumSqB - (sumB * sumB) / nB) / (nB - 1)) : 0;
            const pctA = nnzA / nA;
            const pctB = nnzB / nB;
            const denominator = Math.sqrt((varA / Math.max(1, nA)) + (varB / Math.max(1, nB)) + 1e-12);
            const score = denominator > 0 ? (meanA - meanB) / denominator : (meanA - meanB);
            const log2fc = Math.log2((meanA + eps) / (meanB + eps));
            if (!(log2fc > 0 || score > 0 || pctA > pctB)) return;

            results.push({{
                gene,
                meanA,
                meanB,
                pctA,
                pctB,
                log2fc,
                score,
            }});
        }});

        results.sort((a, b) => {{
            const scoreDiff = (b.score - a.score);
            if (Math.abs(scoreDiff) > 1e-9) return scoreDiff;
            const fcDiff = (b.log2fc - a.log2fc);
            if (Math.abs(fcDiff) > 1e-9) return fcDiff;
            const pctDiff = ((b.pctA - b.pctB) - (a.pctA - a.pctB));
            if (Math.abs(pctDiff) > 1e-9) return pctDiff > 0 ? 1 : -1;
            return a.gene.localeCompare(b.gene);
        }});

        return {{
            available: true,
            reason: null,
            nA: indicesA.length,
            nB: indicesB.length,
            loadedGeneCount: loadedGenes.length,
            totalGeneCount: totalGenes,
            results: results.slice(0, topN),
        }};
    }}

    function getAnnotationDECacheKey(annotationA, annotationB) {{
        if (!annotationA || !annotationB) return '';
        return `${{Number(annotationA.id)}}::${{Number(annotationB.id)}}`;
    }}

    function invalidateAnnotationDEState(clearCache = true) {{
        cancelAnnotationFullDERun();
        if (clearCache) annotationDeFullCache.clear();
    }}

    function cancelAnnotationFullDERun() {{
        annotationDeFullRunToken += 1;
        annotationDeFullRun = null;
    }}

    function resolveViewerRelativeUrl(path) {{
        if (!path) return '';
        if (window.__karospacePackageMode) {{
            return String(path);
        }}
        try {{
            return new URL(String(path), window.location.href).href;
        }} catch (_error) {{
            return String(path);
        }}
    }}

    async function fetchGeneAuxShardForAnalysis(shardUrl) {{
        const resolvedUrl = resolveViewerRelativeUrl(shardUrl);
        const response = await fetch(resolvedUrl, {{ credentials: 'same-origin' }});
        if (!response.ok) {{
            throw new Error(`HTTP ${{response.status}} while loading gene shard`);
        }}
        const payload = await response.json();
        if (!payload || payload.format !== 'karospace-gene-sidecar-shard-v2') {{
            throw new Error('Unsupported gene sidecar shard format');
        }}
        return payload;
    }}

    function computeRegionStatsFromSectionEntry(sectionEntry, selectedIndices, selectedIndexSet, geneMeta = null) {{
        const result = {{ sum: 0, sumSq: 0, nnz: 0 }};
        if (!sectionEntry || !selectedIndices?.length) return result;

        if (typeof sectionEntry.db64 === 'string' || typeof sectionEntry.dq16b64 === 'string' || typeof sectionEntry.dq8b64 === 'string') {{
            const dense = decodeDenseSidecarSection({{
                ...sectionEntry,
                qmin: geneMeta?.vmin,
                qmax: geneMeta?.vmax,
            }});
            if (!dense) return result;
            for (let i = 0; i < selectedIndices.length; i++) {{
                const idx = selectedIndices[i];
                const value = Number(dense[idx]);
                if (!Number.isFinite(value)) continue;
                result.sum += value;
                result.sumSq += value * value;
                if (value > 0) result.nnz += 1;
            }}
            return result;
        }}

        if (Array.isArray(sectionEntry.dense)) {{
            for (let i = 0; i < selectedIndices.length; i++) {{
                const idx = selectedIndices[i];
                const value = Number(sectionEntry.dense[idx]);
                if (!Number.isFinite(value)) continue;
                result.sum += value;
                result.sumSq += value * value;
                if (value > 0) result.nnz += 1;
            }}
            return result;
        }}

        const sparse = sectionEntry.sparse;
        if (!sparse || !selectedIndexSet || !selectedIndexSet.size) return result;

        let idxs = null;
        let vals = null;
        if (typeof sparse.ib64 === 'string' && typeof sparse.vq16b64 === 'string') {{
            idxs = base64ToUint32Array(sparse.ib64);
            vals = decodeQuantizedValues(
                base64ToUint16Array(sparse.vq16b64),
                geneMeta?.vmin,
                geneMeta?.vmax,
                65535
            );
        }} else if (typeof sparse.ib64 === 'string' && typeof sparse.vq8b64 === 'string') {{
            idxs = base64ToUint32Array(sparse.ib64);
            vals = decodeQuantizedValues(
                base64ToUint8Array(sparse.vq8b64),
                geneMeta?.vmin,
                geneMeta?.vmax,
                255
            );
        }} else if (typeof sparse.ib64 === 'string' && typeof sparse.vb64 === 'string') {{
            idxs = base64ToUint32Array(sparse.ib64);
            vals = base64ToFloat32Array(sparse.vb64);
        }} else if (
            (Array.isArray(sparse.i) || ArrayBuffer.isView(sparse.i)) &&
            (Array.isArray(sparse.v) || ArrayBuffer.isView(sparse.v))
        ) {{
            idxs = sparse.i;
            vals = sparse.v;
        }}
        if (!idxs || !vals) return result;

        const m = Math.min(idxs.length, vals.length);
        for (let i = 0; i < m; i++) {{
            const idx = Number(idxs[i]);
            if (!selectedIndexSet.has(idx)) continue;
            const value = Number(vals[i]);
            if (!Number.isFinite(value) || value === 0) continue;
            result.sum += value;
            result.sumSq += value * value;
            result.nnz += 1;
        }}
        return result;
    }}

    function computeRegionAnnotationDEFromSidecarGeneEntry(gene, geneEntry, regionA, regionB) {{
        const sectionEntryA = geneEntry?.sections?.[regionA.sectionId] || null;
        const sectionEntryB = geneEntry?.sections?.[regionB.sectionId] || null;
        const geneMeta = geneAuxManifest?.genes_meta?.[gene] || DATA.genes_meta?.[gene] || null;

        const statsA = computeRegionStatsFromSectionEntry(sectionEntryA, regionA.indices, regionA.indexSet, geneMeta);
        const statsB = computeRegionStatsFromSectionEntry(sectionEntryB, regionB.indices, regionB.indexSet, geneMeta);
        const nA = regionA.indices.length;
        const nB = regionB.indices.length;
        if (!nA || !nB) return null;

        const meanA = statsA.sum / nA;
        const meanB = statsB.sum / nB;
        const varA = nA > 1 ? Math.max(0, (statsA.sumSq - (statsA.sum * statsA.sum) / nA) / (nA - 1)) : 0;
        const varB = nB > 1 ? Math.max(0, (statsB.sumSq - (statsB.sum * statsB.sum) / nB) / (nB - 1)) : 0;
        const pctA = statsA.nnz / nA;
        const pctB = statsB.nnz / nB;
        const denominator = Math.sqrt((varA / Math.max(1, nA)) + (varB / Math.max(1, nB)) + 1e-12);
        const score = denominator > 0 ? (meanA - meanB) / denominator : (meanA - meanB);
        const log2fc = Math.log2((meanA + 1e-6) / (meanB + 1e-6));
        if (!(log2fc > 0 || score > 0 || pctA > pctB)) return null;

        return {{
            gene,
            meanA,
            meanB,
            pctA,
            pctB,
            log2fc,
            score,
        }};
    }}

    async function runFullRegionAnnotationDE(annotationA, annotationB) {{
        if (!annotationA || !annotationB || !DATA.gene_aux_url) return;
        const manifest = await loadGeneAuxManifest();
        if (!manifest) {{
            annotationDeFullRun = {{
                running: false,
                key: getAnnotationDECacheKey(annotationA, annotationB),
                error: 'Failed to load the gene sidecar manifest.',
            }};
            renderAnnotationComparison();
            return;
        }}

        const key = getAnnotationDECacheKey(annotationA, annotationB);
        const token = ++annotationDeFullRunToken;
        const shardEntries = Object.entries(manifest.shards || {{}});
        const sidecarFormat = getGeneAuxSidecarFormat(manifest);
        const totalGenes = shardEntries.reduce((sum, [, genes]) => sum + (Array.isArray(genes) ? genes.length : 0), 0);
        annotationDeFullRun = {{
            token,
            key,
            running: true,
            completedShards: 0,
            totalShards: shardEntries.length,
            completedGenes: 0,
            totalGenes,
            error: null,
        }};
        renderAnnotationComparison();

        const regionA = {{
            sectionId: annotationA.sectionId,
            indices: Array.from(new Set((annotationA.localCellIndices || []).filter((idx) => Number.isInteger(idx) && idx >= 0))),
        }};
        const regionB = {{
            sectionId: annotationB.sectionId,
            indices: Array.from(new Set((annotationB.localCellIndices || []).filter((idx) => Number.isInteger(idx) && idx >= 0))),
        }};
        regionA.indexSet = new Set(regionA.indices);
        regionB.indexSet = new Set(regionB.indices);

        const results = [];
        try {{
            for (let shardIdx = 0; shardIdx < shardEntries.length; shardIdx++) {{
                if (annotationDeFullRunToken !== token) return;
                const [shardUrl, shardGenes] = shardEntries[shardIdx];
                if (sidecarFormat === 'binary-v1') {{
                    const indexInfo = await loadBinaryShardIndex(shardUrl);
                    if (annotationDeFullRunToken !== token) return;
                    const genesInShard = Array.isArray(shardGenes) ? shardGenes : Object.keys(indexInfo?.entries || {{}});
                    for (let geneIdx = 0; geneIdx < genesInShard.length; geneIdx++) {{
                        const gene = genesInShard[geneIdx];
                        const payloadBuffer = await readBinaryGenePayload(shardUrl, gene, indexInfo);
                        if (!payloadBuffer) continue;
                        const geneEntry = parseBinaryGenePayload(
                            payloadBuffer,
                            manifest.section_order || [],
                            manifest?.genes_meta?.[gene] || null
                        );
                        const result = computeRegionAnnotationDEFromSidecarGeneEntry(gene, geneEntry, regionA, regionB);
                        if (result) results.push(result);
                    }}
                }} else {{
                    const shardPayload = await fetchGeneAuxShardForAnalysis(shardUrl);
                    if (annotationDeFullRunToken !== token) return;

                    const genesPayload = shardPayload?.genes || {{}};
                    Object.entries(genesPayload).forEach(([gene, geneEntry]) => {{
                        const result = computeRegionAnnotationDEFromSidecarGeneEntry(gene, geneEntry, regionA, regionB);
                        if (result) results.push(result);
                    }});
                }}

                if (annotationDeFullRunToken !== token) return;
                annotationDeFullRun = {{
                    ...annotationDeFullRun,
                    token,
                    key,
                    running: true,
                    completedShards: shardIdx + 1,
                    totalShards: shardEntries.length,
                    completedGenes: Math.min(totalGenes, (annotationDeFullRun?.completedGenes || 0) + (Array.isArray(shardGenes) ? shardGenes.length : 0)),
                    totalGenes,
                    error: null,
                }};
                renderAnnotationComparison();
                await new Promise((resolve) => window.setTimeout(resolve, 0));
            }}

            if (annotationDeFullRunToken !== token) return;
            results.sort((a, b) => {{
                const scoreDiff = (b.score - a.score);
                if (Math.abs(scoreDiff) > 1e-9) return scoreDiff;
                const fcDiff = (b.log2fc - a.log2fc);
                if (Math.abs(fcDiff) > 1e-9) return fcDiff;
                const pctDiff = ((b.pctA - b.pctB) - (a.pctA - a.pctB));
                if (Math.abs(pctDiff) > 1e-9) return pctDiff > 0 ? 1 : -1;
                return a.gene.localeCompare(b.gene);
            }});

            annotationDeFullCache.set(key, {{
                available: true,
                sourceId: Number(annotationA.id),
                referenceId: Number(annotationB.id),
                nA: regionA.indices.length,
                nB: regionB.indices.length,
                loadedGeneCount: totalGenes,
                totalGeneCount: totalGenes,
                results,
                completedAt: Date.now(),
                mode: 'full-sidecar',
            }});
            annotationDeFullRun = null;
            renderAnnotationComparison();
        }} catch (error) {{
            if (annotationDeFullRunToken !== token) return;
            annotationDeFullRun = {{
                token,
                key,
                running: false,
                completedShards: annotationDeFullRun?.completedShards || 0,
                totalShards: shardEntries.length,
                completedGenes: annotationDeFullRun?.completedGenes || 0,
                totalGenes,
                error: error?.message || 'Unknown error',
            }};
            renderAnnotationComparison();
        }}
    }}

    function getAnnotationDEExportState(annotationA, annotationB, quickResult = null) {{
        if (!annotationA || !annotationB) return null;
        const pairKey = getAnnotationDECacheKey(annotationA, annotationB);
        const fullCached = pairKey ? annotationDeFullCache.get(pairKey) : null;
        if (fullCached?.available) {{
            return {{ mode: 'full-sidecar', result: fullCached }};
        }}
        if (quickResult?.available) {{
            return {{ mode: 'quick-loaded-genes', result: quickResult }};
        }}
        return null;
    }}

    function buildAnnotationDEReport(annotationA, annotationB, exportState) {{
        if (!annotationA || !annotationB || !exportState?.result) return null;
        const result = exportState.result;
        const topN = Math.max(1, Number(annotationDeTopN) || ANNOTATION_DE_TOP_N);
        const exportedGenes = (result.results || []).slice(0, topN).map((entry, index) => ({{
            rank: index + 1,
            gene: entry.gene,
            log2fc_a_vs_b: entry.log2fc,
            score: entry.score,
            pct_expr_a: entry.pctA,
            pct_expr_b: entry.pctB,
            mean_a: entry.meanA,
            mean_b: entry.meanB,
        }}));
        return {{
            format: 'karospace-region-de-report-v1',
            exported_at: new Date().toISOString(),
            result_mode: exportState.mode,
            color_column: currentColor || null,
            top_genes_requested: topN,
            region_a: {{
                id: Number(annotationA.id),
                label: annotationA.label || `Annotation ${{annotationA.id}}`,
                section_id: annotationA.sectionId || null,
                n_cells: Number(result.nA || annotationA.localCellIndices?.length || 0),
                color: annotationA.color || getAnnotationColorById(annotationA.id),
                created_at: annotationA.createdAt || null,
            }},
            region_b: {{
                id: Number(annotationB.id),
                label: annotationB.label || `Annotation ${{annotationB.id}}`,
                section_id: annotationB.sectionId || null,
                n_cells: Number(result.nB || annotationB.localCellIndices?.length || 0),
                color: annotationB.color || getAnnotationColorById(annotationB.id),
                created_at: annotationB.createdAt || null,
            }},
            stats: {{
                loaded_gene_count: Number(result.loadedGeneCount || 0),
                total_gene_count: Number(result.totalGeneCount || result.loadedGeneCount || 0),
                total_hits_available: Array.isArray(result.results) ? result.results.length : 0,
                genes_exported: exportedGenes.length,
            }},
            genes: exportedGenes,
        }};
    }}

    function exportAnnotationDEReport(annotationA, annotationB, exportState) {{
        const report = buildAnnotationDEReport(annotationA, annotationB, exportState);
        if (!report) {{
            alert('No region DE result is available to export yet.');
            return;
        }}
        const sourceLabel = sanitizeFilenamePart(annotationA.label || `annotation-${{annotationA.id}}`);
        const referenceLabel = sanitizeFilenamePart(annotationB.label || `annotation-${{annotationB.id}}`);
        const filename = `karospace-region-de-${{sourceLabel}}-vs-${{referenceLabel}}-${{getScreenshotTimestamp()}}.json`;
        downloadJsonFile(report, filename);
    }}

    function escapeCsvCell(value) {{
        if (value === null || value === undefined) return '';
        const text = String(value);
        if (/[",\\n\\r]/.test(text)) {{
            return `"${{text.replace(/"/g, '""')}}"`;
        }}
        return text;
    }}

    function buildAnnotationDECsv(annotationA, annotationB, exportState) {{
        const report = buildAnnotationDEReport(annotationA, annotationB, exportState);
        if (!report) return '';
        const rows = [];
        rows.push([
            'rank',
            'gene',
            'log2fc_a_vs_b',
            'score',
            'pct_expr_a',
            'pct_expr_b',
            'mean_a',
            'mean_b',
            'result_mode',
            'color_column',
            'region_a_id',
            'region_a_label',
            'region_a_section_id',
            'region_a_n_cells',
            'region_b_id',
            'region_b_label',
            'region_b_section_id',
            'region_b_n_cells',
            'loaded_gene_count',
            'total_gene_count',
            'exported_at',
        ]);
        (report.genes || []).forEach((entry) => {{
            rows.push([
                entry.rank,
                entry.gene,
                entry.log2fc_a_vs_b,
                entry.score,
                entry.pct_expr_a,
                entry.pct_expr_b,
                entry.mean_a,
                entry.mean_b,
                report.result_mode,
                report.color_column,
                report.region_a?.id,
                report.region_a?.label,
                report.region_a?.section_id,
                report.region_a?.n_cells,
                report.region_b?.id,
                report.region_b?.label,
                report.region_b?.section_id,
                report.region_b?.n_cells,
                report.stats?.loaded_gene_count,
                report.stats?.total_gene_count,
                report.exported_at,
            ]);
        }});
        return rows
            .map((row) => row.map((value) => escapeCsvCell(value)).join(','))
            .join('\\n');
    }}

    function exportAnnotationDECsv(annotationA, annotationB, exportState) {{
        const csvText = buildAnnotationDECsv(annotationA, annotationB, exportState);
        if (!csvText) {{
            alert('No region DE result is available to export yet.');
            return;
        }}
        const sourceLabel = sanitizeFilenamePart(annotationA.label || `annotation-${{annotationA.id}}`);
        const referenceLabel = sanitizeFilenamePart(annotationB.label || `annotation-${{annotationB.id}}`);
        const filename = `karospace-region-de-${{sourceLabel}}-vs-${{referenceLabel}}-${{getScreenshotTimestamp()}}.csv`;
        downloadTextFile(csvText, filename, 'text/csv;charset=utf-8');
    }}

    function renderSelectionComparisonHtml(summaryA) {{
        const summaryB = computeSelectionSummary(selectedCellsB);
        let html = `<div class="selection-summary-compare-header">
            <span class="selection-summary-compare-label region-a">Region A (${{summaryA.total.toLocaleString()}} cells)</span>
            <span class="selection-summary-compare-label region-b">Region B (${{summaryB.total.toLocaleString()}} cells)</span>
        </div>`;

        // Cell type breakdown side-by-side
        if (summaryA.typeColumn && (summaryA.types.length || summaryB.types.length)) {{
            html += `<div class="selection-summary-title">Cell types: A vs B</div>`;
            const allTypes = new Set([
                ...summaryA.types.map(([t]) => t),
                ...summaryB.types.map(([t]) => t),
            ]);
            const mapA = new Map(summaryA.types);
            const mapB = new Map(summaryB.types);
            for (const label of allTypes) {{
                const cA = mapA.get(label) || 0;
                const cB = mapB.get(label) || 0;
                const pA = summaryA.total > 0 ? Math.round(100 * cA / summaryA.total) : 0;
                const pB = summaryB.total > 0 ? Math.round(100 * cB / summaryB.total) : 0;
                html += `<div class="selection-summary-compare-row">
                    <span class="selection-summary-compare-type">${{escapeHtml(label)}}</span>
                    <span class="selection-summary-compare-a">${{cA.toLocaleString()}} (${{pA}}%)</span>
                    <span class="selection-summary-compare-sep">|</span>
                    <span class="selection-summary-compare-b">${{cB.toLocaleString()}} (${{pB}}%)</span>
                </div>`;
            }}
        }}

        // Gene expression A vs B
        const exprA = computeSelectionMeanExpression(summaryA);
        const exprB = computeSelectionMeanExpression(summaryB);
        if (exprA.length && exprB.length) {{
            const mapB = new Map(exprB.map(e => [e.gene, e.meanSel]));
            const top = exprA.slice(0, 6);
            const vmax = Math.max(top[0].meanSel || 1, ...top.map(e => mapB.get(e.gene) || 0));
            html += '<div class="selection-summary-expr">';
            html += '<div class="selection-summary-title">Gene expression: A (pink) vs B (blue)</div>';
            top.forEach(({{gene, meanSel}}) => {{
                const meanB = mapB.get(gene) || 0;
                const pA = Math.round(100 * meanSel / vmax);
                const pB = Math.round(100 * meanB / vmax);
                const fc = meanB > 0 ? (meanSel / meanB).toFixed(1) + 'x' : '—';
                html += `<div class="selection-summary-expr-row">
                    <span class="selection-summary-expr-gene" title="${{escapeHtml(gene)}}">${{escapeHtml(gene)}}</span>
                    <div class="selection-summary-expr-bars">
                        <div class="selection-summary-expr-bar sel" style="width:${{pA}}%"></div>
                        <div class="selection-summary-expr-bar region-b" style="width:${{pB}}%"></div>
                    </div>
                    <span class="selection-summary-expr-fc">${{fc}}</span>
                </div>`;
            }});
            html += '</div>';
        }}

        html += '<button class="selection-summary-compare-btn" type="button" id="lasso-clear-b-btn">Clear Region B</button>';
        return html;
    }}

    function renderSelectionSummaryHtml(summary, options = {{}}) {{
        if (!summary || summary.total === 0) {{
            return '<div class="selection-summary-meta">Draw a region with Magic Wand to inspect selected cells.</div>';
        }}

        // Comparison mode: two regions drawn
        if (selectedCellsB.size > 0) {{
            return renderSelectionComparisonHtml(summary);
        }}

        const topSections = summary.sections.slice(0, 3).map(([sid, count]) => `${{sid}} (${{count.toLocaleString()}})`);
        const extraSections = summary.sections.length > 3 ? ` +${{(summary.sections.length - 3).toLocaleString()}} more` : '';

        let html = '';
        if (topSections.length) {{
            html += `<div class="selection-summary-meta">Sections: ${{topSections.join(' • ')}}${{extraSections}}</div>`;
        }}

        if (summary.typeColumn && summary.types.length) {{
            const visibleLimit = 8;
            const topTypes = selectionSummaryExpanded ? summary.types : summary.types.slice(0, visibleLimit);
            html += `<div class="selection-summary-title">Counts by ${{formatMetadataLabel(summary.typeColumn)}}</div>`;
            topTypes.forEach(([label, count]) => {{
                html += `<div class="selection-summary-row"><span class="selection-summary-label">${{label}}</span><span class="selection-summary-count">${{count.toLocaleString()}}</span></div>`;
            }});
            if (summary.types.length > visibleLimit) {{
                if (!selectionSummaryExpanded) {{
                    html += `<button class="selection-summary-toggle" type="button">+${{(summary.types.length - visibleLimit).toLocaleString()}} more categories</button>`;
                }} else {{
                    html += '<button class="selection-summary-toggle" type="button">Show fewer categories</button>';
                }}
            }}
            if (summary.missingTypeValues > 0) {{
                html += `<div class="selection-summary-meta">${{summary.missingTypeValues.toLocaleString()}} cells missing this annotation.</div>`;
            }}
        }} else {{
            html += '<div class="selection-summary-meta">No categorical annotation is available for type counts.</div>';
        }}

        const exprData = computeSelectionMeanExpression(summary);
        if (exprData.length) {{
            const top = exprData.slice(0, 6);
            const vmax = top[0].meanSel || 1;
            html += '<div class="selection-summary-expr">';
            html += '<div class="selection-summary-title">Mean expression — selection (pink) vs rest (grey)</div>';
            top.forEach(({{gene, meanSel, meanRest}}) => {{
                const pctSel = Math.round(100 * meanSel / vmax);
                const pctRest = Math.round(100 * meanRest / vmax);
                const fc = meanRest > 0 ? (meanSel / meanRest).toFixed(1) + 'x' : '—';
                html += `<div class="selection-summary-expr-row">
                    <span class="selection-summary-expr-gene" title="${{escapeHtml(gene)}}">${{escapeHtml(gene)}}</span>
                    <div class="selection-summary-expr-bars">
                        <div class="selection-summary-expr-bar sel" style="width:${{pctSel}}%"></div>
                        <div class="selection-summary-expr-bar rest" style="width:${{pctRest}}%"></div>
                    </div>
                    <span class="selection-summary-expr-fc">${{fc}}</span>
                </div>`;
            }});
            html += '</div>';
        }}

        const actionButtons = [];
        const allowFocusSubview = !!options.allowFocusSubview && !modalSubview;
        const canFocusSubview = allowFocusSubview && modalSection && !!buildModalSubviewStateFromSelection(modalSection);
        if (canFocusSubview) {{
            actionButtons.push('<button class="selection-summary-compare-btn" type="button" id="selection-focus-subview-btn">Open focused view</button>');
        }}

        // Compare region button (only when no B yet and lasso mode is available)
        if (!lassoModeB && !modalSubview) {{
            actionButtons.push('<button class="selection-summary-compare-btn" type="button" id="lasso-compare-region-btn">Compare region</button>');
        }} else if (lassoModeB) {{
            html += '<div class="selection-summary-meta">Draw Region B with Magic Wand…</div>';
        }}
        if (modalSubview) {{
            html += '<div class="selection-summary-meta">Focused view active. Use Back view to return to the full section.</div>';
        }}
        if (actionButtons.length) {{
            html += `<div class="selection-summary-actions">${{actionButtons.join('')}}</div>`;
        }}

        return html;
    }}

    function bindSelectionSummaryInteractions(container) {{
        if (!container) return;
        if (!container.dataset.scrollBound) {{
            // Prevent wheel/touch events from bubbling to canvas zoom handlers.
            container.addEventListener('wheel', (e) => {{
                e.stopPropagation();
            }});
            container.addEventListener('touchmove', (e) => {{
                e.stopPropagation();
            }});
            container.dataset.scrollBound = '1';
        }}
        container.querySelectorAll('.selection-summary-toggle').forEach((btn) => {{
            btn.addEventListener('click', (e) => {{
                e.preventDefault();
                e.stopPropagation();
                selectionSummaryExpanded = !selectionSummaryExpanded;
                updateSelectionInfo();
            }});
        }});
        const compareBtn = container.querySelector('#lasso-compare-region-btn');
        if (compareBtn) {{
            compareBtn.addEventListener('click', (e) => {{
                e.preventDefault();
                e.stopPropagation();
                lassoModeB = true;
                updateSelectionInfo();
            }});
        }}
        const focusSubviewBtn = container.querySelector('#selection-focus-subview-btn');
        if (focusSubviewBtn) {{
            focusSubviewBtn.addEventListener('click', (e) => {{
                e.preventDefault();
                e.stopPropagation();
                activateModalSubviewFromSelection();
            }});
        }}
        const clearBBtn = container.querySelector('#lasso-clear-b-btn');
        if (clearBBtn) {{
            clearBBtn.addEventListener('click', (e) => {{
                e.preventDefault();
                e.stopPropagation();
                selectedCellsB.clear();
                lassoModeB = false;
                updateSelectionInfo();
            }});
        }}
    }}

    // UMAP rendering
    function renderUMAP() {{
        if (!DATA.has_umap || !umapVisible) return;

        const canvas = document.getElementById('umap-canvas');
        const ctx = canvas.getContext('2d');
        const dpr = getRenderDpr();
        const container = document.getElementById('umap-canvas-container');
        const rect = container.getBoundingClientRect();
        canvas.width = rect.width * dpr;
        canvas.height = rect.height * dpr;
        ctx.scale(dpr, dpr);

        const width = rect.width, height = rect.height;
        ctx.fillStyle = getPanelBg();
        ctx.fillRect(0, 0, width, height);

        const bounds = DATA.umap_bounds;
        if (!bounds) return;

        const dataWidth = bounds.xmax - bounds.xmin;
        const dataHeight = bounds.ymax - bounds.ymin;
        const padding = 20;
        const baseScale = Math.min((width - 2*padding) / dataWidth, (height - 2*padding) / dataHeight);
        const scale = baseScale * umapZoom;

        const centerX = width / 2 + umapPanX;
        const centerY = height / 2 + umapPanY;
        const dataCenterX = (bounds.xmin + bounds.xmax) / 2;
        const dataCenterY = (bounds.ymin + bounds.ymax) / 2;

        const config = getColorConfig();
        const adjustedSpotSize = Math.max(0.25, umapSpotSize * umapZoom * 0.5);
        const activeSpotlight = getLinkedSpotlightCategory(config);
        const hasSpotlight = !!activeSpotlight;

        // Check if any categories are hidden
        const hasHidden = hiddenCategories.size > 0 && !config.is_continuous;

        // First pass: draw hidden categories as grey (if any are hidden)
        if (hasHidden) {{
            ctx.fillStyle = '#888888';
            ctx.globalAlpha = 0.2;
            DATA.sections.forEach(section => {{
                ensureSectionUMAP(section);
                if (!section.umap_x || !section.umap_y) return;
                const values = getSectionValues(section);

                for (let i = 0; i < section.umap_x.length; i++) {{
                    const val = values[i];
                    if (val === null || val === undefined) continue;

                    const catIdx = Math.round(val);
                    const catName = config.categories[catIdx];
                    if (!hiddenCategories.has(catName)) continue; // Only draw hidden cells in first pass

                    const x = centerX + (section.umap_x[i] - dataCenterX) * scale;
                    const y = centerY - (section.umap_y[i] - dataCenterY) * scale;

                    if (x < -adjustedSpotSize || x > width + adjustedSpotSize ||
                        y < -adjustedSpotSize || y > height + adjustedSpotSize) continue;

                    ctx.beginPath();
                    ctx.arc(x, y, adjustedSpotSize, 0, Math.PI * 2);
                    ctx.fill();
                }}
            }});
            ctx.globalAlpha = 1;
        }}

        // Second pass: draw visible categories with full color
        DATA.sections.forEach(section => {{
            ensureSectionUMAP(section);
            if (!section.umap_x || !section.umap_y) return;

            const values = getSectionValues(section);

            for (let i = 0; i < section.umap_x.length; i++) {{
                const val = values[i];
                if (val === null || val === undefined) continue;

                // Skip hidden categories (they were drawn in first pass)
                if (!config.is_continuous) {{
                    const catIdx = Math.round(val);
                    const catName = config.categories[catIdx];
                    if (hiddenCategories.has(catName)) continue;
                }}

                const x = centerX + (section.umap_x[i] - dataCenterX) * scale;
                const y = centerY - (section.umap_y[i] - dataCenterY) * scale;

                // Skip if outside canvas
                if (x < -adjustedSpotSize || x > width + adjustedSpotSize ||
                    y < -adjustedSpotSize || y > height + adjustedSpotSize) continue;

                let color;
                let isSpotlightCategory = false;
                if (config.is_continuous) {{
                    const t = (val - config.vmin) / (config.vmax - config.vmin);
                    color = magma(Math.max(0, Math.min(1, t)));
                }} else {{
                    const catIdx = Math.round(val);
                    const catName = config.categories[catIdx];
                    isSpotlightCategory = hasSpotlight && catName === activeSpotlight;
                    color = getCategoryColor(catIdx);
                }}

                if (hasSpotlight && !isSpotlightCategory) {{
                    ctx.fillStyle = '#bbbbbb';
                    ctx.globalAlpha = 0.12;
                }} else {{
                    ctx.fillStyle = color;
                    ctx.globalAlpha = 1;
                }}
                ctx.beginPath();
                ctx.arc(x, y, adjustedSpotSize, 0, Math.PI * 2);
                ctx.fill();

                // Draw selection highlight
                if (isCellSelected(section.id, i)) {{
                    ctx.strokeStyle = '#ffd700';
                    ctx.lineWidth = 2;
                    ctx.stroke();
                }}
            }}
        }});
        ctx.globalAlpha = 1;

        // Draw lasso path if currently drawing
        if (isDrawingLasso && lassoPath.length > 1) {{
            ctx.strokeStyle = 'rgba(255, 255, 255, 0.8)';
            ctx.lineWidth = 2;
            ctx.setLineDash([5, 5]);
            ctx.beginPath();
            ctx.moveTo(lassoPath[0].x, lassoPath[0].y);
            for (let i = 1; i < lassoPath.length; i++) {{
                ctx.lineTo(lassoPath[i].x, lassoPath[i].y);
            }}
            ctx.stroke();
            ctx.setLineDash([]);
        }}
    }}

    // Perform lasso selection on UMAP
    function performLassoSelection() {{
        if (lassoPath.length < 3) return;

        const canvas = document.getElementById('umap-canvas');
        const container = document.getElementById('umap-canvas-container');
        const rect = container.getBoundingClientRect();
        const width = rect.width, height = rect.height;

        const bounds = DATA.umap_bounds;
        if (!bounds) return;

        const dataWidth = bounds.xmax - bounds.xmin;
        const dataHeight = bounds.ymax - bounds.ymin;
        const padding = 20;
        const baseScale = Math.min((width - 2*padding) / dataWidth, (height - 2*padding) / dataHeight);
        const scale = baseScale * umapZoom;

        const centerX = width / 2 + umapPanX;
        const centerY = height / 2 + umapPanY;
        const dataCenterX = (bounds.xmin + bounds.xmax) / 2;
        const dataCenterY = (bounds.ymin + bounds.ymax) / 2;

        const config = getColorConfig();

        // Collect cells inside the lasso
        const newCells = new Set();

        // Check all cells in all sections
        DATA.sections.forEach(section => {{
            ensureSectionUMAP(section);
            if (!section.umap_x || !section.umap_y) return;

            const values = getSectionValues(section);

            for (let i = 0; i < section.umap_x.length; i++) {{
                const val = values[i];
                if (val === null || val === undefined) continue;

                // Skip hidden categories
                if (!config.is_continuous) {{
                    const catIdx = Math.round(val);
                    const catName = config.categories[catIdx];
                    if (hiddenCategories.has(catName)) continue;
                }}

                const x = centerX + (section.umap_x[i] - dataCenterX) * scale;
                const y = centerY - (section.umap_y[i] - dataCenterY) * scale;

                if (pointInPolygon(x, y, lassoPath)) {{
                    newCells.add(`${{section.id}}:${{i}}`);
                }}
            }}
        }});

        if (lassoModeB) {{
            selectedCellsB = newCells;
            lassoModeB = false;
        }} else {{
            selectedCells = newCells;
            selectedCellsB.clear();
        }}

        selectionSummaryExpanded = false;
        updateSelectionInfo();
        renderUMAP();
        renderAllSections();
        if (modalSection) renderModalSection();
    }}

    function escapeHtml(text) {{
        return String(text ?? '')
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#39;');
    }}

    function getViewerScopedStorageKey(name) {{
        const scope = `${{window.location.origin || 'file://'}}${{window.location.pathname || ''}}`;
        return `karospace:${{name}}:${{scope}}`;
    }}

    const GENE_RECENTS_STORAGE_KEY = getViewerScopedStorageKey('gene-recents');
    const GENE_PANELS_STORAGE_KEY = getViewerScopedStorageKey('gene-panels');

    function readViewerJsonStorage(key, fallback) {{
        try {{
            const raw = VIEWER_STORAGE.getItem(key);
            if (!raw) return fallback;
            return JSON.parse(raw);
        }} catch (error) {{
            console.warn('Failed to read local storage key', key, error);
            return fallback;
        }}
    }}

    function writeViewerJsonStorage(key, value) {{
        try {{
            VIEWER_STORAGE.setItem(key, JSON.stringify(value));
            return true;
        }} catch (error) {{
            console.warn('Failed to write local storage key', key, error);
            return false;
        }}
    }}

    function resolveCanonicalGeneName(token) {{
        const text = String(token || '').trim();
        if (!text) return '';
        if (AVAILABLE_GENE_SET.has(text)) return text;
        return GENE_NAME_BY_LOWER.get(text.toLowerCase()) || '';
    }}

    function loadRecentGenes() {{
        const stored = readViewerJsonStorage(GENE_RECENTS_STORAGE_KEY, []);
        if (!Array.isArray(stored)) return [];
        return stored
            .map(resolveCanonicalGeneName)
            .filter(Boolean)
            .slice(0, GENE_DISCOVERY_RECENT_LIMIT);
    }}

    function persistRecentGenes() {{
        writeViewerJsonStorage(GENE_RECENTS_STORAGE_KEY, recentGenes);
    }}

    function recordRecentGene(gene) {{
        const token = resolveCanonicalGeneName(gene);
        if (!token) return;
        recentGenes = [token, ...recentGenes.filter(item => item !== token)].slice(0, GENE_DISCOVERY_RECENT_LIMIT);
        persistRecentGenes();
    }}

    function loadSavedGenePanels() {{
        const stored = readViewerJsonStorage(GENE_PANELS_STORAGE_KEY, []);
        if (!Array.isArray(stored)) return [];
        return stored
            .map((entry) => {{
                const name = String(entry?.name || '').trim();
                const genes = Array.isArray(entry?.genes)
                    ? entry.genes.map(resolveCanonicalGeneName).filter(Boolean)
                    : [];
                if (!name) return null;
                return {{ name, genes: [...new Set(genes)] }};
            }})
            .filter(Boolean);
    }}

    function persistSavedGenePanels() {{
        writeViewerJsonStorage(GENE_PANELS_STORAGE_KEY, savedGenePanels);
    }}

    function getGenePanelSeedToken() {{
        const fromCurrent = resolveCanonicalGeneName(currentGene);
        if (fromCurrent) return fromCurrent;
        const geneInput = document.getElementById('gene-input');
        return resolveCanonicalGeneName(geneInput?.value || '');
    }}

    function upsertSavedGenePanel(panelName, gene = '') {{
        const normalizedName = String(panelName || '').trim();
        if (!normalizedName) return false;
        const geneToken = resolveCanonicalGeneName(gene);
        const existingIndex = savedGenePanels.findIndex(
            panel => panel.name.toLowerCase() === normalizedName.toLowerCase()
        );
        if (existingIndex >= 0) {{
            if (geneToken && !savedGenePanels[existingIndex].genes.includes(geneToken)) {{
                savedGenePanels[existingIndex].genes.push(geneToken);
            }}
            savedGenePanels[existingIndex].name = normalizedName;
        }} else {{
            savedGenePanels.push({{
                name: normalizedName,
                genes: geneToken ? [geneToken] : [],
            }});
        }}
        persistSavedGenePanels();
        return true;
    }}

    function deleteSavedGenePanel(panelName) {{
        const normalizedName = String(panelName || '').trim().toLowerCase();
        if (!normalizedName) return;
        savedGenePanels = savedGenePanels.filter(panel => panel.name.toLowerCase() !== normalizedName);
        persistSavedGenePanels();
    }}

    function fuzzyGeneMatchScore(candidate, query) {{
        const haystack = String(candidate || '').toLowerCase();
        const needle = String(query || '').toLowerCase();
        if (!needle) return null;
        if (haystack === needle) return {{ bucket: 0, score: 0 }};
        if (haystack.startsWith(needle)) return {{ bucket: 1, score: haystack.length - needle.length }};
        const substringIndex = haystack.indexOf(needle);
        if (substringIndex >= 0) {{
            return {{ bucket: 2, score: substringIndex * 100 + (haystack.length - needle.length) }};
        }}

        let cursor = 0;
        let penalty = 0;
        for (let i = 0; i < needle.length; i++) {{
            const next = haystack.indexOf(needle[i], cursor);
            if (next < 0) return null;
            penalty += Math.max(0, next - cursor);
            cursor = next + 1;
        }}
        return {{ bucket: 3, score: penalty + Math.max(0, haystack.length - needle.length) }};
    }}

    function getGeneSearchResults(query, limit = GENE_DISCOVERY_MAX_RESULTS) {{
        const token = String(query || '').trim();
        if (!token) return [];
        return (DATA.available_genes || [])
            .map((gene) => {{
                const match = fuzzyGeneMatchScore(gene, token);
                if (!match) return null;
                return {{
                    gene,
                    bucket: match.bucket,
                    score: match.score,
                    loaded: !!DATA.genes_meta?.[gene],
                }};
            }})
            .filter(Boolean)
            .sort((a, b) => {{
                if (a.bucket !== b.bucket) return a.bucket - b.bucket;
                if (a.score !== b.score) return a.score - b.score;
                if (a.loaded !== b.loaded) return a.loaded ? -1 : 1;
                if (a.gene.length !== b.gene.length) return a.gene.length - b.gene.length;
                return a.gene.localeCompare(b.gene);
            }})
            .slice(0, Math.max(1, Number(limit) || GENE_DISCOVERY_MAX_RESULTS))
            .map(entry => entry.gene);
    }}

    function getAvailableMarkerGeneColors() {{
        return Object.entries(DATA.marker_genes || {{}})
            .filter(([, groups]) => groups && typeof groups === 'object' && Object.keys(groups).length > 0)
            .map(([color]) => color)
            .sort((a, b) => a.localeCompare(b));
    }}

    function getMarkerGenesForColorCategory(colorCol, category) {{
        if (!colorCol || category === null || category === undefined || category === BLEND_ALL_CATEGORIES) return [];
        const byColor = (DATA.marker_genes || {{}})[colorCol];
        if (!byColor || typeof byColor !== 'object') return [];
        if (Array.isArray(byColor[category])) return byColor[category];
        const catKey = String(category);
        if (Array.isArray(byColor[catKey])) return byColor[catKey];
        const matchedKey = Object.keys(byColor).find(k => String(k).toLowerCase() === catKey.toLowerCase());
        if (matchedKey && Array.isArray(byColor[matchedKey])) return byColor[matchedKey];
        return [];
    }}

    function getAvailableComparisonColors() {{
        const withDE = new Set(Object.keys(DATA.cluster_de || {{}}));
        return getCategoricalColorColumns().sort((a, b) => {{
            const aHas = withDE.has(a), bHas = withDE.has(b);
            if (aHas !== bHas) return bHas - aHas;
            return a.localeCompare(b);
        }});
    }}

    function getClusterDECategories(colorCol) {{
        if (!colorCol) return [];
        const fromMeta = getCategoriesForColorColumn(colorCol).map(value => String(value));
        const fromData = Object.keys((DATA.cluster_de || {{}})[colorCol] || {{}});
        return Array.from(new Set([...fromMeta, ...fromData]));
    }}

    function normalizeGeneEntries(genes, limit = 0) {{
        const seen = new Set();
        const entries = [];
        (Array.isArray(genes) ? genes : []).forEach((gene) => {{
            const raw = String(gene || '').trim();
            if (!raw) return;
            const canonical = resolveCanonicalGeneName(raw);
            const key = (canonical || raw).toLowerCase();
            if (seen.has(key)) return;
            seen.add(key);
            entries.push({{ raw, canonical }});
        }});
        return limit > 0 ? entries.slice(0, limit) : entries;
    }}

    function getMarkerGeneEntries(colorCol, category, limit = 0) {{
        return normalizeGeneEntries(getMarkerGenesForColorCategory(colorCol, category), limit);
    }}

    function getMarkerOverlapEntries(colorCol, sourceCategory, referenceCategory, limit = 0) {{
        const sourceEntries = getMarkerGeneEntries(colorCol, sourceCategory, 0);
        const sourceMap = new Map();
        sourceEntries.forEach((entry) => {{
            sourceMap.set((entry.canonical || entry.raw).toLowerCase(), entry);
        }});
        const overlap = [];
        const seen = new Set();
        getMarkerGeneEntries(colorCol, referenceCategory, 0).forEach((entry) => {{
            const key = (entry.canonical || entry.raw).toLowerCase();
            if (!sourceMap.has(key) || seen.has(key)) return;
            seen.add(key);
            overlap.push(sourceMap.get(key));
        }});
        return limit > 0 ? overlap.slice(0, limit) : overlap;
    }}

    function getCategoryCountsForColor(colorCol) {{
        const key = String(colorCol || '').trim();
        if (!key) return null;
        if (comparisonCountsCache.has(key)) return comparisonCountsCache.get(key);

        const categories = getCategoriesForColorColumn(key).map(value => String(value));
        if (!categories.length) {{
            comparisonCountsCache.set(key, null);
            return null;
        }}

        const counts = Object.fromEntries(categories.map((category) => [category, 0]));
        let total = 0;
        DATA.sections.forEach((section) => {{
            const values = getSectionColorValues(section, key);
            if (!values) return;
            for (let i = 0; i < values.length; i++) {{
                const value = values[i];
                if (value === null || value === undefined || Number.isNaN(value)) continue;
                const idx = Math.round(value);
                const category = categories[idx];
                if (category === undefined) continue;
                counts[category] = (counts[category] || 0) + 1;
                total += 1;
            }}
        }});

        const result = {{ categories, counts, total }};
        comparisonCountsCache.set(key, result);
        return result;
    }}

    function getCategoryColorForValue(colorCol, category) {{
        const categories = getCategoriesForColorColumn(colorCol).map(value => String(value));
        const idx = categories.indexOf(String(category));
        return idx >= 0 ? getCategoryColor(idx) : '#999';
    }}

    function getNeighborPairSummary(colorCol, sourceCategory, targetCategory) {{
        const stats = (DATA.neighbor_stats || {{}})[colorCol];
        if (!stats || !stats.counts || !stats.categories) return null;
        const categories = (stats.categories || []).map(category => String(category));
        const sourceIdx = categories.indexOf(String(sourceCategory));
        const targetIdx = categories.indexOf(String(targetCategory));
        if (sourceIdx < 0 || targetIdx < 0) return null;

        const row = Array.isArray(stats.counts[sourceIdx]) ? stats.counts[sourceIdx] : [];
        const total = row.reduce((sum, value) => sum + (Number.isFinite(value) ? value : 0), 0);
        const count = Number(row[targetIdx] ?? 0);
        const z = (stats.zscore && stats.zscore[sourceIdx] && Number.isFinite(stats.zscore[sourceIdx][targetIdx]))
            ? Number(stats.zscore[sourceIdx][targetIdx])
            : null;
        return {{
            sourceCategory: String(sourceCategory),
            targetCategory: String(targetCategory),
            count,
            total,
            pct: total > 0 ? (count / total) * 100 : 0,
            z,
            nCells: Number(stats.n_cells?.[sourceIdx] ?? NaN),
            meanDegree: Number(stats.mean_degree?.[sourceIdx] ?? NaN),
            permN: Number(stats.perm_n ?? 0),
        }};
    }}

    function getInteractionPairSummary(colorCol, sourceCategory, targetCategory) {{
        const colorData = (DATA.interaction_markers || {{}})[colorCol] || {{}};
        const result = (colorData[String(sourceCategory)] || {{}})[String(targetCategory)];
        if (!result) return null;
        return {{
            sourceCategory: String(sourceCategory),
            targetCategory: String(targetCategory),
            result,
            genes: normalizeGeneEntries(result.genes || [], 4),
        }};
    }}

    function getPairwiseClusterDEResult(colorCol, sourceCategory, referenceCategory) {{
        return (((DATA.cluster_de || {{}})[colorCol] || {{}})[sourceCategory] || {{}})[referenceCategory] || null;
    }}

    function getGeneSuggestionGroups() {{
        const config = getColorConfig();
        if (!config || config.is_continuous) {{
            return {{
                title: 'Suggestions unavailable',
                subtitle: 'Switch to a categorical color to use marker-based suggestions.',
                groups: [],
                hiddenCount: 0,
            }};
        }}

        const markersByColor = (DATA.marker_genes || {{}})[currentColor];
        if (!markersByColor || typeof markersByColor !== 'object') {{
            const availableColors = getAvailableMarkerGeneColors();
            const availableLabel = availableColors.length
                ? ` Available for: ${{availableColors.join(', ')}}.`
                : ' No marker genes are embedded in this viewer.';
            return {{
                title: `Suggested from ${{formatMetadataLabel(currentColor)}}`,
                subtitle: `No marker genes are available for the active color.${{availableLabel}}`,
                groups: [],
                hiddenCount: 0,
            }};
        }}

        const categoryOrder = (config.categories || Object.keys(markersByColor)).map(cat => String(cat));
        const groups = categoryOrder
            .map((category) => {{
                const genes = getMarkerGenesForColorCategory(currentColor, category)
                    .map(resolveCanonicalGeneName)
                    .filter(Boolean)
                    .slice(0, GENE_DISCOVERY_SUGGESTION_GENES_PER_GROUP);
                if (!genes.length) return null;
                return {{ category, genes }};
            }})
            .filter(Boolean);

        return {{
            title: `Suggested from ${{formatMetadataLabel(currentColor)}}`,
            subtitle: groups.length ? '' : 'No marker genes are available for the active color.',
            groups: groups.slice(0, GENE_DISCOVERY_SUGGESTION_GROUP_LIMIT),
            hiddenCount: Math.max(0, groups.length - GENE_DISCOVERY_SUGGESTION_GROUP_LIMIT),
        }};
    }}

    function getModalViewTransform(rect) {{
        if (!modalSection || !rect) return null;
        return createSectionViewTransform(modalSection, {{
            width: rect.width,
            height: rect.height,
            padding: 20,
            zoom: modalZoom,
            panX: modalPanX,
            panY: modalPanY,
            isModal: true,
            boundsOverride: getModalActiveBounds(),
        }});
    }}

    function updateModalViewportInfo(transform = null) {{
        const zoomInfo = document.getElementById('zoom-info');
        if (zoomInfo) zoomInfo.textContent = `${{Math.round(modalZoom * 100)}}%`;
        const rotationInfo = document.getElementById('modal-rotation-info');
        if (rotationInfo) {{
            const angleDeg = transform ? transform.angleDeg : (modalSection ? getSectionRotationDeg(modalSection) : 0);
            rotationInfo.textContent = formatRotationLabel(angleDeg);
        }}
    }}

    function clearModalInteractionCommitTimer() {{
        if (modalInteractionCommitTimer !== null) {{
            window.clearTimeout(modalInteractionCommitTimer);
            modalInteractionCommitTimer = null;
        }}
    }}

    function scheduleModalBlendRender() {{
        if (!modalSection) return;
        if (modalBlendRenderFrame !== null) return;
        modalBlendRenderFrame = window.requestAnimationFrame(() => {{
            modalBlendRenderFrame = null;
            if (!renderModalBlendFromCache()) {{
                renderModalSection();
            }}
        }});
    }}

    function updateModalBlendLabels() {{
        if (modalBlendMixLabel) {{
            const pctA = Math.round(modalBlendMix * 100);
            const pctB = 100 - pctA;
            modalBlendMixLabel.textContent = `A left ${{pctA}}% • B right ${{pctB}}%`;
        }}
        if (modalBlendMeta) {{
            const summary = `${{getModalBlendVariableLabel(modalBlendSpec.a)}} \u2194 ${{getModalBlendVariableLabel(modalBlendSpec.b)}}`;
            const aMarkers = getModalBlendMarkerHtml('a');
            const bMarkers = getModalBlendMarkerHtml('b');
            const markerBlocks = [aMarkers, bMarkers].filter(Boolean);
            const hasMarkers = markerBlocks.length > 0;
            const toggleLabel = modalBlendMarkersCollapsed ? 'Show marker genes' : 'Hide marker genes';
            modalBlendMeta.innerHTML = [
                `<div class="modal-blend-summary">${{escapeHtml(summary)}}</div>`,
                hasMarkers
                    ? `<button type="button" class="modal-blend-markers-toggle" id="modal-blend-markers-toggle" aria-expanded="${{modalBlendMarkersCollapsed ? 'false' : 'true'}}">${{toggleLabel}}</button>`
                    : '',
                hasMarkers
                    ? `<div class="modal-blend-markers ${{modalBlendMarkersCollapsed ? 'collapsed' : ''}}">${{markerBlocks.join('')}}</div>`
                    : '',
            ].join('');
            const markersToggle = document.getElementById('modal-blend-markers-toggle');
            markersToggle?.addEventListener('click', () => {{
                modalBlendMarkersCollapsed = !modalBlendMarkersCollapsed;
                updateModalBlendLabels();
            }});
        }}
    }}

    function clearModalInteractionPreview() {{
        const canvas = document.getElementById('modal-canvas');
        if (!canvas) return;
        canvas.style.transform = '';
        canvas.style.transformOrigin = '0 0';
    }}

    function cacheModalRenderedView(rect, transform) {{
        if (!modalSection || !rect || !transform) {{
            modalRenderedView = null;
            return;
        }}
        modalRenderedView = {{
            sectionId: modalSection.id,
            width: rect.width,
            height: rect.height,
            dpr: getRenderDpr(),
            centerX: transform.centerX,
            centerY: transform.centerY,
            scale: transform.scale,
            angleDeg: transform.angleDeg,
            subviewToken: getModalSubviewRenderToken(),
        }};
    }}

    function invalidateModalRenderedView() {{
        clearModalInteractionCommitTimer();
        if (modalBlendRenderFrame !== null) {{
            window.cancelAnimationFrame(modalBlendRenderFrame);
            modalBlendRenderFrame = null;
        }}
        clearModalInteractionPreview();
        modalRenderedView = null;
        modalBlendRenderCache = null;
    }}

    function applyModalInteractionPreview() {{
        if (!modalSection || !modalRenderedView) return false;
        const canvas = document.getElementById('modal-canvas');
        const container = document.getElementById('modal-canvas-container');
        if (!canvas || !container) return false;
        const rect = container.getBoundingClientRect();
        if (rect.width <= 0 || rect.height <= 0) return false;

        const rendered = modalRenderedView;
        if (rendered.sectionId !== modalSection.id) return false;
        if ((rendered.subviewToken || '') !== getModalSubviewRenderToken()) return false;
        if (Math.abs(rendered.width - rect.width) > 0.5 || Math.abs(rendered.height - rect.height) > 0.5) return false;
        if (Math.abs(rendered.dpr - getRenderDpr()) > 0.01) return false;

        const transform = getModalViewTransform(rect);
        if (!transform) return false;
        if (Math.abs(rendered.angleDeg - transform.angleDeg) > ROTATION_EPSILON) return false;
        if (!Number.isFinite(rendered.scale) || rendered.scale <= 0) return false;

        const scaleRatio = transform.scale / rendered.scale;
        if (!Number.isFinite(scaleRatio) || scaleRatio <= 0) return false;

        const translateX = transform.centerX - scaleRatio * rendered.centerX;
        const translateY = transform.centerY - scaleRatio * rendered.centerY;
        canvas.style.transformOrigin = '0 0';
        canvas.style.transform = `matrix(${{scaleRatio}}, 0, 0, ${{scaleRatio}}, ${{translateX}}, ${{translateY}})`;
        updateModalViewportInfo(transform);
        return true;
    }}

    function scheduleModalInteractionCommit(delayMs = 120) {{
        clearModalInteractionCommitTimer();
        if (!modalSection) return;
        modalInteractionCommitTimer = window.setTimeout(() => {{
            modalInteractionCommitTimer = null;
            renderModalSection();
        }}, Math.max(0, Number(delayMs) || 0));
    }}

    function screenPointToModalData(x, y, transform) {{
        return transform.screenToData(x, y);
    }}

    function modalDataPointToScreen(x, y, transform) {{
        return transform.dataToScreen(x, y);
    }}

    function getAnnotationColorById(id) {{
        const index = Math.max(0, (Number(id) || 1) - 1);
        return MODAL_ANNOTATION_COLORS[index % MODAL_ANNOTATION_COLORS.length];
    }}

    function getSectionAnnotations(sectionId) {{
        return modalAnnotations.filter(annotation => annotation.sectionId === sectionId);
    }}

    function getModalAnnotationById(annotationId) {{
        const target = Number(annotationId);
        if (!Number.isFinite(target)) return null;
        return modalAnnotations.find(annotation => annotation.id === target) || null;
    }}

    function computeCellsInsideDataPolygon(section, polygonData) {{
        ensureSectionXY(section);
        ensureSectionObsIndices(section);

        const polygonBounds = getBoundsFromPoints(polygonData);
        const candidateIndices = polygonBounds ? querySectionSpatialIndex(section, polygonBounds) : null;
        const localIndices = [];
        const globalIndices = [];
        const nCandidates = candidateIndices ? candidateIndices.length : section.x.length;
        for (let k = 0; k < nCandidates; k++) {{
            const i = candidateIndices ? candidateIndices[k] : k;
            const x = section.x[i];
            const y = section.y[i];
            if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
            if (pointInPolygon(x, y, polygonData)) {{
                localIndices.push(i);
                const rawGlobal = section.obs_idx?.[i];
                if (Number.isFinite(rawGlobal)) {{
                    globalIndices.push(Number(rawGlobal));
                }}
            }}
        }}
        return {{ localIndices, globalIndices }};
    }}

    function createModalAnnotationFromPath() {{
        if (!modalSection || modalAnnotationPath.length < 3) return false;
        const container = document.getElementById('modal-canvas-container');
        if (!container) return false;
        const rect = container.getBoundingClientRect();
        const transform = getModalViewTransform(rect);
        if (!transform) return false;

        const polygonData = modalAnnotationPath.map((point) => screenPointToModalData(point.x, point.y, transform));
        const cells = computeCellsInsideDataPolygon(modalSection, polygonData);
        const annotationId = modalNextAnnotationId++;
        const annotation = {{
            id: annotationId,
            sectionId: modalSection.id,
            label: `Annotation ${{annotationId}}`,
            color: getAnnotationColorById(annotationId),
            createdAt: new Date().toISOString(),
            vertices: polygonData,
            localCellIndices: cells.localIndices,
            globalCellIndices: cells.globalIndices,
        }};
        modalAnnotations.push(annotation);
        renderModalAnnotationPanel();
        updateAnnotationComparisonTabVisibility();
        if (document.getElementById('color-tab-regions')?.classList.contains('active')) renderAnnotationComparison();
        return true;
    }}

    function removeModalAnnotation(annotationId) {{
        const target = Number(annotationId);
        if (!Number.isFinite(target)) return;
        invalidateAnnotationDEState();
        modalAnnotations = modalAnnotations.filter(annotation => annotation.id !== target);
        renderModalAnnotationPanel();
        updateAnnotationComparisonTabVisibility();
        if (document.getElementById('color-tab-regions')?.classList.contains('active')) renderAnnotationComparison();
        if (modalSection) renderModalSection();
    }}

    function clearModalAnnotationsForSection(sectionId) {{
        invalidateAnnotationDEState();
        modalAnnotations = modalAnnotations.filter(annotation => annotation.sectionId !== sectionId);
        renderModalAnnotationPanel();
        updateAnnotationComparisonTabVisibility();
        if (document.getElementById('color-tab-regions')?.classList.contains('active')) renderAnnotationComparison();
        if (modalSection) renderModalSection();
    }}

    function clearAllModalAnnotations() {{
        invalidateAnnotationDEState();
        modalAnnotations = [];
        renderModalAnnotationPanel();
        updateAnnotationComparisonTabVisibility();
        if (document.getElementById('color-tab-regions')?.classList.contains('active')) renderAnnotationComparison();
        if (modalSection) renderModalSection();
    }}

    function selectCellsFromAnnotation(annotationId) {{
        const annotation = getModalAnnotationById(annotationId);
        if (!annotation) return;
        selectedCells.clear();
        (annotation.localCellIndices || []).forEach((cellIdx) => {{
            if (Number.isInteger(cellIdx) && cellIdx >= 0) {{
                selectedCells.add(`${{annotation.sectionId}}:${{cellIdx}}`);
            }}
        }});
        selectionSummaryExpanded = false;
        updateSelectionInfo();
        renderAllSections();
        if (umapVisible) renderUMAP();
        if (modalSection) renderModalSection();
    }}

    function buildModalAnnotationExport() {{
        const polygons = modalAnnotations.map((annotation) => {{
            const vertices = (annotation.vertices || []).map((point) => {{
                const x = Number(point.x);
                const y = Number(point.y);
                return {{
                    x: Number.isFinite(x) ? Number(x.toFixed(6)) : null,
                    y: Number.isFinite(y) ? Number(y.toFixed(6)) : null,
                }};
            }});
            return {{
                id: annotation.id,
                label: annotation.label || `Annotation ${{annotation.id}}`,
                color: annotation.color || getAnnotationColorById(annotation.id),
                created_at: annotation.createdAt || null,
                section_id: annotation.sectionId,
                n_vertices: vertices.length,
                vertices,
                n_cells: (annotation.localCellIndices || []).length,
                cell_local_indices: (annotation.localCellIndices || []).slice(),
                cell_global_indices: (annotation.globalCellIndices || []).slice(),
            }};
        }});
        return {{
            format: 'karospace-polygon-annotations-v1',
            created_at: new Date().toISOString(),
            groupby: DATA.groupby || null,
            n_polygons: polygons.length,
            polygons,
        }};
    }}

    function exportModalAnnotations() {{
        if (!modalAnnotations.length) {{
            alert('No polygon annotations to export yet.');
            return;
        }}
        const payload = buildModalAnnotationExport();
        downloadJsonFile(payload, `karospace-annotations-${{getScreenshotTimestamp()}}.json`);
    }}

    function renderModalAnnotationPanel() {{
        const listEl = document.getElementById('modal-annotation-list');
        if (!listEl) return;
        if (!modalSection) {{
            listEl.innerHTML = '<div class="modal-annotation-empty">Open a section to annotate.</div>';
            layoutModalAnnotationPanel();
            return;
        }}

        const sectionAnnotations = getSectionAnnotations(modalSection.id);
        if (!sectionAnnotations.length) {{
            listEl.innerHTML = '<div class="modal-annotation-empty">Enable Annotate and draw a region to create polygon annotations.</div>';
            layoutModalAnnotationPanel();
            return;
        }}

        listEl.innerHTML = sectionAnnotations.map((annotation) => {{
            const count = Number(annotation.localCellIndices?.length || 0).toLocaleString();
            const label = escapeHtml(annotation.label || `Annotation ${{annotation.id}}`);
            return `
                <div class="modal-annotation-row" data-annotation-id="${{annotation.id}}">
                    <div class="modal-annotation-row-main">
                        <span class="modal-annotation-color" style="background: ${{annotation.color || getAnnotationColorById(annotation.id)}}"></span>
                        <input class="modal-annotation-label" type="text" value="${{label}}" data-annotation-label="${{annotation.id}}">
                        <span class="modal-annotation-count">${{count}} cells</span>
                    </div>
                    <div class="modal-annotation-row-actions">
                        <button type="button" data-annotation-select="${{annotation.id}}">Select cells</button>
                        <button type="button" data-annotation-delete="${{annotation.id}}">Delete</button>
                    </div>
                </div>
            `;
        }}).join('');

        listEl.querySelectorAll('[data-annotation-label]').forEach((input) => {{
            input.addEventListener('change', () => {{
                const annotation = getModalAnnotationById(input.dataset.annotationLabel);
                if (!annotation) return;
                const nextLabel = String(input.value || '').trim();
                annotation.label = nextLabel || `Annotation ${{annotation.id}}`;
                renderModalAnnotationPanel();
                if (modalSection) renderModalSection();
            }});
        }});
        listEl.querySelectorAll('[data-annotation-select]').forEach((btn) => {{
            btn.addEventListener('click', () => {{
                selectCellsFromAnnotation(btn.dataset.annotationSelect);
            }});
        }});
        listEl.querySelectorAll('[data-annotation-delete]').forEach((btn) => {{
            btn.addEventListener('click', () => {{
                removeModalAnnotation(btn.dataset.annotationDelete);
            }});
        }});
        layoutModalAnnotationPanel();
    }}

    // Perform lasso selection on the active modal section
    function performModalLassoSelection() {{
        if (!modalSection || modalLassoPath.length < 3) return;
        ensureSectionXY(modalSection);

        const container = document.getElementById('modal-canvas-container');
        if (!container) return;
        const rect = container.getBoundingClientRect();
        const transform = getModalViewTransform(rect);
        if (!transform) return;

        const config = getColorConfig();
        const values = getSectionValues(modalSection);
        const polygonData = modalLassoPath.map((point) => screenPointToModalData(point.x, point.y, transform));
        const polygonBounds = getBoundsFromPoints(polygonData);
        const candidateIndices = filterModalCandidateIndices(
            modalSection,
            polygonBounds ? querySectionSpatialIndex(modalSection, polygonBounds) : null,
        );
        selectedCells.clear();

        const nCandidates = candidateIndices ? candidateIndices.length : modalSection.x.length;
        for (let k = 0; k < nCandidates; k++) {{
            const i = candidateIndices ? candidateIndices[k] : k;
            const val = values[i];
            if (val === null || val === undefined) continue;

            if (!config.is_continuous) {{
                const catIdx = Math.round(val);
                const catName = config.categories[catIdx];
                if (hiddenCategories.has(catName)) continue;
            }}

            const x = modalSection.x[i];
            const y = modalSection.y[i];
            if (pointInPolygon(x, y, polygonData)) {{
                selectedCells.add(`${{modalSection.id}}:${{i}}`);
            }}
        }}

        selectionSummaryExpanded = false;
        updateSelectionInfo();
        renderAllSections();
        if (umapVisible) renderUMAP();
        if (modalSection) renderModalSection();
    }}

    // Update selection info display
    function updateSelectionInfo() {{
        if (modalSection && modalSubview) {{
            const nextSubview = buildModalSubviewStateFromSelection(modalSection);
            if (nextSubview) {{
                modalSubview = nextSubview;
            }} else {{
                modalSubview = null;
                modalZoom = 1;
                modalPanX = 0;
                modalPanY = 0;
                invalidateModalRenderedView();
            }}
            updateModalHeader();
        }}
        const countText = selectedCells.size === 0
            ? 'No cells selected'
            : `${{selectedCells.size.toLocaleString()}} cells selected`;

        const umapInfo = document.getElementById('umap-selection-info');
        if (umapInfo) umapInfo.textContent = countText;
        const modalInfo = document.getElementById('modal-selection-info');
        if (modalInfo) modalInfo.textContent = countText;

        const summary = computeSelectionSummary();
        const umapSummary = document.getElementById('umap-selection-summary');
        if (umapSummary) {{
            umapSummary.classList.toggle('expanded', selectionSummaryExpanded);
            umapSummary.innerHTML = renderSelectionSummaryHtml(summary);
            bindSelectionSummaryInteractions(umapSummary);
        }}
        const modalSummary = document.getElementById('modal-selection-summary');
        if (modalSummary) {{
            modalSummary.classList.toggle('expanded', selectionSummaryExpanded);
            modalSummary.innerHTML = renderSelectionSummaryHtml(summary, {{ allowFocusSubview: true }});
            bindSelectionSummaryInteractions(modalSummary);
        }}
        updateModalToolbarState();
    }}

    // Clear selection
    function clearSelection() {{
        selectedCells.clear();
        selectedCellsB.clear();
        lassoModeB = false;
        selectionSummaryExpanded = false;
        updateSelectionInfo();
        renderUMAP();
        renderAllSections();
        if (modalSection) renderModalSection();
    }}

    function clampUMAPPanelSize(size) {{
        const viewportLimit = Math.floor(Math.min(window.innerWidth, window.innerHeight) * 0.8);
        const maxSize = Math.max(UMAP_PANEL_MIN_SIZE, Math.min(UMAP_PANEL_MAX_SIZE, viewportLimit));
        if (!Number.isFinite(size)) return Math.min(320, maxSize);
        return Math.max(UMAP_PANEL_MIN_SIZE, Math.min(maxSize, Math.round(size)));
    }}

    function getUMAPDockLabel(dock) {{
        if (dock === 'top-right') return 'TR';
        if (dock === 'top-left') return 'TL';
        if (dock === 'bottom-right') return 'BR';
        if (dock === 'bottom-left') return 'BL';
        return 'TR';
    }}

    function updateUMAPGridReserve() {{
        const contentColumn = document.getElementById('content-column');
        if (!contentColumn) return;
        let reserveLeft = 0;
        let reserveRight = 0;
        if (umapVisible) {{
            const reserve = Math.max(0, Math.round(umapPanelSize) + UMAP_GRID_RESERVE_MARGIN);
            if (umapPanelDock.endsWith('left')) reserveLeft = reserve;
            else reserveRight = reserve;
        }}
        contentColumn.style.setProperty('--umap-reserve-left', `${{reserveLeft}}px`);
        contentColumn.style.setProperty('--umap-reserve-right', `${{reserveRight}}px`);
    }}

    function loadUMAPPanelState() {{
        try {{
            const saved = VIEWER_STORAGE.getItem(UMAP_PANEL_STORAGE_KEY);
            if (!saved) return;
            const parsed = JSON.parse(saved);
            if (parsed && typeof parsed === 'object') {{
                if (UMAP_PANEL_DOCKS.includes(parsed.dock)) umapPanelDock = parsed.dock;
                if (Number.isFinite(parsed.size)) umapPanelSize = clampUMAPPanelSize(parsed.size);
            }}
        }} catch (err) {{
            // Ignore malformed localStorage values.
        }}
    }}

    function saveUMAPPanelState() {{
        try {{
            VIEWER_STORAGE.setItem(
                UMAP_PANEL_STORAGE_KEY,
                JSON.stringify({{ dock: umapPanelDock, size: umapPanelSize }})
            );
        }} catch (err) {{
            // Ignore storage failures (private mode, quota, etc.).
        }}
    }}

    function applyUMAPPanelState() {{
        const panel = document.getElementById('umap-panel');
        if (!panel) return;
        umapPanelSize = clampUMAPPanelSize(umapPanelSize);
        panel.classList.remove('dock-top-right', 'dock-top-left', 'dock-bottom-right', 'dock-bottom-left');
        panel.classList.add(`dock-${{umapPanelDock}}`);
        panel.style.width = `${{umapPanelSize}}px`;
        updateUMAPGridReserve();
        const dockBtn = document.getElementById('umap-dock-btn');
        if (dockBtn) dockBtn.textContent = getUMAPDockLabel(umapPanelDock);
        saveUMAPPanelState();
    }}

    function adjustUMAPPanelSize(delta) {{
        umapPanelSize = clampUMAPPanelSize(umapPanelSize + delta);
        applyUMAPPanelState();
        if (umapVisible) renderUMAP();
    }}

    // Toggle UMAP panel
    function toggleUMAP() {{
        umapVisible = !umapVisible;
        const panel = document.getElementById('umap-panel');
        const btn = document.getElementById('umap-toggle');
        panel.classList.toggle('visible', umapVisible);
        btn.classList.toggle('active', umapVisible);
        updateUMAPGridReserve();
        // Re-render after layout change to fix grid sizing
        requestAnimationFrame(() => {{
            renderAllSections();
            if (umapVisible) renderUMAP();
        }});
    }}

    // Initialize UMAP panel
    function initUMAP() {{
        if (!DATA.has_umap) return;

        // Show UMAP toggle button
        document.getElementById('umap-toggle').style.display = 'inline-block';
        document.getElementById('umap-toggle').addEventListener('click', toggleUMAP);
        loadUMAPPanelState();
        applyUMAPPanelState();

        const dockBtn = document.getElementById('umap-dock-btn');
        if (dockBtn) {{
            dockBtn.addEventListener('click', () => {{
                const idx = UMAP_PANEL_DOCKS.indexOf(umapPanelDock);
                umapPanelDock = UMAP_PANEL_DOCKS[(idx + 1 + UMAP_PANEL_DOCKS.length) % UMAP_PANEL_DOCKS.length];
                applyUMAPPanelState();
                if (umapVisible) renderUMAP();
            }});
        }}
        document.getElementById('umap-panel-smaller')?.addEventListener('click', () => {{
            adjustUMAPPanelSize(-UMAP_PANEL_SIZE_STEP);
        }});
        document.getElementById('umap-panel-larger')?.addEventListener('click', () => {{
            adjustUMAPPanelSize(UMAP_PANEL_SIZE_STEP);
        }});

        // Magic wand button
        document.getElementById('magic-wand-btn').addEventListener('click', () => {{
            magicWandActive = !magicWandActive;
            const btn = document.getElementById('magic-wand-btn');
            btn.classList.toggle('active', magicWandActive);
            const canvas = document.getElementById('umap-canvas');
            canvas.style.cursor = magicWandActive ? 'crosshair' : 'grab';
        }});

        // Clear selection button
        document.getElementById('clear-selection-btn').addEventListener('click', clearSelection);

        // UMAP spot size slider
        document.getElementById('umap-spot-size').addEventListener('input', (e) => {{
            umapSpotSize = parseFloat(e.target.value);
            document.getElementById('umap-spot-size-label').textContent = umapSpotSize.toFixed(2);
            renderUMAP();
        }});

        // UMAP canvas events
        const canvas = document.getElementById('umap-canvas');
        const container = document.getElementById('umap-canvas-container');

        canvas.addEventListener('mousedown', (e) => {{
            const rect = container.getBoundingClientRect();
            const x = e.clientX - rect.left;
            const y = e.clientY - rect.top;

            if (magicWandActive) {{
                // Start lasso drawing
                isDrawingLasso = true;
                lassoPath = [{{ x, y }}];
            }} else {{
                // Start panning
                isUmapDragging = true;
                umapDragStartX = e.clientX;
                umapDragStartY = e.clientY;
                umapLastPanX = umapPanX;
                umapLastPanY = umapPanY;
                canvas.style.cursor = 'grabbing';
            }}
        }});

        canvas.addEventListener('mousemove', (e) => {{
            const rect = container.getBoundingClientRect();
            const x = e.clientX - rect.left;
            const y = e.clientY - rect.top;

            if (isDrawingLasso) {{
                lassoPath.push({{ x, y }});
                renderUMAP();
            }} else if (isUmapDragging) {{
                umapPanX = umapLastPanX + (e.clientX - umapDragStartX);
                umapPanY = umapLastPanY + (e.clientY - umapDragStartY);
                renderUMAP();
            }}
        }});

        document.addEventListener('mouseup', (e) => {{
            if (isDrawingLasso) {{
                isDrawingLasso = false;
                performLassoSelection();
                lassoPath = [];
            }}
            if (isUmapDragging) {{
                isUmapDragging = false;
                canvas.style.cursor = magicWandActive ? 'crosshair' : 'grab';
            }}
        }});

        // Zoom with mouse wheel
        container.addEventListener('wheel', (e) => {{
            e.preventDefault();
            const rect = container.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;
            const bounds = DATA.umap_bounds;
            if (!bounds) return;

            const dataWidth = bounds.xmax - bounds.xmin;
            const dataHeight = bounds.ymax - bounds.ymin;
            const padding = 20;
            const baseScale = Math.min((rect.width - 2 * padding) / dataWidth, (rect.height - 2 * padding) / dataHeight);
            const oldScale = baseScale * umapZoom;
            const nextZoom = Math.max(0.1, Math.min(20, umapZoom * (e.deltaY > 0 ? 0.9 : 1.1)));
            const newScale = baseScale * nextZoom;

            const dataCenterX = (bounds.xmin + bounds.xmax) / 2;
            const dataCenterY = (bounds.ymin + bounds.ymax) / 2;
            const centerX = rect.width / 2 + umapPanX;
            const centerY = rect.height / 2 + umapPanY;

            const dataX = dataCenterX + (mouseX - centerX) / oldScale;
            const dataY = dataCenterY - (mouseY - centerY) / oldScale;

            const newCenterX = mouseX - (dataX - dataCenterX) * newScale;
            const newCenterY = mouseY + (dataY - dataCenterY) * newScale;
            umapPanX = newCenterX - rect.width / 2;
            umapPanY = newCenterY - rect.height / 2;
            umapZoom = nextZoom;

            renderUMAP();
        }});

        canvas.style.cursor = 'grab';
    }}

    // Rendering
    let renderAllJobId = 0;

    function hideLoader() {{
        const loader = document.getElementById('loading-overlay');
        if (loader) loader.style.display = 'none';
        if (typeof window.__karospaceRemoveLoadingWarning === 'function') {{
            window.__karospaceRemoveLoadingWarning();
        }}
    }}

    function renderSection(section, canvas) {{
        ensureSectionXY(section);
        const ctx = canvas.getContext('2d');
        const dpr = getRenderDpr();
        const rect = canvas.getBoundingClientRect();
        canvas.width = rect.width * dpr;
        canvas.height = rect.height * dpr;
        ctx.scale(dpr, dpr);

        const width = rect.width, height = rect.height, padding = 8;
        ctx.fillStyle = getPanelBg();
        ctx.fillRect(0, 0, width, height);

        if (section.x.length === 0) return;
        const transform = createSectionViewTransform(section, {{
            width,
            height,
            padding,
            isModal: false,
        }});
        if (!transform) return;

        const config = getColorConfig();
        const values = getSectionValues(section);

        const edges = getSectionEdgesPacked(section);
        if (showGraph && edges && edges.length) {{
            const graphColor = getGraphColor();
            ctx.strokeStyle = graphColor;
            ctx.lineWidth = Math.max(0.3, spotSize * 0.15);
            ctx.beginPath();
            if (Array.isArray(edges)) {{
                for (let e = 0; e < edges.length; e++) {{
                    const edge = edges[e];
                    const i = edge[0];
                    const j = edge[1];
                    if (i >= section.x.length || j >= section.x.length) continue;
                    const p1 = transform.dataToScreen(section.x[i], section.y[i]);
                    const p2 = transform.dataToScreen(section.x[j], section.y[j]);
                    const x1 = p1.x;
                    const y1 = p1.y;
                    const x2 = p2.x;
                    const y2 = p2.y;
                    ctx.moveTo(x1, y1);
                    ctx.lineTo(x2, y2);
                }}
            }} else if (edges instanceof Uint32Array) {{
                for (let e = 0; e + 1 < edges.length; e += 2) {{
                    const i = edges[e];
                    const j = edges[e + 1];
                    if (i >= section.x.length || j >= section.x.length) continue;
                    const p1 = transform.dataToScreen(section.x[i], section.y[i]);
                    const p2 = transform.dataToScreen(section.x[j], section.y[j]);
                    const x1 = p1.x;
                    const y1 = p1.y;
                    const x2 = p2.x;
                    const y2 = p2.y;
                    ctx.moveTo(x1, y1);
                    ctx.lineTo(x2, y2);
                }}
            }}
            ctx.stroke();
        }}

        // First pass: draw grey background for hidden categories (if any are hidden)
        const hasHidden = hiddenCategories.size > 0 && !config.is_continuous;
        if (hasHidden) {{
            ctx.fillStyle = '#cccccc';
            ctx.globalAlpha = 0.2;
            for (let i = 0; i < section.x.length; i++) {{
                const val = values[i];
                if (val === null || val === undefined) continue;
                const catIdx = Math.round(val);
                const catName = config.categories[catIdx];
                if (!hiddenCategories.has(catName)) continue;  // Only draw hidden ones

                const point = transform.dataToScreen(section.x[i], section.y[i]);
                const x = point.x;
                const y = point.y;
                ctx.beginPath();
                ctx.arc(x, y, spotSize, 0, Math.PI * 2);
                ctx.fill();
            }}
            ctx.globalAlpha = 1;
        }}

        // Second pass: draw visible categories on top (with optional selected-category focus)
        const activeSpotlight = getLinkedSpotlightCategory(config);
        const focusCategory = activeSpotlight || modalSelectedCategory;
        const hasTypeFocus = !config.is_continuous && focusCategory;
        for (let i = 0; i < section.x.length; i++) {{
            const val = values[i];
            if (val === null || val === undefined) continue;

            let color;
            let isSelectedCat = false;
            if (config.is_continuous) {{
                const t = (val - config.vmin) / (config.vmax - config.vmin);
                color = magma(Math.max(0, Math.min(1, t)));
            }} else {{
                const catIdx = Math.round(val);
                const catName = config.categories[catIdx];
                if (hiddenCategories.has(catName)) continue;  // Skip hidden, already drawn as grey
                isSelectedCat = focusCategory && catName === focusCategory;
                color = getCategoryColor(catIdx);
            }}

            const point = transform.dataToScreen(section.x[i], section.y[i]);
            const x = point.x;
            const y = point.y;
            if (hasTypeFocus && !isSelectedCat) {{
                ctx.fillStyle = '#bbbbbb';
                ctx.globalAlpha = 0.15;
            }} else {{
                ctx.fillStyle = color;
                ctx.globalAlpha = 1;
            }}
            ctx.beginPath();
            ctx.arc(x, y, spotSize, 0, Math.PI * 2);
            ctx.fill();
        }}
        ctx.globalAlpha = 1;

        // Third pass: draw selection highlights
        if (selectedCells.size > 0) {{
            ctx.strokeStyle = '#ffd700';
            ctx.lineWidth = 2;
            for (let i = 0; i < section.x.length; i++) {{
                if (!isCellSelected(section.id, i)) continue;

                const val = values[i];
                if (val === null || val === undefined) continue;

                // Skip hidden categories
                if (!config.is_continuous) {{
                    const catIdx = Math.round(val);
                    const catName = config.categories[catIdx];
                    if (hiddenCategories.has(catName)) continue;
                }}

                const point = transform.dataToScreen(section.x[i], section.y[i]);
                const x = point.x;
                const y = point.y;
                ctx.beginPath();
                ctx.arc(x, y, spotSize + 1, 0, Math.PI * 2);
                ctx.stroke();
            }}
        }}
    }}

    function renderAllSections() {{
        renderAllJobId += 1;
        const jobId = renderAllJobId;
        const panels = document.querySelectorAll('.section-panel');
        const grid = document.getElementById('grid');
        const gridRect = grid ? grid.getBoundingClientRect() : null;
        const isInView = (panel) => {{
            if (!gridRect) return true;
            const r = panel.getBoundingClientRect();
            const margin = 200;
            return (
                r.bottom >= gridRect.top - margin &&
                r.top <= gridRect.bottom + margin &&
                r.right >= gridRect.left - margin &&
                r.left <= gridRect.right + margin
            );
        }};
        let visibleCount = 0;
        let totalCells = 0;
        const drawList = [];

        DATA.sections.forEach((section, idx) => {{
            const panel = panels[idx];
            if (!panel) return;

            const passes = sectionPassesFilter(section);
            panel.classList.toggle('filtered-out', !passes);

            if (passes) {{
                visibleCount++;
                totalCells += section.n_cells;
                const canvas = panel.querySelector('canvas');
                if (canvas && isInView(panel)) drawList.push({{ section, canvas }});
            }}
        }});

        // Update stats
        const colorLabel = currentGene || currentColor;
        document.getElementById('stats-text').textContent =
            `${{visibleCount}}/${{DATA.n_sections}} sections | ${{totalCells.toLocaleString()}} cells | ${{colorLabel}}`;

        // When filters collapse to a single visible sample, avoid full-width stretching.
        if (grid) {{
            const useSingleVisibleFilterLayout = visibleCount === 1 && (DATA.n_sections || 0) > 1;
            grid.classList.toggle('single-visible-filter-layout', useSingleVisibleFilterLayout);
        }}

        // Show no results message
        let noResults = document.querySelector('.no-results');
        if (visibleCount === 0) {{
            if (!noResults) {{
                noResults = document.createElement('div');
                noResults.className = 'no-results';
                noResults.textContent = 'No sections match the current filters';
                document.getElementById('grid').appendChild(noResults);
            }}
        }} else if (noResults) {{
            noResults.remove();
        }}

        // Draw visible sections incrementally to keep the UI responsive.
        let i = 0;
        const step = () => {{
            if (jobId !== renderAllJobId) return;
            const start = performance.now();
            while (i < drawList.length && (performance.now() - start) < 10) {{
                const item = drawList[i++];
                try {{
                    renderSection(item.section, item.canvas);
                }} catch (e) {{
                    console.error('renderSection failed', e);
                }}
            }}
            if (i < drawList.length) {{
                requestAnimationFrame(step);
            }}
        }};
        requestAnimationFrame(step);
    }}

    // Tooltip functionality
    const tooltip = document.getElementById('cell-tooltip');
    let tooltipTimeout = null;

    function showTooltip(x, y, content) {{
        tooltip.innerHTML = content;
        tooltip.classList.add('visible');
        // Position tooltip, keeping it on screen
        const rect = tooltip.getBoundingClientRect();
        const tooltipX = Math.min(x + 15, window.innerWidth - rect.width - 10);
        const tooltipY = Math.min(y + 15, window.innerHeight - rect.height - 10);
        tooltip.style.left = tooltipX + 'px';
        tooltip.style.top = tooltipY + 'px';
    }}

    function hideTooltip() {{
        tooltip.classList.remove('visible');
    }}

    function initHelpTooltips() {{
        const helpTooltip = document.getElementById('help-tooltip');
        if (!helpTooltip) return;
        let activeTarget = null;

        const positionHelpTooltip = (target) => {{
            if (!target) return;
            helpTooltip.style.left = '-9999px';
            helpTooltip.style.top = '-9999px';
            const tipRect = helpTooltip.getBoundingClientRect();
            const rect = target.getBoundingClientRect();
            let left = rect.left + (rect.width - tipRect.width) / 2;
            left = Math.max(8, Math.min(window.innerWidth - tipRect.width - 8, left));
            let top = rect.bottom + 8;
            if (top + tipRect.height > window.innerHeight - 8) {{
                top = rect.top - tipRect.height - 8;
            }}
            top = Math.max(8, top);
            helpTooltip.style.left = `${{left}}px`;
            helpTooltip.style.top = `${{top}}px`;
        }};

        const hideHelpTooltip = () => {{
            helpTooltip.classList.remove('visible');
            helpTooltip.setAttribute('aria-hidden', 'true');
            activeTarget = null;
        }};

        const showHelpTooltip = (target) => {{
            const text = (target?.dataset?.help || '').trim();
            if (!text) {{
                hideHelpTooltip();
                return;
            }}
            activeTarget = target;
            helpTooltip.textContent = text;
            helpTooltip.classList.add('visible');
            helpTooltip.setAttribute('aria-hidden', 'false');
            positionHelpTooltip(target);
        }};

        document.querySelectorAll('[data-help]').forEach((el) => {{
            if (el.hasAttribute('title')) {{
                el.setAttribute('data-native-title', el.getAttribute('title') || '');
                el.removeAttribute('title');
            }}
            el.addEventListener('mouseenter', () => showHelpTooltip(el));
            el.addEventListener('mouseleave', () => hideHelpTooltip());
            el.addEventListener('focus', () => showHelpTooltip(el));
            el.addEventListener('blur', () => hideHelpTooltip());
        }});

        document.addEventListener('pointerdown', (event) => {{
            const target = event.target;
            if (target && target.closest && target.closest('[data-help]')) return;
            hideHelpTooltip();
        }});
        window.addEventListener('scroll', () => {{
            if (!activeTarget || !helpTooltip.classList.contains('visible')) return;
            if (!document.body.contains(activeTarget)) {{
                hideHelpTooltip();
                return;
            }}
            positionHelpTooltip(activeTarget);
        }}, true);
        window.addEventListener('resize', () => {{
            if (!activeTarget || !helpTooltip.classList.contains('visible')) return;
            positionHelpTooltip(activeTarget);
        }});
    }}

    function findNearestCell(section, mouseX, mouseY, canvasRect, transform, options = {{}}) {{
        ensureSectionXY(section);
        const config = getColorConfig();
        const values = getSectionValues(section);
        const ignoreMissing = !!options.ignoreMissing;
        const ignoreHidden = !!options.ignoreHidden;
        const activeIndexSet = transform.isModal ? getModalActiveCellIndexSet(section?.id) : null;
        const searchRadius = transform.isModal ? modalSpotSize * modalZoom * 2 : spotSize * 3;
        const candidateBounds = getDataBoundsFromScreenRect(
            transform,
            mouseX - searchRadius,
            mouseY - searchRadius,
            mouseX + searchRadius,
            mouseY + searchRadius,
        );
        const candidateIndices = (transform.isModal && candidateBounds)
            ? filterModalCandidateIndices(section, querySectionSpatialIndex(section, candidateBounds))
            : null;

        let nearestIdx = -1;
        let nearestDist = Infinity;

        const nCandidates = candidateIndices ? candidateIndices.length : section.x.length;
        for (let k = 0; k < nCandidates; k++) {{
            const i = candidateIndices ? candidateIndices[k] : k;
            if (activeIndexSet && !activeIndexSet.has(i)) continue;
            const val = values[i];
            if (!ignoreMissing && (val === null || val === undefined)) continue;

            // Skip hidden categories
            if (!ignoreHidden && !config.is_continuous && Number.isFinite(val)) {{
                const catIdx = Math.round(val);
                const catName = config.categories[catIdx];
                if (hiddenCategories.has(catName)) continue;
            }}

            const point = transform.dataToScreen(section.x[i], section.y[i]);
            const screenX = point.x;
            const screenY = point.y;

            const dist = Math.sqrt((mouseX - screenX) ** 2 + (mouseY - screenY) ** 2);
            if (dist < nearestDist && dist < searchRadius) {{
                nearestDist = dist;
                nearestIdx = i;
            }}
        }}

        return nearestIdx;
    }}

    function getCellTooltipContent(section, cellIdx) {{
        const config = getColorConfig();
        const values = getSectionValues(section);
        const val = values[cellIdx];
        const colorLabel = currentGene || currentColor;

        if (config.is_continuous) {{
            const t = (val - config.vmin) / (config.vmax - config.vmin);
            const color = magma(Math.max(0, Math.min(1, t)));
            return `<span class="cell-tooltip-color" style="background: ${{color}}"></span>
                    <span class="cell-tooltip-label">${{colorLabel}}:</span>
                    <span class="cell-tooltip-value">${{val.toFixed(3)}}</span>`;
        }} else {{
            const catIdx = Math.round(val);
            const catName = config.categories[catIdx];
            const color = getCategoryColor(catIdx);
            return `<span class="cell-tooltip-color" style="background: ${{color}}"></span>
                    <span class="cell-tooltip-label">${{catName}}</span>`;
        }}
    }}

    function getModalSplitTooltipContent(section, cellIdx, transform) {{
        const blendRuntimes = getModalBlendRuntimes(section);
        if (!blendRuntimes) return null;

        const cellX = transform.dataToScreen(section.x[cellIdx], section.y[cellIdx]).x;
        const splitX = transform.width * modalBlendMix;
        const side = cellX <= splitX ? 'a' : 'b';
        const sideLabel = side === 'a' ? 'A (left)' : 'B (right)';
        const spec = modalBlendSpec[side];
        const runtime = blendRuntimes[side];
        if (!runtime || !spec) return null;

        const raw = runtime.values?.[cellIdx];
        const color = rgbToCss(getModalBlendCellRgb(runtime, cellIdx));
        const baseLabel = `${{sideLabel}} • ${{getModalBlendVariableLabel(spec)}}`;

        if (runtime.kind === 'gene') {{
            const valueText = Number.isFinite(raw) ? raw.toFixed(3) : 'n/a';
            return `<span class="cell-tooltip-color" style="background: ${{color}}"></span>
                    <span class="cell-tooltip-label">${{baseLabel}}:</span>
                    <span class="cell-tooltip-value">${{valueText}}</span>`;
        }}

        const categories = getCategoriesForColorColumn(spec.color);
        const catIdx = Number.isFinite(raw) ? Math.round(raw) : -1;
        const catName = (catIdx >= 0 && catIdx < categories.length) ? categories[catIdx] : 'unknown';
        if (runtime.kind === 'cell-all') {{
            return `<span class="cell-tooltip-color" style="background: ${{color}}"></span>
                    <span class="cell-tooltip-label">${{baseLabel}}:</span>
                    <span class="cell-tooltip-value">${{catName}}</span>`;
        }}

        const isMatch = Number.isFinite(raw) && Math.round(raw) === runtime.catIdx;
        const valueText = Number.isFinite(raw) ? `${{catName}} (${{isMatch ? 'match' : 'other'}})` : 'n/a';
        return `<span class="cell-tooltip-color" style="background: ${{color}}"></span>
                <span class="cell-tooltip-label">${{baseLabel}}:</span>
                <span class="cell-tooltip-value">${{valueText}}</span>`;
    }}

    function drawModalAnnotations(ctx, transform) {{
        if (!modalSection) return;
        const sectionAnnotations = getSectionAnnotations(modalSection.id);
        if (!sectionAnnotations.length) return;

        sectionAnnotations.forEach((annotation) => {{
            const vertices = annotation.vertices || [];
            if (vertices.length < 2) return;
            ctx.save();
            ctx.lineWidth = 2;
            ctx.strokeStyle = annotation.color || getAnnotationColorById(annotation.id);
            ctx.fillStyle = annotation.color || getAnnotationColorById(annotation.id);
            ctx.globalAlpha = 0.18;
            const first = modalDataPointToScreen(vertices[0].x, vertices[0].y, transform);
            ctx.beginPath();
            ctx.moveTo(first.x, first.y);
            for (let i = 1; i < vertices.length; i++) {{
                const point = modalDataPointToScreen(vertices[i].x, vertices[i].y, transform);
                ctx.lineTo(point.x, point.y);
            }}
            ctx.closePath();
            ctx.fill();
            ctx.globalAlpha = 0.95;
            ctx.stroke();

            let cx = 0;
            let cy = 0;
            for (let i = 0; i < vertices.length; i++) {{
                const point = modalDataPointToScreen(vertices[i].x, vertices[i].y, transform);
                cx += point.x;
                cy += point.y;
            }}
            cx /= vertices.length;
            cy /= vertices.length;
            const count = Number(annotation.localCellIndices?.length || 0).toLocaleString();
            ctx.fillStyle = annotation.color || getAnnotationColorById(annotation.id);
            ctx.globalAlpha = 1;
            ctx.font = '11px sans-serif';
            ctx.fillText(`${{annotation.label || `Annotation ${{annotation.id}}`}} (${{count}})`, cx + 6, cy - 6);
            ctx.restore();
        }});
    }}

    function createModalBlendLayerCanvas(width, height, dpr) {{
        const canvas = document.createElement('canvas');
        canvas.width = Math.max(1, Math.round(width * dpr));
        canvas.height = Math.max(1, Math.round(height * dpr));
        const ctx = canvas.getContext('2d');
        ctx.scale(dpr, dpr);
        return {{ canvas, ctx }};
    }}

    function getModalBlendSpecCacheKey(spec, runtime) {{
        const parts = [
            spec?.kind || '',
            spec?.color || '',
            spec?.category || '',
            spec?.gene || '',
            runtime?.kind || '',
        ];
        if (runtime?.kind === 'gene') {{
            parts.push(
                Number(runtime.vmin || 0).toFixed(6),
                Number(runtime.vmax || 0).toFixed(6),
            );
        }} else if (Number.isFinite(runtime?.catIdx)) {{
            parts.push(String(runtime.catIdx));
        }}
        return parts.join('|');
    }}

    function getModalBlendRenderCacheKey(section, transform, width, height, dpr, adjustedSpotSize, runtimes) {{
        return [
            section?.id || '',
            width,
            height,
            dpr.toFixed(3),
            transform.centerX.toFixed(3),
            transform.centerY.toFixed(3),
            transform.scale.toFixed(6),
            transform.angleDeg.toFixed(3),
            adjustedSpotSize.toFixed(4),
            getModalSubviewRenderToken(),
            getModalBlendSpecCacheKey(modalBlendSpec.a, runtimes?.a),
            getModalBlendSpecCacheKey(modalBlendSpec.b, runtimes?.b),
        ].join('::');
    }}

    function renderModalBlendRuntimeLayer(ctx, section, runtime, transform, adjustedSpotSize, candidateIndices = null) {{
        const nCandidates = candidateIndices ? candidateIndices.length : section.x.length;
        for (let k = 0; k < nCandidates; k++) {{
            const i = candidateIndices ? candidateIndices[k] : k;
            const point = transform.dataToScreen(section.x[i], section.y[i]);
            const x = point.x;
            const y = point.y;
            if (!transform.isPointVisible(x, y, adjustedSpotSize)) continue;
            ctx.fillStyle = rgbToCss(getModalBlendCellRgb(runtime, i));
            ctx.beginPath();
            ctx.arc(x, y, adjustedSpotSize, 0, Math.PI * 2);
            ctx.fill();
        }}
    }}

    function ensureModalBlendRenderCache(section, transform, width, height, dpr, adjustedSpotSize, runtimes, candidateIndices = null) {{
        if (!section || !transform || !runtimes) return null;
        const cacheKey = getModalBlendRenderCacheKey(section, transform, width, height, dpr, adjustedSpotSize, runtimes);
        if (modalBlendRenderCache && modalBlendRenderCache.key === cacheKey) {{
            return modalBlendRenderCache;
        }}

        const layerA = createModalBlendLayerCanvas(width, height, dpr);
        const layerB = createModalBlendLayerCanvas(width, height, dpr);
        renderModalBlendRuntimeLayer(layerA.ctx, section, runtimes.a, transform, adjustedSpotSize, candidateIndices);
        renderModalBlendRuntimeLayer(layerB.ctx, section, runtimes.b, transform, adjustedSpotSize, candidateIndices);

        modalBlendRenderCache = {{
            key: cacheKey,
            width,
            height,
            dpr,
            canvasA: layerA.canvas,
            canvasB: layerB.canvas,
        }};
        return modalBlendRenderCache;
    }}

    function drawModalBlendCachedLayers(ctx, width, height, cache, splitX) {{
        if (!ctx || !cache) return false;
        const clampedSplitX = Math.max(0, Math.min(width, splitX));

        if (clampedSplitX > 0) {{
            ctx.save();
            ctx.beginPath();
            ctx.rect(0, 0, clampedSplitX, height);
            ctx.clip();
            ctx.drawImage(cache.canvasA, 0, 0, width, height);
            ctx.restore();
        }}
        if (clampedSplitX < width) {{
            ctx.save();
            ctx.beginPath();
            ctx.rect(clampedSplitX, 0, width - clampedSplitX, height);
            ctx.clip();
            ctx.drawImage(cache.canvasB, 0, 0, width, height);
            ctx.restore();
        }}

        ctx.strokeStyle = 'rgba(255, 255, 255, 0.8)';
        ctx.lineWidth = 1.5;
        ctx.setLineDash([5, 3]);
        ctx.beginPath();
        ctx.moveTo(clampedSplitX, 0);
        ctx.lineTo(clampedSplitX, height);
        ctx.stroke();
        ctx.setLineDash([]);
        return true;
    }}

    function renderModalBlendFromCache() {{
        if (!modalSection || !modalBlendEnabled) return false;
        if (showGraph || hoverNeighbors || isDrawingModalLasso || isDrawingModalAnnotation) return false;
        if (selectedCells.size > 0 && !getModalActiveCellIndexSet(modalSection.id)) return false;

        const canvas = document.getElementById('modal-canvas');
        const container = document.getElementById('modal-canvas-container');
        if (!canvas || !container) return false;

        clearModalInteractionPreview();
        const ctx = canvas.getContext('2d');
        const dpr = getRenderDpr();
        const rect = container.getBoundingClientRect();
        const width = rect.width;
        const height = rect.height;
        if (width <= 0 || height <= 0) return false;

        const transform = getModalViewTransform(rect);
        const runtimes = getModalBlendRuntimes(modalSection);
        if (!transform || !runtimes || !modalBlendRenderCache) return false;

        const cacheKey = getModalBlendRenderCacheKey(
            modalSection,
            transform,
            width,
            height,
            dpr,
            Math.max(0.25, modalSpotSize * modalZoom * 0.8),
            runtimes,
        );
        if (modalBlendRenderCache.key !== cacheKey) return false;

        canvas.width = width * dpr;
        canvas.height = height * dpr;
        ctx.scale(dpr, dpr);
        ctx.fillStyle = getPanelBg();
        ctx.fillRect(0, 0, width, height);
        drawModalBlendCachedLayers(ctx, width, height, modalBlendRenderCache, width * modalBlendMix);
        cacheModalRenderedView(rect, transform);
        updateModalViewportInfo(transform);
        return true;
    }}

    // Modal rendering
    function renderModalSectionExact() {{
        if (!modalSection) return;
        ensureSectionXY(modalSection);

        const canvas = document.getElementById('modal-canvas');
        clearModalInteractionPreview();
        const ctx = canvas.getContext('2d');
        const dpr = getRenderDpr();
        const container = document.getElementById('modal-canvas-container');
        const rect = container.getBoundingClientRect();
        canvas.width = rect.width * dpr;
        canvas.height = rect.height * dpr;
        ctx.scale(dpr, dpr);

        const width = rect.width, height = rect.height;
        ctx.fillStyle = getPanelBg();
        ctx.fillRect(0, 0, width, height);

        if (modalSection.x.length === 0) {{
            modalRenderedView = null;
            updateModalViewportInfo();
            return;
        }}
        const transform = getModalViewTransform(rect);
        if (!transform) {{
            modalRenderedView = null;
            updateModalViewportInfo();
            return;
        }}
        const adjustedSpotSize = Math.max(0.25, modalSpotSize * modalZoom * 0.8);

        const config = getColorConfig();
        const values = getSectionValues(modalSection);
        const blendRuntimes = getModalBlendRuntimes(modalSection);
        const blendActive = !!blendRuntimes;
        const candidateIndices = filterModalCandidateIndices(
            modalSection,
            getSectionVisibleScreenCandidates(modalSection, transform, adjustedSpotSize),
        );
        const nCandidates = candidateIndices ? candidateIndices.length : modalSection.x.length;
        const focusedSubviewActive = !!getModalActiveCellIndexSet(modalSection.id);

        const modalEdges = getSectionEdgesPacked(modalSection);
        if (showGraph && modalEdges && modalEdges.length) {{
            const graphColor = getGraphColor();
            ctx.strokeStyle = graphColor;
            ctx.lineWidth = Math.max(0.3, Math.min(2.2, modalSpotSize * modalZoom * 0.12));
            ctx.beginPath();
            const n = modalSection.x.length;
            const visibleIndexSet = candidateIndices ? new Set(candidateIndices) : null;
            const drawEdge = (i, j) => {{
                if (i < 0 || j < 0 || i >= n || j >= n) return;
                if (focusedSubviewActive && (!modalSubview.indexSet.has(i) || !modalSubview.indexSet.has(j))) return;
                if (visibleIndexSet && !visibleIndexSet.has(i) && !visibleIndexSet.has(j)) return;
                const p1 = transform.dataToScreen(modalSection.x[i], modalSection.y[i]);
                const p2 = transform.dataToScreen(modalSection.x[j], modalSection.y[j]);
                const x1 = p1.x;
                const y1 = p1.y;
                const x2 = p2.x;
                const y2 = p2.y;
                const minX = Math.min(x1, x2);
                const maxX = Math.max(x1, x2);
                const minY = Math.min(y1, y2);
                const maxY = Math.max(y1, y2);
                if (maxX < -adjustedSpotSize || minX > width + adjustedSpotSize ||
                    maxY < -adjustedSpotSize || minY > height + adjustedSpotSize) return;
                ctx.moveTo(x1, y1);
                ctx.lineTo(x2, y2);
            }};
            if (Array.isArray(modalEdges)) {{
                for (let e = 0; e < modalEdges.length; e++) {{
                    const edge = modalEdges[e];
                    drawEdge(edge[0], edge[1]);
                }}
            }} else if (modalEdges instanceof Uint32Array) {{
                for (let e = 0; e + 1 < modalEdges.length; e += 2) {{
                    drawEdge(modalEdges[e], modalEdges[e + 1]);
                }}
            }}
            ctx.stroke();
        }}

        // First pass: draw grey background for hidden categories
        const hasHidden = !blendActive && hiddenCategories.size > 0 && !config.is_continuous;
        if (hasHidden) {{
            ctx.fillStyle = '#cccccc';
            ctx.globalAlpha = 0.2;
            for (let k = 0; k < nCandidates; k++) {{
                const i = candidateIndices ? candidateIndices[k] : k;
                const val = values[i];
                if (val === null || val === undefined) continue;
                const catIdx = Math.round(val);
                const catName = config.categories[catIdx];
                if (!hiddenCategories.has(catName)) continue;

                const point = transform.dataToScreen(modalSection.x[i], modalSection.y[i]);
                const x = point.x;
                const y = point.y;
                if (!transform.isPointVisible(x, y, adjustedSpotSize)) continue;

                ctx.beginPath();
                ctx.arc(x, y, adjustedSpotSize, 0, Math.PI * 2);
                ctx.fill();
            }}
            ctx.globalAlpha = 1;
        }}

        // Second pass: draw visible categories on top (with optional selected-category focus)
        if (blendActive) {{
            const splitX = width * modalBlendMix;
            for (let k = 0; k < nCandidates; k++) {{
                const i = candidateIndices ? candidateIndices[k] : k;
                const point = transform.dataToScreen(modalSection.x[i], modalSection.y[i]);
                const x = point.x;
                const y = point.y;
                if (!transform.isPointVisible(x, y, adjustedSpotSize)) continue;
                const runtime = x <= splitX ? blendRuntimes.a : blendRuntimes.b;
                ctx.fillStyle = rgbToCss(getModalBlendCellRgb(runtime, i));
                ctx.globalAlpha = 1;
                ctx.beginPath();
                ctx.arc(x, y, adjustedSpotSize, 0, Math.PI * 2);
                ctx.fill();
            }}
            ctx.strokeStyle = 'rgba(255, 255, 255, 0.8)';
            ctx.lineWidth = 1.5;
            ctx.setLineDash([5, 3]);
            ctx.beginPath();
            ctx.moveTo(splitX, 0);
            ctx.lineTo(splitX, height);
            ctx.stroke();
            ctx.setLineDash([]);
            ensureModalBlendRenderCache(modalSection, transform, width, height, dpr, adjustedSpotSize, blendRuntimes, candidateIndices);
        }} else {{
            const activeSpotlight = getLinkedSpotlightCategory(config);
            const focusCategory = activeSpotlight || modalSelectedCategory;
            const hasTypeFocus = !config.is_continuous && focusCategory;
            for (let k = 0; k < nCandidates; k++) {{
                const i = candidateIndices ? candidateIndices[k] : k;
                const val = values[i];
                if (val === null || val === undefined) continue;

                let color;
                let isSelectedCat = false;
                if (config.is_continuous) {{
                    const t = (val - config.vmin) / (config.vmax - config.vmin);
                    color = magma(Math.max(0, Math.min(1, t)));
                }} else {{
                    const catIdx = Math.round(val);
                    const catName = config.categories[catIdx];
                    if (hiddenCategories.has(catName)) continue;
                    isSelectedCat = focusCategory && catName === focusCategory;
                    color = getCategoryColor(catIdx);
                }}

                const point = transform.dataToScreen(modalSection.x[i], modalSection.y[i]);
                const x = point.x;
                const y = point.y;
                if (!transform.isPointVisible(x, y, adjustedSpotSize)) continue;

                if (hasTypeFocus && !isSelectedCat) {{
                    ctx.fillStyle = '#bbbbbb';
                    ctx.globalAlpha = 0.15;
                }} else {{
                    ctx.fillStyle = color;
                    ctx.globalAlpha = 1;
                }}
                ctx.beginPath();
                ctx.arc(x, y, adjustedSpotSize, 0, Math.PI * 2);
                ctx.fill();
            }}
        }}
        ctx.globalAlpha = 1;

        // Third pass: draw selection highlights
        if (selectedCells.size > 0 && !focusedSubviewActive) {{
            ctx.strokeStyle = '#ffd700';
            ctx.lineWidth = 3;
            for (let k = 0; k < nCandidates; k++) {{
                const i = candidateIndices ? candidateIndices[k] : k;
                if (!isCellSelected(modalSection.id, i)) continue;

                const val = values[i];
                if (!blendActive && (val === null || val === undefined)) continue;

                // Skip hidden categories
                if (!blendActive && !config.is_continuous) {{
                    const catIdx = Math.round(val);
                    const catName = config.categories[catIdx];
                    if (hiddenCategories.has(catName)) continue;
                }}

                const point = transform.dataToScreen(modalSection.x[i], modalSection.y[i]);
                const x = point.x;
                const y = point.y;
                if (!transform.isPointVisible(x, y, adjustedSpotSize)) continue;

                ctx.beginPath();
                ctx.arc(x, y, adjustedSpotSize + 2, 0, Math.PI * 2);
                ctx.stroke();
            }}
        }}

        // No extra highlight needed: non-selected categories are greyed out above

        drawNeighborHighlights(ctx, modalSection, transform, adjustedSpotSize);
        drawModalAnnotations(ctx, transform);

        if (isDrawingModalLasso && modalLassoPath.length > 1) {{
            ctx.strokeStyle = 'rgba(255, 255, 255, 0.85)';
            ctx.lineWidth = 2;
            ctx.setLineDash([6, 4]);
            ctx.beginPath();
            ctx.moveTo(modalLassoPath[0].x, modalLassoPath[0].y);
            for (let i = 1; i < modalLassoPath.length; i++) {{
                ctx.lineTo(modalLassoPath[i].x, modalLassoPath[i].y);
            }}
            ctx.stroke();
            ctx.setLineDash([]);
        }}

        if (isDrawingModalAnnotation && modalAnnotationPath.length > 1) {{
            ctx.strokeStyle = 'rgba(255, 127, 80, 0.95)';
            ctx.lineWidth = 2;
            ctx.setLineDash([4, 3]);
            ctx.beginPath();
            ctx.moveTo(modalAnnotationPath[0].x, modalAnnotationPath[0].y);
            for (let i = 1; i < modalAnnotationPath.length; i++) {{
                ctx.lineTo(modalAnnotationPath[i].x, modalAnnotationPath[i].y);
            }}
            ctx.stroke();
            ctx.setLineDash([]);
        }}

        cacheModalRenderedView(rect, transform);
        updateModalViewportInfo(transform);
    }}

    function renderModalSection() {{
        clearModalInteractionCommitTimer();
        renderModalSectionExact();
    }}

    // Legend
    function renderLegend(targetId = 'legend') {{
        const legend = document.getElementById(targetId);
        if (!legend) return;
        const config = getColorConfig();
        const colorLabel = currentGene || currentColor;
        const splitLegendActive = targetId === 'modal-legend' && !!getModalBlendRuntimes(modalSection);
        if (targetId === 'modal-legend') {{
            legend.classList.toggle('split-expanded', splitLegendActive);
        }}

        if (splitLegendActive) {{
            const splitEntries = ['a', 'b']
                .map(side => getModalSplitLegendEntry(side))
                .filter(Boolean);
            let html = `
                <div class="legend-title">Split legend</div>
                <div class="split-legend-list">
            `;
            splitEntries.forEach((entry) => {{
                const categoriesHtml = Array.isArray(entry.categories) && entry.categories.length
                    ? `<div class="split-legend-categories">${{entry.categories.map((cat) => `
                        <div class="split-legend-cat">
                            <span class="split-legend-swatch split-legend-swatch-dot" style="background: ${{cat.color}};"></span>
                            <span class="split-legend-cat-label">${{cat.label}}</span>
                        </div>
                    `).join('')}}</div>`
                    : '';
                html += `
                    <div class="split-legend-item">
                        <div class="split-legend-row">
                            <span class="split-legend-side split-legend-side-${{entry.side}}">${{entry.sideLabel}}</span>
                            <div class="split-legend-label">${{entry.label}}</div>
                        </div>
                        <div class="split-legend-key">
                            <span class="split-legend-swatch ${{entry.swatchClass}}" style="background: ${{entry.swatchBackground}};"></span>
                            <span>${{entry.detail}}</span>
                        </div>
                        ${{categoriesHtml}}
                    </div>
                `;
            }});
            html += '</div>';
            legend.innerHTML = html;
            return;
        }}

        if (config.is_continuous) {{
            legend.innerHTML = `
                <div class="legend-title">${{colorLabel}}</div>
                <div class="colorbar-container">
                    <canvas class="colorbar" id="${{targetId}}-colorbar"></canvas>
                    <div class="colorbar-labels">
                        <span>${{config.vmax.toFixed(2)}}</span>
                        <span>${{((config.vmax + config.vmin) / 2).toFixed(2)}}</span>
                        <span>${{config.vmin.toFixed(2)}}</span>
                    </div>
                </div>
            `;
            const colorbar = document.getElementById(`${{targetId}}-colorbar`);
            const ctx = colorbar.getContext('2d');
            const dpr = getRenderDpr();
            colorbar.width = 16 * dpr;
            colorbar.height = 150 * dpr;
            ctx.scale(dpr, dpr);
            for (let i = 0; i < 150; i++) {{
                ctx.fillStyle = magma(1 - i / 149);
                ctx.fillRect(0, i, 16, 1);
            }}
        }} else {{
            const activeSpotlight = getLinkedSpotlightCategory(config);
            let html = `
                <div class="legend-title">${{colorLabel}}</div>
                <div class="legend-actions">
                    <button class="legend-btn" id="${{targetId}}-show-all">Show All</button>
                    <button class="legend-btn" id="${{targetId}}-hide-all">Hide All</button>
                    <button class="legend-btn ${{linkedSpotlightEnabled ? 'active' : ''}}" id="${{targetId}}-spotlight-toggle" title="Hover or click a category to spotlight">Spotlight</button>
                </div>
            `;
            (config.categories || []).forEach((cat, idx) => {{
                const hiddenClass = hiddenCategories.has(cat) ? 'hidden' : '';
                const selectedClass = modalSelectedCategory && cat === modalSelectedCategory ? 'selected' : '';
                const spotlightClass = activeSpotlight && cat === activeSpotlight ? 'spotlight' : '';
                const dimmedClass = activeSpotlight && cat !== activeSpotlight ? 'dimmed' : '';
                html += `<div class="legend-item ${{hiddenClass}} ${{selectedClass}} ${{spotlightClass}} ${{dimmedClass}}" data-category="${{cat}}">
                    <div class="legend-color" style="background: ${{getCategoryColor(idx)}}"></div>
                    <span>${{cat}}</span>
                </div>`;
            }});
            legend.innerHTML = html;

            document.getElementById(`${{targetId}}-show-all`)?.addEventListener('click', () => {{
                hiddenCategories.clear();
                if (modalSelectedCategory && config.categories?.includes(modalSelectedCategory)) {{
                    // Keep selected category visible; focus is handled by rendering
                    hiddenCategories.delete(modalSelectedCategory);
                }}
                renderLegend('legend');
                renderLegend('modal-legend');
                renderAllSections();
                if (modalSection) renderModalSection();
                if (umapVisible) renderUMAP();
            }});

            document.getElementById(`${{targetId}}-hide-all`)?.addEventListener('click', () => {{
                (config.categories || []).forEach(cat => hiddenCategories.add(cat));
                if (modalSelectedCategory && config.categories?.includes(modalSelectedCategory)) {{
                    hiddenCategories.delete(modalSelectedCategory);
                }}
                renderLegend('legend');
                renderLegend('modal-legend');
                renderAllSections();
                if (modalSection) renderModalSection();
                if (umapVisible) renderUMAP();
            }});

            document.getElementById(`${{targetId}}-spotlight-toggle`)?.addEventListener('click', () => {{
                linkedSpotlightEnabled = !linkedSpotlightEnabled;
                if (!linkedSpotlightEnabled) {{
                    spotlightPinnedCategory = null;
                    spotlightHoverCategory = null;
                }}
                updateAllLegendSpotlightClasses();
                rerenderForSpotlightChange();
            }});

            legend.querySelectorAll('.legend-item').forEach(item => {{
                item.addEventListener('mouseenter', () => {{
                    if (!linkedSpotlightEnabled) return;
                    const cat = item.dataset.category;
                    if (!cat || spotlightHoverCategory === cat) return;
                    spotlightHoverCategory = cat;
                    updateAllLegendSpotlightClasses();
                    rerenderForSpotlightChange();
                }});
                item.addEventListener('mouseleave', () => {{
                    if (!linkedSpotlightEnabled) return;
                    const cat = item.dataset.category;
                    if (!cat || spotlightHoverCategory !== cat) return;
                    spotlightHoverCategory = null;
                    updateAllLegendSpotlightClasses();
                    rerenderForSpotlightChange();
                }});
                item.addEventListener('click', () => {{
                    const cat = item.dataset.category;
                    if (linkedSpotlightEnabled) {{
                        if (hiddenCategories.has(cat)) hiddenCategories.delete(cat);
                        if (spotlightPinnedCategory === cat) spotlightPinnedCategory = null;
                        else spotlightPinnedCategory = cat;
                        spotlightHoverCategory = null;
                        renderLegend('legend');
                        renderLegend('modal-legend');
                        rerenderForSpotlightChange();
                        return;
                    }}
                    if (hiddenCategories.has(cat)) hiddenCategories.delete(cat);
                    else hiddenCategories.add(cat);
                    renderLegend('legend');
                    renderLegend('modal-legend');
                    renderAllSections();
                    if (modalSection) renderModalSection();
                    if (umapVisible) renderUMAP();
                }});
            }});
            updateLegendSpotlightClasses(targetId);
        }}
    }}

    function updateAnnotationComparisonTabVisibility() {{
        const tab = document.getElementById('color-tab-regions');
        if (!tab) return;
        const hasAnnotations = modalAnnotations.length > 0;
        tab.style.display = hasAnnotations ? '' : 'none';
        // If tab was active but no annotations remain, switch to Stats
        if (!hasAnnotations && tab.classList.contains('active')) {{
            tab.classList.remove('active');
            document.getElementById('color-tab-regions-content')?.classList.remove('active');
            document.getElementById('color-tab-aggregate')?.classList.add('active');
            document.getElementById('color-tab-aggregate-content')?.classList.add('active');
        }}
    }}

    function renderAnnotationRegionDESection(annotations) {{
        const {{ annotations: availableAnnotations, source, reference }} = syncAnnotationDEState(annotations);
        let html = '<div class="annot-de-controls">';
        html += '<div class="selection-summary-title" style="margin-top:10px">Region-to-Region DE</div>';
        html += '<div class="agg-group-meta">Quick DE uses genes already loaded in the viewer. Full region DE can scan the sidecar in the background when available.</div>';

        if (availableAnnotations.length < 2) {{
            html += '<div class="agg-group-meta">Draw at least two annotations to compare region-level DE.</div>';
            html += '</div>';
            return html;
        }}

        const options = availableAnnotations.map((annotation) => {{
            const id = Number(annotation.id);
            const label = escapeHtml(getAnnotationDisplayLabel(annotation));
            return `<option value="${{id}}">${{label}}</option>`;
        }}).join('');

        html += `
            <div class="cluster-de-controls">
                <div class="cluster-de-select-row">
                    <div>
                        <label>Region A</label>
                        <select id="annotation-de-source">${{options}}</select>
                    </div>
                    <div>
                        <label>Region B</label>
                        <select id="annotation-de-reference">${{options}}</select>
                    </div>
                    <div>
                        <label>Top Genes</label>
                        <input id="annotation-de-topn" type="number" min="3" max="100" step="1" value="${{Math.max(3, Number(annotationDeTopN) || ANNOTATION_DE_TOP_N)}}">
                    </div>
                </div>
                <div style="display: flex; justify-content: flex-end;">
                    <button class="legend-btn" id="annotation-de-swap" type="button">Swap A/B</button>
                </div>
            </div>
        `;

        const deResult = computeRegionAnnotationDE(source, reference, {{ topN: annotationDeTopN }});
        const exportState = getAnnotationDEExportState(source, reference, deResult);
        const pairKey = getAnnotationDECacheKey(source, reference);
        const fullCached = pairKey ? annotationDeFullCache.get(pairKey) : null;
        const fullRun = (annotationDeFullRun && annotationDeFullRun.key === pairKey) ? annotationDeFullRun : null;
        const sidecarAvailable = !!DATA.gene_aux_url;
        const renderCards = (result) => {{
            const topN = Math.max(1, Number(annotationDeTopN) || ANNOTATION_DE_TOP_N);
            return (result.results || []).slice(0, topN).map((entry) => {{
                return `
                    <div class="comparison-card">
                        <div class="comparison-card-title">
                            ${{renderGeneTokenButton(entry.gene, {{
                                isActive: entry.gene === currentGene,
                                showMeta: false,
                                title: 'Load region DE gene into the viewer',
                            }})}}
                            ${{renderGeneGoogleSearchButton(entry.gene, {{
                                title: 'Search Google for this gene',
                            }})}}
                        </div>
                        <div class="comparison-metric-grid">
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">log2FC A/B</span>
                                <span class="comparison-metric-value">${{formatScaleNumber(entry.log2fc)}}</span>
                            </div>
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">Score</span>
                                <span class="comparison-metric-value">${{formatScaleNumber(entry.score)}}</span>
                            </div>
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">% expr A</span>
                                <span class="comparison-metric-value">${{formatClusterDEPct(entry.pctA)}}</span>
                            </div>
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">% expr B</span>
                                <span class="comparison-metric-value">${{formatClusterDEPct(entry.pctB)}}</span>
                            </div>
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">Mean A</span>
                                <span class="comparison-metric-value">${{formatScaleNumber(entry.meanA)}}</span>
                            </div>
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">Mean B</span>
                                <span class="comparison-metric-value">${{formatScaleNumber(entry.meanB)}}</span>
                            </div>
                        </div>
                    </div>
                `;
            }}).join('');
        }};
        const sourceLabel = source ? escapeHtml(getAnnotationDisplayLabel(source)) : '';
        const referenceLabel = reference ? escapeHtml(getAnnotationDisplayLabel(reference)) : '';
        const topNLabel = Math.max(1, Number(annotationDeTopN) || ANNOTATION_DE_TOP_N).toLocaleString();
        const quickSummaryHtml = (deResult.available && deResult.results.length)
            ? `
                <div class="agg-group-meta">
                    Region A: <strong>${{sourceLabel}}</strong> (${{Number(deResult.nA || 0).toLocaleString()}} cells)
                    vs Region B: <strong>${{referenceLabel}}</strong> (${{Number(deResult.nB || 0).toLocaleString()}} cells).
                </div>
                <div class="agg-group-meta">
                    ${{
                        Number(deResult.loadedGeneCount || 0) < Number(deResult.totalGeneCount || 0)
                            ? `Using ${{Number(deResult.loadedGeneCount || 0).toLocaleString()}} of ${{Number(deResult.totalGeneCount || 0).toLocaleString()}} genes currently loaded in the viewer.`
                            : `Using ${{Number(deResult.loadedGeneCount || 0).toLocaleString()}} loaded genes.`
                    }}
                    Showing top ${{topNLabel}} genes. Positive scores indicate enrichment in Region A.
                </div>
                <div class="comparison-stack">${{renderCards(deResult)}}</div>
            `
            : '';
        const cachedTimestamp = fullCached?.completedAt
            ? new Date(fullCached.completedAt).toLocaleTimeString([], {{ hour: '2-digit', minute: '2-digit', second: '2-digit' }})
            : '';
        const cachedSummaryHtml = fullCached?.available
            ? (fullCached.results || []).length
                ? `
                    <div class="agg-group-meta">
                        Region A: <strong>${{sourceLabel}}</strong> (${{Number(fullCached.nA || 0).toLocaleString()}} cells)
                        vs Region B: <strong>${{referenceLabel}}</strong> (${{Number(fullCached.nB || 0).toLocaleString()}} cells).
                    </div>
                    <div class="agg-group-meta">
                        Cached full sidecar DE${{cachedTimestamp ? ` from ${{cachedTimestamp}}` : ''}}
                        across ${{Number(fullCached.totalGeneCount || 0).toLocaleString()}} genes.
                        Showing top ${{topNLabel}} genes. Positive scores indicate enrichment in Region A.
                    </div>
                    <div class="comparison-stack">${{renderCards(fullCached)}}</div>
                `
                : `
                    <div class="agg-group-meta">
                        Region A: <strong>${{sourceLabel}}</strong> (${{Number(fullCached.nA || 0).toLocaleString()}} cells)
                        vs Region B: <strong>${{referenceLabel}}</strong> (${{Number(fullCached.nB || 0).toLocaleString()}} cells).
                    </div>
                    <div class="agg-group-meta">
                        Cached full sidecar DE${{cachedTimestamp ? ` from ${{cachedTimestamp}}` : ''}} found no enriched genes across ${{Number(fullCached.totalGeneCount || 0).toLocaleString()}} genes.
                    </div>
                `
            : '';
        let resultHtml = '';
        if (!source || !reference) {{
            resultHtml = '<div class="agg-group-meta">Choose two different annotations to compare.</div>';
        }} else if (fullRun?.running) {{
            const totalShards = Number(fullRun.totalShards || 0);
            const completedShards = Number(fullRun.completedShards || 0);
            const totalGenes = Number(fullRun.totalGenes || 0);
            const completedGenes = Number(fullRun.completedGenes || 0);
            const progressMax = Math.max(1, totalGenes || totalShards || 1);
            const progressValue = Math.min(progressMax, totalGenes ? completedGenes : completedShards);
            const progressPct = progressMax > 0 ? Math.round((progressValue / progressMax) * 100) : 0;
            resultHtml = `
                <div class="agg-group-meta">Scanning the full sidecar in the background. The viewer stays interactive while this runs.</div>
                <progress value="${{progressValue}}" max="${{progressMax}}" style="width:100%; height:12px;"></progress>
                <div class="agg-group-meta">
                    ${{progressPct}}% complete
                    · ${{completedGenes.toLocaleString()}} / ${{totalGenes.toLocaleString()}} genes
                    · ${{completedShards.toLocaleString()}} / ${{totalShards.toLocaleString()}} shards
                </div>
                <div style="display:flex; justify-content:flex-end;">
                    <button class="legend-btn" id="annotation-de-cancel" type="button">Cancel Full DE</button>
                </div>
            `;
            if (cachedSummaryHtml) {{
                resultHtml += '<div class="agg-group-meta">Showing the cached full result while the refresh runs:</div>';
                resultHtml += cachedSummaryHtml;
            }} else if (quickSummaryHtml) {{
                resultHtml += '<div class="agg-group-meta">Preview from the genes already loaded in the viewer:</div>';
                resultHtml += quickSummaryHtml;
            }}
        }} else if (fullRun?.error) {{
            resultHtml = `
                <div class="agg-group-meta">Full sidecar DE failed: ${{escapeHtml(fullRun.error)}}</div>
                <div style="display:flex; justify-content:flex-end;">
                    <button class="legend-btn" id="annotation-de-run-full" type="button">Retry Full Region DE</button>
                </div>
            `;
            if (cachedSummaryHtml) {{
                resultHtml += '<div class="agg-group-meta">Showing the previous cached full result:</div>';
                resultHtml += cachedSummaryHtml;
            }} else if (quickSummaryHtml) {{
                resultHtml += '<div class="agg-group-meta">Quick DE preview from currently loaded genes:</div>';
                resultHtml += quickSummaryHtml;
            }}
        }} else if (fullCached?.available) {{
            resultHtml = cachedSummaryHtml;
            resultHtml += '<div style="display:flex; justify-content:flex-end;"><button class="legend-btn" id="annotation-de-refresh-full" type="button">Refresh Full DE</button></div>';
        }} else if (!deResult.available && deResult.reason === 'empty_region') {{
            resultHtml = '<div class="agg-group-meta">One of the selected annotations has no cells.</div>';
        }} else if (!deResult.available && deResult.reason === 'no_loaded_genes' && sidecarAvailable) {{
            const totalGenes = Number(deResult.totalGeneCount || 0).toLocaleString();
            resultHtml = `
                <div class="agg-group-meta">
                    No genes are currently loaded for quick region DE.
                    ${{totalGenes !== '0' ? `This viewer knows about ${{totalGenes}} sidecar genes that can be scanned on demand.` : ''}}
                </div>
                <div style="display:flex; justify-content:flex-end;">
                    <button class="legend-btn" id="annotation-de-run-full" type="button">Run Full Region DE</button>
                </div>
            `;
        }} else if (!deResult.available && deResult.reason === 'no_loaded_genes') {{
            resultHtml = '<div class="agg-group-meta">No genes are currently loaded for region DE. Load genes in the Genes tab or click marker genes first.</div>';
        }} else if (!deResult.available) {{
            resultHtml = '<div class="agg-group-meta">Choose two different annotations to compare.</div>';
        }} else if (!deResult.results.length) {{
            if (sidecarAvailable && Number(deResult.loadedGeneCount || 0) < Number(deResult.totalGeneCount || 0)) {{
                resultHtml = `
                    <div class="agg-group-meta">No enriched genes were found among the currently loaded genes.</div>
                    <div style="display:flex; justify-content:flex-end;">
                        <button class="legend-btn" id="annotation-de-run-full" type="button">Run Full Region DE</button>
                    </div>
                `;
            }} else {{
                resultHtml = '<div class="agg-group-meta">No enriched genes were found among the genes currently loaded in the viewer.</div>';
            }}
        }} else {{
            resultHtml = quickSummaryHtml;
            if (sidecarAvailable && Number(deResult.loadedGeneCount || 0) < Number(deResult.totalGeneCount || 0)) {{
                resultHtml += '<div style="display:flex; justify-content:flex-end;"><button class="legend-btn" id="annotation-de-run-full" type="button">Run Full Region DE</button></div>';
            }}
        }}

        if (exportState) {{
            resultHtml += `
                <div style="display:flex; justify-content:flex-end; gap:6px; margin-top:6px;">
                    <button class="legend-btn" id="annotation-de-export-json" type="button">Export JSON</button>
                    <button class="legend-btn" id="annotation-de-export-csv" type="button">Export CSV</button>
                </div>
            `;
        }}

        html += `<div id="annotation-de-results">${{resultHtml}}</div>`;
        html += '</div>';
        return html;
    }}

    function renderAnnotationComparison() {{
        const container = document.getElementById('annotation-comparison');
        if (!container) return;

        if (!modalAnnotations.length) {{
            container.innerHTML = '<div class="agg-group-meta">Open a section, draw annotations, then switch to this tab to compare them.</div>';
            return;
        }}

        // Build per-annotation summaries
        const annotSummaries = modalAnnotations.map((annotation) => {{
            const cellSet = new Set(
                (annotation.localCellIndices || [])
                    .filter(idx => Number.isInteger(idx) && idx >= 0)
                    .map(idx => `${{annotation.sectionId}}:${{idx}}`)
            );
            const summary = computeSelectionSummary(cellSet);
            return {{ annotation, summary }};
        }});

        let html = '';

        // Composition section: stacked bars per annotation
        const typeColumn = annotSummaries[0]?.summary?.typeColumn;
        if (typeColumn) {{
            // Collect all type names
            const allTypes = new Set();
            annotSummaries.forEach((as) => as.summary.types.forEach(([t]) => allTypes.add(t)));
            const typeList = Array.from(allTypes);

            html += `<div class="selection-summary-title">Cell Composition by Annotation</div>`;
            annotSummaries.forEach((as) => {{
                const annotation = as.annotation;
                const summary = as.summary;
                const color = annotation.color || '#888';
                const label = escapeHtml(annotation.label || `Annotation ${{annotation.id}}`);
                const total = summary.total;
                const typeMap = new Map(summary.types);

                html += `<div class="annot-comp-row">`;
                html += `<div class="annot-comp-header"><span class="annot-comp-dot" style="background:${{color}}"></span><span>${{label}}</span><span class="annot-comp-count">${{total.toLocaleString()}} cells</span></div>`;
                if (total > 0) {{
                    html += `<div class="annot-comp-bar-track">`;
                    let offset = 0;
                    typeList.forEach((type, ti) => {{
                        const cnt = typeMap.get(type) || 0;
                        const pct = Math.round(100 * cnt / total);
                        if (pct <= 0) return;
                        const bgColor = getBarColor(ti);
                        html += `<div class="annot-comp-segment" title="${{escapeHtml(type)}}: ${{cnt}} (${{pct}}%)" style="width:${{pct}}%;background:${{bgColor}}"></div>`;
                    }});
                    html += `</div>`;
                }}
                html += `</div>`;
            }});

            // Legend
            if (typeList.length) {{
                html += `<div class="annot-comp-legend">`;
                typeList.slice(0, 8).forEach((type, ti) => {{
                    html += `<span class="annot-comp-legend-item"><span class="annot-comp-dot" style="background:${{getBarColor(ti)}}"></span>${{escapeHtml(type)}}</span>`;
                }});
                if (typeList.length > 8) {{
                    html += `<span class="annot-comp-legend-item agg-group-meta">+${{typeList.length - 8}} more</span>`;
                }}
                html += `</div>`;
            }}
        }}

        // Expression section: top 6 genes, one row per gene
        const exprArrays = annotSummaries.map((as) => computeSelectionMeanExpression(as.summary));
        const firstNonEmpty = exprArrays.find(arr => arr.length > 0);
        if (firstNonEmpty) {{
            const topGenes = firstNonEmpty.slice(0, 6).map(e => e.gene);
            const vmax = Math.max(1, ...exprArrays.flatMap(arr =>
                topGenes.map(g => (arr.find(e => e.gene === g)?.meanSel || 0))
            ));

            html += `<div class="selection-summary-title" style="margin-top:10px">Gene Expression by Annotation</div>`;
            topGenes.forEach((gene) => {{
                html += `<div class="annot-expr-row"><span class="selection-summary-expr-gene">${{escapeHtml(gene)}}</span>`;
                html += `<div class="annot-expr-bars">`;
                exprArrays.forEach((arr, ai) => {{
                    const val = arr.find(e => e.gene === gene)?.meanSel || 0;
                    const pct = Math.round(100 * val / vmax);
                    const color = annotSummaries[ai].annotation.color || '#888';
                    html += `<div class="annot-expr-bar" style="width:${{pct}}%;background:${{color}}" title="${{escapeHtml(annotSummaries[ai].annotation.label || '')}}: ${{val.toFixed(2)}}"></div>`;
                }});
                html += `</div></div>`;
            }});
        }}

        html += renderAnnotationRegionDESection(modalAnnotations);

        container.innerHTML = html || '<div class="agg-group-meta">No data available.</div>';
        const {{ source, reference }} = syncAnnotationDEState(modalAnnotations);
        const sourceSelect = container.querySelector('#annotation-de-source');
        const referenceSelect = container.querySelector('#annotation-de-reference');
        const topNInput = container.querySelector('#annotation-de-topn');
        const swapBtn = container.querySelector('#annotation-de-swap');
        if (sourceSelect) {{
            sourceSelect.value = annotationDeSourceId ? String(annotationDeSourceId) : '';
            sourceSelect.addEventListener('change', () => {{
                cancelAnnotationFullDERun();
                annotationDeSourceId = Number(sourceSelect.value) || null;
                renderAnnotationComparison();
            }});
        }}
        if (referenceSelect) {{
            referenceSelect.value = annotationDeReferenceId ? String(annotationDeReferenceId) : '';
            referenceSelect.addEventListener('change', () => {{
                cancelAnnotationFullDERun();
                annotationDeReferenceId = Number(referenceSelect.value) || null;
                renderAnnotationComparison();
            }});
        }}
        if (topNInput) {{
            topNInput.value = String(Math.max(3, Number(annotationDeTopN) || ANNOTATION_DE_TOP_N));
            topNInput.addEventListener('change', () => {{
                const next = Math.min(100, Math.max(3, Number(topNInput.value) || ANNOTATION_DE_TOP_N));
                annotationDeTopN = next;
                renderAnnotationComparison();
            }});
        }}
        if (swapBtn) {{
            swapBtn.addEventListener('click', () => {{
                cancelAnnotationFullDERun();
                const sourceId = annotationDeSourceId;
                annotationDeSourceId = annotationDeReferenceId;
                annotationDeReferenceId = sourceId;
                renderAnnotationComparison();
            }});
        }}
        const runFullBtn = container.querySelector('#annotation-de-run-full');
        if (runFullBtn && source && reference) {{
            runFullBtn.addEventListener('click', () => {{
                runFullRegionAnnotationDE(source, reference);
            }});
        }}
        const refreshFullBtn = container.querySelector('#annotation-de-refresh-full');
        if (refreshFullBtn && source && reference) {{
            refreshFullBtn.addEventListener('click', () => {{
                runFullRegionAnnotationDE(source, reference);
            }});
        }}
        const cancelBtn = container.querySelector('#annotation-de-cancel');
        if (cancelBtn) {{
            cancelBtn.addEventListener('click', () => {{
                cancelAnnotationFullDERun();
                renderAnnotationComparison();
            }});
        }}
        const exportJsonBtn = container.querySelector('#annotation-de-export-json');
        if (exportJsonBtn && source && reference) {{
            exportJsonBtn.addEventListener('click', () => {{
                const quickResult = computeRegionAnnotationDE(source, reference, {{ topN: annotationDeTopN }});
                const exportState = getAnnotationDEExportState(source, reference, quickResult);
                exportAnnotationDEReport(source, reference, exportState);
            }});
        }}
        const exportCsvBtn = container.querySelector('#annotation-de-export-csv');
        if (exportCsvBtn && source && reference) {{
            exportCsvBtn.addEventListener('click', () => {{
                const quickResult = computeRegionAnnotationDE(source, reference, {{ topN: annotationDeTopN }});
                const exportState = getAnnotationDEExportState(source, reference, quickResult);
                exportAnnotationDECsv(source, reference, exportState);
            }});
        }}
        bindGeneActivateButtons(container, renderAnnotationComparison);
        bindGeneGoogleSearchButtons(container);
    }}

    function getBarColor(index) {{
        const palette = ['#4cc9f0','#f4a261','#2a9d8f','#e63946','#06d6a0','#ffd166','#118ab2','#ef476f','#43aa8b','#ff7f50'];
        return palette[index % palette.length];
    }}

    function buildColorPanel() {{
        const panel = document.getElementById('color-panel');
        if (!panel) return;
        const metadataKeys = Object.keys(DATA.metadata_filters || {{}});
        const hasMetadata = metadataKeys.length > 0;

        const options = ['<option value="">None</option>']
            .concat(metadataKeys.map(key => `<option value="${{key}}">${{formatMetadataLabel(key)}}</option>`))
            .join('');

        panel.innerHTML = `
            <div class="color-panel-header">
                <div class="color-panel-title">Insights</div>
            </div>
            <div class="color-panel-section">
                <label>Search Colors</label>
                <input class="color-search" id="color-search" type="text" placeholder="Type to filter...">
            </div>
            <div class="color-panel-section">
                <label>Available Colors</label>
                <div class="color-list" id="color-list"></div>
            </div>
            <div class="color-panel-section">
                <label>Details</label>
                <div class="color-tabs">
                    <button class="color-tab active" id="color-tab-aggregate" type="button">Summary</button>
                    <button class="color-tab" id="color-tab-compare" type="button">Compare</button>
                    <button class="color-tab" id="color-tab-genes" type="button">Genes</button>
                    <button class="color-tab" id="color-tab-neighbors" type="button">Neighborhood</button>
                    <button class="color-tab" id="color-tab-regions" type="button" style="display:none">Regions</button>
                </div>
                <div class="color-tab-content active" id="color-tab-aggregate-content">
                    <div>
                        <label>Aggregate By</label>
                        <select id="color-groupby" ${{!hasMetadata ? 'disabled' : ''}}>
                            ${{options}}
                        </select>
                    </div>
                    <div style="display: flex; justify-content: flex-end;">
                        <button class="legend-btn" id="color-aggregation-toggle" type="button">Collapse Stats</button>
                    </div>
                    <div class="color-aggregation" id="color-aggregation">
                        <div class="agg-group-meta">${{hasMetadata ? 'Pick a metadata column to summarize.' : 'No metadata columns available for aggregation.'}}</div>
                    </div>
                    <div>
                        <label>Cell Type Trend</label>
                        <input class="color-search" id="celltype-search" type="text" placeholder="Search Cell Type...">
                    </div>
                    <div class="color-aggregation" id="celltype-trend">
                        <div class="agg-group-meta">${{hasMetadata ? 'Search for a category to see counts across the selected metadata.' : 'No metadata columns available.'}}</div>
                    </div>
                </div>
                <div class="color-tab-content" id="color-tab-compare-content">
                    <div class="cluster-de-controls">
                        <div>
                            <label>Annotation</label>
                            <select id="cluster-de-groupby"></select>
                        </div>
                        <div class="cluster-de-select-row">
                            <div>
                                <label>Category A</label>
                                <select id="cluster-de-source"></select>
                            </div>
                            <div>
                                <label>Category B</label>
                                <select id="cluster-de-reference"></select>
                            </div>
                        </div>
                        <div style="display: flex; justify-content: flex-end;">
                            <button class="legend-btn" id="cluster-de-swap" type="button">Swap A/B</button>
                        </div>
                    </div>
                    <div class="color-aggregation" id="cluster-de-results">
                        <div class="agg-group-meta">Choose two categories to compare.</div>
                    </div>
                </div>
                <div class="color-tab-content" id="color-tab-neighbors-content">
                    <div>
                        <label>Search Cell Type</label>
                        <input class="color-search" id="neighbor-search" type="text" placeholder="Search Cell Type...">
                    </div>
                    <div style="display: flex; justify-content: flex-end;">
                        <button class="legend-btn" id="neighbor-stats-toggle" type="button">Collapse Neighbor Stats</button>
                    </div>
                    <div class="color-aggregation" id="neighbor-stats">
                        <div class="agg-group-meta">Select a categorical color to view neighbor stats.</div>
                    </div>
                    <div>
                        <label>Interaction Source</label>
                        <select id="interaction-source"></select>
                    </div>
                    <div>
                        <label>Target Filter</label>
                        <input class="color-search" id="interaction-search" type="text" placeholder="Filter Target Cell Types...">
                    </div>
                    <div class="color-aggregation" id="interaction-browser">
                        <div class="agg-group-meta">Select a source cell type to browse interactions.</div>
                    </div>
                </div>
                <div class="color-tab-content" id="color-tab-genes-content">
                    <div class="color-tabs">
                        <button class="color-tab active" id="genes-tab-dotplot" type="button">Dotplot</button>
                        <button class="color-tab" id="genes-tab-markers" type="button">Markers</button>
                    </div>
                    <div class="color-tab-content active" id="genes-tab-dotplot-content">
                        <div class="dotplot-controls">
                            <div>
                                <label>Group By (Categorical Color)</label>
                                <select id="dotplot-groupby"></select>
                            </div>
                            <div>
                                <label>Aggregate By (Metadata)</label>
                                <select id="dotplot-aggregate-by" ${{!hasMetadata ? 'disabled' : ''}}>
                                    ${{options}}
                                </select>
                            </div>
                            <div id="dotplot-aggregate-value-wrap" style="display: none;">
                                <label>Aggregate Value</label>
                                <select id="dotplot-aggregate-value"></select>
                            </div>
                            <div>
                                <label>Genes (Comma-Separated)</label>
                                <input class="color-search" id="dotplot-genes" type="text" placeholder="e.g. Cd4, Cd8a, Gfap" list="gene-list">
                                <div class="scale-hint">Dot size = % expressing; dot color = mean expression. Press Tab to autocomplete the current gene token.</div>
                            </div>
                            <div style="display: flex; gap: 6px; justify-content: flex-end;">
                                <button class="legend-btn" id="dotplot-use-hvgs" type="button">Suggest Genes</button>
                                <button class="legend-btn" id="dotplot-run" type="button">Update</button>
                            </div>
                            <div class="agg-group-meta" id="dotplot-status">Pick a categorical color + genes to compute a dotplot.</div>
                            <div class="dotplot-grid" id="dotplot-grid"></div>
                        </div>
                    </div>
                    <div class="color-tab-content" id="genes-tab-markers-content">
                        <input class="marker-search" id="marker-gene-search" type="text" placeholder="Search marker genes...">
                        <div class="marker-genes" id="marker-genes"></div>
                    </div>
                </div>
                <div class="color-tab-content" id="color-tab-regions-content">
                    <div class="color-aggregation" id="annotation-comparison">
                        <div class="agg-group-meta">Open a section, draw annotations, then switch to this tab to compare them.</div>
                    </div>
                </div>
            </div>
        `;

        const search = document.getElementById('color-search');
        search.addEventListener('input', () => {{
            renderColorList(search.value);
        }});

        const groupBy = document.getElementById('color-groupby');
        groupBy.addEventListener('change', () => {{
            renderColorAggregation();
            renderCellTypeTrend();
            const dpAgg = document.getElementById('dotplot-aggregate-by');
            if (dpAgg) dpAgg.value = groupBy.value;
            updateDotplotAggregateValueOptions();
            const genesTop = document.getElementById('color-tab-genes');
            const genesDot = document.getElementById('genes-tab-dotplot');
            if (genesTop?.classList.contains('active') && genesDot?.classList.contains('active')) {{
                renderDotplot();
            }}
        }});

        const aggregationToggle = document.getElementById('color-aggregation-toggle');
        const aggregationContainer = document.getElementById('color-aggregation');
        aggregationToggle.addEventListener('click', () => {{
            const isCollapsed = aggregationContainer.classList.toggle('collapsed');
            aggregationToggle.textContent = isCollapsed ? 'Show stats' : 'Collapse stats';
        }});

        const aggregateTab = document.getElementById('color-tab-aggregate');
        const compareTab = document.getElementById('color-tab-compare');
        const genesTab = document.getElementById('color-tab-genes');
        const neighborTab = document.getElementById('color-tab-neighbors');
        const regionsTab = document.getElementById('color-tab-regions');
        const aggregateContent = document.getElementById('color-tab-aggregate-content');
        const compareContent = document.getElementById('color-tab-compare-content');
        const genesContent = document.getElementById('color-tab-genes-content');
        const neighborContent = document.getElementById('color-tab-neighbors-content');
        const regionsContent = document.getElementById('color-tab-regions-content');

        const genesDotTab = document.getElementById('genes-tab-dotplot');
        const genesMarkersTab = document.getElementById('genes-tab-markers');
        const genesDotContent = document.getElementById('genes-tab-dotplot-content');
        const genesMarkersContent = document.getElementById('genes-tab-markers-content');
        const topLevelTabs = [
            [aggregateTab, aggregateContent],
            [compareTab, compareContent],
            [genesTab, genesContent],
            [neighborTab, neighborContent],
            [regionsTab, regionsContent],
        ];
        const activateTopLevelInsightsTab = (activeTab, activeContent) => {{
            topLevelTabs.forEach(([tab, content]) => {{
                tab?.classList.toggle('active', tab === activeTab);
                content?.classList.toggle('active', content === activeContent);
            }});
        }};

        aggregateTab.addEventListener('click', () => {{
            activateTopLevelInsightsTab(aggregateTab, aggregateContent);
            renderColorAggregation();
            renderCellTypeTrend();
        }});
        compareTab.addEventListener('click', () => {{
            activateTopLevelInsightsTab(compareTab, compareContent);
            renderClusterDE();
        }});
        genesTab.addEventListener('click', () => {{
            activateTopLevelInsightsTab(genesTab, genesContent);
            if (!genesDotTab.classList.contains('active') && !genesMarkersTab.classList.contains('active')) {{
                genesDotTab.classList.add('active');
                genesMarkersTab.classList.remove('active');
                genesDotContent.classList.add('active');
                genesMarkersContent.classList.remove('active');
            }}
            if (genesDotTab.classList.contains('active')) renderDotplot();
            else renderMarkerGenes();
        }});
        neighborTab.addEventListener('click', () => {{
            activateTopLevelInsightsTab(neighborTab, neighborContent);
            renderNeighborStats();
            renderInteractionBrowser();
        }});

        genesDotTab.addEventListener('click', () => {{
            genesDotTab.classList.add('active');
            genesMarkersTab.classList.remove('active');
            genesDotContent.classList.add('active');
            genesMarkersContent.classList.remove('active');
            renderDotplot();
        }});
        genesMarkersTab.addEventListener('click', () => {{
            genesMarkersTab.classList.add('active');
            genesDotTab.classList.remove('active');
            genesMarkersContent.classList.add('active');
            genesDotContent.classList.remove('active');
            renderMarkerGenes();
        }});
        regionsTab?.addEventListener('click', () => {{
            activateTopLevelInsightsTab(regionsTab, regionsContent);
            renderAnnotationComparison();
        }});

        const markerSearch = document.getElementById('marker-gene-search');
        markerSearch.addEventListener('input', () => {{
            if (genesTab.classList.contains('active') && genesMarkersTab.classList.contains('active')) {{
                renderMarkerGenes();
            }}
        }});

        const celltypeSearch = document.getElementById('celltype-search');
        celltypeSearch.addEventListener('input', () => {{
            renderCellTypeTrend();
        }});

        const neighborSearch = document.getElementById('neighbor-search');
        neighborSearch.addEventListener('input', () => {{
            renderNeighborStats();
        }});
        const neighborStatsToggle = document.getElementById('neighbor-stats-toggle');
        const neighborStatsContainer = document.getElementById('neighbor-stats');
        neighborStatsToggle.addEventListener('click', () => {{
            const isCollapsed = neighborStatsContainer.classList.toggle('collapsed');
            neighborStatsToggle.textContent = isCollapsed ? 'Show neighbor stats' : 'Collapse neighbor stats';
        }});
        const interactionSource = document.getElementById('interaction-source');
        interactionSource.addEventListener('change', () => {{
            interactionSourceCategory = interactionSource.value || null;
            renderInteractionBrowser();
        }});
        const interactionSearch = document.getElementById('interaction-search');
        interactionSearch.addEventListener('input', () => {{
            renderInteractionBrowser();
        }});

        // Dotplot setup
        const dotplotGroupby = document.getElementById('dotplot-groupby');
        const dotplotGenes = document.getElementById('dotplot-genes');
        const dotplotRun = document.getElementById('dotplot-run');
        const dotplotUseHvgs = document.getElementById('dotplot-use-hvgs');
        const dotplotAgg = document.getElementById('dotplot-aggregate-by');
        const dotplotAggValue = document.getElementById('dotplot-aggregate-value');
        const clusterDeGroupbySelect = document.getElementById('cluster-de-groupby');
        const clusterDeSourceSelect = document.getElementById('cluster-de-source');
        const clusterDeReferenceSelect = document.getElementById('cluster-de-reference');
        const clusterDeSwap = document.getElementById('cluster-de-swap');

        const catCols = getCategoricalColorColumns();
        if (dotplotGroupby) {{
            dotplotGroupby.innerHTML = catCols.map(c => `<option value="${{c}}">${{c}}</option>`).join('');
            if (catCols.includes(currentColor)) dotplotGroupby.value = currentColor;
            dotplotGroupby.addEventListener('change', () => renderDotplot());
        }}
        if (dotplotAgg) {{
            dotplotAgg.value = groupBy.value;
            dotplotAgg.addEventListener('change', () => {{
                groupBy.value = dotplotAgg.value;
                renderColorAggregation();
                renderCellTypeTrend();
                updateDotplotAggregateValueOptions();
                renderDotplot();
            }});
        }}
        updateDotplotAggregateValueOptions();
        dotplotAggValue?.addEventListener('change', () => renderDotplot());
        dotplotRun?.addEventListener('click', () => renderDotplot());
        dotplotGenes?.addEventListener('keydown', (e) => {{
            if (e.key !== 'Tab') return;
            if (e.altKey || e.ctrlKey || e.metaKey) return;
            if (autocompleteDotplotGeneToken(dotplotGenes)) e.preventDefault();
        }});
        dotplotGenes?.addEventListener('change', () => renderDotplot());
        dotplotUseHvgs?.addEventListener('click', () => {{
            const hvgs = (DATA.available_genes || []).slice(0, 8);
            if (dotplotGenes) dotplotGenes.value = hvgs.join(', ');
            renderDotplot();
        }});
        syncClusterDEControls();
        clusterDeGroupbySelect?.addEventListener('change', () => {{
            clusterDeGroupby = clusterDeGroupbySelect.value || null;
            clusterDeSourceCategory = null;
            clusterDeReferenceCategory = null;
            renderClusterDE();
        }});
        clusterDeSourceSelect?.addEventListener('change', () => {{
            clusterDeSourceCategory = clusterDeSourceSelect.value || null;
            renderClusterDE();
        }});
        clusterDeReferenceSelect?.addEventListener('change', () => {{
            clusterDeReferenceCategory = clusterDeReferenceSelect.value || null;
            renderClusterDE();
        }});
        clusterDeSwap?.addEventListener('click', () => {{
            const source = clusterDeSourceCategory;
            clusterDeSourceCategory = clusterDeReferenceCategory;
            clusterDeReferenceCategory = source;
            renderClusterDE();
        }});

        // Keep initial load fast: only render the list. Other Insights content is computed on-demand
        // when the panel/tab is opened.
        renderColorList('');
        updateExpressionScaleUI();

        const exprVmin = document.getElementById('expr-vmin');
        const exprVmax = document.getElementById('expr-vmax');
        const exprAuto = document.getElementById('expr-auto');
        if (exprVmin && exprVmax) {{
            const applyExpressionScale = () => {{
                if (!currentGene) return;
                const vmin = parseFloat(exprVmin.value);
                const vmax = parseFloat(exprVmax.value);
                if (!Number.isFinite(vmin) || !Number.isFinite(vmax)) return;
                let adjMin = vmin;
                let adjMax = vmax;
                if (adjMin === adjMax) {{
                    adjMin -= 1e-6;
                    adjMax += 1e-6;
                }} else if (adjMin > adjMax) {{
                    const tmp = adjMin;
                    adjMin = adjMax;
                    adjMax = tmp;
                }}
                geneScaleOverrides[currentGene] = {{ vmin: adjMin, vmax: adjMax }};
                rerenderForExpressionScaleChange();
            }};
            exprVmin.addEventListener('change', applyExpressionScale);
            exprVmax.addEventListener('change', applyExpressionScale);
        }}
        if (exprAuto) {{
            exprAuto.addEventListener('click', () => {{
                if (!currentGene) return;
                const autoScale = computeGenePercentiles(currentGene);
                if (autoScale) {{
                    geneScaleAuto[currentGene] = autoScale;
                    delete geneScaleOverrides[currentGene];
                    rerenderForExpressionScaleChange();
                }}
            }});
        }}
    }}

    function renderColorList(query) {{
        const list = document.getElementById('color-list');
        if (!list) return;
        const q = (query || '').trim().toLowerCase();
        const items = DATA.available_colors.filter(col => col.toLowerCase().includes(q));
        if (items.length === 0) {{
            list.innerHTML = `<div class="agg-group-meta">No matches.</div>`;
            return;
        }}
        list.innerHTML = items.map(col => `
            <div class="color-item ${{col === currentColor && !currentGene ? 'active' : ''}}" data-color="${{col}}">
                ${{col}}
            </div>
        `).join('');

        list.querySelectorAll('.color-item').forEach(item => {{
            item.addEventListener('click', () => {{
                const col = item.dataset.color;
                if (!col) return;
                currentColor = col;
                currentGene = null;
                modalSelectedCategory = null;
                modalTypeSelectEnabled = false;
                document.getElementById('color-select').value = col;
                document.getElementById('gene-input').value = '';
                hiddenCategories.clear();
                if (DATA.colors_meta?.[col] && !DATA.colors_meta[col].is_continuous) {{
                    selectionSummaryColor = col;
                }}
                updateExpressionScaleUI();
                renderLegend('legend');
                renderLegend('modal-legend');
                renderAllSections();
                if (modalSection) renderModalSection();
                if (umapVisible) renderUMAP();
                updateSelectionInfo();
                renderGeneDiscoveryPanel();
                renderColorList(document.getElementById('color-search').value);
                renderColorAggregation();
                renderCellTypeTrend();
                renderNeighborStats();
                renderInteractionBrowser();
                renderMarkerGenes();
                if (document.getElementById('color-tab-compare')?.classList.contains('active')) {{
                    renderClusterDE();
                }}
            }});
        }});
    }}

    function renderColorAggregation() {{
        const container = document.getElementById('color-aggregation');
        const groupBy = document.getElementById('color-groupby');
        if (!container || !groupBy) return;
        const groupKey = groupBy.value;
        if (!groupKey) {{
            container.innerHTML = '<div class="agg-group-meta">Pick a metadata column to summarize.</div>';
            return;
        }}

        if (currentGene) {{
            container.innerHTML = '<div class="agg-group-meta">Aggregation is disabled while a gene is active. Clear the gene input to aggregate by categorical colors.</div>';
            return;
        }}

        const config = getColorConfig();
        const label = currentGene || currentColor;
        const groups = new Map();

        DATA.sections.forEach(section => {{
            const groupVal = section.metadata?.[groupKey] || 'unknown';
            if (!groups.has(groupVal)) {{
                groups.set(groupVal, {{ total: 0, counts: {{}}, sum: 0, min: null, max: null }});
            }}
            const group = groups.get(groupVal);
            const values = getSectionValues(section);
            values.forEach(val => {{
                if (val === null || val === undefined || Number.isNaN(val)) return;
                group.total += 1;
                if (config.is_continuous) {{
                    group.sum += val;
                    if (group.min === null || val < group.min) group.min = val;
                    if (group.max === null || val > group.max) group.max = val;
                }} else {{
                    const catIdx = Math.round(val);
                    const catName = config.categories?.[catIdx] || 'unknown';
                    group.counts[catName] = (group.counts[catName] || 0) + 1;
                }}
            }});
        }});

        if (groups.size === 0) {{
            container.innerHTML = '<div class="agg-group-meta">No data to summarize.</div>';
            return;
        }}

        const entries = Array.from(groups.entries());
        entries.sort((a, b) => a[0].localeCompare(b[0]));

        if (config.is_continuous) {{
            container.innerHTML = entries.map(([groupVal, stats]) => {{
                const mean = stats.total > 0 ? (stats.sum / stats.total) : 0;
                const min = stats.min !== null ? stats.min.toFixed(2) : 'n/a';
                const max = stats.max !== null ? stats.max.toFixed(2) : 'n/a';
                return `
                    <div class="agg-group">
                        <div class="agg-group-title">${{formatMetadataLabel(groupKey)}}: ${{groupVal}}</div>
                        <div class="agg-group-meta">${{label}} · n=${{stats.total}}</div>
                        <div class="agg-row"><span class="agg-label">Mean</span><span class="agg-value">${{mean.toFixed(2)}}</span></div>
                        <div class="agg-row"><span class="agg-label">Min</span><span class="agg-value">${{min}}</span></div>
                        <div class="agg-row"><span class="agg-label">Max</span><span class="agg-value">${{max}}</span></div>
                    </div>
                `;
            }}).join('');
            return;
        }}

        container.innerHTML = entries.map(([groupVal, stats]) => {{
            const total = stats.total || 0;
            const counts = Object.entries(stats.counts);
            counts.sort((a, b) => b[1] - a[1]);
            const isExpanded = expandedAggGroups.has(groupVal);
            const top = isExpanded ? counts : counts.slice(0, 6);
            const shownTotal = top.reduce((sum, [, c]) => sum + c, 0);
            const other = total - shownTotal;
            const toggleLabel = isExpanded ? 'Show top 6' : 'Show all';

            const rows = top.map(([cat, count]) => {{
                const pct = total > 0 ? Math.round((count / total) * 100) : 0;
                const catIdx = config.categories?.indexOf(cat) ?? -1;
                const color = catIdx >= 0 ? getCategoryColor(catIdx) : '#999';
                return `
                    <div class="agg-row">
                        <span class="agg-dot" style="background: ${{color}}"></span>
                        <span class="agg-label">${{cat}}</span>
                        <span class="agg-value">${{pct}}% (${{count}})</span>
                    </div>
                `;
            }}).join('');

            const otherRow = other > 0 ? `
                <div class="agg-row">
                    <span class="agg-dot" style="background: #bbb"></span>
                    <span class="agg-label">Other</span>
                    <span class="agg-value">${{Math.round((other / total) * 100)}}% (${{other}})</span>
                </div>
            ` : '';

            return `
                <div class="agg-group">
                    <div class="agg-group-title">${{formatMetadataLabel(groupKey)}}: ${{groupVal}}</div>
                    <div class="agg-group-meta">${{label}} · n=${{total}}</div>
                    <button class="legend-btn" data-agg-toggle="${{groupVal}}">${{toggleLabel}}</button>
                    ${{rows}}
                    ${{otherRow}}
                </div>
            `;
        }}).join('');

        container.querySelectorAll('[data-agg-toggle]').forEach(btn => {{
            btn.addEventListener('click', () => {{
                const key = btn.getAttribute('data-agg-toggle');
                if (!key) return;
                if (expandedAggGroups.has(key)) expandedAggGroups.delete(key);
                else expandedAggGroups.add(key);
                renderColorAggregation();
            }});
        }});
    }}

    function renderMarkerGenes() {{
        const container = document.getElementById('marker-genes');
        if (!container) return;
        const searchInput = document.getElementById('marker-gene-search');
        const query = (searchInput?.value || '').trim().toLowerCase();
        const markers = DATA.marker_genes || {{}};

        const colorMeta = DATA.colors_meta?.[currentColor];
        if (!colorMeta || colorMeta.is_continuous) {{
            container.innerHTML = '<div class="agg-group-meta">Marker genes are available for categorical colors only.</div>';
            return;
        }}

        const groupMarkers = markers[currentColor];
        if (!groupMarkers || Object.keys(groupMarkers).length === 0) {{
            const availableColors = getAvailableMarkerGeneColors();
            const extra = availableColors.length
                ? ` Available for: ${{availableColors.join(', ')}}.`
                : ' No marker genes were embedded in this viewer.';
            container.innerHTML = `<div class="agg-group-meta">No marker genes available for this color.${{escapeHtml(extra)}}</div>`;
            return;
        }}

        const categories = colorMeta.categories || Object.keys(groupMarkers);
        const rows = categories.map(cat => {{
            const key = String(cat);
            const genes = getMarkerGenesForColorCategory(currentColor, key)
                .map((gene) => {{
                    const raw = String(gene || '').trim();
                    if (!raw) return null;
                    return {{
                        raw,
                        canonical: resolveCanonicalGeneName(raw),
                    }};
                }})
                .filter(Boolean);
            if (query) {{
                const hasMatch = genes.some((entry) => {{
                    const rawMatch = entry.raw.toLowerCase().includes(query);
                    const canonicalMatch = entry.canonical
                        ? entry.canonical.toLowerCase().includes(query)
                        : false;
                    return rawMatch || canonicalMatch;
                }});
                if (!hasMatch) return '';
            }}
            const geneButtons = genes.length
                ? genes.map((entry) => renderGeneTokenButton(entry.raw, {{
                    allowUnknown: true,
                    isActive: !!entry.canonical && entry.canonical === currentGene,
                    showMeta: false,
                    title: entry.canonical
                        ? 'Load marker gene into the viewer'
                        : 'Marker gene name only; this gene was not embedded in the viewer',
                }})).join('')
                : '<div class="agg-group-meta">No marker genes found.</div>';
            return `
                <div class="marker-group">
                    <div class="marker-group-title">${{key}}</div>
                    <div class="gene-token-grid">${{geneButtons}}</div>
                </div>
            `;
        }}).filter(Boolean);

        if (rows.length === 0) {{
            container.innerHTML = '<div class="agg-group-meta">No marker genes match your search.</div>';
            return;
        }}

        const loadedGeneNote = currentGene
            ? `<div class="agg-group-meta">Loaded gene: <strong>${{escapeHtml(currentGene)}}</strong>. Click another marker gene to switch.</div>`
            : '<div class="agg-group-meta">Click a marker gene to load it in the viewer. Genes not embedded in this viewer are shown as names only.</div>';

        container.innerHTML = loadedGeneNote + rows.join('');
        bindGeneActivateButtons(container, renderMarkerGenes);
    }}

    function formatClusterDEPct(value) {{
        if (!Number.isFinite(value)) return 'n/a';
        return `${{(100 * value).toFixed(1)}}%`;
    }}

    function renderComparisonGeneTokenGrid(entries, emptyMessage, activateTitle = 'Load gene into the viewer') {{
        if (!entries || !entries.length) {{
            return `<div class="agg-group-meta">${{escapeHtml(emptyMessage)}}</div>`;
        }}
        return `
            <div class="gene-token-grid">
                ${{entries.map((entry) => renderGeneTokenButton(entry.raw, {{
                    allowUnknown: true,
                    isActive: !!entry.canonical && entry.canonical === currentGene,
                    showMeta: false,
                    title: entry.canonical
                        ? activateTitle
                        : 'Gene name only; this gene was not embedded in the viewer',
                }})).join('')}}
            </div>
        `;
    }}

    function renderComparisonCountSummary(colorCol, sourceCategory, referenceCategory) {{
        const countData = getCategoryCountsForColor(colorCol);
        if (!countData) {{
            return `
                <div class="agg-group">
                    <div class="agg-group-title">Cell Counts</div>
                    <div class="agg-group-meta">No categorical counts are available for this annotation.</div>
                </div>
            `;
        }}

        const total = Number(countData.total || 0);
        const sourceCount = Number(countData.counts?.[sourceCategory] || 0);
        const referenceCount = Number(countData.counts?.[referenceCategory] || 0);
        const sourcePct = total > 0 ? (sourceCount / total) * 100 : 0;
        const referencePct = total > 0 ? (referenceCount / total) * 100 : 0;
        const sourceColor = getCategoryColorForValue(colorCol, sourceCategory);
        const referenceColor = getCategoryColorForValue(colorCol, referenceCategory);

        return `
            <div class="agg-group">
                <div class="agg-group-title">Cell Counts</div>
                <div class="agg-group-meta">Cells assigned in ${{formatMetadataLabel(colorCol)}}.</div>
                <div class="comparison-stack">
                    <div class="comparison-card">
                        <div class="comparison-card-title">
                            <span class="agg-dot" style="background: ${{sourceColor}}"></span>
                            <span>${{escapeHtml(sourceCategory)}}</span>
                        </div>
                        <div class="comparison-metric-grid">
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">Cells</span>
                                <span class="comparison-metric-value">${{sourceCount.toLocaleString()}}</span>
                            </div>
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">Share</span>
                                <span class="comparison-metric-value">${{sourcePct.toFixed(1)}}%</span>
                            </div>
                        </div>
                    </div>
                    <div class="comparison-card">
                        <div class="comparison-card-title">
                            <span class="agg-dot" style="background: ${{referenceColor}}"></span>
                            <span>${{escapeHtml(referenceCategory)}}</span>
                        </div>
                        <div class="comparison-metric-grid">
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">Cells</span>
                                <span class="comparison-metric-value">${{referenceCount.toLocaleString()}}</span>
                            </div>
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">Share</span>
                                <span class="comparison-metric-value">${{referencePct.toFixed(1)}}%</span>
                            </div>
                        </div>
                    </div>
                </div>
                <div class="agg-group-meta">Total Labeled Cells: ${{total.toLocaleString()}}</div>
            </div>
        `;
    }}

    function renderComparisonNeighborSummary(colorCol, sourceCategory, referenceCategory) {{
        if (!DATA.has_neighbors) {{
            return `
                <div class="agg-group">
                    <div class="agg-group-title">Neighbor Relationship</div>
                    <div class="agg-group-meta">No neighbor graph was found in this dataset.</div>
                </div>
            `;
        }}

        const sourceToReference = getNeighborPairSummary(colorCol, sourceCategory, referenceCategory);
        const referenceToSource = getNeighborPairSummary(colorCol, referenceCategory, sourceCategory);
        if (!sourceToReference && !referenceToSource) {{
            return `
                <div class="agg-group">
                    <div class="agg-group-title">Neighbor Relationship</div>
                    <div class="agg-group-meta">Neighbor stats were not precomputed for this annotation.</div>
                </div>
            `;
        }}

        const cards = [sourceToReference, referenceToSource]
            .filter(Boolean)
            .map((entry) => {{
                const targetColor = getCategoryColorForValue(colorCol, entry.targetCategory);
                const zLabel = entry.z === null ? 'n/a' : entry.z.toFixed(2);
                const degreeLabel = Number.isFinite(entry.meanDegree) ? entry.meanDegree.toFixed(2) : 'n/a';
                return `
                    <div class="comparison-card">
                        <div class="comparison-card-title">
                            <span class="agg-dot" style="background: ${{targetColor}}"></span>
                            <span>${{escapeHtml(entry.sourceCategory)}} → ${{escapeHtml(entry.targetCategory)}}</span>
                        </div>
                        <div class="comparison-metric-grid">
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">Share</span>
                                <span class="comparison-metric-value">${{entry.pct.toFixed(1)}}%</span>
                            </div>
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">Edges</span>
                                <span class="comparison-metric-value">${{formatNeighborCount(entry.count)}}</span>
                            </div>
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">Z</span>
                                <span class="comparison-metric-value">${{zLabel}}</span>
                            </div>
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">Mean degree</span>
                                <span class="comparison-metric-value">${{degreeLabel}}</span>
                            </div>
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">n cells</span>
                                <span class="comparison-metric-value">${{Number.isFinite(entry.nCells) ? Math.round(entry.nCells).toLocaleString() : 'n/a'}}</span>
                            </div>
                            <div class="comparison-metric">
                                <span class="comparison-metric-label">Permutations</span>
                                <span class="comparison-metric-value">${{entry.permN > 0 ? entry.permN.toLocaleString() : '0'}}</span>
                            </div>
                        </div>
                    </div>
                `;
            }})
            .join('');

        return `
            <div class="agg-group">
                <div class="agg-group-title">Neighbor Relationship</div>
                <div class="agg-group-meta">Directional neighbor enrichment between the selected categories.</div>
                <div class="comparison-stack">${{cards}}</div>
            </div>
        `;
    }}

    function renderComparisonMarkerSummary(colorCol, sourceCategory, referenceCategory) {{
        const sourceMarkers = getMarkerGeneEntries(colorCol, sourceCategory, 6);
        const referenceMarkers = getMarkerGeneEntries(colorCol, referenceCategory, 6);
        const overlap = getMarkerOverlapEntries(colorCol, sourceCategory, referenceCategory, 6);

        if (!sourceMarkers.length && !referenceMarkers.length) {{
            return `
                <div class="agg-group">
                    <div class="agg-group-title">Marker Summary</div>
                    <div class="agg-group-meta">Marker genes were not precomputed for this annotation.</div>
                </div>
            `;
        }}

        return `
            <div class="agg-group">
                <div class="agg-group-title">Marker Summary</div>
                <div class="agg-group-meta">Top markers for each category. Click a gene to load it when available.</div>
            </div>
            <div class="agg-group">
                <div class="agg-group-title">${{escapeHtml(sourceCategory)}} markers</div>
                ${{renderComparisonGeneTokenGrid(sourceMarkers, 'No marker genes available for Category A.', 'Load Category A marker gene into the viewer')}}
            </div>
            <div class="agg-group">
                <div class="agg-group-title">${{escapeHtml(referenceCategory)}} markers</div>
                ${{renderComparisonGeneTokenGrid(referenceMarkers, 'No marker genes available for Category B.', 'Load Category B marker gene into the viewer')}}
            </div>
            <div class="agg-group">
                <div class="agg-group-title">Marker Overlap</div>
                <div class="agg-group-meta">Shared names across the top marker lists.</div>
                ${{renderComparisonGeneTokenGrid(overlap, 'No shared top markers found.', 'Load shared marker gene into the viewer')}}
            </div>
        `;
    }}

    function renderInteractionComparisonCard(colorCol, sourceCategory, targetCategory) {{
        const summary = getInteractionPairSummary(colorCol, sourceCategory, targetCategory);
        const title = `${{escapeHtml(sourceCategory)}} → ${{escapeHtml(targetCategory)}}`;
        if (!summary) {{
            return `
                <div class="agg-group">
                    <div class="agg-group-title">${{title}}</div>
                    <div class="agg-group-meta">No contact-conditioned marker result for this direction.</div>
                </div>
            `;
        }}

        const result = summary.result || {{}};
        const nContact = Number(result.n_contact ?? 0).toLocaleString();
        const nNonContact = Number(result.n_non_contact ?? 0).toLocaleString();
        const available = result.available !== false;
        const meta = available
            ? `n+ ${{nContact}} | n- ${{nNonContact}}`
            : `Unavailable: need >= ${{Number(result.min_cells_required ?? 0).toLocaleString()}} cells per group (got n+ ${{nContact}} | n- ${{nNonContact}})`;

        return `
            <div class="agg-group">
                <div class="agg-group-title">${{title}}</div>
                <div class="agg-group-meta">${{meta}}</div>
                ${{renderComparisonGeneTokenGrid(summary.genes, available ? 'No contact-conditioned genes returned.' : 'Contact-conditioned markers unavailable for this direction.', 'Load contact-conditioned marker gene into the viewer')}}
            </div>
        `;
    }}

    function renderComparisonInteractionSummary(colorCol, sourceCategory, referenceCategory) {{
        const colorData = (DATA.interaction_markers || {{}})[colorCol] || {{}};
        if (!Object.keys(colorData).length) {{
            return `
                <div class="agg-group">
                    <div class="agg-group-title">Contact-Conditioned Markers</div>
                    <div class="agg-group-meta">Interaction markers were not precomputed for this annotation.</div>
                </div>
            `;
        }}

        return `
            <div class="agg-group">
                <div class="agg-group-title">Contact-Conditioned Markers</div>
                <div class="agg-group-meta">Directional DE between contact-positive and contact-negative source cells.</div>
            </div>
            ${{renderInteractionComparisonCard(colorCol, sourceCategory, referenceCategory)}}
            ${{renderInteractionComparisonCard(colorCol, referenceCategory, sourceCategory)}}
        `;
    }}

    function renderClusterDEResultSection(colorCol, sourceCategory, referenceCategory) {{
        const result = getPairwiseClusterDEResult(colorCol, sourceCategory, referenceCategory);
        if (!result) {{
            return `
                <div class="agg-group">
                    <div class="agg-group-title">Pairwise DE</div>
                    <div class="agg-group-meta">No pairwise DE result is available for this comparison.</div>
                </div>
            `;
        }}

        if (result.available === false) {{
            const nSource = Number(result.n_source ?? 0).toLocaleString();
            const nReference = Number(result.n_reference ?? 0).toLocaleString();
            const minCells = Number(result.min_cells_required ?? 0).toLocaleString();
            let detail = `n=${{nSource}} vs n=${{nReference}}.`;
            if (result.reason === 'insufficient_cells') {{
                detail = `Need >= ${{minCells}} cells per category (got ${{nSource}} vs ${{nReference}}).`;
            }}
            return `
                <div class="agg-group">
                    <div class="agg-group-title">Pairwise DE</div>
                    <div class="agg-group-meta">${{escapeHtml(detail)}}</div>
                </div>
            `;
        }}

        const genes = Array.isArray(result.genes) ? result.genes : [];
        const logfc = Array.isArray(result.logfoldchanges) ? result.logfoldchanges : [];
        const pvalsAdj = Array.isArray(result.pvals_adj) ? result.pvals_adj : [];
        const scores = Array.isArray(result.scores) ? result.scores : [];
        const pctSource = Array.isArray(result.pct_source) ? result.pct_source : [];
        const pctReference = Array.isArray(result.pct_reference) ? result.pct_reference : [];
        if (!genes.length) {{
            return `
                <div class="agg-group">
                    <div class="agg-group-title">Pairwise DE</div>
                    <div class="agg-group-meta">No DE genes were returned for this comparison.</div>
                </div>
            `;
        }}

        const rows = genes.map((gene, idx) => {{
            const logfcValue = logfc[idx];
            const pvalAdjValue = pvalsAdj[idx];
            const scoreValue = scores[idx];
            const pctSourceValue = pctSource[idx];
            const pctReferenceValue = pctReference[idx];
            return `
                <tr>
                    <td><button type="button" class="cluster-de-gene-btn" data-cluster-de-gene="${{escapeHtml(gene)}}">${{escapeHtml(gene)}}</button></td>
                    <td>${{formatScaleNumber(Number.isFinite(logfcValue) ? logfcValue : NaN)}}</td>
                    <td>${{formatScaleNumber(Number.isFinite(pvalAdjValue) ? pvalAdjValue : NaN)}}</td>
                    <td>${{formatScaleNumber(Number.isFinite(scoreValue) ? scoreValue : NaN)}}</td>
                    <td>${{formatClusterDEPct(Number.isFinite(pctSourceValue) ? pctSourceValue : NaN)}}</td>
                    <td>${{formatClusterDEPct(Number.isFinite(pctReferenceValue) ? pctReferenceValue : NaN)}}</td>
                </tr>
            `;
        }}).join('');

        return `
            <div class="cluster-de-result">
                <div class="cluster-de-meta">
                    Pairwise DE for <strong>${{escapeHtml(sourceCategory)}}</strong> vs <strong>${{escapeHtml(referenceCategory)}}</strong>.
                    Click a gene to load it in the viewer.
                </div>
                <table class="cluster-de-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>logFC</th>
                            <th>adj. p</th>
                            <th>Score</th>
                            <th>% A</th>
                            <th>% B</th>
                        </tr>
                    </thead>
                    <tbody>${{rows}}</tbody>
                </table>
            </div>
        `;
    }}

    function syncClusterDEControls() {{
        const groupbySelect = document.getElementById('cluster-de-groupby');
        const sourceSelect = document.getElementById('cluster-de-source');
        const referenceSelect = document.getElementById('cluster-de-reference');
        const availableGroupbys = getAvailableComparisonColors();

        if (groupbySelect) {{
            groupbySelect.innerHTML = availableGroupbys.map(col => `<option value="${{col}}">${{col}}</option>`).join('');
        }}

        if (!availableGroupbys.length) {{
            clusterDeGroupby = null;
            clusterDeSourceCategory = null;
            clusterDeReferenceCategory = null;
            if (sourceSelect) sourceSelect.innerHTML = '';
            if (referenceSelect) referenceSelect.innerHTML = '';
            return {{ availableGroupbys, categories: [] }};
        }}

        if (!clusterDeGroupby || !availableGroupbys.includes(clusterDeGroupby)) {{
            clusterDeGroupby = availableGroupbys.includes(currentColor) ? currentColor : availableGroupbys[0];
        }}
        if (groupbySelect) groupbySelect.value = clusterDeGroupby;

        const categories = getClusterDECategories(clusterDeGroupby);
        const options = categories.map(category => `<option value="${{escapeHtml(category)}}">${{escapeHtml(category)}}</option>`).join('');
        if (sourceSelect) sourceSelect.innerHTML = options;
        if (referenceSelect) referenceSelect.innerHTML = options;

        if (!categories.length) {{
            clusterDeSourceCategory = null;
            clusterDeReferenceCategory = null;
            return {{ availableGroupbys, categories }};
        }}

        if (!clusterDeSourceCategory || !categories.includes(clusterDeSourceCategory)) {{
            clusterDeSourceCategory = categories[0];
        }}
        if (!clusterDeReferenceCategory || !categories.includes(clusterDeReferenceCategory) || clusterDeReferenceCategory === clusterDeSourceCategory) {{
            clusterDeReferenceCategory = categories.find(category => category !== clusterDeSourceCategory) || null;
        }}

        if (sourceSelect) sourceSelect.value = clusterDeSourceCategory || '';
        if (referenceSelect) {{
            referenceSelect.value = clusterDeReferenceCategory || '';
            referenceSelect.disabled = categories.length < 2;
        }}

        return {{ availableGroupbys, categories }};
    }}

    function renderClusterDE() {{
        const container = document.getElementById('cluster-de-results');
        if (!container) return;

        const {{ availableGroupbys, categories }} = syncClusterDEControls();
        if (!availableGroupbys.length) {{
            container.innerHTML = '<div class="agg-group-meta">Category comparison is available for categorical colors only.</div>';
            return;
        }}
        if (!clusterDeGroupby) {{
            container.innerHTML = '<div class="agg-group-meta">Choose a categorical annotation to compare.</div>';
            return;
        }}
        if (!categories.length) {{
            container.innerHTML = '<div class="agg-group-meta">No categories are available for this annotation.</div>';
            return;
        }}
        if (!clusterDeSourceCategory || !clusterDeReferenceCategory) {{
            container.innerHTML = '<div class="agg-group-meta">Choose two different categories to compare.</div>';
            return;
        }}
        if (clusterDeSourceCategory === clusterDeReferenceCategory) {{
            container.innerHTML = '<div class="agg-group-meta">Choose two different categories to compare.</div>';
            return;
        }}

        container.innerHTML = `
            <div class="cluster-de-meta">
                Comparing <strong>${{escapeHtml(clusterDeSourceCategory)}}</strong> vs <strong>${{escapeHtml(clusterDeReferenceCategory)}}</strong>
                in <strong>${{escapeHtml(formatMetadataLabel(clusterDeGroupby))}}</strong>.
            </div>
            ${{renderComparisonCountSummary(clusterDeGroupby, clusterDeSourceCategory, clusterDeReferenceCategory)}}
            ${{renderComparisonNeighborSummary(clusterDeGroupby, clusterDeSourceCategory, clusterDeReferenceCategory)}}
            ${{renderComparisonMarkerSummary(clusterDeGroupby, clusterDeSourceCategory, clusterDeReferenceCategory)}}
            ${{renderComparisonInteractionSummary(clusterDeGroupby, clusterDeSourceCategory, clusterDeReferenceCategory)}}
            ${{renderClusterDEResultSection(clusterDeGroupby, clusterDeSourceCategory, clusterDeReferenceCategory)}}
        `;

        bindGeneActivateButtons(container, renderClusterDE);
        container.querySelectorAll('[data-cluster-de-gene]').forEach((btn) => {{
            btn.addEventListener('click', async () => {{
                const gene = btn.getAttribute('data-cluster-de-gene') || '';
                if (!gene) return;
                const ok = await activateViewerGene(gene, {{ showErrors: true }});
                if (ok) renderClusterDE();
            }});
        }});
    }}

    function renderCellTypeTrend() {{
        const container = document.getElementById('celltype-trend');
        const groupBy = document.getElementById('color-groupby');
        if (!container || !groupBy) return;
        const groupKey = groupBy.value;
        if (!groupKey) {{
            container.innerHTML = '<div class="agg-group-meta">Pick a metadata column to summarize.</div>';
            return;
        }}

        if (currentGene) {{
            container.innerHTML = '<div class="agg-group-meta">Clear the gene input to view categorical trends.</div>';
            return;
        }}

        const config = getColorConfig();
        if (config.is_continuous) {{
            container.innerHTML = '<div class="agg-group-meta">Trends are available for categorical colors only.</div>';
            return;
        }}

        const input = document.getElementById('celltype-search');
        const query = (input?.value || '').trim().toLowerCase();
        if (!query) {{
            container.innerHTML = '<div class="agg-group-meta">Search for a category to see counts across metadata groups.</div>';
            return;
        }}

        const categories = config.categories || [];
        const matches = categories.filter(cat => String(cat).toLowerCase().includes(query));
        if (matches.length === 0) {{
            container.innerHTML = '<div class="agg-group-meta">No matching categories.</div>';
            return;
        }}

        if (!celltypeTrendTarget || !matches.includes(celltypeTrendTarget)) {{
            celltypeTrendTarget = matches[0];
        }}
        const target = celltypeTrendTarget;
        const groups = new Map();

        DATA.sections.forEach(section => {{
            const groupVal = section.metadata?.[groupKey] || 'unknown';
            if (!groups.has(groupVal)) {{
                groups.set(groupVal, {{ total: 0, count: 0 }});
            }}
            const group = groups.get(groupVal);
            const values = getSectionValues(section);
            values.forEach(val => {{
                if (val === null || val === undefined || Number.isNaN(val)) return;
                group.total += 1;
                const catIdx = Math.round(val);
                const catName = config.categories?.[catIdx];
                if (catName === target) group.count += 1;
            }});
        }});

        const entries = Array.from(groups.entries()).sort((a, b) => String(a[0]).localeCompare(String(b[0])));
        const rows = entries.map(([groupVal, stats]) => {{
            const pct = stats.total > 0 ? (stats.count / stats.total) * 100 : 0;
            return `
                <tr>
                    <td>${{groupVal}}</td>
                    <td>${{stats.count}}</td>
                    <td>${{stats.total}}</td>
                    <td>${{pct.toFixed(1)}}%</td>
                </tr>
            `;
        }}).join('');

        const matchBadges = matches.length > 1
            ? `<div class="agg-group-meta" style="display:flex;flex-wrap:wrap;gap:4px;align-items:center;">
                Matches:
                ${{matches.map(cat => `<button type="button" class="legend-btn${{cat === target ? ' active' : ''}}" data-trend-target="${{escapeHtml(cat)}}">${{escapeHtml(cat)}}</button>`).join('')}}
               </div>`
            : `<div class="agg-group-meta">Showing: <strong>${{escapeHtml(target)}}</strong></div>`;

        container.innerHTML = `
            ${{matchBadges}}
            <table class="trend-table">
                <thead>
                    <tr>
                        <th>${{formatMetadataLabel(groupKey)}}</th>
                        <th>Count</th>
                        <th>Total</th>
                        <th>%</th>
                    </tr>
                </thead>
                <tbody>
                    ${{rows}}
                </tbody>
            </table>
        `;
        container.querySelectorAll('[data-trend-target]').forEach(btn => {{
            btn.addEventListener('click', () => {{
                celltypeTrendTarget = btn.getAttribute('data-trend-target');
                renderCellTypeTrend();
            }});
        }});
    }}

    function renderNeighborStats() {{
        const container = document.getElementById('neighbor-stats');
        if (!container) return;

        if (!DATA.has_neighbors) {{
            container.innerHTML = '<div class="agg-group-meta">No neighbor graph was found in this dataset.</div>';
            return;
        }}

        if (currentGene) {{
            container.innerHTML = '<div class="agg-group-meta">Clear the gene input to view neighbor stats.</div>';
            return;
        }}

        const config = getColorConfig();
        if (config.is_continuous) {{
            container.innerHTML = '<div class="agg-group-meta">Neighbor stats are available for categorical colors only.</div>';
            return;
        }}

        const stats = (DATA.neighbor_stats || {{}})[currentColor];
        if (!stats || !stats.counts || !stats.categories) {{
            container.innerHTML = '<div class="agg-group-meta">No neighbor stats available for this color.</div>';
            return;
        }}

        const categories = stats.categories || [];
        const counts = stats.counts || [];
        const nCells = stats.n_cells || [];
        const meanDegree = stats.mean_degree || [];
        const zscores = stats.zscore || null;
        const permN = stats.perm_n || 0;
        if (categories.length === 0 || counts.length === 0) {{
            container.innerHTML = '<div class="agg-group-meta">Neighbor stats are empty for this color.</div>';
            return;
        }}

        const input = document.getElementById('neighbor-search');
        const query = (input?.value || '').trim().toLowerCase();
        const matches = categories
            .map((cat, idx) => (query && !String(cat).toLowerCase().includes(query)) ? null : idx)
            .filter(idx => idx !== null);

        if (matches.length === 0) {{
            container.innerHTML = '<div class="agg-group-meta">No matching cell types.</div>';
            return;
        }}

        const rows = matches.map((idx) => {{
            const source = String(categories[idx]);
            const row = counts[idx] || [];
            const total = row.reduce((sum, val) => sum + (Number.isFinite(val) ? val : 0), 0);
            const entries = row
                .map((val, j) => (Number.isFinite(val) && val > 0) ? [j, val] : null)
                .filter(Boolean)
                .sort((a, b) => b[1] - a[1]);
            const key = `${{currentColor}}::${{idx}}`;
            const isExpanded = expandedNeighborGroups.has(key);
            const top = isExpanded ? entries : entries.slice(0, 6);
            const shownTotal = top.reduce((sum, [, v]) => sum + v, 0);
            const other = total - shownTotal;
            const toggleLabel = isExpanded ? 'Show top 6' : 'Show all';
            const formatCount = (value) => {{
                if (!Number.isFinite(value)) return '0';
                if (Math.abs(value - Math.round(value)) < 1e-6) return Math.round(value).toLocaleString();
                return value.toFixed(2);
            }};

            const rowsHtml = top.map(([j, val]) => {{
                const pct = total > 0 ? (val / total) * 100 : 0;
                const target = String(categories[j] ?? 'unknown');
                const color = j >= 0 ? getCategoryColor(j) : '#999';
                let zLabel = '';
                if (zscores && zscores[idx] && Number.isFinite(zscores[idx][j])) {{
                    zLabel = ` z=${{zscores[idx][j].toFixed(2)}}`;
                }}
                return `
                    <div class="agg-row">
                        <span class="agg-dot" style="background: ${{color}}"></span>
                        <span class="agg-label">${{target}}</span>
                        <span class="agg-value">${{pct.toFixed(1)}}% (${{formatCount(val)}})${{zLabel}}</span>
                    </div>
                `;
            }}).join('');

            const otherRow = other > 0 ? `
                <div class="agg-row">
                    <span class="agg-dot" style="background: #bbb"></span>
                    <span class="agg-label">Other</span>
                    <span class="agg-value">${{((other / total) * 100).toFixed(1)}}% (${{formatCount(other)}})</span>
                </div>
            ` : '';

            const totalLabel = formatCount(total);
            const nLabel = (nCells[idx] ?? 0).toLocaleString();
            const degreeLabel = Number.isFinite(meanDegree[idx]) ? meanDegree[idx].toFixed(2) : '0.00';
            const permLabel = permN ? ` | perms=${{permN}}` : '';

            if (entries.length === 0) {{
                return `
                    <div class="agg-group">
                        <div class="agg-group-title">${{source}}</div>
                        <div class="agg-group-meta">n=${{nLabel}} | mean degree=${{degreeLabel}}</div>
                        <div class="agg-group-meta">No neighbors found for this cell type.</div>
                    </div>
                `;
            }}

            return `
                <div class="agg-group">
                    <div class="agg-group-title">${{source}}</div>
                    <div class="agg-group-meta">n=${{nLabel}} | mean degree=${{degreeLabel}} | neighbor edges=${{totalLabel}}${{permLabel}}</div>
                    <button class="legend-btn" data-neighbor-toggle="${{idx}}">${{toggleLabel}}</button>
                    ${{rowsHtml}}
                    ${{otherRow}}
                </div>
            `;
        }}).join('');

        container.innerHTML = rows;

        container.querySelectorAll('[data-neighbor-toggle]').forEach(btn => {{
            btn.addEventListener('click', () => {{
                const idx = btn.getAttribute('data-neighbor-toggle');
                if (idx === null) return;
                const key = `${{currentColor}}::${{idx}}`;
                if (expandedNeighborGroups.has(key)) expandedNeighborGroups.delete(key);
                else expandedNeighborGroups.add(key);
                renderNeighborStats();
            }});
        }});
    }}

    function renderInteractionBrowser() {{
        const container = document.getElementById('interaction-browser');
        const sourceSelect = document.getElementById('interaction-source');
        if (!container || !sourceSelect) return;

        if (!DATA.has_neighbors) {{
            container.innerHTML = '<div class="agg-group-meta">No neighbor graph was found in this dataset.</div>';
            sourceSelect.innerHTML = '';
            return;
        }}
        if (currentGene) {{
            container.innerHTML = '<div class="agg-group-meta">Clear the gene input to browse cell-cell interactions.</div>';
            sourceSelect.innerHTML = '';
            return;
        }}

        const config = getColorConfig();
        if (config.is_continuous) {{
            container.innerHTML = '<div class="agg-group-meta">Interaction browser is available for categorical colors only.</div>';
            sourceSelect.innerHTML = '';
            return;
        }}

        const stats = (DATA.neighbor_stats || {{}})[currentColor];
        if (!stats || !stats.counts || !stats.categories) {{
            container.innerHTML = '<div class="agg-group-meta">No neighbor stats available for this color.</div>';
            sourceSelect.innerHTML = '';
            return;
        }}

        const categories = (stats.categories || []).map(cat => String(cat));
        const counts = stats.counts || [];
        const zscores = stats.zscore || null;
        const nCells = stats.n_cells || [];
        const meanDegree = stats.mean_degree || [];
        const markers = (DATA.marker_genes || {{}})[currentColor] || {{}};
        const interactionMarkersByColor = (DATA.interaction_markers || {{}})[currentColor] || {{}};
        const hasInteractionMarkers = Object.keys(interactionMarkersByColor).length > 0;
        if (categories.length === 0 || counts.length === 0) {{
            container.innerHTML = '<div class="agg-group-meta">Interaction data is empty for this color.</div>';
            sourceSelect.innerHTML = '';
            return;
        }}

        const currentOptions = Array.from(sourceSelect.options).map(opt => opt.value);
        const needsOptionsUpdate = (
            currentOptions.length !== categories.length ||
            currentOptions.some((value, idx) => value !== categories[idx])
        );
        if (needsOptionsUpdate) {{
            sourceSelect.innerHTML = '';
            categories.forEach(cat => {{
                const opt = document.createElement('option');
                opt.value = cat;
                opt.textContent = cat;
                sourceSelect.appendChild(opt);
            }});
        }}

        if (interactionSourceCategory && categories.includes(interactionSourceCategory)) {{
            sourceSelect.value = interactionSourceCategory;
        }} else if (sourceSelect.value && categories.includes(sourceSelect.value)) {{
            interactionSourceCategory = sourceSelect.value;
        }} else {{
            interactionSourceCategory = categories[0];
            sourceSelect.value = interactionSourceCategory;
        }}

        const source = interactionSourceCategory;
        const sourceIdx = categories.indexOf(source);
        if (sourceIdx < 0) {{
            container.innerHTML = '<div class="agg-group-meta">Choose a source cell type.</div>';
            return;
        }}
        const sourceInteractionMarkers = interactionMarkersByColor[source] || {{}};

        const row = counts[sourceIdx] || [];
        const total = row.reduce((sum, value) => sum + (Number.isFinite(value) ? value : 0), 0);
        const targetQuery = (document.getElementById('interaction-search')?.value || '').trim().toLowerCase();
        const sourceMarkers = (markers[source] || []).slice(0, 6);
        const sortedEntries = categories
            .map((target, targetIdx) => {{
                const count = Number(row[targetIdx] ?? 0);
                const pct = total > 0 ? (count / total) * 100 : 0;
                const z = (zscores && zscores[sourceIdx] && Number.isFinite(zscores[sourceIdx][targetIdx]))
                    ? Number(zscores[sourceIdx][targetIdx])
                    : null;
                const targetMarkers = (markers[target] || []).slice(0, 4);
                const contact = sourceInteractionMarkers[target] || null;
                const contactMarkers = contact && Array.isArray(contact.genes)
                    ? contact.genes.slice(0, 4).filter(Boolean)
                    : [];
                return {{ target, targetIdx, count, pct, z, targetMarkers, contact, contactMarkers }};
            }})
            .filter(entry => !targetQuery || entry.target.toLowerCase().includes(targetQuery))
            .filter(entry => entry.count > 0 || (entry.z !== null && entry.z > 0))
            .sort((a, b) => {{
                if (a.z !== null && b.z !== null && a.z !== b.z) return b.z - a.z;
                if (a.count !== b.count) return b.count - a.count;
                return a.target.localeCompare(b.target);
            }});

        if (sortedEntries.length === 0) {{
            container.innerHTML = '<div class="agg-group-meta">No target cell types match the current filter.</div>';
            return;
        }}

        const topEntries = sortedEntries.slice(0, 12);
        const sourceMarkerLabel = sourceMarkers.length ? sourceMarkers.join(', ') : 'No marker genes available.';
        const sourceN = (nCells[sourceIdx] ?? 0).toLocaleString();
        const degreeLabel = Number.isFinite(meanDegree[sourceIdx]) ? meanDegree[sourceIdx].toFixed(2) : '0.00';
        const withContactMarkers = topEntries.filter(entry => !!entry.contact).length;
        const rows = topEntries.map(entry => {{
            const color = getCategoryColor(entry.targetIdx);
            const zLabel = entry.z === null ? 'n/a' : entry.z.toFixed(2);
            const markerLabel = entry.targetMarkers.length ? entry.targetMarkers.join(', ') : '—';
            const contactLabel = entry.contactMarkers.length ? entry.contactMarkers.join(', ') : '—';
            let contactMeta = 'not precomputed';
            if (entry.contact) {{
                const nPos = Number(entry.contact.n_contact ?? 0);
                const nNeg = Number(entry.contact.n_non_contact ?? 0);
                if (entry.contact.available === false) {{
                    const minReq = Number(entry.contact.min_cells_required ?? 0);
                    contactMeta = `n+ ${{nPos}} / n- ${{nNeg}} (need >= ${{minReq}} each)`;
                }} else {{
                    contactMeta = `n+ ${{nPos}} / n- ${{nNeg}}`;
                }}
            }}
            return `
                <tr>
                    <td>
                        <div class="interaction-target">
                            <span class="agg-dot" style="background: ${{color}}"></span>
                            <span>${{entry.target}}</span>
                        </div>
                    </td>
                    <td>${{entry.pct.toFixed(1)}}%</td>
                    <td>${{formatNeighborCount(entry.count)}}</td>
                    <td>${{zLabel}}</td>
                    <td class="interaction-markers">${{contactLabel}}<br><span style="opacity:0.75;">${{contactMeta}}</span></td>
                    <td class="interaction-markers">${{markerLabel}}</td>
                </tr>
            `;
        }}).join('');

        container.innerHTML = `
            <div class="agg-group">
                <div class="agg-group-title">${{source}} → targets</div>
                <div class="agg-group-meta">n=${{sourceN}} | mean degree=${{degreeLabel}} | neighbor edges=${{formatNeighborCount(total)}}</div>
                <div class="agg-group-meta">Source markers: ${{sourceMarkerLabel}}</div>
                <div class="agg-group-meta">Contact-conditioned markers available for ${{withContactMarkers}}/${{topEntries.length}} shown targets.</div>
                ${{hasInteractionMarkers ? '' : '<div class="agg-group-meta">Contact markers not precomputed for this color (set interaction_markers_groupby during export).</div>'}}
            </div>
            <table class="trend-table">
                <thead>
                    <tr>
                        <th>Target</th>
                        <th>Share</th>
                        <th>Edges</th>
                        <th>Z</th>
                        <th>Contact markers</th>
                        <th>Type markers</th>
                    </tr>
                </thead>
                <tbody>
                    ${{rows}}
                </tbody>
            </table>
        `;
    }}

    function stepRange(rangeEl, delta) {{
        if (!rangeEl) return;
        const min = parseFloat(rangeEl.min || '0');
        const max = parseFloat(rangeEl.max || '100');
        const step = parseFloat(rangeEl.step || '1');
        const current = parseFloat(rangeEl.value || '0');
        const next = Math.min(max, Math.max(min, current + delta * step));
        rangeEl.value = String(next);
        rangeEl.dispatchEvent(new Event('input', {{ bubbles: true }}));
    }}

    // Filters
    function initFilters() {{
        const filterBar = document.getElementById('filter-bar');
        const filters = DATA.metadata_filters || {{}};

        if (Object.keys(filters).length === 0) {{
            filterBar.style.display = 'none';
            return;
        }}

        let html = '';
        for (const [key, values] of Object.entries(filters)) {{
            activeFilters[key] = new Set();  // Start with all shown
            html += `<div class="filter-group">
                <label>${{formatMetadataLabel(key)}}:</label>
                <div class="filter-chips" data-filter="${{key}}">
                    ${{values.map(v => `<span class="filter-chip" data-value="${{v}}">${{v}}</span>`).join('')}}
                </div>
            </div>`;
        }}

        // Add outline legend if outline metadata exists
        if (OUTLINE_BY && filters[OUTLINE_BY] && filters[OUTLINE_BY].length > 0) {{
            const outlineLabel = formatMetadataLabel(OUTLINE_BY);
            html += `<div class="filter-group" style="margin-left: auto;">
                <label>Outline (${{
                    outlineLabel
                }}):</label>
                <div style="display: flex; gap: 8px; align-items: center;">
                    ${{filters[OUTLINE_BY].map(v => `<span style="display: flex; align-items: center; gap: 3px; font-size: 10px;">
                        <span style="width: 12px; height: 12px; border: 3px solid ${{getOutlineColor(v)}}; border-radius: 2px;"></span>
                        ${{v}}
                    </span>`).join('')}}
                </div>
            </div>`;
        }}

        html += `<div class="filter-group">
            <button class="filter-reset-btn" id="filter-reset-btn" type="button" title="Clear all active metadata filters" data-help="Clear all metadata-based sample filters and show every section again.">Reset filters</button>
        </div>`;

        filterBar.innerHTML = html;

        const updateFilterResetState = () => {{
            const btn = document.getElementById('filter-reset-btn');
            if (!btn) return;
            const hasActive = Object.values(activeFilters).some(set => set && set.size > 0);
            btn.disabled = !hasActive;
        }};

        const resetAllFilters = () => {{
            Object.values(activeFilters).forEach(set => {{
                if (set && typeof set.clear === 'function') set.clear();
            }});
            filterBar.querySelectorAll('.filter-chip').forEach(chip => {{
                chip.classList.remove('active');
                chip.classList.remove('inactive');
            }});
            updateFilterResetState();
            renderAllSections();
        }};

        document.getElementById('filter-reset-btn')?.addEventListener('click', resetAllFilters);

        filterBar.querySelectorAll('.filter-chip').forEach(chip => {{
            chip.addEventListener('click', () => {{
                const filterKey = chip.parentElement.dataset.filter;
                const value = chip.dataset.value;
                const filterSet = activeFilters[filterKey];

                if (chip.classList.contains('active')) {{
                    chip.classList.remove('active');
                    filterSet.delete(value);
                }} else {{
                    chip.classList.add('active');
                    filterSet.add(value);
                }}

                // Update chip states
                const chips = chip.parentElement.querySelectorAll('.filter-chip');
                const anyActive = filterSet.size > 0;
                chips.forEach(c => {{
                    c.classList.toggle('inactive', anyActive && !c.classList.contains('active'));
                }});

                updateFilterResetState();
                renderAllSections();
            }});
        }});
        updateFilterResetState();
    }}

    // Modal
    function openModal(sectionId) {{
        modalSection = DATA.sections.find(s => s.id === sectionId);
        if (!modalSection) return;
        modalSubview = null;
        invalidateModalRenderedView();
        modalZoom = 1; modalPanX = 0; modalPanY = 0;
        modalPointerMoved = false;
        modalMagicWandActive = false;
        isDrawingModalLasso = false;
        modalLassoPath = [];
        document.getElementById('modal-magic-wand-btn')?.classList.remove('active');
        modalAnnotationModeActive = false;
        isDrawingModalAnnotation = false;
        modalAnnotationPath = [];
        document.getElementById('modal-annotate-btn')?.classList.remove('active');

        updateModalHeader();
        document.getElementById('modal').classList.add('active');
        updateSelectionInfo();
        updateSectionRotationIndicators(sectionId);
        renderLegend('modal-legend');
        requestAnimationFrame(() => {{
            const canvas = document.getElementById('modal-canvas');
            if (canvas) canvas.style.cursor = 'grab';
            setModalControlsCollapsed(modalControlsCollapsed);
            updateModalToolbarState();
            layoutModalControls();
            layoutModalAnnotationPanel();
            renderModalAnnotationPanel();
            renderModalSection();
        }});
    }}

    function closeModal() {{
        invalidateModalRenderedView();
        document.getElementById('modal').classList.remove('active');
        modalMagicWandActive = false;
        isDrawingModalLasso = false;
        modalLassoPath = [];
        modalAnnotationModeActive = false;
        isDrawingModalAnnotation = false;
        modalAnnotationPath = [];
        modalPointerMoved = false;
        document.getElementById('modal-magic-wand-btn')?.classList.remove('active');
        document.getElementById('modal-annotate-btn')?.classList.remove('active');
        document.querySelector('#modal .modal-controls')?.classList.remove('dragging');
        modalSubview = null;
        modalSection = null;
        updateModalHeader();
        updateSectionRotationIndicators();
        renderModalAnnotationPanel();
        hideTooltip();
    }}

    function getModalAnnotationPanelBounds() {{
        const container = document.getElementById('modal-canvas-container');
        const panel = document.getElementById('modal-annotation-panel');
        if (!container || !panel) return null;
        const containerRect = container.getBoundingClientRect();
        if (containerRect.width <= 0 || containerRect.height <= 0) return null;

        const margin = 8;
        const controls = container.querySelector('.modal-controls');
        let controlsTop = containerRect.height - margin;
        if (controls && !controls.classList.contains('hidden')) {{
            const controlsRect = controls.getBoundingClientRect();
            controlsTop = Math.max(
                margin,
                Math.min(controlsTop, controlsRect.top - containerRect.top - margin)
            );
        }}

        const panelWidth = panel.offsetWidth || Math.min(360, Math.max(120, containerRect.width - margin * 2));
        const panelHeight = panel.offsetHeight || 200;
        const maxHeight = Math.max(120, controlsTop - margin);

        const minX = margin;
        const maxX = Math.max(minX, containerRect.width - panelWidth - margin);
        const minY = margin;
        const maxY = Math.max(minY, controlsTop - panelHeight - margin);

        return {{ minX, maxX, minY, maxY, maxHeight }};
    }}

    function clampModalAnnotationPanelPosition(x, y) {{
        const bounds = getModalAnnotationPanelBounds();
        if (!bounds) return {{ x, y }};
        return {{
            x: Math.min(bounds.maxX, Math.max(bounds.minX, x)),
            y: Math.min(bounds.maxY, Math.max(bounds.minY, y)),
        }};
    }}

    function applyModalAnnotationPanelPosition(position) {{
        const panel = document.getElementById('modal-annotation-panel');
        if (!panel) return;
        if (!position) {{
            panel.style.left = '';
            panel.style.top = '';
            panel.style.right = '8px';
            panel.style.bottom = '82px';
            return;
        }}
        panel.style.left = `${{position.x}}px`;
        panel.style.top = `${{position.y}}px`;
        panel.style.right = 'auto';
        panel.style.bottom = 'auto';
    }}

    function layoutModalAnnotationPanel(forceReset = false) {{
        const panel = document.getElementById('modal-annotation-panel');
        if (!panel) return;
        if (forceReset) modalAnnotationPanelCustomPosition = null;

        const initialBounds = getModalAnnotationPanelBounds();
        if (!initialBounds) return;
        panel.style.maxHeight = `${{Math.floor(initialBounds.maxHeight)}}px`;

        const bounds = getModalAnnotationPanelBounds();
        if (!bounds) return;

        if (modalAnnotationPanelCustomPosition) {{
            modalAnnotationPanelCustomPosition = clampModalAnnotationPanelPosition(
                modalAnnotationPanelCustomPosition.x,
                modalAnnotationPanelCustomPosition.y
            );
            applyModalAnnotationPanelPosition(modalAnnotationPanelCustomPosition);
            return;
        }}

        applyModalAnnotationPanelPosition({{ x: bounds.maxX, y: bounds.maxY }});
    }}

    function initModalAnnotationPanelDragging() {{
        const panel = document.getElementById('modal-annotation-panel');
        const handle = panel?.querySelector('.modal-annotation-title');
        const container = document.getElementById('modal-canvas-container');
        if (!panel || !handle || !container) return;

        const stopDragging = () => {{
            if (!isDraggingModalAnnotationPanel) return;
            isDraggingModalAnnotationPanel = false;
            panel.classList.remove('dragging');
        }};

        handle.addEventListener('mousedown', (e) => {{
            if (e.button !== 0) return;
            const containerRect = container.getBoundingClientRect();
            const panelRect = panel.getBoundingClientRect();
            const next = clampModalAnnotationPanelPosition(
                panelRect.left - containerRect.left,
                panelRect.top - containerRect.top
            );
            modalAnnotationPanelCustomPosition = next;
            applyModalAnnotationPanelPosition(next);
            modalAnnotationPanelDragOffsetX = e.clientX - panelRect.left;
            modalAnnotationPanelDragOffsetY = e.clientY - panelRect.top;
            isDraggingModalAnnotationPanel = true;
            panel.classList.add('dragging');
            e.preventDefault();
        }});

        document.addEventListener('mousemove', (e) => {{
            if (!isDraggingModalAnnotationPanel) return;
            const containerRect = container.getBoundingClientRect();
            const next = clampModalAnnotationPanelPosition(
                e.clientX - containerRect.left - modalAnnotationPanelDragOffsetX,
                e.clientY - containerRect.top - modalAnnotationPanelDragOffsetY
            );
            modalAnnotationPanelCustomPosition = next;
            applyModalAnnotationPanelPosition(next);
        }});

        document.addEventListener('mouseup', stopDragging);
        window.addEventListener('blur', stopDragging);
    }}

    // Grid
    function initGrid() {{
        const grid = document.getElementById('grid');
        grid.innerHTML = '';
        grid.classList.toggle('single-section-layout', (DATA.n_sections || 0) === 1);

        DATA.sections.forEach(section => {{
            const panel = document.createElement('div');
            panel.className = 'section-panel';
            panel.dataset.sectionId = section.id;

            // Apply outline color
            const outlineValue = OUTLINE_BY ? section.metadata?.[OUTLINE_BY] : null;
            const borderColor = getOutlineColor(outlineValue);
            if (borderColor) {{
                panel.style.borderColor = borderColor;
                panel.style.borderWidth = '3px';
            }}

            const metaParts = Object.entries(section.metadata || {{}})
                .map(([k, v]) => `${{formatMetadataLabel(k)}}: ${{v}}`).join(' | ');
            const metaHtml = metaParts ? `<div class="section-meta">${{metaParts}}</div>` : '';
            const rotationLabel = formatRotationLabel(getSectionRotationDeg(section));

            panel.innerHTML = `
                <div class="section-header">
                    <div class="section-header-main">${{section.id}}${{metaHtml}}</div>
                    <div class="section-header-actions">
                        <span class="section-rotation-label" data-section-rotation-label>${{rotationLabel}}</span>
                        <div class="section-rotate-group">
                            <button class="section-rotate-btn" type="button" data-rotate-step="-45" title="Rotate section -45 degrees">-45</button>
                            <button class="section-rotate-btn" type="button" data-rotate-step="45" title="Rotate section +45 degrees">+45</button>
                            <button class="section-rotate-btn" type="button" data-rotate-reset title="Reset section rotation">0</button>
                        </div>
                        <span class="expand-icon">&#x26F6;</span>
                    </div>
                </div>
                <canvas class="section-canvas"></canvas>
            `;
            panel.querySelectorAll('[data-rotate-step], [data-rotate-reset]').forEach((btn) => {{
                btn.addEventListener('click', (event) => {{
                    event.preventDefault();
                    event.stopPropagation();
                    if (btn.hasAttribute('data-rotate-reset')) {{
                        resetSectionRotation(section.id);
                        return;
                    }}
                    rotateSectionBy(section.id, Number(btn.dataset.rotateStep || 0));
                }});
                btn.addEventListener('mousedown', (event) => event.stopPropagation());
            }});
            panel.addEventListener('click', () => openModal(section.id));
            grid.appendChild(panel);
        }});

        // Lazy render thumbnails while scrolling (skip offscreen panels).
        grid.addEventListener('scroll', () => {{
            requestAnimationFrame(renderAllSections);
        }});
    }}

    // Controls
    function initControls() {{
        const colorSelect = document.getElementById('color-select');
        DATA.available_colors.forEach(col => {{
            const opt = document.createElement('option');
            opt.value = col;
            opt.textContent = col;
            opt.selected = col === currentColor;
            colorSelect.appendChild(opt);
        }});

        const isInsightsVisible = () => {{
            const panel = document.getElementById('color-panel');
            return panel && !panel.classList.contains('collapsed');
        }};

        const refreshInsights = () => {{
            if (!isInsightsVisible()) return;
            renderColorList(document.getElementById('color-search')?.value || '');
            const isSummary = document.getElementById('color-tab-aggregate')?.classList.contains('active');
            const isCompare = document.getElementById('color-tab-compare')?.classList.contains('active');
            const isGenes = document.getElementById('color-tab-genes')?.classList.contains('active');
            const isNeighborhood = document.getElementById('color-tab-neighbors')?.classList.contains('active');
            const isRegions = document.getElementById('color-tab-regions')?.classList.contains('active');
            if (isSummary) {{
                renderColorAggregation();
                renderCellTypeTrend();
            }} else if (isCompare) {{
                renderClusterDE();
            }} else if (isNeighborhood) {{
                renderNeighborStats();
                renderInteractionBrowser();
            }} else if (isGenes) {{
                const isDot = document.getElementById('genes-tab-dotplot')?.classList.contains('active');
                const isMarkers = document.getElementById('genes-tab-markers')?.classList.contains('active');
                if (isDot) renderDotplot();
                if (isMarkers) renderMarkerGenes();
            }} else if (isRegions) {{
                renderAnnotationComparison();
            }}
        }};

        colorSelect.addEventListener('change', (e) => {{
            currentColor = e.target.value;
            currentGene = null;
            celltypeTrendTarget = null;
            modalSelectedCategory = null;
            modalTypeSelectEnabled = false;
            document.getElementById('gene-input').value = '';
            (DATA.sections || []).forEach(s => {{ if (s && s._colorCache) s._colorCache = {{}}; }});
            hiddenCategories.clear();
            if (DATA.colors_meta?.[currentColor] && !DATA.colors_meta[currentColor].is_continuous) {{
                selectionSummaryColor = currentColor;
            }}
            updateExpressionScaleUI();
            renderLegend('legend');
            renderLegend('modal-legend');
            renderAllSections();
            if (modalSection) renderModalSection();
            if (umapVisible) renderUMAP();
            updateSelectionInfo();
            renderGeneDiscoveryPanel();
            refreshInsights();
        }});

        const geneList = document.getElementById('gene-list');
        (DATA.available_genes || []).forEach(gene => {{
            const opt = document.createElement('option');
            opt.value = gene;
            geneList.appendChild(opt);
        }});

        const geneInput = document.getElementById('gene-input');
        const geneInputShell = document.getElementById('gene-input-shell');
        const geneDiscoveryPanel = document.getElementById('gene-discovery-panel');
        const genePanelNew = document.getElementById('gene-panel-new');
        recentGenes = loadRecentGenes();
        savedGenePanels = loadSavedGenePanels();
        renderGeneDiscoveryPanel();

        geneInput?.addEventListener('focus', () => {{
            setGeneDiscoveryOpen(true);
            renderGeneDiscoveryPanel();
        }});
        geneInput?.addEventListener('input', () => {{
            setGeneDiscoveryOpen(true);
            geneDiscoveryResults = getGeneSearchResults(geneInput.value);
            geneDiscoveryActiveIndex = geneDiscoveryResults.length ? 0 : -1;
            renderGeneDiscoveryPanel();
        }});
        geneInput?.addEventListener('keydown', async (e) => {{
            if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {{
                setGeneDiscoveryOpen(true);
                geneDiscoveryResults = getGeneSearchResults(geneInput.value);
                if (!geneDiscoveryResults.length) {{
                    geneDiscoveryActiveIndex = -1;
                    renderGeneDiscoveryPanel();
                    return;
                }}
                e.preventDefault();
                const delta = e.key === 'ArrowDown' ? 1 : -1;
                geneDiscoveryActiveIndex = (geneDiscoveryActiveIndex + delta + geneDiscoveryResults.length) % geneDiscoveryResults.length;
                renderGeneDiscoveryPanel();
                return;
            }}
            if (e.key === 'Enter') {{
                const highlightedGene = geneDiscoveryResults[geneDiscoveryActiveIndex];
                const exactGene = resolveCanonicalGeneName(geneInput.value);
                if (!highlightedGene && !exactGene && geneInput.value.trim()) return;
                e.preventDefault();
                await activateViewerGene(highlightedGene || exactGene || geneInput.value, {{ showErrors: false }});
                refreshInsights();
                return;
            }}
            if (e.key === 'Escape') {{
                setGeneDiscoveryOpen(false);
            }}
        }});
        geneInput?.addEventListener('change', async () => {{
            const raw = geneInput.value.trim();
            if (!raw) {{
                await activateViewerGene('', {{ showErrors: false }});
                refreshInsights();
                return;
            }}
            const exactGene = resolveCanonicalGeneName(raw);
            if (exactGene) {{
                await activateViewerGene(exactGene, {{ showErrors: false }});
                refreshInsights();
            }} else {{
                renderGeneDiscoveryPanel();
            }}
        }});
        genePanelNew?.addEventListener('click', () => {{
            const suggestedName = currentGene ? `${{currentGene}} panel` : '';
            const panelName = prompt('Panel name', suggestedName);
            if (!panelName) return;
            upsertSavedGenePanel(panelName, getGenePanelSeedToken());
            renderGeneDiscoveryPanel();
            setGeneDiscoveryOpen(true);
        }});
        geneDiscoveryPanel?.addEventListener('click', async (event) => {{
            const target = event.target;
            const activateBtn = target?.closest?.('[data-gene-activate]');
            if (activateBtn) {{
                event.preventDefault();
                await activateViewerGene(activateBtn.getAttribute('data-gene-activate') || '', {{ showErrors: false }});
                refreshInsights();
                return;
            }}

            const addBtn = target?.closest?.('[data-gene-panel-add]');
            if (addBtn) {{
                event.preventDefault();
                const panelName = addBtn.getAttribute('data-gene-panel-add') || '';
                const geneToken = getGenePanelSeedToken();
                if (!geneToken) {{
                    alert('Load or type an exact gene first, then add it to a panel.');
                    return;
                }}
                upsertSavedGenePanel(panelName, geneToken);
                renderGeneDiscoveryPanel();
                setGeneDiscoveryOpen(true);
                return;
            }}

            const deleteBtn = target?.closest?.('[data-gene-panel-delete]');
            if (deleteBtn) {{
                event.preventDefault();
                const panelName = deleteBtn.getAttribute('data-gene-panel-delete') || '';
                if (!panelName) return;
                if (!confirm(`Delete saved panel "${{panelName}}"?`)) return;
                deleteSavedGenePanel(panelName);
                renderGeneDiscoveryPanel();
                setGeneDiscoveryOpen(true);
            }}
        }});
        document.addEventListener('mousedown', (event) => {{
            if (!geneInputShell || geneInputShell.contains(event.target)) return;
            setGeneDiscoveryOpen(false);
        }});

        const spotRange = document.getElementById('spot-size');
        if (spotRange) {{
            const min = parseFloat(spotRange.min || '0');
            const max = parseFloat(spotRange.max || '100');
            spotSize = Math.min(max, Math.max(min, spotSize));
            spotRange.value = String(spotSize);
            spotRange.addEventListener('input', (e) => {{
                spotSize = parseFloat(e.target.value);
                renderAllSections();
                if (modalSection) renderModalSection();
            }});
        }}
        document.getElementById('spot-size-dec')?.addEventListener('click', () => stepRange(spotRange, -1));
        document.getElementById('spot-size-inc')?.addEventListener('click', () => stepRange(spotRange, 1));

        const umapRange = document.getElementById('umap-spot-size');
        if (umapRange) {{
            umapRange.addEventListener('input', (e) => {{
                umapSpotSize = parseFloat(e.target.value);
                document.getElementById('umap-spot-size-label').textContent = umapSpotSize.toFixed(2);
                if (umapVisible) renderUMAP();
            }});
        }}
        document.getElementById('umap-spot-size-dec')?.addEventListener('click', () => stepRange(umapRange, -1));
        document.getElementById('umap-spot-size-inc')?.addEventListener('click', () => stepRange(umapRange, 1));

        document.getElementById('screenshot-btn').addEventListener('click', screenshotFullPage);

        // Legend toggle
        document.getElementById('legend-toggle').addEventListener('click', toggleLegendPanel);

        // Insights panel toggle
        buildColorPanel();
        const colorToggle = document.getElementById('color-toggle');
        colorToggle.addEventListener('click', toggleInsightsPanel);

        const infoTrigger = document.getElementById('info-trigger');
        if (infoTrigger) {{
            const infoPopover = document.getElementById('info-popover');
            infoTrigger.addEventListener('click', (event) => {{
                event.stopPropagation();
                if (!infoPopover) return;
                const isActive = infoPopover.classList.toggle('active');
                infoPopover.setAttribute('aria-hidden', isActive ? 'false' : 'true');
            }});
            document.addEventListener('click', (event) => {{
                if (!infoPopover || !infoPopover.classList.contains('active')) return;
                if (infoPopover.contains(event.target) || event.target === infoTrigger) return;
                infoPopover.classList.remove('active');
                infoPopover.setAttribute('aria-hidden', 'true');
            }});
            document.addEventListener('keydown', (event) => {{
                if (!infoPopover || !infoPopover.classList.contains('active')) return;
                if (event.key === 'Escape') {{
                    infoPopover.classList.remove('active');
                    infoPopover.setAttribute('aria-hidden', 'true');
                }}
            }});
        }}

    }}

    function initModal() {{
        document.getElementById('modal-close').addEventListener('click', closeModal);
        document.getElementById('modal-controls-toggle')?.addEventListener('click', () => {{
            setModalControlsCollapsed(!modalControlsCollapsed);
        }});
        document.getElementById('modal').addEventListener('click', (e) => {{
            if (e.target.id === 'modal') closeModal();
        }});

        document.getElementById('zoom-in').addEventListener('click', () => {{
            modalZoom = Math.min(modalZoom * 1.5, 20);
            renderModalSection();
        }});
        document.getElementById('zoom-out').addEventListener('click', () => {{
            modalZoom = Math.max(modalZoom / 1.5, 0.1);
            renderModalSection();
        }});
        document.getElementById('zoom-reset').addEventListener('click', zoomModalToFit);
        document.getElementById('modal-rotate-left')?.addEventListener('click', () => {{
            if (!modalSection) return;
            rotateSectionBy(modalSection.id, -ROTATION_STEP_DEG);
        }});
        document.getElementById('modal-rotate-right')?.addEventListener('click', () => {{
            if (!modalSection) return;
            rotateSectionBy(modalSection.id, ROTATION_STEP_DEG);
        }});
        document.getElementById('modal-rotate-reset')?.addEventListener('click', () => {{
            if (!modalSection) return;
            resetSectionRotation(modalSection.id);
        }});
        document.getElementById('modal-screenshot-btn')?.addEventListener('click', screenshotModalView);
        document.getElementById('modal-exit-subview-btn')?.addEventListener('click', () => {{
            exitModalSubview();
        }});

        const modalRange = document.getElementById('modal-spot-size');
        if (modalRange) {{
            const min = parseFloat(modalRange.min || '0');
            const max = parseFloat(modalRange.max || '100');
            modalSpotSize = Math.min(max, Math.max(min, modalSpotSize));
            modalRange.value = String(modalSpotSize);
            document.getElementById('modal-spot-size-label').textContent = modalSpotSize;
            modalRange.addEventListener('input', (e) => {{
                modalSpotSize = parseFloat(e.target.value);
                document.getElementById('modal-spot-size-label').textContent = modalSpotSize;
                renderModalSection();
            }});
        }}
        document.getElementById('modal-spot-size-dec')?.addEventListener('click', () => stepRange(modalRange, -1));
        document.getElementById('modal-spot-size-inc')?.addEventListener('click', () => stepRange(modalRange, 1));

        const container = document.getElementById('modal-canvas-container');
        container.addEventListener('wheel', (e) => {{
            if (!modalSection) return;
            e.preventDefault();
            const rect = container.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;
            const currentTransform = getModalViewTransform(rect);
            if (!currentTransform) return;
            const anchor = currentTransform.screenToData(mouseX, mouseY);
            const nextZoom = Math.max(0.1, Math.min(20, modalZoom * (e.deltaY > 0 ? 0.9 : 1.1)));
            modalZoom = nextZoom;
            const nextTransform = getModalViewTransform(rect);
            if (nextTransform) {{
                const nextAnchorScreen = nextTransform.dataToScreen(anchor.x, anchor.y);
                modalPanX += mouseX - nextAnchorScreen.x;
                modalPanY += mouseY - nextAnchorScreen.y;
            }}
            if (!applyModalInteractionPreview()) {{
                renderModalSection();
                return;
            }}
            scheduleModalInteractionCommit(120);
        }});

        const canvas = document.getElementById('modal-canvas');
        const typeToggleBtn = document.getElementById('modal-type-toggle');
        const typeClearBtn = document.getElementById('modal-type-clear');
        const modalBlendToggleBtn = document.getElementById('modal-blend-toggle');
        const modalBlendPanel = document.getElementById('modal-blend-panel');
        const modalBlendMixRange = document.getElementById('modal-blend-mix');
        const modalBlendMixLabel = document.getElementById('modal-blend-mix-label');
        const modalBlendMeta = document.getElementById('modal-blend-meta');
        const modalBlendControls = {{
            a: {{
                kind: document.getElementById('modal-blend-a-kind'),
                color: document.getElementById('modal-blend-a-color'),
                category: document.getElementById('modal-blend-a-category'),
                gene: document.getElementById('modal-blend-a-gene'),
                scaleRow: document.getElementById('modal-blend-a-scale-row'),
                scaleMin: document.getElementById('modal-blend-a-vmin'),
                scaleMax: document.getElementById('modal-blend-a-vmax'),
                scaleAuto: document.getElementById('modal-blend-a-auto'),
            }},
            b: {{
                kind: document.getElementById('modal-blend-b-kind'),
                color: document.getElementById('modal-blend-b-color'),
                category: document.getElementById('modal-blend-b-category'),
                gene: document.getElementById('modal-blend-b-gene'),
                scaleRow: document.getElementById('modal-blend-b-scale-row'),
                scaleMin: document.getElementById('modal-blend-b-vmin'),
                scaleMax: document.getElementById('modal-blend-b-vmax'),
                scaleAuto: document.getElementById('modal-blend-b-auto'),
            }},
        }};
        const modalWandBtn = document.getElementById('modal-magic-wand-btn');
        const modalAnnotateBtn = document.getElementById('modal-annotate-btn');
        const modalClearSelectionBtn = document.getElementById('modal-clear-selection-btn');
        const modalAnnotationPanel = document.getElementById('modal-annotation-panel');
        const modalAnnotationExportBtn = document.getElementById('modal-annotations-export');
        const modalAnnotationClearSectionBtn = document.getElementById('modal-annotations-clear-section');
        const modalAnnotationClearAllBtn = document.getElementById('modal-annotations-clear-all');
        initModalAnnotationPanelDragging();
        initModalControlsDragging();

        function updateModalCanvasCursor() {{
            if (!canvas) return;
            if (isSplitDragging) {{
                canvas.style.cursor = 'ew-resize';
            }} else if (isDragging) {{
                canvas.style.cursor = 'grabbing';
            }} else if (modalMagicWandActive || modalAnnotationModeActive) {{
                canvas.style.cursor = 'crosshair';
            }} else {{
                canvas.style.cursor = 'grab';
            }}
        }}

        function setSelectOptions(selectEl, values, selectedValue) {{
            if (!selectEl) return;
            selectEl.innerHTML = '';
            values.forEach((entry) => {{
                const value = (entry && typeof entry === 'object') ? entry.value : entry;
                const label = (entry && typeof entry === 'object') ? (entry.label ?? entry.value) : entry;
                const opt = document.createElement('option');
                opt.value = value;
                opt.textContent = label;
                selectEl.appendChild(opt);
            }});
            if (selectedValue !== undefined && selectedValue !== null) {{
                selectEl.value = selectedValue;
            }}
        }}

        function ensureModalBlendDefaults() {{
            const catCols = getCategoricalColorColumns();
            const preferredCol = catCols.includes(currentColor) ? currentColor : (catCols[0] || null);
            const genes = DATA.available_genes || [];
            const firstGene = genes.length ? genes[0] : '';

            const normalize = (entry, preferSecondCategory = false) => {{
                if (!entry.kind) entry.kind = catCols.length ? 'cell' : 'gene';
                if (entry.kind === 'cell') {{
                    if (!catCols.length) {{
                        entry.kind = 'gene';
                    }} else {{
                        if (!entry.color || !catCols.includes(entry.color)) entry.color = preferredCol;
                        const cats = getCategoriesForColorColumn(entry.color);
                        const hasSelectedCategory = entry.category === BLEND_ALL_CATEGORIES || cats.includes(entry.category);
                        if (!entry.category || !hasSelectedCategory) {{
                            if (preferSecondCategory && cats.length > 1) entry.category = cats[1];
                            else entry.category = BLEND_ALL_CATEGORIES;
                        }}
                    }}
                }}
                if (entry.kind === 'gene') {{
                    if (!firstGene) {{
                        if (catCols.length) entry.kind = 'cell';
                    }} else if (!entry.gene || !DATA.genes_meta?.[entry.gene]) {{
                        if (preferSecondCategory && genes.length > 1) entry.gene = genes[1];
                        else entry.gene = firstGene;
                    }}
                }}
                if (entry.kind === 'cell') {{
                    const cats = getCategoriesForColorColumn(entry.color);
                    if (!cats.length && firstGene) {{
                        entry.kind = 'gene';
                        entry.gene = firstGene;
                    }}
                }}
            }};

            normalize(modalBlendSpec.a, false);
            normalize(modalBlendSpec.b, true);
        }}

        function getModalBlendMarkerHtml(side) {{
            const spec = modalBlendSpec?.[side];
            if (!spec || spec.kind !== 'cell') return '';
            const sideLabel = side === 'a' ? 'A' : 'B';
            const colorLabel = spec.color ? formatMetadataLabel(spec.color) : 'Cell type';
            const byColor = (DATA.marker_genes || {{}})[spec.color];
            if (!byColor || typeof byColor !== 'object') {{
                return `<div class="modal-blend-marker-header">${{escapeHtml(`${{sideLabel}} markers (${{colorLabel}}): unavailable`)}}</div>`;
            }}
            const category = spec.category;
            if (category && category !== BLEND_ALL_CATEGORIES) {{
                const genes = getMarkerGenesForColorCategory(spec.color, category);
                if (!genes.length) {{
                    return `<div class="modal-blend-marker-header">${{escapeHtml(`${{sideLabel}} markers (${{colorLabel}} · ${{category}}): unavailable`)}}</div>`;
                }}
                return [
                    `<div class="modal-blend-marker-header">${{escapeHtml(`${{sideLabel}} markers (${{colorLabel}} · ${{category}})`)}}</div>`,
                    `<div class="modal-blend-marker-row">${{escapeHtml(genes.join(' '))}}</div>`,
                ].join('');
            }}

            const categories = getCategoriesForColorColumn(spec.color);
            const categoryOrder = categories.length ? categories : Object.keys(byColor);
            const rows = categoryOrder.map(cat => {{
                const genes = getMarkerGenesForColorCategory(spec.color, cat);
                if (!genes.length) return '';
                return `<div class="modal-blend-marker-row"><span class="modal-blend-marker-name">${{escapeHtml(String(cat))}}</span>: ${{escapeHtml(genes.join(' '))}}</div>`;
            }}).filter(Boolean);
            if (!rows.length) {{
                return `<div class="modal-blend-marker-header">${{escapeHtml(`${{sideLabel}} markers (${{colorLabel}}): unavailable`)}}</div>`;
            }}
            return [
                `<div class="modal-blend-marker-header">${{escapeHtml(`${{sideLabel}} markers (${{colorLabel}} · all categories)`)}}</div>`,
                ...rows,
            ].join('');
        }}

        function updateModalBlendLabels() {{
            if (modalBlendMixLabel) {{
                const pctA = Math.round(modalBlendMix * 100);
                const pctB = 100 - pctA;
                modalBlendMixLabel.textContent = `A left ${{pctA}}% • B right ${{pctB}}%`;
            }}
            if (modalBlendMeta) {{
                const summary = `${{getModalBlendVariableLabel(modalBlendSpec.a)}} \u2194 ${{getModalBlendVariableLabel(modalBlendSpec.b)}}`;
                const aMarkers = getModalBlendMarkerHtml('a');
                const bMarkers = getModalBlendMarkerHtml('b');
                const markerBlocks = [aMarkers, bMarkers].filter(Boolean);
                const hasMarkers = markerBlocks.length > 0;
                const toggleLabel = modalBlendMarkersCollapsed ? 'Show marker genes' : 'Hide marker genes';
                modalBlendMeta.innerHTML = [
                    `<div class="modal-blend-summary">${{escapeHtml(summary)}}</div>`,
                    hasMarkers
                        ? `<button type="button" class="modal-blend-markers-toggle" id="modal-blend-markers-toggle" aria-expanded="${{modalBlendMarkersCollapsed ? 'false' : 'true'}}">${{toggleLabel}}</button>`
                        : '',
                    hasMarkers
                        ? `<div class="modal-blend-markers ${{modalBlendMarkersCollapsed ? 'collapsed' : ''}}">${{markerBlocks.join('')}}</div>`
                        : '',
                ].join('');
                const markersToggle = document.getElementById('modal-blend-markers-toggle');
                markersToggle?.addEventListener('click', () => {{
                    modalBlendMarkersCollapsed = !modalBlendMarkersCollapsed;
                    updateModalBlendLabels();
                }});
            }}
        }}

        function syncModalBlendSide(side) {{
            const controls = modalBlendControls[side];
            const spec = modalBlendSpec[side];
            if (!controls || !spec) return;
            controls.kind.value = spec.kind;
            const isCell = spec.kind === 'cell';
            controls.color.style.display = isCell ? '' : 'none';
            controls.category.style.display = isCell ? '' : 'none';
            controls.gene.style.display = isCell ? 'none' : '';
            if (controls.scaleRow) controls.scaleRow.style.display = isCell ? 'none' : '';
            if (isCell) {{
                const cols = getCategoricalColorColumns();
                setSelectOptions(controls.color, cols, spec.color);
                const cats = getCategoriesForColorColumn(spec.color);
                const catOptions = [{{ value: BLEND_ALL_CATEGORIES, label: 'All categories' }}]
                    .concat(cats.map(cat => ({{ value: cat, label: cat }})));
                setSelectOptions(controls.category, catOptions, spec.category);
            }} else {{
                controls.gene.value = spec.gene || '';
                const gene = (spec.gene || '').trim();
                const hasGene = !!(gene && DATA.genes_meta?.[gene]);
                if (hasGene) ensureGeneAutoScale(gene);
                const scale = hasGene ? getGeneScaleRange(gene) : null;
                if (controls.scaleMin) {{
                    controls.scaleMin.value = scale ? scale.vmin.toFixed(3) : '';
                    controls.scaleMin.disabled = !hasGene;
                }}
                if (controls.scaleMax) {{
                    controls.scaleMax.value = scale ? scale.vmax.toFixed(3) : '';
                    controls.scaleMax.disabled = !hasGene;
                }}
                if (controls.scaleAuto) controls.scaleAuto.disabled = !hasGene;
            }}
        }}

        function syncModalBlendUI() {{
            ensureModalBlendDefaults();
            syncModalBlendSide('a');
            syncModalBlendSide('b');
            if (modalBlendPanel) modalBlendPanel.classList.toggle('visible', modalBlendEnabled);
            if (modalBlendToggleBtn) modalBlendToggleBtn.classList.toggle('active', modalBlendEnabled);
            if (modalBlendMixRange) modalBlendMixRange.value = String(Math.round(modalBlendMix * 100));
            updateModalBlendLabels();
            updateModalToolbarState();
        }}

        function applyModalBlendControlChange() {{
            syncModalBlendUI();
            renderLegend('modal-legend');
            if (modalSection) renderModalSection();
        }}

        function applyModalBlendGeneScale(side, useAuto = false) {{
            const controls = modalBlendControls[side];
            const spec = modalBlendSpec[side];
            if (!controls || !spec || spec.kind !== 'gene') return;
            const gene = (spec.gene || '').trim();
            if (!gene || !DATA.genes_meta?.[gene]) return;

            if (useAuto) {{
                const autoScale = computeGenePercentiles(gene);
                if (!autoScale) return;
                geneScaleAuto[gene] = autoScale;
                delete geneScaleOverrides[gene];
            }} else {{
                const vmin = parseFloat(controls.scaleMin?.value);
                const vmax = parseFloat(controls.scaleMax?.value);
                if (!Number.isFinite(vmin) || !Number.isFinite(vmax)) return;
                let adjMin = vmin;
                let adjMax = vmax;
                if (adjMin === adjMax) {{
                    adjMin -= 1e-6;
                    adjMax += 1e-6;
                }} else if (adjMin > adjMax) {{
                    const tmp = adjMin;
                    adjMin = adjMax;
                    adjMax = tmp;
                }}
                geneScaleOverrides[gene] = {{ vmin: adjMin, vmax: adjMax }};
            }}

            syncModalBlendUI();
            rerenderForExpressionScaleChange();
        }}

        modalBlendToggleBtn?.addEventListener('click', () => {{
            modalBlendEnabled = !modalBlendEnabled;
            syncModalBlendUI();
            renderLegend('modal-legend');
            if (modalSection) renderModalSection();
        }});
        modalBlendMixRange?.addEventListener('input', (e) => {{
            modalBlendMix = clamp01(parseFloat(e.target.value) / 100);
            updateModalBlendLabels();
            scheduleModalBlendRender();
        }});
        ['a', 'b'].forEach((side) => {{
            const controls = modalBlendControls[side];
            if (!controls) return;
            controls.kind?.addEventListener('change', () => {{
                modalBlendSpec[side].kind = controls.kind.value === 'gene' ? 'gene' : 'cell';
                applyModalBlendControlChange();
            }});
            controls.color?.addEventListener('change', () => {{
                modalBlendSpec[side].color = controls.color.value;
                modalBlendSpec[side].category = BLEND_ALL_CATEGORIES;
                applyModalBlendControlChange();
            }});
            controls.category?.addEventListener('change', () => {{
                modalBlendSpec[side].category = controls.category.value;
                applyModalBlendControlChange();
            }});
            controls.gene?.addEventListener('change', async () => {{
                await runAsyncUIAction(`Modal gene selection (${{side.toUpperCase()}})`, async () => {{
                    const gene = controls.gene.value.trim();
                    if (!gene) {{
                        modalBlendSpec[side].gene = '';
                        applyModalBlendControlChange();
                        return;
                    }}
                    if (!await ensureGeneAvailable(gene)) {{
                        controls.gene.value = modalBlendSpec[side].gene || '';
                        return;
                    }}
                    modalBlendSpec[side].gene = gene;
                    ensureGeneAutoScale(gene);
                    applyModalBlendControlChange();
                }});
            }});
            controls.scaleMin?.addEventListener('change', () => applyModalBlendGeneScale(side));
            controls.scaleMax?.addEventListener('change', () => applyModalBlendGeneScale(side));
            controls.scaleAuto?.addEventListener('click', () => applyModalBlendGeneScale(side, true));
        }});
        modalBlendPanel?.addEventListener('wheel', (e) => e.stopPropagation());
        modalBlendPanel?.addEventListener('touchmove', (e) => e.stopPropagation());
        modalAnnotationPanel?.addEventListener('wheel', (e) => e.stopPropagation());
        modalAnnotationPanel?.addEventListener('touchmove', (e) => e.stopPropagation());
        syncModalBlendUI();
        renderModalAnnotationPanel();
        layoutModalAnnotationPanel();
        updateModalToolbarState();

        modalAnnotationExportBtn?.addEventListener('click', () => {{
            exportModalAnnotations();
        }});
        modalAnnotationClearSectionBtn?.addEventListener('click', () => {{
            if (!modalSection) return;
            clearModalAnnotationsForSection(modalSection.id);
        }});
        modalAnnotationClearAllBtn?.addEventListener('click', () => {{
            clearAllModalAnnotations();
        }});

        typeToggleBtn.addEventListener('click', () => {{
            const config = getColorConfig();
            if (config.is_continuous) return;
            modalTypeSelectEnabled = !modalTypeSelectEnabled;
            updateModalToolbarState();
        }});
        typeClearBtn.addEventListener('click', () => {{
            modalSelectedCategory = null;
            updateModalToolbarState();
            renderAllSections();
            renderModalSection();
        }});
        modalWandBtn?.addEventListener('click', () => {{
            modalMagicWandActive = !modalMagicWandActive;
            modalWandBtn.classList.toggle('active', modalMagicWandActive);
            if (modalMagicWandActive) {{
                modalAnnotationModeActive = false;
                isDrawingModalAnnotation = false;
                modalAnnotationPath = [];
                modalAnnotateBtn?.classList.remove('active');
            }}
            if (!modalMagicWandActive) {{
                isDrawingModalLasso = false;
                modalLassoPath = [];
            }}
            updateModalToolbarState();
            updateModalCanvasCursor();
            renderModalSection();
        }});
        modalAnnotateBtn?.addEventListener('click', () => {{
            modalAnnotationModeActive = !modalAnnotationModeActive;
            modalAnnotateBtn.classList.toggle('active', modalAnnotationModeActive);
            if (modalAnnotationModeActive) {{
                modalMagicWandActive = false;
                isDrawingModalLasso = false;
                modalLassoPath = [];
                modalWandBtn?.classList.remove('active');
            }} else {{
                isDrawingModalAnnotation = false;
                modalAnnotationPath = [];
            }}
            updateModalToolbarState();
            updateModalCanvasCursor();
            renderModalSection();
        }});
        modalClearSelectionBtn?.addEventListener('click', clearSelection);
        canvas.addEventListener('mousedown', (e) => {{
            modalPointerMoved = false;
            hideTooltip();
            if (modalAnnotationModeActive) {{
                if (!modalSection) return;
                const rect = container.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const y = e.clientY - rect.top;
                isDrawingModalAnnotation = true;
                modalAnnotationPath = [{{ x, y }}];
                renderModalSection();
                return;
            }}
            if (modalMagicWandActive) {{
                if (!modalSection) return;
                const rect = container.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const y = e.clientY - rect.top;
                isDrawingModalLasso = true;
                modalLassoPath = [{{ x, y }}];
                renderModalSection();
                return;
            }}
            if (isSplitDragging) return;
            isDragging = true;
            dragStartX = e.clientX; dragStartY = e.clientY;
            lastPanX = modalPanX; lastPanY = modalPanY;
            updateModalCanvasCursor();
        }});
        canvas.addEventListener('pointerdown', (e) => {{
            if (modalAnnotationModeActive || modalMagicWandActive) return;
            const blendActiveNow = !!(modalSection && getModalBlendRuntimes(modalSection));
            if (!blendActiveNow) return;
            const rect = container.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const splitX = rect.width * modalBlendMix;
            if (Math.abs(mouseX - splitX) > 12) return;
            canvas.setPointerCapture(e.pointerId);
            isSplitDragging = true;
            modalPointerMoved = true;
            updateModalCanvasCursor();
            e.preventDefault();
        }});
        canvas.addEventListener('pointermove', (e) => {{
            if (!isSplitDragging) return;
            const rect = container.getBoundingClientRect();
            modalBlendMix = clamp01((e.clientX - rect.left) / rect.width);
            if (modalBlendMixRange) modalBlendMixRange.value = String(Math.round(modalBlendMix * 100));
            updateModalBlendLabels();
            scheduleModalBlendRender();
            e.preventDefault();
        }});
        canvas.addEventListener('pointerup', (e) => {{
            if (!isSplitDragging) return;
            isSplitDragging = false;
            updateModalCanvasCursor();
            e.preventDefault();
        }});
        canvas.addEventListener('click', (e) => {{
            if (!modalSection || isDragging || isDrawingModalLasso || isDrawingModalAnnotation || modalPointerMoved || modalMagicWandActive || modalAnnotationModeActive) return;
            const config = getColorConfig();
            if (config.is_continuous || !modalTypeSelectEnabled) return;

            const rect = container.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;
            const transform = getModalViewTransform(rect);
            if (!transform) return;

            const cellIdx = findNearestCell(modalSection, mouseX, mouseY, rect, transform);
            if (cellIdx >= 0) {{
                const values = getSectionValues(modalSection);
                const val = values[cellIdx];
                if (val !== null && val !== undefined) {{
                    const catIdx = Math.round(val);
                    const catName = config.categories[catIdx];
                    modalSelectedCategory = catName || null;
                    updateModalToolbarState();
                    renderAllSections();
                    renderModalSection();
                }}
            }}
        }});
        document.addEventListener('mousemove', (e) => {{
            if (isDrawingModalAnnotation) {{
                const rect = container.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const y = e.clientY - rect.top;
                const last = modalAnnotationPath[modalAnnotationPath.length - 1];
                if (!last || Math.abs(x - last.x) + Math.abs(y - last.y) >= 1) {{
                    modalAnnotationPath.push({{ x, y }});
                    modalPointerMoved = true;
                    renderModalSection();
                }}
                return;
            }}
            if (isDrawingModalLasso) {{
                const rect = container.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const y = e.clientY - rect.top;
                const last = modalLassoPath[modalLassoPath.length - 1];
                if (!last || Math.abs(x - last.x) + Math.abs(y - last.y) >= 1) {{
                    modalLassoPath.push({{ x, y }});
                    modalPointerMoved = true;
                    renderModalSection();
                }}
                return;
            }}
            if (!isDragging) return;
            modalPanX = lastPanX + (e.clientX - dragStartX);
            modalPanY = lastPanY + (e.clientY - dragStartY);
            modalPointerMoved = true;
            if (!applyModalInteractionPreview()) {{
                renderModalSection();
                return;
            }}
            scheduleModalInteractionCommit(90);
        }});
        document.addEventListener('mouseup', () => {{
            if (isDrawingModalAnnotation) {{
                isDrawingModalAnnotation = false;
                createModalAnnotationFromPath();
                modalAnnotationPath = [];
            }}
            if (isDrawingModalLasso) {{
                isDrawingModalLasso = false;
                performModalLassoSelection();
                modalLassoPath = [];
            }}
            if (isDragging) {{
                isDragging = false;
                renderModalSection();
            }}
            updateModalCanvasCursor();
        }});
        updateModalCanvasCursor();

        // Tooltip on hover in modal
        canvas.addEventListener('mousemove', (e) => {{
            if (isDragging || isDrawingModalLasso || isDrawingModalAnnotation || modalMagicWandActive || modalAnnotationModeActive || !modalSection) return;

            const rect = container.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;
            const transform = getModalViewTransform(rect);
            if (!transform) return;

            const blendActive = !!getModalBlendRuntimes(modalSection);
            const cellIdx = findNearestCell(
                modalSection,
                mouseX,
                mouseY,
                rect,
                transform,
                blendActive ? {{ ignoreMissing: true, ignoreHidden: true }} : {{}}
            );
            if (cellIdx >= 0) {{
                const content = blendActive
                    ? (getModalSplitTooltipContent(modalSection, cellIdx, transform) || getCellTooltipContent(modalSection, cellIdx))
                    : getCellTooltipContent(modalSection, cellIdx);
                showTooltip(e.clientX, e.clientY, content);
                const changed = updateHoverNeighbors(modalSection, cellIdx);
                if (changed) renderModalSection();
            }} else {{
                hideTooltip();
                if (hoverNeighbors) {{
                    hoverNeighbors = null;
                    renderModalSection();
                }}
            }}
        }});

        canvas.addEventListener('mouseleave', () => {{
            hideTooltip();
            if (hoverNeighbors) {{
                hoverNeighbors = null;
                renderModalSection();
            }}
        }});

        if (DATA.has_neighbors) {{
            const modalGraphBtn = document.getElementById('modal-graph-toggle');
            modalGraphBtn.classList.toggle('active', showGraph);
            modalGraphBtn.addEventListener('click', () => {{
                showGraph = !showGraph;
                modalGraphBtn.classList.toggle('active', showGraph);
                updateModalToolbarState();
                renderAllSections();
                if (modalSection) renderModalSection();
            }});

            const modalNeighborBtn = document.getElementById('modal-neighbor-hover-toggle');
            modalNeighborBtn.classList.toggle('active', neighborHoverEnabled);
            modalNeighborBtn.addEventListener('click', () => {{
                neighborHoverEnabled = !neighborHoverEnabled;
                modalNeighborBtn.classList.toggle('active', neighborHoverEnabled);
                updateModalToolbarState();
                if (!neighborHoverEnabled) {{
                    hoverNeighbors = null;
                    if (modalSection) renderModalSection();
                }}
            }});

            const modalHopSelect = document.getElementById('modal-neighbor-hop-select');
            modalHopSelect.value = neighborHopMode;
            modalHopSelect.addEventListener('change', () => {{
                neighborHopMode = modalHopSelect.value;
                if (modalSection) renderModalSection();
            }});
        }}

        updateModalToolbarState();
        layoutModalControls();
        layoutModalAnnotationPanel();
    }}

    // Initialize (don't wait for external resources)
    document.addEventListener('DOMContentLoaded', () => {{
        let initSucceeded = false;
        try {{
            hydratePackedSections();
            initTheme();
            initGrid();
            initControls();
            initFilters();
            initModal();
            initUMAP();
            initKeyboardShortcuts();
            initHelpTooltips();
            updateSectionRotationIndicators();
            renderLegend('legend');
            updateSelectionInfo();
            initSucceeded = true;
        }} catch (e) {{
            console.error('KaroSpace initialization failed', e);
            if (typeof window.__karospaceShowStartupError === 'function') {{
                window.__karospaceShowStartupError(`KaroSpace failed to initialize: ${{e.message || 'Unknown error'}} (open DevTools console).`);
            }}
        }} finally {{
            // Always clear loader/warning, even if an init step throws.
            hideLoader();
        }}
        if (initSucceeded) requestAnimationFrame(renderAllSections);
    }});
    window.addEventListener('resize', () => {{
        if (DATA.has_umap) applyUMAPPanelState();
        renderAllSections();
        if (modalSection) {{
            layoutModalControls();
            layoutModalAnnotationPanel();
            renderModalSection();
        }}
        if (umapVisible) renderUMAP();
    }});
    </script>
    {footer_logo}
</body>
</html>
'''


def export_to_html(
    dataset: SpatialDataset,
    output_path: str,
    color: str = "leiden",
    title: str = "Spatial Viewer",
    min_panel_size: int = 150,
    spot_size: Union[float, str, None] = "auto",
    downsample: Optional[int] = None,
    theme: str = "light",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    outline_by: Optional[str] = "course",
    viewer_info_html: Optional[str] = None,
    additional_colors: Optional[List[str]] = None,
    genes: Optional[List[str]] = None,
    gene_encoding: str = "auto",
    gene_value_encoding: str = "float32",
    gene_sidecar_format: str = GENE_SIDECAR_FORMAT_JSON_V2,
    gene_storage: str = "embedded",
    gene_aux_path: Optional[str] = None,
    gene_sidecar_shard_size: int = GENE_SIDECAR_SHARD_SIZE,
    gene_sparse_zero_threshold: float = 0.8,
    pack_arrays: bool = True,
    pack_arrays_min_len: int = 1024,
    hvg_limit: int = 20,
    marker_genes_groupby: Optional[List[str]] = None,
    marker_genes_top_n: int = 30,
    use_hvgs: bool = True,
    cluster_de_groupby: Optional[List[str]] = None,
    cluster_de_top_n: int = 20,
    cluster_de_method: str = "wilcoxon",
    cluster_de_layer: Optional[str] = "normalized",
    cluster_de_min_cells: int = 20,
    neighbor_stats_groupby: Optional[List[str]] = None,
    neighbor_stats_permutations: Union[int, None] = None,
    neighbor_stats_seed: int = 0,
    interaction_markers_groupby: Optional[List[str]] = None,
    interaction_markers_top_targets: int = 8,
    interaction_markers_top_genes: int = 20,
    interaction_markers_min_cells: int = 30,
    interaction_markers_min_neighbors: int = 1,
    interaction_markers_method: str = "wilcoxon",
    interaction_markers_layer: Optional[str] = "normalized",
    section_rotations: Optional[Mapping[str, Union[int, float]]] = None,
    gene_correlation_top_n: int = 10,
    cluster_means_n_genes: int = 500,
    spatial_variable_genes_n: int = 200,
) -> str:
    """
    Export spatial dataset to a standalone HTML file.

    Parameters
    ----------
    dataset : SpatialDataset
        Dataset to export
    output_path : str
        Path for output HTML file, or a `.karospace` package when
        `gene_storage="sidecar"`.
    color : str
        Initial color column or gene name
    title : str
        Page title
    min_panel_size : int
        Minimum width of each section panel in pixels (default 150).
        The grid auto-adjusts columns based on screen width.
    spot_size : float | str | None
        Default spot size. Use a positive number for a fixed value, or
        "auto"/"adaptive"/"density" (or None) to estimate from cell density.
    downsample : int, optional
        Downsample cells per section (for large datasets)
    theme : str
        'light' or 'dark'
    vmin, vmax : float, optional
        Min/max for continuous color scale
    outline_by : str, optional
        Metadata column used to color panel outlines (default: "course")
    viewer_info_html : str, optional
        HTML string shown in the Info tab of the color panel.
    additional_colors : list, optional
        Additional obs columns to include for color switching
    genes : list, optional
        Gene names to include for expression visualization
    gene_encoding : str
        "dense", "sparse", or "auto" (default). "auto" uses sparse encoding for
        zero-inflated genes to reduce HTML size.
    gene_value_encoding : str
        Only used for sidecar/package genes. "float32" (default) stores full-
        precision values; "uint16" and "uint8" store linearly quantized values
        using each gene's exported vmin/vmax range.
    gene_sidecar_format : str
        Sidecar storage format. "json-v2" (default) writes JSON shard files.
        "binary-v1" writes indexed binary `.bin` shards for faster random access.
    gene_storage : str
        "embedded" (default) keeps gene data in the HTML. "sidecar" embeds only
        the selected/preloaded genes and writes remaining genes to an auxiliary
        JSON file for downstream lazy loading. `.karospace` package export
        requires `"sidecar"`.
    gene_aux_path : str, optional
        Output path for the sidecar gene JSON when gene_storage="sidecar".
    gene_sidecar_shard_size : int
        Number of genes to write per sidecar shard when gene_storage="sidecar".
        Reduce this for very large datasets to avoid oversized shard JSON files.
    gene_sparse_zero_threshold : float
        Only used when gene_encoding="auto". Use sparse encoding when the
        fraction of zeros is >= this threshold (default: 0.8).
    pack_arrays : bool
        Pack large per-section numeric arrays (coordinates, colors, UMAP) as base64 typed arrays
        for smaller HTML and faster load. Default: True.
    pack_arrays_min_len : int
        Only pack per-section arrays when section cell count is >= this value. Default: 1024.
    hvg_limit : int
        Max number of highly variable genes to include (default 20)
    gene_correlation_top_n : int
        Number of top correlated genes to show per embedded gene (default 10). Set to 0 to disable.
    cluster_means_n_genes : int
        Number of top variable genes to precompute per-cluster means for (default 500).
        These means power the mean-expression panel in the lasso selection summary for
        all genes including sidecar ones. Set to 0 to disable.
    spatial_variable_genes_n : int
        Number of top variable genes to score with Moran's I spatial autocorrelation
        (default 200). Requires a spatial weight matrix in adata.obsp. Set to 0 to disable.
    marker_genes_groupby : list, optional
        Obs columns to compute marker genes for (categorical only). If None/empty,
        marker genes are not computed.
    marker_genes_top_n : int
        Number of top marker genes to keep per group
    cluster_de_groupby : list, optional
        Obs columns to compute pairwise cluster-vs-cluster DE for (categorical only).
        If None/empty, cluster DE is not computed.
    cluster_de_top_n : int
        Number of top DE genes to keep per cluster pair.
    cluster_de_method : str
        Method passed to scanpy.tl.rank_genes_groups for cluster DE
        (e.g. "wilcoxon", "t-test").
    cluster_de_layer : str, optional
        Layer used for cluster DE (default: "normalized" if present).
    cluster_de_min_cells : int
        Minimum cells required in both clusters to report pairwise DE.
    neighbor_stats_groupby : list, optional
        Obs columns to compute neighbor composition stats for (categorical only).
        If None/empty, neighbor stats are not computed.
    neighbor_stats_permutations : int
        Number of permutations for neighbor enrichment z-scores (0 disables)
    neighbor_stats_seed : int
        Random seed used for neighbor permutations
    interaction_markers_groupby : list, optional
        Obs columns to compute contact-conditioned interaction markers for.
    interaction_markers_top_targets : int
        Number of target cell types to evaluate per source.
    interaction_markers_top_genes : int
        Number of top DE genes to keep per source-target interaction.
    interaction_markers_min_cells : int
        Minimum cells required in both contact+ and contact- groups.
    interaction_markers_min_neighbors : int
        Minimum target neighbors to classify source cells as contact+.
    interaction_markers_method : str
        Method passed to scanpy.tl.rank_genes_groups (e.g. "wilcoxon", "t-test").
    interaction_markers_layer : str, optional
        Layer used for interaction DE (default: "normalized" if present).
    section_rotations : mapping, optional
        Optional mapping of section_id -> initial rotation angle in degrees.
        Angles are stored exactly (normalized modulo 360) for the interactive viewer.

    Returns
    -------
    str
        Path to created HTML file
    """
    requested_output_path = Path(output_path).expanduser()
    package_mode = requested_output_path.suffix.lower() == ".karospace"
    package_output_path_obj: Optional[Path] = None
    package_loader_output_path: Optional[Path] = None
    package_temp_dir: Optional[tempfile.TemporaryDirectory] = None
    package_entry_html = "index.html"
    package_gene_manifest_name: Optional[str] = None
    if package_mode:
        package_output_path_obj = requested_output_path.resolve()
        package_loader_output_path = package_output_path_obj.with_suffix(".loader.html")
        if str(gene_storage or "embedded").strip().lower() != "sidecar":
            raise ValueError(".karospace export requires gene_storage='sidecar'")
        if gene_aux_path:
            aux_name = Path(gene_aux_path).name
            if aux_name != str(Path(gene_aux_path)):
                raise ValueError("gene_aux_path must be a filename when exporting .karospace")
            if Path(aux_name).suffix.lower() != ".json":
                raise ValueError("gene_aux_path must end with .json when exporting .karospace")
            package_gene_manifest_name = aux_name
        else:
            package_gene_manifest_name = f"{package_output_path_obj.stem}.genes.json"
        package_temp_dir = tempfile.TemporaryDirectory(prefix="karospace-package-")
        requested_output_path = Path(package_temp_dir.name) / package_entry_html
        gene_aux_path = package_gene_manifest_name
        output_path = str(requested_output_path)

    # Theme colors
    if theme == "dark":
        colors = {
            "background": "#1a1a1a",
            "text_color": "#e0e0e0",
            "header_bg": "#2a2a2a",
            "panel_bg": "#2a2a2a",
            "border_color": "#404040",
            "input_bg": "#333333",
            "muted_color": "#888888",
            "hover_bg": "#3a3a3a",
            "graph_color": "rgba(255, 255, 255, 0.12)",
        }
    else:
        colors = {
            "background": "#f5f5f5",
            "text_color": "#1a1a1a",
            "header_bg": "#ffffff",
            "panel_bg": "#ffffff",
            "border_color": "#e0e0e0",
            "input_bg": "#ffffff",
            "muted_color": "#666666",
            "hover_bg": "#f0f0f0",
            "graph_color": "rgba(0, 0, 0, 0.12)",
        }

    gene_storage = str(gene_storage or "embedded").strip().lower()
    if gene_storage not in {"embedded", "sidecar"}:
        raise ValueError("gene_storage must be one of: 'embedded', 'sidecar'")
    gene_sidecar_format = _validate_gene_sidecar_format(gene_sidecar_format)
    gene_sidecar_shard_size = int(gene_sidecar_shard_size)
    if gene_sidecar_shard_size < 1:
        raise ValueError("gene_sidecar_shard_size must be >= 1")
    if gene_sidecar_format == GENE_SIDECAR_FORMAT_BINARY_V1 and gene_value_encoding not in {"uint8", "uint16"}:
        raise ValueError("binary-v1 sidecar format requires gene_value_encoding='uint8' or 'uint16'")
    resolved_section_rotations = _resolve_section_rotations(dataset, section_rotations)

    # Prefer highly variable genes for expression if available; otherwise use provided genes
    hv_genes = None
    if use_hvgs and "highly_variable" in dataset.adata.var.columns:
        hv_mask = dataset.adata.var["highly_variable"].to_numpy(dtype=bool)
        if hv_mask.any():
            hv_genes = dataset.adata.var_names[hv_mask].tolist()[:max(0, int(hvg_limit))]
    if hv_genes is not None:
        genes = hv_genes
    embedded_genes = list(dict.fromkeys([g for g in (genes or []) if g in dataset.adata.var_names]))
    sidecar_genes: List[str] = []
    resolved_gene_aux_path: Optional[Path] = None
    resolved_gene_aux_dir: Optional[Path] = None
    gene_aux_url: Optional[str] = None
    sidecar_export_indices = None
    if gene_storage == "sidecar":
        sidecar_genes = [g for g in dataset.var_names if g not in set(embedded_genes)]
        sidecar_export_indices = dataset._get_export_section_indices(downsample=downsample)
        output_path_obj = requested_output_path
        if gene_aux_path:
            resolved_gene_aux_path = Path(gene_aux_path).expanduser()
        else:
            resolved_gene_aux_path = output_path_obj.with_suffix(".genes.json")
        if not resolved_gene_aux_path.is_absolute():
            resolved_gene_aux_path = (output_path_obj.parent / resolved_gene_aux_path).resolve()
        resolved_gene_aux_dir = resolved_gene_aux_path.with_suffix("")

    if outline_by and outline_by not in dataset.metadata_columns:
        print(f"  Warning: outline_by '{outline_by}' not in metadata columns; no outlines will be shown.")
    if int(cluster_de_top_n) < 1:
        raise ValueError("cluster_de_top_n must be >= 1")
    if int(cluster_de_min_cells) < 1:
        raise ValueError("cluster_de_min_cells must be >= 1")

    if viewer_info_html is None:
        viewer_info_html = (
            '<div class="info-block">'
            '<div class="info-title">Viewer</div>'
            '<div class="info-text">KaroSpace interactive spatial viewer for exploring '
            'sections, cell types, and gene expression.</div>'
            '</div>'
            '<div class="info-block">'
            '<div class="info-title">Contact</div>'
            '<div class="info-list">'
            '<div>Karolinska Institutet, Stockholm</div>'
            '<div><a class="info-link" href="mailto:christoffer.mattsson.langseth@ki.se">'
            'christoffer.mattsson.langseth@ki.se</a></div>'
            '<div><a class="info-link" href="https://ki.se/personer/christoffer-mattsson-langseth" '
            'target="_blank" rel="noopener noreferrer">Profile</a></div>'
            '<div><a class="info-link" href="https://orcid.org/0000-0003-2230-8594" '
            'target="_blank" rel="noopener noreferrer">ORCID</a></div>'
            '<div><a class="info-link" href="https://www.linkedin.com/in/christoffer-mattsson-langseth-76427011a" '
            'target="_blank" rel="noopener noreferrer">LinkedIn</a></div>'
            '</div>'
            '</div>'
        )
    viewer_info_html_safe = viewer_info_html.replace('{', '{{').replace('}', '}}')

    # Get data with multiple color layers and genes
    if neighbor_stats_groupby is None:
        neighbor_stats_groupby = [color]
        if additional_colors:
            neighbor_stats_groupby.extend(additional_colors)

    if neighbor_stats_permutations is None:
        # Auto-tune: permutation z-scores are expensive for very large datasets.
        # Keep the feature (counts + mean degree) always, but skip permutations unless requested.
        neighbor_stats_permutations = 0 if int(dataset.adata.n_obs) >= 200_000 else 20

    marker_genes_groupby = _resolve_marker_genes_groupby(
        dataset=dataset,
        color=color,
        marker_genes_groupby=marker_genes_groupby,
    )
    companion_analytics = dataset.get_companion_analytics()

    data = dataset.to_json_data(
        color,
        downsample=downsample,
        vmin=vmin,
        vmax=vmax,
        additional_colors=additional_colors,
        genes=embedded_genes,
        gene_encoding=gene_encoding,
        gene_sparse_zero_threshold=gene_sparse_zero_threshold,
        section_array_pack=pack_arrays,
        section_array_pack_min_len=pack_arrays_min_len,
        marker_genes_groupby=marker_genes_groupby,
        marker_genes_top_n=marker_genes_top_n,
        cluster_de_groupby=cluster_de_groupby,
        cluster_de_top_n=cluster_de_top_n,
        cluster_de_method=cluster_de_method,
        cluster_de_layer=cluster_de_layer,
        cluster_de_min_cells=cluster_de_min_cells,
        neighbor_stats_groupby=neighbor_stats_groupby,
        neighbor_stats_permutations=neighbor_stats_permutations,
        neighbor_stats_seed=neighbor_stats_seed,
        interaction_markers_groupby=interaction_markers_groupby,
        interaction_markers_top_targets=interaction_markers_top_targets,
        interaction_markers_top_genes=interaction_markers_top_genes,
        interaction_markers_min_cells=interaction_markers_min_cells,
        interaction_markers_min_neighbors=interaction_markers_min_neighbors,
        interaction_markers_method=interaction_markers_method,
        interaction_markers_layer=interaction_markers_layer,
        section_rotations=resolved_section_rotations,
    )
    if gene_storage == "sidecar":
        assert resolved_gene_aux_path is not None
        data["available_genes"] = list(dataset.var_names)
        data["gene_aux_url"] = Path(
            os.path.relpath(resolved_gene_aux_path, start=requested_output_path.resolve().parent)
        ).as_posix()
        gene_aux_url = data["gene_aux_url"]
    else:
        data["gene_aux_url"] = None

    if int(gene_correlation_top_n) > 0 and embedded_genes:
        if "gene_correlations" in companion_analytics:
            print("Using KaroSpaceCompanion gene correlations.")
            data["gene_correlations"] = companion_analytics["gene_correlations"]
        else:
            data["gene_correlations"] = _compute_gene_correlations(
                dataset.adata, embedded_genes, top_n=int(gene_correlation_top_n)
            )
    else:
        data["gene_correlations"] = {}

    if int(spatial_variable_genes_n) > 0:
        if "spatial_variable_genes" in companion_analytics:
            print("Using KaroSpaceCompanion spatially variable genes.")
            data["spatial_variable_genes"] = companion_analytics["spatial_variable_genes"]
        else:
            print(f"  - computing Moran's I for up to {int(spatial_variable_genes_n)} spatially variable genes...")
            data["spatial_variable_genes"] = _compute_morans_i(
                dataset.adata, list(dataset.var_names), n_genes=int(spatial_variable_genes_n)
            )
    else:
        data["spatial_variable_genes"] = []

    if int(cluster_means_n_genes) > 0:
        if "cluster_gene_means" in companion_analytics:
            print("Using KaroSpaceCompanion cluster gene means.")
            data["cluster_gene_means"] = companion_analytics["cluster_gene_means"]
        else:
            cmeans_genes = _select_top_variable_genes(dataset.adata, int(cluster_means_n_genes))
            categorical_cols = [
                col for col in (data.get("available_colors") or [])
                if not (data.get("colors_meta") or {}).get(col, {}).get("is_continuous", True)
            ]
            if categorical_cols and cmeans_genes:
                print(f"  - precomputing cluster gene means ({len(cmeans_genes)} genes × {len(categorical_cols)} columns)...")
                cmeans_columns = {}
                for col in categorical_cols:
                    result = _compute_cluster_gene_means(dataset.adata, cmeans_genes, col)
                    if result:
                        cmeans_columns[col] = result
                data["cluster_gene_means"] = {"genes": cmeans_genes, "columns": cmeans_columns} if cmeans_columns else None
            else:
                data["cluster_gene_means"] = None
    else:
        data["cluster_gene_means"] = None

    resolved_spot_size, used_auto_spot_size = _resolve_spot_size(
        dataset=dataset,
        spot_size=spot_size,
        min_panel_size=min_panel_size,
    )

    # Theme settings
    theme_icon = "☀️" if theme == "dark" else "🌙"
    initial_theme = theme

    # Load logo for favicon
    logo_base64 = _load_logo_base64()
    if logo_base64:
        favicon_link = f'<link rel="icon" type="image/png" href="data:image/png;base64,{logo_base64}">'
    else:
        favicon_link = ""
    footer_logo = (
        '<div class="footer-logo">'
        '<span>KaroSpace</span>'
        '<a class="footer-link" href="https://github.com/christoffermattssonlangseth/karospace" '
        'target="_blank" rel="noopener noreferrer">GitHub</a>'
        '</div>'
    )

    # Generate HTML
    metadata_labels = {
        "last_score": "disease score",
        "last_day": "day of sacrifice",
    }
    max_panel_size = int(min_panel_size * 2)

    data_json_safe = json.dumps(data, separators=(',', ':')).replace("</", "<\\/")

    html = HTML_TEMPLATE.format(
        title=title,
        min_panel_size=min_panel_size,
        max_panel_size=max_panel_size,
        spot_size=resolved_spot_size,
        data_json=data_json_safe,
        palette_json=json.dumps(DEFAULT_CATEGORICAL_PALETTE),
        metadata_labels_json=json.dumps(metadata_labels),
        outline_by_json=json.dumps(outline_by),
        viewer_info_html_json=json.dumps(viewer_info_html),
        viewer_info_html=viewer_info_html_safe,
        theme_icon=theme_icon,
        initial_theme=initial_theme,
        favicon_link=favicon_link,
        footer_logo=footer_logo,
        **colors
    )

    # Write file
    output_path = str(requested_output_path.resolve())
    if resolved_gene_aux_path is not None:
        resolved_gene_aux_path.parent.mkdir(parents=True, exist_ok=True)
        assert resolved_gene_aux_dir is not None
        resolved_gene_aux_dir.mkdir(parents=True, exist_ok=True)
        total_sidecar_genes = len(sidecar_genes)
        shard_groups = _chunked(sidecar_genes, gene_sidecar_shard_size)
        total_shards = len(shard_groups)
        sidecar_t0 = time.perf_counter()
        progress = None
        if total_sidecar_genes:
            print(
                f"Building gene sidecar: {total_sidecar_genes} genes across "
                f"{total_shards} shard{'s' if total_shards != 1 else ''}..."
            )
            if tqdm is not None:
                progress = tqdm(
                    total=total_sidecar_genes,
                    desc="Gene sidecar",
                    unit="gene",
                    dynamic_ncols=True,
                )
        manifest = {
            "format": (
                "karospace-gene-sidecar-manifest-v3"
                if gene_sidecar_format == GENE_SIDECAR_FORMAT_BINARY_V1
                else "karospace-gene-sidecar-manifest-v2"
            ),
            "gene_sidecar_format": gene_sidecar_format,
            "shards": {},
            "genes_meta": {},
            "gene_encodings": {},
            "gene_value_encodings": {},
            "gene_to_shard": {},
        }
        if gene_sidecar_format == GENE_SIDECAR_FORMAT_BINARY_V1:
            manifest["section_order"] = [section.section_id for section in dataset.sections]
        output_parent = Path(output_path).resolve().parent
        section_cell_counts = {
            section_id: int(len(indices))
            for section_id, indices in (sidecar_export_indices or {}).items()
        }
        genes_written = 0
        for shard_idx, shard_genes in enumerate(shard_groups):
            shard_suffix = ".bin" if gene_sidecar_format == GENE_SIDECAR_FORMAT_BINARY_V1 else ".json"
            shard_filename = f"{shard_idx:03d}{shard_suffix}"
            shard_path = resolved_gene_aux_dir / shard_filename
            shard_rel = Path(os.path.relpath(shard_path, start=output_parent)).as_posix()
            shard_start = genes_written + 1
            shard_end = genes_written + len(shard_genes)
            shard_t0 = time.perf_counter()
            if shard_genes and progress is None:
                print(
                    f"  - shard {shard_idx + 1}/{total_shards}: genes {shard_start}-{shard_end} "
                    f"({shard_genes[0]} .. {shard_genes[-1]})"
                )
            elif shard_genes and progress is not None:
                progress.set_postfix_str(
                    f"shard {shard_idx + 1}/{total_shards} {shard_genes[0]}..{shard_genes[-1]}"
                )
            shard_data = dataset.to_gene_sidecar_data(
                genes=shard_genes,
                downsample=downsample,
                export_indices=sidecar_export_indices,
                gene_encoding=gene_encoding,
                gene_value_encoding=gene_value_encoding,
                gene_sparse_zero_threshold=gene_sparse_zero_threshold,
            )
            manifest["shards"][shard_rel] = shard_genes
            for gene in shard_genes:
                manifest["gene_to_shard"][gene] = shard_rel
                if gene in shard_data.get("genes_meta", {}):
                    manifest["genes_meta"][gene] = shard_data["genes_meta"][gene]
                if gene in shard_data.get("gene_encodings", {}):
                    manifest["gene_encodings"][gene] = shard_data["gene_encodings"][gene]
                if gene in shard_data.get("gene_value_encodings", {}):
                    manifest["gene_value_encodings"][gene] = shard_data["gene_value_encodings"][gene]
            if gene_sidecar_format == GENE_SIDECAR_FORMAT_BINARY_V1:
                _write_binary_gene_shard(
                    shard_path=shard_path,
                    shard_genes=shard_genes,
                    shard_data=shard_data,
                    section_order=manifest["section_order"],
                    section_cell_counts=section_cell_counts,
                )
            else:
                shard_data["format"] = "karospace-gene-sidecar-shard-v2"
                shard_data.pop("genes_meta", None)
                shard_data.pop("gene_encodings", None)
                shard_data.pop("gene_value_encodings", None)
                with open(shard_path, "w", encoding="utf-8") as f:
                    json.dump(shard_data, f, separators=(",", ":"))
            genes_written += len(shard_genes)
            shard_elapsed = time.perf_counter() - shard_t0
            if progress is not None:
                progress.update(len(shard_genes))
                progress.set_postfix_str(
                    f"shard {shard_idx + 1}/{total_shards} wrote {shard_filename} in {shard_elapsed:.1f}s"
                )
            else:
                print(
                    f"    wrote {shard_filename} ({genes_written}/{total_sidecar_genes} genes) "
                    f"in {shard_elapsed:.1f}s"
                )
        with open(resolved_gene_aux_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, separators=(",", ":"))
        if total_sidecar_genes:
            total_elapsed = time.perf_counter() - sidecar_t0
            if progress is not None:
                progress.close()
            print(f"Completed gene sidecar build in {total_elapsed:.1f}s")
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)

    if package_mode:
        assert package_output_path_obj is not None
        assert resolved_gene_aux_path is not None
        assert package_gene_manifest_name is not None
        _write_karospace_package(
            package_path=package_output_path_obj,
            source_root=Path(output_path).parent,
            entry_html=package_entry_html,
            gene_manifest_path=package_gene_manifest_name,
            gene_shard_dir=Path(package_gene_manifest_name).with_suffix("").as_posix(),
            title=title,
            n_sections=data["n_sections"],
            total_cells=data["total_cells"],
        )
        written_loader_path = None
        if package_loader_output_path is not None:
            written_loader_path = _write_karospace_package_loader(
                loader_path=package_loader_output_path,
            )
        print(f"Exported KaroSpace package to: {package_output_path_obj}")
        print(f"  - entry HTML: {package_entry_html}")
        print(f"  - packaged gene manifest: {package_gene_manifest_name}")
        print(f"  - packaged shard directory: {Path(package_gene_manifest_name).with_suffix('').as_posix()}")
        if written_loader_path is not None:
            print(f"  - local package loader: {written_loader_path}")
        else:
            print(f"  - local package loader: unavailable (template not found)")
    else:
        print(f"Exported HTML viewer to: {output_path}")
        if resolved_gene_aux_path is not None:
            print(f"Exported gene sidecar to: {resolved_gene_aux_path}")
            print(f"  - shard directory: {resolved_gene_aux_dir}")
    print(f"  - {data['n_sections']} sections")
    print(f"  - {data['total_cells']:,} cells")
    if used_auto_spot_size:
        print(f"  - spot size auto-resolved to {resolved_spot_size:.2f}")
    else:
        print(f"  - spot size {resolved_spot_size:.2f}")
    print(f"  - {len(data['available_colors'])} color options")
    if embedded_genes:
        print(f"  - {len(data['genes_meta'])} genes loaded")
        enc = data.get("gene_encodings") or {}
        if enc:
            n_sparse = sum(1 for v in enc.values() if v == "sparse")
            n_dense = sum(1 for v in enc.values() if v == "dense")
            print(f"  - gene encoding: {n_sparse} sparse, {n_dense} dense")
    if gene_aux_url:
        print(f"  - sidecar genes available via {gene_aux_url}")
        print(f"  - sidecar value encoding: {gene_value_encoding}")
    if package_temp_dir is not None:
        package_temp_dir.cleanup()
    return str(package_output_path_obj) if package_output_path_obj is not None else output_path
