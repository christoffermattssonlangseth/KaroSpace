"""Unit tests for the pure-Python analytics and data-loading helpers.

These cover the numeric category ordering, spatial-key resolution, and the
analytics functions (Moran's I, gene correlations, cluster means, variable-gene
selection) with small deterministic fixtures.
"""
from collections import Counter

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from karospace.data_loader import (
    _numeric_category_perm,
    _resolve_spatial_key,
    load_spatial_data,
)
from karospace.exporter import (
    _compute_cluster_gene_means,
    _compute_gene_correlations,
    _compute_morans_i,
    _select_top_variable_genes,
)


# --------------------------------------------------------------------------- #
# _numeric_category_perm
# --------------------------------------------------------------------------- #
def test_numeric_perm_orders_numeric_strings():
    cats = ["0", "1", "10", "11", "2"]
    order = _numeric_category_perm(cats)
    assert order is not None
    assert [cats[i] for i in order] == ["0", "1", "2", "10", "11"]


def test_numeric_perm_none_for_non_numeric():
    assert _numeric_category_perm(["low", "high", "med"]) is None


def test_numeric_perm_none_when_already_sorted():
    assert _numeric_category_perm(["0", "1", "2"]) is None


def test_numeric_perm_none_for_single_or_empty():
    assert _numeric_category_perm([]) is None
    assert _numeric_category_perm(["7"]) is None


def test_numeric_perm_handles_floats():
    cats = ["1.5", "0.2", "10"]
    order = _numeric_category_perm(cats)
    assert order is not None
    assert [cats[i] for i in order] == ["0.2", "1.5", "10"]


# --------------------------------------------------------------------------- #
# Numeric category ordering — end-to-end through to_json_data (the bug fix)
# --------------------------------------------------------------------------- #
def _write_leiden_h5ad(tmp_path, labels, palette=None, spatial_key="spatial"):
    n = len(labels)
    X = np.random.RandomState(0).poisson(1.0, size=(n, 4)).astype(float)
    obs = pd.DataFrame(
        {
            "leiden": pd.Categorical([str(v) for v in labels]),
            "sample": pd.Categorical(["s1"] * n),
        },
        index=[f"c{i}" for i in range(n)],
    )
    adata = AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(4)]))
    adata.obsm[spatial_key] = np.random.RandomState(1).rand(n, 2)
    if palette is not None:
        adata.uns["leiden_colors"] = palette
    path = tmp_path / "leiden.h5ad"
    adata.write(path)
    return path


def test_to_json_data_numeric_category_order_and_palette(tmp_path):
    labels = [str(i % 13) for i in range(26)]  # "0".."12"
    # Palette aligned to pandas' (lexical) categorical order; each color encodes
    # its integer label in the red byte so we can check alignment after reorder.
    panda_cats = sorted(set(labels))  # lexical: 0,1,10,11,12,2,...
    palette = [f"#{int(c):02x}0000" for c in panda_cats]

    path = _write_leiden_h5ad(tmp_path, labels, palette=palette)
    dataset = load_spatial_data(str(path), groupby="sample")
    data = dataset.to_json_data("leiden")

    meta = data["colors_meta"]["leiden"]
    cats = [str(c) for c in meta["categories"]]

    # Categories are numerically ordered, not lexical.
    assert cats == [str(i) for i in range(13)]

    # Palette followed the same reorder: color i still encodes label cats[i].
    pal = meta["palette"]
    assert all(int(pal[i][1:3], 16) == int(cats[i]) for i in range(len(cats)))

    # Per-cell codes still resolve to the correct labels (remap stayed consistent).
    embedded = Counter()
    for section in data["sections"]:
        for code in section["colors"].get("leiden") or []:
            if code is not None:
                embedded[cats[int(round(code))]] += 1
    assert dict(embedded) == dict(Counter(str(v) for v in labels))


# --------------------------------------------------------------------------- #
# _resolve_spatial_key
# --------------------------------------------------------------------------- #
def _adata_with_obsm(**obsm):
    adata = AnnData(X=np.zeros((3, 2), dtype=float))
    for key, arr in obsm.items():
        adata.obsm[key] = arr
    return adata


def test_resolve_spatial_key_present():
    adata = _adata_with_obsm(spatial=np.zeros((3, 2)))
    assert _resolve_spatial_key(adata, "spatial") == "spatial"


def test_resolve_spatial_key_falls_back_to_common_name():
    adata = _adata_with_obsm(X_spatial=np.zeros((3, 2)))
    assert _resolve_spatial_key(adata, "spatial") == "X_spatial"


def test_resolve_spatial_key_missing_raises_and_lists_keys():
    adata = _adata_with_obsm(X_umap=np.zeros((3, 2)))
    with pytest.raises(ValueError) as excinfo:
        _resolve_spatial_key(adata, "spatial")
    assert "X_umap" in str(excinfo.value)


def test_resolve_spatial_key_rejects_bad_shape():
    adata = _adata_with_obsm(spatial=np.zeros((3, 1)))  # only one column
    with pytest.raises(ValueError):
        _resolve_spatial_key(adata, "spatial")


def test_load_spatial_data_uses_fallback_key_end_to_end(tmp_path):
    path = _write_leiden_h5ad(tmp_path, [0, 1, 0, 1], spatial_key="X_spatial")
    dataset = load_spatial_data(str(path), groupby="sample")
    assert dataset.spatial_key == "X_spatial"
    # to_json_data must honor the resolved key (previously hardcoded "spatial").
    data = dataset.to_json_data("leiden")
    assert data["sections"]


# --------------------------------------------------------------------------- #
# _select_top_variable_genes
# --------------------------------------------------------------------------- #
def test_select_top_variable_genes_orders_by_variance():
    X = np.array(
        [[0.0, 5.0, 1.0], [10.0, 5.0, 1.0], [0.0, 5.0, 1.0], [10.0, 5.0, 1.0]]
    )
    adata = AnnData(X=X, var=pd.DataFrame(index=["G_hi", "G_const", "G_mid"]))
    assert _select_top_variable_genes(adata, 1) == ["G_hi"]


def test_select_top_variable_genes_sparse_matches_dense():
    dense = np.array([[0.0, 1.0], [10.0, 1.0], [0.0, 1.0], [10.0, 1.0]])
    adata = AnnData(X=sp.csr_matrix(dense), var=pd.DataFrame(index=["A", "B"]))
    assert _select_top_variable_genes(adata, 1) == ["A"]


# --------------------------------------------------------------------------- #
# _compute_cluster_gene_means
# --------------------------------------------------------------------------- #
def test_cluster_gene_means_values():
    obs = pd.DataFrame(
        {"ct": pd.Categorical(["a", "a", "b", "b"])},
        index=[f"c{i}" for i in range(4)],
    )
    X = np.array([[2.0, 0.0], [4.0, 0.0], [10.0, 1.0], [20.0, 1.0]])
    adata = AnnData(X=X, obs=obs, var=pd.DataFrame(index=["G1", "G2"]))
    res = _compute_cluster_gene_means(adata, ["G1", "G2"], "ct")
    assert res["categories"] == ["a", "b"]
    assert res["means"]["a"] == [3.0, 0.0]
    assert res["means"]["b"] == [15.0, 1.0]
    assert res["background"] == [9.0, 0.5]


def test_cluster_gene_means_missing_column_returns_none():
    adata = AnnData(X=np.zeros((2, 2)), var=pd.DataFrame(index=["A", "B"]))
    assert _compute_cluster_gene_means(adata, ["A"], "nope") is None


# --------------------------------------------------------------------------- #
# _compute_gene_correlations
# --------------------------------------------------------------------------- #
def test_gene_correlations_positive_only():
    base = np.random.RandomState(0).rand(50)
    X = np.column_stack([base, base * 2 + 1, -base])  # G2 +corr, G3 -corr
    adata = AnnData(X=X, var=pd.DataFrame(index=["G1", "G2", "G3"]))
    res = _compute_gene_correlations(adata, ["G1", "G2", "G3"], top_n=2)
    partners = {d["gene"] for d in res["G1"]}
    assert "G2" in partners
    assert "G3" not in partners  # negative correlation excluded
    g2_r = next(d["r"] for d in res["G1"] if d["gene"] == "G2")
    assert g2_r > 0.99


def test_gene_correlations_single_gene_returns_empty():
    adata = AnnData(X=np.zeros((5, 1)), var=pd.DataFrame(index=["G1"]))
    assert _compute_gene_correlations(adata, ["G1"]) == {"G1": []}


# --------------------------------------------------------------------------- #
# _compute_morans_i
# --------------------------------------------------------------------------- #
def test_morans_i_without_graph_returns_empty():
    adata = AnnData(X=np.zeros((4, 2)), var=pd.DataFrame(index=["A", "B"]))
    assert _compute_morans_i(adata, ["A", "B"]) == []


def test_morans_i_detects_spatial_structure():
    n = 20
    rows, cols = [], []
    for i in range(n - 1):  # 1D chain: i <-> i+1
        rows += [i, i + 1]
        cols += [i + 1, i]
    W = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))

    smooth = np.arange(n, dtype=float)  # gradient -> strong positive autocorr
    alt = np.array([i % 2 for i in range(n)], dtype=float)  # checkerboard -> negative
    X = np.column_stack([smooth, alt])
    adata = AnnData(X=X, var=pd.DataFrame(index=["smooth", "alt"]))
    adata.obsp["spatial_connectivities"] = W

    res = _compute_morans_i(adata, ["smooth", "alt"], n_genes=2)
    by_gene = {r["gene"]: r["I"] for r in res}
    assert by_gene["smooth"] > 0.5
    assert by_gene["alt"] < 0
    assert all(-1.0 <= v <= 1.0 for v in by_gene.values())


# --------------------------------------------------------------------------- #
# CLI integration
# --------------------------------------------------------------------------- #
def test_cli_accepts_spatial_key(monkeypatch, tmp_path):
    import karospace.cli as cli_module

    input_path = tmp_path / "input.h5ad"
    input_path.write_text("placeholder", encoding="utf-8")
    captured = {}

    def fake_load(path, groupby="sample_id", spatial_key="spatial"):
        captured["load"] = {"path": path, "groupby": groupby, "spatial_key": spatial_key}
        class MockDataset:
            metadata_columns = []
            def to_json_data(self, *args, **kwargs):
                return {}
        return MockDataset()

    def fake_export(dataset, **kwargs):
        captured["export_kwargs"] = kwargs
        return str(tmp_path / "viewer.html")

    monkeypatch.setattr("karospace.data_loader.load_spatial_data", fake_load)
    monkeypatch.setattr("karospace.exporter.export_to_html", fake_export)
    monkeypatch.setattr(
        cli_module.sys,
        "argv",
        [
            "karospace",
            str(input_path),
            "--spatial-key",
            "X_my_spatial",
        ],
    )

    cli_module.main()

    assert captured["load"]["spatial_key"] == "X_my_spatial"

