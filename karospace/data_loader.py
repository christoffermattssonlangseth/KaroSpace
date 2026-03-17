"""
Data loading utilities for spatial transcriptomics data.

Handles loading h5ad files with scanpy and extracting spatial coordinates,
gene expression, and metadata for visualization.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from pandas.api.types import CategoricalDtype
from scipy.sparse import issparse
import scipy.sparse as sp
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
import json
import os
import shutil
import tempfile

def _strip_null_encoded_h5ad_entries(src_path: str) -> Tuple[str, List[str]]:
    """Copy an h5ad file and remove null-encoded datasets unsupported by older anndata builds."""
    import h5py

    fd, tmp_path = tempfile.mkstemp(suffix=".h5ad")
    os.close(fd)
    shutil.copy2(src_path, tmp_path)

    removed_paths: List[str] = []
    with h5py.File(tmp_path, "r+") as handle:
        def _walk(group, prefix: str = "") -> None:
            for key in list(group.keys()):
                obj = group[key]
                path = f"{prefix}/{key}" if prefix else f"/{key}"
                if obj.attrs.get("encoding-type") == "null":
                    removed_paths.append(path)
                    del group[key]
                    continue
                if isinstance(obj, h5py.Group):
                    _walk(obj, path)

        _walk(handle)

    return tmp_path, removed_paths


def _read_h5ad_with_fallback(path: str) -> sc.AnnData:
    """Read h5ad, retrying with null-encoded uns entries stripped if needed."""
    try:
        return sc.read_h5ad(path)
    except Exception as exc:
        if "encoding_type='null'" not in str(exc):
            raise

        print("  Detected unsupported null-encoded H5AD fields; retrying with a sanitized temporary copy...")
        temp_path = None
        try:
            temp_path, removed_paths = _strip_null_encoded_h5ad_entries(path)
            if removed_paths:
                print(f"  Removed {len(removed_paths)} null field(s): {', '.join(removed_paths)}")
            return sc.read_h5ad(temp_path)
        finally:
            if temp_path and os.path.exists(temp_path):
                os.remove(temp_path)


@dataclass
class SectionData:
    """Data for a single tissue section."""
    section_id: str
    coordinates: np.ndarray  # (n_cells, 2) array of x, y coordinates
    metadata: Dict[str, str] = field(default_factory=dict)

    @property
    def n_cells(self) -> int:
        return self.coordinates.shape[0]

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """Return (xmin, xmax, ymin, ymax)."""
        return (
            float(self.coordinates[:, 0].min()),
            float(self.coordinates[:, 0].max()),
            float(self.coordinates[:, 1].min()),
            float(self.coordinates[:, 1].max()),
        )


@dataclass
class SpatialDataset:
    """Container for spatial transcriptomics dataset."""
    adata: sc.AnnData
    sections: List[SectionData]
    groupby: str
    obs_columns: List[str]
    var_names: List[str]
    metadata_columns: List[str]
    metadata_value_order: Optional[Dict[str, List[str]]] = None

    @property
    def n_sections(self) -> int:
        return len(self.sections)

    @property
    def n_cells(self) -> int:
        return self.adata.n_obs

    @property
    def has_umap(self) -> bool:
        """Check if UMAP coordinates are available."""
        return "X_umap" in self.adata.obsm

    def get_color_data(
        self,
        color: str,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None
    ) -> Tuple[np.ndarray, bool, Optional[List[str]]]:
        """
        Get color values for all cells.

        Parameters
        ----------
        color : str
            Column in obs or gene name
        vmin, vmax : float, optional
            Min/max for continuous data

        Returns
        -------
        values : np.ndarray
            Numeric values for each cell
        is_continuous : bool
            Whether data is continuous
        categories : list or None
            Category names if categorical, else None
        """
        if color in self.adata.obs.columns:
            col = self.adata.obs[color]
            if isinstance(col.dtype, CategoricalDtype):
                categories = list(col.cat.categories)
                values = col.cat.codes.to_numpy().astype(float)
                # Handle NaN codes (-1)
                values[values < 0] = np.nan
                return values, False, categories
            elif pd.api.types.is_numeric_dtype(col):
                values = col.to_numpy(dtype=float)
                return values, True, None
            else:
                # Convert to categorical
                cat = col.astype("category")
                categories = list(cat.cat.categories)
                values = cat.cat.codes.to_numpy().astype(float)
                values[values < 0] = np.nan
                return values, False, categories
        elif color in self.adata.var_names:
            # Gene expression (prefer normalized layer when available)
            gene_idx = self.adata.var_names.get_loc(color)
            expr_layer = None
            if "normalized" in self.adata.layers:
                expr_layer = self.adata.layers["normalized"]
            x = expr_layer[:, gene_idx] if expr_layer is not None else self.adata.X[:, gene_idx]
            if issparse(x):
                values = np.asarray(x.toarray()).ravel()
            else:
                values = np.asarray(x).ravel()
            return values, True, None
        else:
            raise KeyError(f"{color!r} not found in obs columns or var_names")

    def get_section_indices(self) -> Dict[str, np.ndarray]:
        """Get cell indices for each section."""
        indices = {}
        gvals = self.adata.obs[self.groupby].astype(str).to_numpy()
        for section in self.sections:
            indices[section.section_id] = np.flatnonzero(gvals == section.section_id)
        return indices

    def get_metadata_filters(self) -> Dict[str, List[str]]:
        """Get unique values for filterable metadata columns."""
        filters = {}
        for col in self.metadata_columns:
            if col in self.adata.obs.columns:
                unique_vals = list(self.adata.obs[col].dropna().astype(str).unique())
                custom_order = None
                if self.metadata_value_order and col in self.metadata_value_order:
                    custom_order = [str(v) for v in self.metadata_value_order[col]]
                if custom_order:
                    custom_set = set(custom_order)
                    ordered = [v for v in custom_order if v in unique_vals]
                    remaining = [v for v in unique_vals if v not in custom_set]
                    if col == "last_day":
                        def _sort_key(v):
                            try:
                                return (0, float(v))
                            except ValueError:
                                return (1, v)
                        remaining = sorted(remaining, key=_sort_key)
                    else:
                        remaining = sorted(remaining)
                    filters[col] = ordered + remaining
                elif col == "last_day":
                    def _sort_key(v):
                        try:
                            return (0, float(v))
                        except ValueError:
                            return (1, v)
                    filters[col] = sorted(unique_vals, key=_sort_key)
                else:
                    filters[col] = sorted(unique_vals)
        return filters

    def _get_export_section_indices(self, downsample: Optional[int] = None) -> Dict[str, np.ndarray]:
        """Return per-section obs indices used for export, including deterministic downsampling."""
        section_indices = self.get_section_indices()
        export_indices: Dict[str, np.ndarray] = {}
        for section in self.sections:
            idx = section_indices[section.section_id]
            if downsample and len(idx) > downsample:
                rng = np.random.default_rng(42)
                idx = rng.choice(idx, size=downsample, replace=False)
                idx = np.sort(idx)
            export_indices[section.section_id] = np.asarray(idx)
        return export_indices

    def _collect_gene_data(self, genes: Optional[List[str]] = None) -> Dict[str, Dict[str, Union[np.ndarray, float]]]:
        """Collect gene vectors and ranges for export."""
        gene_data: Dict[str, Dict[str, Union[np.ndarray, float]]] = {}
        for gene in genes or []:
            if gene not in self.adata.var_names:
                continue
            try:
                vals, _, _ = self.get_color_data(gene)
                finite = np.isfinite(vals)
                gene_vmin = float(np.nanmin(vals[finite])) if finite.any() else 0.0
                gene_vmax = float(np.nanmax(vals[finite])) if finite.any() else 1.0
                gene_data[gene] = {
                    "values": vals,
                    "vmin": gene_vmin,
                    "vmax": gene_vmax,
                }
            except Exception as e:
                print(f"  Warning: Could not load gene '{gene}': {e}")
        return gene_data

    @staticmethod
    def _resolve_gene_encodings(
        gene_data: Dict[str, Dict[str, Union[np.ndarray, float]]],
        gene_encoding: str,
        gene_sparse_zero_threshold: float,
    ) -> Dict[str, str]:
        gene_encodings: Dict[str, str] = {}
        for gene, gdata in gene_data.items():
            if gene_encoding == "dense":
                gene_encodings[gene] = "dense"
            elif gene_encoding == "sparse":
                gene_encodings[gene] = "sparse"
            else:
                vals = np.asarray(gdata["values"])
                finite = np.isfinite(vals)
                nonzero = finite & (vals != 0)
                zero_frac = 1.0
                if vals.size:
                    zero_frac = 1.0 - (float(np.count_nonzero(nonzero)) / float(vals.size))
                gene_encodings[gene] = "sparse" if zero_frac >= float(gene_sparse_zero_threshold) else "dense"
        return gene_encodings

    @staticmethod
    def _validate_gene_export_options(
        gene_encoding: str,
        gene_sparse_zero_threshold: float,
        gene_sparse_pack_min_nnz: int,
    ) -> None:
        gene_encoding = str(gene_encoding or "auto").lower()
        if gene_encoding not in {"auto", "dense", "sparse"}:
            raise ValueError("gene_encoding must be one of: 'auto', 'dense', 'sparse'")
        if not (0.0 <= float(gene_sparse_zero_threshold) <= 1.0):
            raise ValueError("gene_sparse_zero_threshold must be between 0 and 1")
        if int(gene_sparse_pack_min_nnz) < 0:
            raise ValueError("gene_sparse_pack_min_nnz must be >= 0")

    @staticmethod
    def _serialize_gene_section_values(
        section_vals: np.ndarray,
        mode: str,
        gene_sparse_pack: bool,
        gene_sparse_pack_min_nnz: int,
        b64_encoder,
    ) -> Dict:
        if mode == "sparse":
            finite = np.isfinite(section_vals)
            nonzero = finite & (section_vals != 0)
            nz_idx = np.flatnonzero(nonzero).astype(np.uint32)
            nz_vals = np.asarray(section_vals[nonzero], dtype=np.float32)
            if bool(gene_sparse_pack) and int(nz_idx.size) >= int(gene_sparse_pack_min_nnz):
                sparse_entry = {
                    "ib64": b64_encoder(np.asarray(nz_idx, dtype="<u4")),
                    "vb64": b64_encoder(np.asarray(nz_vals, dtype="<f4")),
                }
            else:
                sparse_entry = {
                    "i": nz_idx.astype(int).tolist(),
                    "v": nz_vals.astype(float).tolist(),
                }
            nan_idx = np.flatnonzero(np.isnan(section_vals)).astype(int)
            if nan_idx.size:
                sparse_entry["nan"] = nan_idx.tolist()
            return {"sparse": sparse_entry}

        return {
            "dense": [float(v) if np.isfinite(v) else None for v in section_vals]
        }

    def to_gene_sidecar_data(
        self,
        genes: Optional[List[str]] = None,
        downsample: Optional[int] = None,
        export_indices: Optional[Dict[str, np.ndarray]] = None,
        gene_encoding: str = "auto",
        gene_sparse_zero_threshold: float = 0.8,
        gene_sparse_pack: bool = True,
        gene_sparse_pack_min_nnz: int = 256,
    ) -> Dict:
        """Export a gene-major sidecar payload for downstream viewer loading."""
        self._validate_gene_export_options(gene_encoding, gene_sparse_zero_threshold, gene_sparse_pack_min_nnz)
        export_indices = export_indices or self._get_export_section_indices(downsample=downsample)
        unique_genes = []
        seen = set()
        for gene in genes or []:
            if gene in seen or gene not in self.adata.var_names:
                continue
            seen.add(gene)
            unique_genes.append(gene)

        def _b64(arr: np.ndarray) -> str:
            import base64

            carr = np.ascontiguousarray(arr)
            return base64.b64encode(carr.tobytes(order="C")).decode("ascii")

        genes_payload = {}
        genes_meta = {}
        gene_encodings: Dict[str, str] = {}
        if unique_genes:
            gene_indices = [self.adata.var_names.get_loc(gene) for gene in unique_genes]
            expr_layer = self.adata.layers["normalized"] if "normalized" in self.adata.layers else self.adata.X
            batch = expr_layer[:, gene_indices]
            resolved_gene_encoding = str(gene_encoding or "auto").lower()

            if issparse(batch):
                batch_csr = batch.tocsr()
                batch_csc = batch.tocsc()
                section_batches = {
                    section.section_id: batch_csr[export_indices[section.section_id], :].tocsc()
                    for section in self.sections
                }
                n_obs = int(batch.shape[0])
                for gene_pos, gene in enumerate(unique_genes):
                    col = batch_csc.getcol(gene_pos)
                    raw_vals = np.asarray(col.data, dtype=float).ravel()
                    finite_vals = raw_vals[np.isfinite(raw_vals)]
                    if finite_vals.size:
                        raw_min = float(np.min(finite_vals))
                        raw_max = float(np.max(finite_vals))
                        has_implicit_zeros = int(np.count_nonzero(col.data != 0)) < n_obs
                        gene_vmin = min(0.0, raw_min) if has_implicit_zeros else raw_min
                        gene_vmax = max(0.0, raw_max) if has_implicit_zeros else raw_max
                    else:
                        gene_vmin = 0.0
                        gene_vmax = 0.0

                    if resolved_gene_encoding == "dense":
                        mode = "dense"
                    elif resolved_gene_encoding == "sparse":
                        mode = "sparse"
                    else:
                        nonzero = np.isfinite(raw_vals) & (raw_vals != 0)
                        zero_frac = 1.0
                        if n_obs:
                            zero_frac = 1.0 - (float(np.count_nonzero(nonzero)) / float(n_obs))
                        mode = "sparse" if zero_frac >= float(gene_sparse_zero_threshold) else "dense"
                    gene_encodings[gene] = mode

                    sections_payload = {}
                    for section in self.sections:
                        sec_col = section_batches[section.section_id].getcol(gene_pos)
                        if mode == "sparse":
                            section_vals = np.asarray(sec_col.data, dtype=float).ravel()
                            section_idx = np.asarray(sec_col.indices, dtype=np.int64).ravel()
                            finite = np.isfinite(section_vals)
                            nonzero = finite & (section_vals != 0)
                            nz_idx = section_idx[nonzero].astype(np.uint32, copy=False)
                            nz_vals = section_vals[nonzero].astype(np.float32, copy=False)
                            if bool(gene_sparse_pack) and int(nz_idx.size) >= int(gene_sparse_pack_min_nnz):
                                sparse_entry = {
                                    "ib64": _b64(np.asarray(nz_idx, dtype="<u4")),
                                    "vb64": _b64(np.asarray(nz_vals, dtype="<f4")),
                                }
                            else:
                                sparse_entry = {
                                    "i": nz_idx.astype(int).tolist(),
                                    "v": nz_vals.astype(float).tolist(),
                                }
                            nan_idx = section_idx[np.isnan(section_vals)].astype(int, copy=False)
                            if nan_idx.size:
                                sparse_entry["nan"] = nan_idx.tolist()
                            sections_payload[section.section_id] = {"sparse": sparse_entry}
                        else:
                            dense_vals = np.asarray(sec_col.toarray()).ravel()
                            sections_payload[section.section_id] = {
                                "dense": [
                                    float(v) if np.isfinite(v) else None for v in dense_vals
                                ]
                            }

                    genes_payload[gene] = {"sections": sections_payload}
                    genes_meta[gene] = {"vmin": gene_vmin, "vmax": gene_vmax}
            else:
                batch_dense = np.asarray(batch)
                if batch_dense.ndim == 1:
                    batch_dense = batch_dense.reshape(-1, 1)
                for gene_pos, gene in enumerate(unique_genes):
                    vals = np.asarray(batch_dense[:, gene_pos], dtype=float).ravel()
                    finite = np.isfinite(vals)
                    gene_vmin = float(np.nanmin(vals[finite])) if finite.any() else 0.0
                    gene_vmax = float(np.nanmax(vals[finite])) if finite.any() else 0.0
                    if resolved_gene_encoding == "dense":
                        mode = "dense"
                    elif resolved_gene_encoding == "sparse":
                        mode = "sparse"
                    else:
                        nonzero = finite & (vals != 0)
                        zero_frac = 1.0
                        if vals.size:
                            zero_frac = 1.0 - (float(np.count_nonzero(nonzero)) / float(vals.size))
                        mode = "sparse" if zero_frac >= float(gene_sparse_zero_threshold) else "dense"
                    gene_encodings[gene] = mode

                    sections_payload = {}
                    for section in self.sections:
                        idx = export_indices[section.section_id]
                        section_vals = batch_dense[idx, gene_pos]
                        sections_payload[section.section_id] = self._serialize_gene_section_values(
                            section_vals=np.asarray(section_vals),
                            mode=mode,
                            gene_sparse_pack=gene_sparse_pack,
                            gene_sparse_pack_min_nnz=gene_sparse_pack_min_nnz,
                            b64_encoder=_b64,
                        )
                    genes_payload[gene] = {"sections": sections_payload}
                    genes_meta[gene] = {"vmin": gene_vmin, "vmax": gene_vmax}

        return {
            "format": "karospace-gene-sidecar-v1",
            "genes_meta": genes_meta,
            "gene_encodings": gene_encodings,
            "genes": genes_payload,
        }

    def to_json_data(
        self,
        color: str,
        downsample: Optional[int] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        additional_colors: Optional[List[str]] = None,
        genes: Optional[List[str]] = None,
        gene_encoding: str = "auto",
        gene_sparse_zero_threshold: float = 0.8,
        gene_sparse_pack: bool = True,
        gene_sparse_pack_min_nnz: int = 256,
        section_array_pack: bool = True,
        section_array_pack_min_len: int = 1024,
        marker_genes_groupby: Optional[List[str]] = None,
        marker_genes_top_n: int = 30,
        neighbor_stats_groupby: Optional[List[str]] = None,
        neighbor_stats_permutations: int = 0,
        neighbor_stats_seed: int = 0,
        interaction_markers_groupby: Optional[List[str]] = None,
        interaction_markers_top_targets: int = 8,
        interaction_markers_top_genes: int = 20,
        interaction_markers_min_cells: int = 30,
        interaction_markers_min_neighbors: int = 1,
        interaction_markers_method: str = "wilcoxon",
        interaction_markers_layer: Optional[str] = "normalized",
        section_rotations: Optional[Dict[str, float]] = None,
    ) -> Dict:
        """
        Export dataset to JSON-serializable format for the HTML viewer.

        Parameters
        ----------
        color : str
            Initial color column or gene
        downsample : int, optional
            If set, randomly downsample to this many cells per section
        vmin, vmax : float, optional
            Min/max for continuous color scale
        additional_colors : list, optional
            Additional obs columns to include for color switching
        genes : list, optional
            Gene names to include for expression visualization
        gene_encoding : str
            "dense", "sparse", or "auto" (default: "auto"). When "sparse" (or when
            "auto" decides to use sparse), per-section gene vectors are stored as
            (index, value) pairs for non-zero entries to reduce HTML size for
            zero-inflated expression matrices.
        gene_sparse_zero_threshold : float
            Only used when gene_encoding="auto". Use sparse encoding when the
            fraction of zeros is >= this threshold (default: 0.8).
        gene_sparse_pack : bool
            When using sparse gene encoding, store indices/values as base64 typed arrays
            (smaller + faster JSON parse for large datasets). Default: True.
        gene_sparse_pack_min_nnz : int
            Only pack sparse arrays when non-zero entries in a section are >= this value.
            Default: 256.
        section_array_pack : bool
            Pack large per-section numeric arrays (coordinates, colors, UMAP) as base64 typed arrays
            for smaller HTML and faster JSON parse. Default: True.
        section_array_pack_min_len : int
            Only pack per-section arrays when section cell count is >= this value. Default: 1024.
        marker_genes_groupby : list, optional
            Obs columns to compute marker genes for (categorical only)
        marker_genes_top_n : int
            Number of top marker genes to keep per group
        neighbor_stats_groupby : list, optional
            Obs columns to compute neighbor composition stats for (categorical only)
        neighbor_stats_permutations : int
            Number of permutations for neighbor enrichment z-scores (0 disables)
        neighbor_stats_seed : int
            Random seed used for neighbor permutations
        interaction_markers_groupby : list, optional
            Obs columns to compute contact-conditioned interaction markers for.
            For each source/target pair: source cells contacting target vs source cells not contacting target.
        interaction_markers_top_targets : int
            Number of target cell types to evaluate per source (ranked by z-score or edge count).
        interaction_markers_top_genes : int
            Number of top genes to keep per source-target interaction.
        interaction_markers_min_cells : int
            Minimum cells required in both contact+ and contact- groups.
        interaction_markers_min_neighbors : int
            Minimum number of target neighbors for a source cell to be labeled contact+.
        interaction_markers_method : str
            Differential expression method for scanpy.tl.rank_genes_groups (e.g. "wilcoxon", "t-test").
        interaction_markers_layer : str, optional
            AnnData layer to use for DE, if available (default: "normalized").

        Returns
        -------
        dict
            JSON-serializable data structure
        """
        coords = np.asarray(self.adata.obsm["spatial"])[:, :2]
        export_section_indices = self._get_export_section_indices(downsample=downsample)

        # Get UMAP coordinates if available
        umap_coords = None
        umap_bounds = None
        if self.has_umap:
            umap_coords = np.asarray(self.adata.obsm["X_umap"])[:, :2]
            # Compute global UMAP bounds for consistent scaling across all sections
            umap_bounds = {
                "xmin": float(umap_coords[:, 0].min()),
                "xmax": float(umap_coords[:, 0].max()),
                "ymin": float(umap_coords[:, 1].min()),
                "ymax": float(umap_coords[:, 1].max()),
            }

        # Get neighborhood graph if available
        neighbor_graph = None
        neighbor_graph_key = None
        for key in ("spatial_connectivities", "connectivities", "neighbors", "neighbor_graph"):
            if key in self.adata.obsp:
                neighbor_graph = self.adata.obsp[key]
                neighbor_graph_key = key
                break
        if neighbor_graph is not None:
            if not issparse(neighbor_graph):
                neighbor_graph = sp.csr_matrix(neighbor_graph)
            else:
                neighbor_graph = neighbor_graph.tocsr()

        # Get initial color data
        values, is_continuous, categories = self.get_color_data(color, vmin, vmax)

        # Compute global bounds for initial color
        if is_continuous:
            finite_mask = np.isfinite(values)
            if finite_mask.any():
                global_vmin = float(np.nanmin(values[finite_mask])) if vmin is None else vmin
                global_vmax = float(np.nanmax(values[finite_mask])) if vmax is None else vmax
            else:
                global_vmin, global_vmax = 0.0, 1.0
        else:
            global_vmin, global_vmax = None, None

        # Build list of all colors to export
        all_colors = [color]
        if additional_colors:
            all_colors.extend([c for c in additional_colors if c != color and c in self.obs_columns])

        # Pre-compute all color data
        color_data = {}
        for col in all_colors:
            try:
                vals, is_cont, cats = self.get_color_data(col)
                if is_cont:
                    finite = np.isfinite(vals)
                    col_vmin = float(np.nanmin(vals[finite])) if finite.any() else 0.0
                    col_vmax = float(np.nanmax(vals[finite])) if finite.any() else 1.0
                else:
                    col_vmin, col_vmax = None, None
                color_data[col] = {
                    "values": vals,
                    "is_continuous": is_cont,
                    "categories": cats,
                    "vmin": col_vmin,
                    "vmax": col_vmax,
                }
            except Exception as e:
                print(f"  Warning: Could not load color '{col}': {e}")

        # Pre-compute gene expression data
        gene_data = self._collect_gene_data(genes)

        gene_encoding = str(gene_encoding or "auto").lower()
        self._validate_gene_export_options(gene_encoding, gene_sparse_zero_threshold, gene_sparse_pack_min_nnz)
        if int(section_array_pack_min_len) < 0:
            raise ValueError("section_array_pack_min_len must be >= 0")
        if int(interaction_markers_top_targets) < 1:
            raise ValueError("interaction_markers_top_targets must be >= 1")
        if int(interaction_markers_top_genes) < 1:
            raise ValueError("interaction_markers_top_genes must be >= 1")
        if int(interaction_markers_min_cells) < 1:
            raise ValueError("interaction_markers_min_cells must be >= 1")
        if int(interaction_markers_min_neighbors) < 1:
            raise ValueError("interaction_markers_min_neighbors must be >= 1")

        gene_encodings = self._resolve_gene_encodings(gene_data, gene_encoding, gene_sparse_zero_threshold)

        def _b64(arr: np.ndarray) -> str:
            import base64
            carr = np.ascontiguousarray(arr)
            return base64.b64encode(carr.tobytes(order="C")).decode("ascii")

        # Prepare float32 views to avoid per-section dtype conversions when packing arrays.
        coords_f4 = np.asarray(coords, dtype=np.float32, order="C")
        umap_f4 = None
        if umap_coords is not None:
            umap_f4 = np.asarray(umap_coords, dtype=np.float32, order="C")

        if bool(section_array_pack):
            for col, cdata in color_data.items():
                if cdata.get("is_continuous"):
                    cdata["_values_f4"] = np.asarray(cdata["values"], dtype=np.float32, order="C")
                else:
                    # Categorical codes are already numeric; keep float32 for compatibility (NaN for missing).
                    cdata["_values_f4"] = np.asarray(cdata["values"], dtype=np.float32, order="C")

        # Get metadata filters
        metadata_filters = self.get_metadata_filters()

        # Compute marker genes for requested groupby columns
        marker_genes = {}
        if marker_genes_groupby:
            print(
                f"Computing marker genes for {len(marker_genes_groupby)} "
                f"groupby column{'s' if len(marker_genes_groupby) != 1 else ''}..."
            )
            for groupby in marker_genes_groupby:
                print(f"  - marker genes: {groupby}")
                if groupby not in self.adata.obs.columns:
                    print(f"  Warning: marker_genes groupby '{groupby}' not found in obs.")
                    continue
                col = self.adata.obs[groupby]
                if not isinstance(col.dtype, CategoricalDtype):
                    self.adata.obs[groupby] = col.astype("category")
                key_added = f"rank_genes_groups_{groupby}"
                alt_key_added = f"rank_genes_groups__{groupby}"
                existing_key = None
                if key_added in self.adata.uns:
                    existing_key = key_added
                elif alt_key_added in self.adata.uns:
                    existing_key = alt_key_added

                if existing_key is None:
                    try:
                        sc.tl.rank_genes_groups(
                            self.adata,
                            groupby=groupby,
                            reference="rest",
                            method="t-test",
                            pts=False,
                            key_added=key_added,
                        )
                    except Exception as e:
                        print(f"  Warning: Could not compute marker genes for '{groupby}': {e}")
                        continue

                rg = self.adata.uns.get(existing_key or key_added)
                if not rg:
                    print(f"  Warning: marker genes not found for '{groupby}'.")
                    continue

                names = rg.get("names")
                if names is None:
                    print(f"  Warning: marker genes missing names for '{groupby}'.")
                    continue

                if isinstance(names, pd.DataFrame):
                    marker_genes[groupby] = {
                        col_name: names[col_name].astype(str).tolist()[:marker_genes_top_n]
                        for col_name in names.columns
                    }
                elif isinstance(names, np.ndarray) and names.dtype.names:
                    marker_genes[groupby] = {
                        group: [str(x) for x in names[group][:marker_genes_top_n]]
                        for group in names.dtype.names
                    }
                else:
                    print(f"  Warning: Unrecognized marker gene format for '{groupby}'.")

        # Compute neighbor composition stats
        neighbor_stats = {}
        neighbor_stats_context = {}
        if neighbor_graph is not None and neighbor_stats_groupby:
            print(
                f"Computing neighbor stats for {len(neighbor_stats_groupby)} "
                f"groupby column{'s' if len(neighbor_stats_groupby) != 1 else ''}..."
            )
            for groupby in neighbor_stats_groupby:
                print(f"  - neighbor stats: {groupby}")
                if groupby not in self.adata.obs.columns:
                    print(f"  Warning: neighbor stats groupby '{groupby}' not found in obs.")
                    continue
                col = self.adata.obs[groupby]
                if pd.api.types.is_numeric_dtype(col):
                    print(f"  Warning: neighbor stats '{groupby}' is numeric; skipping.")
                    continue
                if not isinstance(col.dtype, CategoricalDtype):
                    col = col.astype("category")
                categories = list(col.cat.categories)
                codes = col.cat.codes.to_numpy()
                valid_mask = codes >= 0
                if not valid_mask.any():
                    print(f"  Warning: neighbor stats '{groupby}' has no valid categories.")
                    continue

                if valid_mask.all():
                    graph = neighbor_graph
                    labels = codes
                else:
                    valid_idx = np.flatnonzero(valid_mask)
                    graph = neighbor_graph[valid_idx][:, valid_idx]
                    labels = codes[valid_mask]

                n_cells = np.bincount(labels, minlength=len(categories)).astype(int)
                if graph is None or graph.shape[0] == 0:
                    print(f"  Warning: neighbor stats '{groupby}' has empty graph.")
                    continue

                onehot = sp.csr_matrix(
                    (np.ones(len(labels), dtype=float), (np.arange(len(labels)), labels)),
                    shape=(len(labels), len(categories)),
                )
                counts = onehot.T.dot(graph).dot(onehot)
                if issparse(counts):
                    counts = counts.toarray()
                counts = np.asarray(counts, dtype=float)
                row_sums = counts.sum(axis=1)
                mean_degree = np.zeros(len(categories), dtype=float)
                valid_cells = n_cells > 0
                mean_degree[valid_cells] = row_sums[valid_cells] / n_cells[valid_cells]

                zscore = None
                entry = {
                    "categories": categories,
                    "counts": counts.tolist(),
                    "n_cells": n_cells.tolist(),
                    "mean_degree": mean_degree.tolist(),
                }
                if neighbor_stats_permutations and neighbor_stats_permutations > 0:
                    rng = np.random.default_rng(int(neighbor_stats_seed))
                    perm_mean = np.zeros_like(counts, dtype=float)
                    perm_m2 = np.zeros_like(counts, dtype=float)
                    for i in range(int(neighbor_stats_permutations)):
                        perm_labels = rng.permutation(labels)
                        perm_onehot = sp.csr_matrix(
                            (np.ones(len(perm_labels), dtype=float), (np.arange(len(perm_labels)), perm_labels)),
                            shape=(len(perm_labels), len(categories)),
                        )
                        perm_counts = perm_onehot.T.dot(graph).dot(perm_onehot)
                        if issparse(perm_counts):
                            perm_counts = perm_counts.toarray()
                        perm_counts = np.asarray(perm_counts, dtype=float)
                        delta = perm_counts - perm_mean
                        perm_mean += delta / (i + 1)
                        perm_m2 += delta * (perm_counts - perm_mean)
                    if neighbor_stats_permutations > 1:
                        perm_var = perm_m2 / (neighbor_stats_permutations - 1)
                    else:
                        perm_var = np.zeros_like(counts, dtype=float)
                    perm_std = np.sqrt(perm_var)
                    zscore = np.zeros_like(counts, dtype=float)
                    valid_std = perm_std > 0
                    zscore[valid_std] = (counts[valid_std] - perm_mean[valid_std]) / perm_std[valid_std]
                    entry["perm_n"] = int(neighbor_stats_permutations)
                    entry["zscore"] = zscore.tolist()
                neighbor_stats[groupby] = entry
                neighbor_stats_context[groupby] = {
                    "categories": categories,
                    "labels": labels.astype(np.int32, copy=False),
                    "graph": graph.tocsr(),
                    "obs_idx": (
                        np.arange(self.adata.n_obs, dtype=np.int64)
                        if valid_mask.all()
                        else np.flatnonzero(valid_mask).astype(np.int64)
                    ),
                    "counts": counts,
                    "zscore": zscore,
                    "n_cells": n_cells,
                    "mean_degree": mean_degree,
                }

        # Compute contact-conditioned interaction markers:
        # for source S and target T, compare source cells contacting T vs source cells not contacting T.
        interaction_markers = {}
        if neighbor_graph is not None and interaction_markers_groupby:
            print(
                f"Computing interaction markers for {len(interaction_markers_groupby)} "
                f"groupby column{'s' if len(interaction_markers_groupby) != 1 else ''}..."
            )
            method = str(interaction_markers_method or "wilcoxon")
            top_targets = int(interaction_markers_top_targets)
            top_genes = int(interaction_markers_top_genes)
            min_cells = int(interaction_markers_min_cells)
            min_neighbors = int(interaction_markers_min_neighbors)
            de_layer = None
            if interaction_markers_layer:
                if interaction_markers_layer in self.adata.layers:
                    de_layer = str(interaction_markers_layer)
                else:
                    print(
                        f"  Warning: interaction markers layer '{interaction_markers_layer}' not found; "
                        "using adata.X."
                    )

            def _extract_group_values(
                obj: Union[pd.DataFrame, np.ndarray, List],
                group: str,
                n: int,
                cast=None,
            ) -> List:
                vals = []
                if obj is None:
                    return vals
                if isinstance(obj, pd.DataFrame):
                    if group in obj.columns:
                        vals = obj[group].tolist()
                    elif obj.shape[1] > 0:
                        vals = obj.iloc[:, 0].tolist()
                elif isinstance(obj, np.ndarray) and obj.dtype.names:
                    g = group if group in obj.dtype.names else obj.dtype.names[0]
                    vals = list(obj[g])
                else:
                    vals = list(np.asarray(obj).ravel())
                vals = vals[:n]
                if cast is None:
                    return vals
                out = []
                for v in vals:
                    try:
                        out.append(cast(v))
                    except Exception:
                        out.append(None)
                return out

            for groupby in interaction_markers_groupby:
                print(f"  - interaction markers: {groupby}")
                ctx = neighbor_stats_context.get(groupby)
                if ctx is None:
                    print(
                        f"  Warning: interaction markers '{groupby}' unavailable "
                        "(missing neighbor stats for this groupby)."
                    )
                    continue

                categories = [str(c) for c in ctx["categories"]]
                labels = np.asarray(ctx["labels"], dtype=np.int32)
                graph = ctx["graph"]
                obs_idx = np.asarray(ctx["obs_idx"], dtype=np.int64)
                counts = np.asarray(ctx["counts"], dtype=float)
                zscore = ctx.get("zscore")
                n_cells = np.asarray(ctx["n_cells"], dtype=int)

                if graph.shape[0] != len(labels):
                    print(
                        f"  Warning: interaction markers '{groupby}' graph/label size mismatch; skipping."
                    )
                    continue

                group_interactions = {}
                for source_idx, source_name in enumerate(categories):
                    if source_idx >= len(n_cells) or int(n_cells[source_idx]) <= 0:
                        continue
                    source_mask = labels == source_idx
                    if not source_mask.any():
                        continue
                    row = counts[source_idx] if source_idx < counts.shape[0] else None
                    if row is None:
                        continue

                    candidate_targets = [
                        t for t in range(len(categories))
                        if t != source_idx and t < len(row) and float(row[t]) > 0
                    ]
                    if not candidate_targets:
                        continue

                    def _target_sort_key(tidx: int):
                        zval = None
                        if isinstance(zscore, np.ndarray):
                            if source_idx < zscore.shape[0] and tidx < zscore.shape[1]:
                                zval = zscore[source_idx, tidx]
                        if zval is not None and np.isfinite(zval):
                            return (0, -float(zval), -float(row[tidx]), categories[tidx])
                        return (1, 0.0, -float(row[tidx]), categories[tidx])

                    ranked_targets = sorted(candidate_targets, key=_target_sort_key)[:top_targets]
                    source_result = {}

                    for target_idx in ranked_targets:
                        target_name = categories[target_idx]
                        target_vec = (labels == target_idx).astype(np.float32, copy=False)
                        target_neighbor_counts = np.asarray(graph.dot(target_vec)).ravel()

                        pos_mask = source_mask & (target_neighbor_counts >= min_neighbors)
                        neg_mask = source_mask & (target_neighbor_counts == 0)
                        n_pos = int(np.count_nonzero(pos_mask))
                        n_neg = int(np.count_nonzero(neg_mask))
                        if n_pos < min_cells or n_neg < min_cells:
                            source_result[target_name] = {
                                "available": False,
                                "reason": "insufficient_cells",
                                "genes": [],
                                "logfoldchanges": [],
                                "pvals_adj": [],
                                "n_contact": n_pos,
                                "n_non_contact": n_neg,
                                "min_cells_required": min_cells,
                                "pct_contact": float((100.0 * n_pos) / max(1, n_pos + n_neg)),
                                "mean_target_neighbors_contact": float(
                                    np.mean(target_neighbor_counts[pos_mask]) if n_pos > 0 else 0.0
                                ),
                                "mean_target_neighbors_non_contact": float(
                                    np.mean(target_neighbor_counts[neg_mask]) if n_neg > 0 else 0.0
                                ),
                                "target_edge_count": float(row[target_idx]),
                                "target_zscore": (
                                    float(zscore[source_idx, target_idx])
                                    if isinstance(zscore, np.ndarray)
                                    and source_idx < zscore.shape[0]
                                    and target_idx < zscore.shape[1]
                                    and np.isfinite(zscore[source_idx, target_idx])
                                    else None
                                ),
                            }
                            continue

                        selected_mask = pos_mask | neg_mask
                        selected_idx = np.flatnonzero(selected_mask)
                        if selected_idx.size == 0:
                            continue

                        adata_idx = obs_idx[selected_idx]
                        try:
                            pair_adata = self.adata[adata_idx].copy()
                            contact_labels = np.where(pos_mask[selected_idx], "contact+", "contact-")
                            pair_adata.obs["_karospace_contact_group"] = pd.Categorical(
                                contact_labels,
                                categories=["contact+", "contact-"],
                            )

                            rg_kwargs = {
                                "groupby": "_karospace_contact_group",
                                "groups": ["contact+"],
                                "reference": "contact-",
                                "method": method,
                                "pts": False,
                                "key_added": "_karospace_interaction_markers",
                                "n_genes": top_genes,
                            }
                            if de_layer is not None:
                                rg_kwargs["layer"] = de_layer
                            sc.tl.rank_genes_groups(pair_adata, **rg_kwargs)
                            rg = pair_adata.uns.get("_karospace_interaction_markers", {})

                            genes = _extract_group_values(rg.get("names"), "contact+", top_genes, cast=str)
                            logfc = _extract_group_values(
                                rg.get("logfoldchanges"),
                                "contact+",
                                top_genes,
                                cast=lambda x: float(x) if np.isfinite(float(x)) else None,
                            )
                            pvals_adj = _extract_group_values(
                                rg.get("pvals_adj"),
                                "contact+",
                                top_genes,
                                cast=lambda x: float(x) if np.isfinite(float(x)) else None,
                            )

                            source_result[target_name] = {
                                "available": True,
                                "genes": genes,
                                "logfoldchanges": logfc,
                                "pvals_adj": pvals_adj,
                                "n_contact": n_pos,
                                "n_non_contact": n_neg,
                                "pct_contact": float((100.0 * n_pos) / max(1, n_pos + n_neg)),
                                "mean_target_neighbors_contact": float(
                                    np.mean(target_neighbor_counts[pos_mask]) if n_pos > 0 else 0.0
                                ),
                                "mean_target_neighbors_non_contact": float(
                                    np.mean(target_neighbor_counts[neg_mask]) if n_neg > 0 else 0.0
                                ),
                                "target_edge_count": float(row[target_idx]),
                                "target_zscore": (
                                    float(zscore[source_idx, target_idx])
                                    if isinstance(zscore, np.ndarray)
                                    and source_idx < zscore.shape[0]
                                    and target_idx < zscore.shape[1]
                                    and np.isfinite(zscore[source_idx, target_idx])
                                    else None
                                ),
                            }
                        except Exception as e:
                            print(
                                f"  Warning: interaction markers failed for '{groupby}' "
                                f"{source_name}->{target_name}: {e}"
                            )

                    if source_result:
                        group_interactions[source_name] = source_result

                if group_interactions:
                    interaction_markers[groupby] = group_interactions

        # Build section data with all color layers
        sections_data = []
        for section in self.sections:
            idx = export_section_indices[section.section_id]
            rotation_deg = float((section_rotations or {}).get(section.section_id, 0.0))

            section_coords = coords_f4[idx]

            # Get UMAP coordinates for this section if available
            section_umap = None
            if umap_f4 is not None:
                section_umap = umap_f4[idx]

            # Build color values for this section
            section_colors = {}
            section_colors_b64 = {}
            for col, cdata in color_data.items():
                if bool(section_array_pack) and int(len(idx)) >= int(section_array_pack_min_len):
                    section_vals = cdata.get("_values_f4", cdata["values"])[idx]
                    section_colors_b64[col] = _b64(section_vals.astype("<f4", copy=False))
                else:
                    section_vals = cdata["values"][idx]
                    # Convert numpy types to native Python types for JSON serialization
                    section_colors[col] = [
                        float(v) if np.isfinite(v) else None for v in section_vals
                    ]

            # Build gene expression values for this section
            section_genes_dense = {}
            section_genes_sparse = {}
            for gene, gdata in gene_data.items():
                section_vals = gdata["values"][idx]
                mode = gene_encodings.get(gene, "dense")
                payload = self._serialize_gene_section_values(
                    section_vals=np.asarray(section_vals),
                    mode=mode,
                    gene_sparse_pack=gene_sparse_pack,
                    gene_sparse_pack_min_nnz=gene_sparse_pack_min_nnz,
                    b64_encoder=_b64,
                )
                if "sparse" in payload:
                    section_genes_sparse[gene] = payload["sparse"]
                else:
                    section_genes_dense[gene] = payload["dense"]

            section_entry = {
                "id": section.section_id,
                "metadata": section.metadata,
                "n_cells": int(len(idx)),
                "rotation_deg": rotation_deg,
                "x": None,
                "y": None,
                "xb64": None,
                "yb64": None,
                "obs_idx": None,
                "obs_idxb64": None,
                "colors": section_colors,
                "colors_b64": section_colors_b64,
                "genes": section_genes_dense,
                "genes_sparse": section_genes_sparse,
                "bounds": {
                    "xmin": float(section_coords[:, 0].min()) if len(idx) > 0 else 0,
                    "xmax": float(section_coords[:, 0].max()) if len(idx) > 0 else 0,
                    "ymin": float(section_coords[:, 1].min()) if len(idx) > 0 else 0,
                    "ymax": float(section_coords[:, 1].max()) if len(idx) > 0 else 0,
                }
            }

            # Add UMAP coordinates if available
            if section_umap is not None:
                section_entry["umap_x"] = None
                section_entry["umap_y"] = None
                section_entry["umap_xb64"] = None
                section_entry["umap_yb64"] = None
            else:
                section_entry["umap_x"] = None
                section_entry["umap_y"] = None
                section_entry["umap_xb64"] = None
                section_entry["umap_yb64"] = None

            # Coordinates (pack when large)
            if bool(section_array_pack) and int(len(idx)) >= int(section_array_pack_min_len):
                section_entry["xb64"] = _b64(section_coords[:, 0].astype("<f4", copy=False))
                section_entry["yb64"] = _b64(section_coords[:, 1].astype("<f4", copy=False))
                section_entry["obs_idxb64"] = _b64(np.asarray(idx, dtype="<u4"))
                if section_umap is not None:
                    section_entry["umap_xb64"] = _b64(section_umap[:, 0].astype("<f4", copy=False))
                    section_entry["umap_yb64"] = _b64(section_umap[:, 1].astype("<f4", copy=False))
            else:
                section_entry["x"] = section_coords[:, 0].tolist()
                section_entry["y"] = section_coords[:, 1].tolist()
                section_entry["obs_idx"] = np.asarray(idx, dtype=int).tolist()
                if section_umap is not None:
                    section_entry["umap_x"] = section_umap[:, 0].tolist()
                    section_entry["umap_y"] = section_umap[:, 1].tolist()

            if neighbor_graph is not None:
                subgraph = neighbor_graph[idx][:, idx]
                if issparse(subgraph) and subgraph.nnz > 0:
                    upper = sp.triu(subgraph, k=1).tocoo()
                    rows = np.asarray(upper.row, dtype=np.uint32)
                    cols = np.asarray(upper.col, dtype=np.uint32)
                    if bool(section_array_pack) and int(len(idx)) >= int(section_array_pack_min_len):
                        pairs = np.empty(rows.size * 2, dtype=np.uint32)
                        pairs[0::2] = rows
                        pairs[1::2] = cols
                        section_entry["edges"] = None
                        section_entry["edges_b64"] = _b64(pairs.astype("<u4", copy=False))
                    else:
                        section_entry["edges"] = list(zip(rows.astype(int).tolist(), cols.astype(int).tolist()))
                        section_entry["edges_b64"] = None
                else:
                    section_entry["edges"] = []
                    section_entry["edges_b64"] = None

            sections_data.append(section_entry)

        # Build color metadata
        colors_meta = {}
        for col, cdata in color_data.items():
            colors_meta[col] = {
                "is_continuous": cdata["is_continuous"],
                "categories": cdata["categories"],
                "vmin": cdata["vmin"],
                "vmax": cdata["vmax"],
            }

        # Build gene metadata
        genes_meta = {}
        for gene, gdata in gene_data.items():
            genes_meta[gene] = {
                "vmin": gdata["vmin"],
                "vmax": gdata["vmax"],
            }

        return {
            "initial_color": color,
            "groupby": self.groupby,
            "colors_meta": colors_meta,
            "genes_meta": genes_meta,
            "gene_encodings": gene_encodings,
            "metadata_filters": metadata_filters,
            "n_sections": len(sections_data),
            "total_cells": sum(s["n_cells"] for s in sections_data),
            "sections": sections_data,
            "available_colors": list(color_data.keys()),
            "available_genes": list(gene_data.keys()),
            "marker_genes": marker_genes,
            "has_umap": umap_coords is not None,
            "umap_bounds": umap_bounds,
            "has_neighbors": neighbor_graph is not None,
            "neighbors_key": neighbor_graph_key,
            "neighbor_stats": neighbor_stats,
            "interaction_markers": interaction_markers,
        }


def load_spatial_data(
    path: str,
    groupby: str = "sample_id",
    spatial_key: str = "spatial",
    group_order: Optional[List[str]] = None,
    metadata_columns: Optional[List[str]] = None,
    metadata_value_order: Optional[Dict[str, List[str]]] = None,
    metadata_max_columns: Optional[int] = None,
) -> SpatialDataset:
    """
    Load spatial transcriptomics data from h5ad file.

    Parameters
    ----------
    path : str
        Path to .h5ad file
    groupby : str
        Column in obs to group sections by
    spatial_key : str
        Key in obsm containing spatial coordinates
    group_order : list, optional
        Custom order for sections
    metadata_columns : list, optional
        Obs columns to use for section metadata and filter chips
    metadata_value_order : dict, optional
        Custom ordering for metadata values per column (e.g. {"course": ["A", "B"]})
        If group_order is not provided, the first key in this dict is used to order sections
        by that metadata column (unknowns last, then section_id sort).
    metadata_max_columns : int, optional
        Limit the number of metadata columns used (order preserved)

    Returns
    -------
    SpatialDataset
        Loaded dataset ready for visualization
    """
    print(f"Loading {path}...")
    adata = _read_h5ad_with_fallback(path)
    print(f"  Loaded {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial coordinates not found in adata.obsm['{spatial_key}']")

    if groupby not in adata.obs.columns:
        raise ValueError(f"Groupby column '{groupby}' not found in adata.obs")

    # Determine section order
    gser = adata.obs[groupby]
    gser_str = gser.astype(str)
    if group_order is not None:
        section_ids = [str(g) for g in group_order if str(g) in gser_str.unique()]
    else:
        order_by_meta = None
        if metadata_value_order:
            order_by_meta = next(iter(metadata_value_order.keys()), None)
        if order_by_meta and order_by_meta in adata.obs.columns:
            desired_order = [str(v) for v in metadata_value_order.get(order_by_meta, [])]
            desired_index = {v: i for i, v in enumerate(desired_order)}
            section_ids = []
            for sid in gser_str.unique():
                mask = gser_str == str(sid)
                vals = adata.obs.loc[mask, order_by_meta].dropna().astype(str).unique()
                meta_value = vals[0] if len(vals) == 1 else "mixed"
                section_ids.append((str(sid), meta_value))
            def _order_key(item):
                sid, meta_value = item
                if meta_value in desired_index:
                    return (0, desired_index[meta_value], sid)
                return (1, meta_value, sid)
            section_ids = [sid for sid, _ in sorted(section_ids, key=_order_key)]
        elif isinstance(gser.dtype, CategoricalDtype) and gser.cat.ordered:
            section_ids = [str(c) for c in gser.cat.categories if str(c) in gser_str.unique()]
        else:
            section_ids = sorted(gser_str.unique())

    print(f"  Found {len(section_ids)} sections")

    # Determine metadata columns
    if metadata_columns is None:
        metadata_columns = ["course", "region", "condition", "timepoint", "last_score", "last_day"]
    if metadata_max_columns is not None:
        if metadata_max_columns < 0:
            raise ValueError("metadata_max_columns must be >= 0")
        metadata_columns = metadata_columns[:metadata_max_columns]

    # Build section data
    coords = np.asarray(adata.obsm[spatial_key])[:, :2]
    gvals = gser.astype(str).to_numpy()

    sections = []
    for sid in section_ids:
        mask = gvals == sid
        section_coords = coords[mask]

        # Extract metadata
        metadata = {}
        for meta_col in metadata_columns:
            if meta_col in adata.obs.columns:
                vals = adata.obs.loc[mask, meta_col].dropna().astype(str).unique()
                if len(vals) == 1:
                    metadata[meta_col] = vals[0]
                elif len(vals) > 1:
                    metadata[meta_col] = "mixed"

        sections.append(SectionData(
            section_id=sid,
            coordinates=section_coords,
            metadata=metadata,
        ))

    # Get available columns for coloring
    obs_columns = [
        col for col in adata.obs.columns
        if isinstance(adata.obs[col].dtype, CategoricalDtype)
        or pd.api.types.is_numeric_dtype(adata.obs[col])
    ]

    return SpatialDataset(
        adata=adata,
        sections=sections,
        groupby=groupby,
        obs_columns=obs_columns,
        var_names=list(adata.var_names),
        metadata_columns=metadata_columns,
        metadata_value_order=metadata_value_order,
    )
