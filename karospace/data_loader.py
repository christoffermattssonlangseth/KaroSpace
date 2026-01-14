"""
Data loading utilities for spatial transcriptomics data.

Handles loading h5ad files with scanpy and extracting spatial coordinates,
gene expression, and metadata for visualization.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
import json


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
            if pd.api.types.is_categorical_dtype(col):
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
            # Gene expression
            gene_idx = self.adata.var_names.get_loc(color)
            x = self.adata.X[:, gene_idx]
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
        for col in ["course", "region", "condition", "timepoint"]:
            if col in self.adata.obs.columns:
                unique_vals = self.adata.obs[col].dropna().astype(str).unique()
                filters[col] = sorted(unique_vals)
        return filters

    def to_json_data(
        self,
        color: str,
        downsample: Optional[int] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        additional_colors: Optional[List[str]] = None,
        genes: Optional[List[str]] = None,
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

        Returns
        -------
        dict
            JSON-serializable data structure
        """
        coords = np.asarray(self.adata.obsm["spatial"])[:, :2]
        section_indices = self.get_section_indices()

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
        gene_data = {}
        if genes:
            for gene in genes:
                if gene in self.adata.var_names:
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

        # Get metadata filters
        metadata_filters = self.get_metadata_filters()

        # Build section data with all color layers
        sections_data = []
        for section in self.sections:
            idx = section_indices[section.section_id]

            if downsample and len(idx) > downsample:
                rng = np.random.default_rng(42)
                idx = rng.choice(idx, size=downsample, replace=False)
                idx = np.sort(idx)

            section_coords = coords[idx]

            # Get UMAP coordinates for this section if available
            section_umap = None
            if umap_coords is not None:
                section_umap = umap_coords[idx]

            # Build color values for this section
            section_colors = {}
            for col, cdata in color_data.items():
                section_vals = cdata["values"][idx]
                # Convert numpy types to native Python types for JSON serialization
                section_colors[col] = [
                    float(v) if np.isfinite(v) else None for v in section_vals
                ]

            # Build gene expression values for this section
            section_genes = {}
            for gene, gdata in gene_data.items():
                section_vals = gdata["values"][idx]
                section_genes[gene] = [
                    float(v) if np.isfinite(v) else None for v in section_vals
                ]

            section_entry = {
                "id": section.section_id,
                "metadata": section.metadata,
                "n_cells": int(len(idx)),
                "x": section_coords[:, 0].tolist(),
                "y": section_coords[:, 1].tolist(),
                "colors": section_colors,
                "genes": section_genes,
                "bounds": {
                    "xmin": float(section_coords[:, 0].min()) if len(idx) > 0 else 0,
                    "xmax": float(section_coords[:, 0].max()) if len(idx) > 0 else 0,
                    "ymin": float(section_coords[:, 1].min()) if len(idx) > 0 else 0,
                    "ymax": float(section_coords[:, 1].max()) if len(idx) > 0 else 0,
                }
            }

            # Add UMAP coordinates if available
            if section_umap is not None:
                section_entry["umap_x"] = section_umap[:, 0].tolist()
                section_entry["umap_y"] = section_umap[:, 1].tolist()

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
            "colors_meta": colors_meta,
            "genes_meta": genes_meta,
            "metadata_filters": metadata_filters,
            "n_sections": len(sections_data),
            "total_cells": sum(s["n_cells"] for s in sections_data),
            "sections": sections_data,
            "available_colors": list(color_data.keys()),
            "available_genes": list(gene_data.keys()) if gene_data else self.var_names[:100],
            "all_genes": self.var_names[:500],  # For autocomplete
            "has_umap": umap_coords is not None,
            "umap_bounds": umap_bounds,
        }


def load_spatial_data(
    path: str,
    groupby: str = "sample_id",
    spatial_key: str = "spatial",
    group_order: Optional[List[str]] = None,
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

    Returns
    -------
    SpatialDataset
        Loaded dataset ready for visualization
    """
    print(f"Loading {path}...")
    adata = sc.read_h5ad(path)
    print(f"  Loaded {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial coordinates not found in adata.obsm['{spatial_key}']")

    if groupby not in adata.obs.columns:
        raise ValueError(f"Groupby column '{groupby}' not found in adata.obs")

    # Determine section order
    gser = adata.obs[groupby]
    if group_order is not None:
        section_ids = [str(g) for g in group_order if str(g) in gser.astype(str).unique()]
    elif pd.api.types.is_categorical_dtype(gser) and gser.cat.ordered:
        section_ids = [str(c) for c in gser.cat.categories if str(c) in gser.astype(str).unique()]
    else:
        section_ids = sorted(gser.astype(str).unique())

    print(f"  Found {len(section_ids)} sections")

    # Build section data
    coords = np.asarray(adata.obsm[spatial_key])[:, :2]
    gvals = gser.astype(str).to_numpy()

    sections = []
    for sid in section_ids:
        mask = gvals == sid
        section_coords = coords[mask]

        # Extract metadata
        metadata = {}
        for meta_col in ["region", "course", "condition", "timepoint"]:
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
        if pd.api.types.is_categorical_dtype(adata.obs[col])
        or pd.api.types.is_numeric_dtype(adata.obs[col])
    ]

    return SpatialDataset(
        adata=adata,
        sections=sections,
        groupby=groupby,
        obs_columns=obs_columns,
        var_names=list(adata.var_names),
    )
