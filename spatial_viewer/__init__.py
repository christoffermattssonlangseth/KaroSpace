"""
Spatial Viewer - HTML viewer for Xenium spatial transcriptomics data.

A tool to create interactive HTML visualizations of spatial transcriptomics
data stored in h5ad format, designed for large datasets with multiple sections.
"""

from .data_loader import load_spatial_data, SpatialDataset
from .exporter import export_to_html

__version__ = "0.1.0"
__all__ = ["load_spatial_data", "SpatialDataset", "export_to_html"]
