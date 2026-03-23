"""
Spatial Viewer - HTML viewer for Xenium spatial transcriptomics data.

A tool to create interactive HTML visualizations of spatial transcriptomics
data stored in h5ad format, designed for large datasets with multiple sections.
"""

from importlib import import_module

__version__ = "0.1.0"
__all__ = [
    "load_spatial_data",
    "SpatialDataset",
    "export_to_html",
    "package_sidecar_viewer",
    "integrate_polygon_annotations",
]


def __getattr__(name):
    if name in {"load_spatial_data", "SpatialDataset"}:
        module = import_module(".data_loader", __name__)
        return getattr(module, name)
    if name in {"export_to_html", "package_sidecar_viewer"}:
        module = import_module(".exporter", __name__)
        return getattr(module, name)
    if name == "integrate_polygon_annotations":
        module = import_module(".annotations", __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))
