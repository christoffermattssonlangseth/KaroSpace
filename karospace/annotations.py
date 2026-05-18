"""Utilities for integrating KaroSpace polygon annotations into AnnData."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, MutableMapping, Sequence

import numpy as np
import pandas as pd


def _load_annotation_payload(
    annotations: str | Path | Mapping[str, Any],
) -> MutableMapping[str, Any]:
    """Load polygon annotation payload from a JSON path or a dict-like object."""
    if isinstance(annotations, (str, Path)):
        with Path(annotations).expanduser().open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    elif isinstance(annotations, Mapping):
        payload = dict(annotations)
    else:
        raise TypeError("annotations must be a path or a mapping")

    polygons = payload.get("polygons")
    if not isinstance(polygons, list):
        raise ValueError("annotation payload must contain a list at key 'polygons'")
    return payload


def _normalize_indices(values: Sequence[Any] | None, n_obs: int) -> list[int]:
    """Normalize a sequence to unique, sorted, in-range integer indices."""
    if not values:
        return []

    resolved: list[int] = []
    for raw in values:
        if raw is None:
            continue
        try:
            idx = int(raw)
        except (TypeError, ValueError):
            continue
        if 0 <= idx < n_obs:
            resolved.append(idx)

    if not resolved:
        return []
    return sorted(set(resolved))


def _resolve_from_local_indices(
    adata,
    *,
    groupby: str,
    section_id: str,
    local_indices: Sequence[Any] | None,
) -> list[int]:
    """Resolve section-local indices to global obs indices using adata.obs[groupby]."""
    if groupby not in adata.obs.columns:
        return []

    groups = adata.obs[groupby].astype(str).to_numpy()
    section_positions = np.flatnonzero(groups == str(section_id))
    if section_positions.size == 0:
        return []

    local_unique = _normalize_indices(local_indices, int(section_positions.size))
    if not local_unique:
        return []

    return sorted(set(int(section_positions[i]) for i in local_unique))


def _columnarize_polygon_metadata(polygons: Sequence[Mapping[str, Any]]) -> dict[str, list[Any]]:
    """Convert polygon records to a writeable columnar structure for AnnData uns."""
    if not polygons:
        return {
            "id": [],
            "label": [],
            "section_id": [],
            "n_cells": [],
            "cell_global_indices": [],
            "vertices": [],
            "created_at": [],
            "color": [],
        }

    column_order = list(polygons[0].keys())
    columns: dict[str, list[Any]] = {key: [] for key in column_order}
    for polygon in polygons:
        for key in column_order:
            value = polygon.get(key)
            if isinstance(value, (list, dict)):
                value = json.dumps(value, ensure_ascii=False, separators=(",", ":"))
            elif value is None:
                value = ""
            columns[key].append(value)
    return columns


def integrate_polygon_annotations(
    adata,
    annotations: str | Path | Mapping[str, Any],
    *,
    label_key: str = "karospace_polygon_labels",
    count_key: str = "karospace_polygon_count",
    uns_key: str = "karospace_polygon_annotations",
    delimiter: str = "|",
) -> Any:
    """
    Integrate exported KaroSpace polygon annotations into an AnnData object.

    The function writes two obs columns:
    - ``label_key``: joined polygon labels per cell (string; NA if not annotated)
    - ``count_key``: number of polygons covering each cell

    It also stores the full payload and resolved polygon metadata in ``adata.uns[uns_key]``.

    Parameters
    ----------
    adata
        AnnData object to annotate.
    annotations
        Path to an exported annotation JSON file or an equivalent mapping.
    label_key
        obs column name for per-cell polygon labels.
    count_key
        obs column name for number of polygons per cell.
    uns_key
        uns key where the payload and resolved metadata are stored.
    delimiter
        Separator used when multiple polygon labels map to one cell.

    Returns
    -------
    AnnData
        The same AnnData object with updates applied.
    """
    payload = _load_annotation_payload(annotations)
    polygons = payload.get("polygons", [])
    groupby = payload.get("groupby")

    labels_by_cell: list[list[str]] = [[] for _ in range(adata.n_obs)]
    resolved_polygons: list[dict[str, Any]] = []

    for i, polygon in enumerate(polygons, start=1):
        if not isinstance(polygon, Mapping):
            continue

        polygon_id = polygon.get("id", i)
        section_id = str(polygon.get("section_id", ""))
        label = str(polygon.get("label") or f"Annotation {polygon_id}")

        global_indices = _normalize_indices(polygon.get("cell_global_indices"), adata.n_obs)

        if not global_indices and section_id and groupby:
            global_indices = _resolve_from_local_indices(
                adata,
                groupby=str(groupby),
                section_id=section_id,
                local_indices=polygon.get("cell_local_indices"),
            )

        for idx in global_indices:
            labels_by_cell[idx].append(label)

        resolved_polygons.append(
            {
                "id": polygon_id,
                "label": label,
                "section_id": section_id,
                "n_cells": len(global_indices),
                "cell_global_indices": global_indices,
                "vertices": polygon.get("vertices", []),
                "created_at": polygon.get("created_at"),
                "color": polygon.get("color"),
            }
        )

    label_values = [delimiter.join(labels) if labels else None for labels in labels_by_cell]
    count_values = [int(len(labels)) for labels in labels_by_cell]

    adata.obs[label_key] = pd.Series(label_values, index=adata.obs_names, dtype="object")
    adata.obs[count_key] = pd.Series(count_values, index=adata.obs_names, dtype="int32")

    adata.uns[uns_key] = {
        "format": payload.get("format"),
        "created_at": payload.get("created_at"),
        "groupby": groupby,
        "n_polygons": len(resolved_polygons),
        "polygons_storage": "columnar-json-v1",
        "polygons": _columnarize_polygon_metadata(resolved_polygons),
    }
    return adata
