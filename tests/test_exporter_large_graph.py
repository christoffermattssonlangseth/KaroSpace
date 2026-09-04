import karospace.exporter as exporter


def test_oversized_section_neighbor_graph_is_omitted():
    original_limit = exporter.EMBEDDED_VIEWER_DATA_SINGLE_SCRIPT_MAX_CHARS
    exporter.EMBEDDED_VIEWER_DATA_SINGLE_SCRIPT_MAX_CHARS = 96
    try:
        primary_json, fragments = exporter._serialize_embedded_viewer_data(
            {
                "sections": [
                    {
                        "id": "section-a",
                        "n_cells": 2,
                        "edges": None,
                        "edges_b64": "x" * 128,
                    }
                ],
                "has_neighbors": True,
            }
        )
    finally:
        exporter.EMBEDDED_VIEWER_DATA_SINGLE_SCRIPT_MAX_CHARS = original_limit

    assert '"sections":[null]' in primary_json
    assert 'data-karospace-section-key="edges_b64">null</script>' in fragments
    assert "x" * 128 not in fragments


def test_oversized_non_graph_section_field_still_raises():
    original_limit = exporter.EMBEDDED_VIEWER_DATA_SINGLE_SCRIPT_MAX_CHARS
    exporter.EMBEDDED_VIEWER_DATA_SINGLE_SCRIPT_MAX_CHARS = 96
    try:
        try:
            exporter._serialize_embedded_viewer_data(
                {
                    "sections": [
                        {
                            "id": "section-a",
                            "xb64": "x" * 128,
                        }
                    ],
                }
            )
        except ValueError as exc:
            assert "section 0, 'xb64'" in str(exc)
        else:
            raise AssertionError("Expected oversized non-graph section field to raise")
    finally:
        exporter.EMBEDDED_VIEWER_DATA_SINGLE_SCRIPT_MAX_CHARS = original_limit
