"""
Export spatial data to standalone HTML viewer.

Creates self-contained HTML files with embedded data and JavaScript
for interactive visualization of spatial transcriptomics data.
"""

import base64
import json
from pathlib import Path
from typing import Optional, List

from .data_loader import SpatialDataset


def _load_logo_base64() -> Optional[str]:
    """Load logo from assets as base64 string."""
    logo_path = Path(__file__).parent.parent / "assets" / "logo.png"
    if logo_path.exists():
        with open(logo_path, "rb") as f:
            return base64.b64encode(f.read()).decode("utf-8")
    return None


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
        .header h1 {{ font-size: 16px; font-weight: 600; }}
        .controls {{ display: flex; align-items: center; gap: 8px; flex-wrap: wrap; }}
        .control-group {{ display: flex; align-items: center; gap: 4px; }}
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
        .stats {{ font-size: 11px; color: var(--muted-color); }}

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

        .main-container {{ display: flex; flex: 1; min-height: 0; }}
        .grid-container {{
            flex: 1;
            padding: 8px;
            overflow: auto;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax({min_panel_size}px, 1fr));
            gap: 8px;
            align-content: start;
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
            font-size: 11px;
        }}
        .agg-group {{
            padding: 6px;
            border-radius: 6px;
            background: rgba(135, 0, 82, 0.06);
            border: 1px solid var(--border-color);
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
        .legend-color {{
            width: 12px;
            height: 12px;
            border-radius: 50%;
            flex-shrink: 0;
            border: 2px solid transparent;
        }}
        .legend-item.hidden .legend-color {{ border-color: var(--muted-color); background: transparent !important; }}
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
        .modal-header h2 {{ font-size: 15px; font-weight: 600; }}
        .modal-header .modal-meta {{ font-size: 11px; color: var(--muted-color); margin-left: 10px; }}
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
        .modal-canvas-container {{ flex: 1; position: relative; overflow: hidden; }}
        .modal-canvas {{ position: absolute; top: 0; left: 0; width: 100%; height: 100%; }}
        .modal-controls {{
            position: absolute;
            bottom: 12px;
            left: 50%;
            transform: translateX(-50%);
            display: flex;
            gap: 6px;
            background: var(--header-bg);
            padding: 6px 10px;
            border-radius: 6px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.2);
            transition: background 0.3s;
        }}
        .modal-controls button {{
            padding: 5px 10px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            font-size: 12px;
            transition: background 0.3s, border-color 0.3s, color 0.3s;
        }}
        .modal-controls button:hover {{ background: var(--hover-bg); }}
        .modal-legend {{ width: 180px; padding: 12px; border-left: 1px solid var(--border-color); overflow-y: auto; transition: border-color 0.3s; }}
        .zoom-info {{ font-size: 10px; color: var(--muted-color); margin-left: 6px; }}

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

        /* UMAP panel styles */
        .umap-panel {{
            width: 300px;
            background: var(--panel-bg);
            border-left: 1px solid var(--border-color);
            display: none;
            flex-direction: column;
            transition: background 0.3s, border-color 0.3s;
        }}
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
        .umap-canvas-container {{
            flex: 1;
            position: relative;
            min-height: 250px;
        }}
        .umap-canvas {{
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
        }}
        .umap-controls {{
            padding: 8px 12px;
            border-top: 1px solid var(--border-color);
            display: flex;
            flex-wrap: wrap;
            gap: 6px;
            align-items: center;
            transition: border-color 0.3s;
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
            font-size: 11px;
            color: var(--muted-color);
            padding: 4px 12px;
            border-top: 1px solid var(--border-color);
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
            opacity: 0.6;
            transition: opacity 0.2s;
            z-index: 100;
        }}
        .footer-logo:hover {{ opacity: 1; }}
        .footer-logo img {{
            height: 40px;
            width: auto;
        }}
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
        <h1>{title}</h1>
        <div class="controls">
            <div class="control-group">
                <label>Color:</label>
                <select id="color-select"></select>
            </div>
            <div class="control-group">
                <label>Gene:</label>
                <input type="text" id="gene-input" placeholder="e.g. Cd4, Gfap..." list="gene-list">
                <datalist id="gene-list"></datalist>
            </div>
            <div class="control-group">
                <label>Size:</label>
                <input type="range" id="spot-size" min="0.1" max="8" step="0.1" value="{spot_size}" style="width:80px">
            </div>
            <button class="umap-toggle" id="umap-toggle" title="Toggle UMAP view" style="display: none;">
                UMAP
            </button>
            <button class="legend-toggle active" id="legend-toggle" title="Toggle legend panel">
                Legend
            </button>
            <button class="color-toggle" id="color-toggle" title="Toggle color explorer">
                Colors
            </button>
            <button class="graph-toggle" id="graph-toggle" title="Toggle neighborhood graph" style="display: none;">
                Graph
            </button>
            <button class="graph-toggle" id="neighbor-hover-toggle" title="Toggle neighbor rings on hover" style="display: none;">
                Neighbors
            </button>
            <select id="neighbor-hop-select" title="Neighbor hop display" style="display: none; min-width: 90px;">
                <option value="1">1-hop</option>
                <option value="2">2-hop</option>
                <option value="3">3-hop</option>
                <option value="all" selected>All hops</option>
            </select>
            <button class="export-btn" id="screenshot-btn" title="Download screenshot">
                Screenshot
            </button>
            <button class="theme-toggle" id="theme-toggle" title="Toggle dark/light mode">
                <span id="theme-icon">{theme_icon}</span>
            </button>
        </div>
        <div class="stats"><span id="stats-text"></span></div>
    </div>

    <div class="filter-bar" id="filter-bar"></div>

    <div class="main-container">
        <div class="grid-container" id="grid"></div>
        <div class="umap-panel" id="umap-panel">
            <div class="umap-header">
                <h3>UMAP</h3>
            </div>
            <div class="umap-canvas-container" id="umap-canvas-container">
                <canvas class="umap-canvas" id="umap-canvas"></canvas>
            </div>
            <div class="umap-controls">
                <button class="umap-btn" id="magic-wand-btn" title="Draw to select cells">Magic Wand</button>
                <button class="umap-btn" id="clear-selection-btn" title="Clear selection">Clear</button>
                <span style="margin-left: 6px; font-size: 11px; color: var(--muted-color);">Size:</span>
                <input type="range" id="umap-spot-size" min="0.1" max="6" step="0.1" value="2" style="width: 60px;">
                <span id="umap-spot-size-label" style="font-size: 11px; min-width: 20px;">2</span>
            </div>
            <div class="umap-selection-info" id="umap-selection-info">No cells selected</div>
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
                </div>
                <button class="modal-close" id="modal-close">&times;</button>
            </div>
            <div class="modal-body">
                <div class="modal-canvas-container" id="modal-canvas-container">
                    <canvas class="modal-canvas" id="modal-canvas"></canvas>
                    <div class="modal-controls">
                        <button id="zoom-in">+ Zoom</button>
                        <button id="zoom-out">- Zoom</button>
                        <button id="zoom-reset">Reset</button>
                        <button class="graph-toggle" id="modal-graph-toggle" title="Toggle neighborhood graph" style="display: none;">Graph</button>
                        <button class="graph-toggle" id="modal-neighbor-hover-toggle" title="Toggle neighbor rings on hover" style="display: none;">Neighbors</button>
                        <select id="modal-neighbor-hop-select" title="Neighbor hop display" style="display: none; min-width: 90px;">
                            <option value="1">1-hop</option>
                            <option value="2">2-hop</option>
                            <option value="3">3-hop</option>
                            <option value="all" selected>All hops</option>
                        </select>
                        <span style="margin-left: 10px; font-size: 11px; color: {muted_color};">Size:</span>
                        <input type="range" id="modal-spot-size" min="0.1" max="12" step="0.1" value="{spot_size}" style="width: 80px;">
                        <span id="modal-spot-size-label" style="font-size: 11px; min-width: 24px;">{spot_size}</span>
                    </div>
                </div>
                <div class="modal-legend" id="modal-legend"></div>
            </div>
        </div>
    </div>

    <div class="cell-tooltip" id="cell-tooltip"></div>

    <script src="https://cdn.jsdelivr.net/npm/html2canvas@1.4.1/dist/html2canvas.min.js"></script>
    <script>
    const DATA = {data_json};
    const PALETTE = {palette_json};
    const METADATA_LABELS = {metadata_labels_json};
    const OUTLINE_BY = {outline_by_json};

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
    let hiddenCategories = new Set();
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

    // Modal state
    let modalSection = null;
    let modalZoom = 1;
    let modalPanX = 0, modalPanY = 0;
    let modalSpotSize = {spot_size};
    let isDragging = false;
    let dragStartX = 0, dragStartY = 0;
    let lastPanX = 0, lastPanY = 0;

    // UMAP state
    let umapVisible = false;
    let umapZoom = 1;
    let umapPanX = 0, umapPanY = 0;
    let umapSpotSize = 2;
    let isUmapDragging = false;
    let umapDragStartX = 0, umapDragStartY = 0;
    let umapLastPanX = 0, umapLastPanY = 0;

    // Selection state
    let magicWandActive = false;
    let isDrawingLasso = false;
    let lassoPath = [];  // Array of {{x, y}} points
    let selectedCells = new Set();  // Set of "sectionId:cellIdx" strings

    // Theme toggle
    function toggleTheme() {{
        currentTheme = currentTheme === 'light' ? 'dark' : 'light';
        document.documentElement.classList.remove('light', 'dark');
        document.documentElement.classList.add(currentTheme);
        document.getElementById('theme-icon').textContent = currentTheme === 'dark' ? '☀️' : '🌙';
        localStorage.setItem('spatial-viewer-theme', currentTheme);
        // Re-render canvases with new background
        renderAllSections();
        if (modalSection) renderModalSection();
        if (umapVisible) renderUMAP();
    }}

    function initTheme() {{
        // Check for saved preference or use initial theme
        const saved = localStorage.getItem('spatial-viewer-theme');
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

    function downloadCanvasImage(canvas, filename) {{
        if (!canvas) return;
        const link = document.createElement('a');
        link.href = canvas.toDataURL('image/png');
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        link.remove();
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
        if (typeof html2canvas !== 'function') {{
            alert('Screenshot library failed to load. Please check your connection and try again.');
            return;
        }}
        html2canvas(document.body, {{
            backgroundColor: null,
            scale: window.devicePixelRatio || 1,
            useCORS: true
        }}).then(canvas => {{
            downloadCanvasImage(canvas, name);
        }}).catch(() => {{
            alert('Screenshot failed to render.');
        }});
    }}

    // Color utilities
    function viridis(t) {{
        const colors = [
            [0.267, 0.005, 0.329], [0.282, 0.141, 0.458], [0.254, 0.265, 0.530],
            [0.207, 0.372, 0.553], [0.164, 0.471, 0.558], [0.128, 0.567, 0.551],
            [0.135, 0.659, 0.518], [0.267, 0.749, 0.441], [0.478, 0.821, 0.318],
            [0.741, 0.873, 0.150], [0.993, 0.906, 0.144]
        ];
        const idx = Math.min(Math.floor(t * 10), 9);
        const frac = (t * 10) - idx;
        const c1 = colors[idx], c2 = colors[idx + 1];
        const r = c1[0] + frac * (c2[0] - c1[0]);
        const g = c1[1] + frac * (c2[1] - c1[1]);
        const b = c1[2] + frac * (c2[2] - c1[2]);
        return `rgb(${{Math.round(r*255)}}, ${{Math.round(g*255)}}, ${{Math.round(b*255)}})`;
    }}

    function getCategoryColor(idx) {{ return PALETTE[idx % PALETTE.length]; }}

    function formatMetadataLabel(key) {{
        return METADATA_LABELS[key] || key.replace(/_/g, ' ');
    }}

    // Get current color config
    function getColorConfig() {{
        if (currentGene && DATA.genes_meta[currentGene]) {{
            return {{
                is_continuous: true,
                categories: null,
                vmin: DATA.genes_meta[currentGene].vmin,
                vmax: DATA.genes_meta[currentGene].vmax
            }};
        }}
        return DATA.colors_meta[currentColor] || {{ is_continuous: false, categories: [], vmin: 0, vmax: 1 }};
    }}

    // Get values for a section
    function getSectionValues(section) {{
        if (currentGene && section.genes[currentGene]) {{
            return section.genes[currentGene];
        }}
        return section.colors[currentColor] || [];
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
        const n = section.x.length;
        const adj = Array.from({{ length: n }}, () => []);
        if (section.edges) {{
            section.edges.forEach(edge => {{
                const i = edge[0];
                const j = edge[1];
                if (i >= 0 && j >= 0 && i < n && j < n) {{
                    adj[i].push(j);
                    adj[j].push(i);
                }}
            }});
        }}
        section._adj = adj;
        return adj;
    }}

    function computeNeighborRings(section, startIdx, maxHops) {{
        if (!section.edges || section.edges.length === 0) return [];
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
        if (!section || !section.edges || section.edges.length === 0) {{
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

        const config = getColorConfig();
        const values = getSectionValues(section);

        if (neighborHopMode !== 'all') {{
            const hopIdx = Math.max(1, Math.min(MAX_HOVER_HOPS, parseInt(neighborHopMode, 10))) - 1;
            rings = rings[hopIdx] ? [rings[hopIdx]] : [];
        }}

        const xCenter = transform.centerX + (section.x[centerIdx] - transform.dataCenterX) * transform.scale;
        const yCenter = transform.centerY - (section.y[centerIdx] - transform.dataCenterY) * transform.scale;

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
                const val = values[cellIdx];
                if (val === null || val === undefined) return;
                if (!config.is_continuous) {{
                    const catIdx = Math.round(val);
                    const catName = config.categories[catIdx];
                    if (hiddenCategories.has(catName)) return;
                }}
                const x = transform.centerX + (section.x[cellIdx] - transform.dataCenterX) * transform.scale;
                const y = transform.centerY - (section.y[cellIdx] - transform.dataCenterY) * transform.scale;
                if (x < -adjustedSpotSize || x > transform.width + adjustedSpotSize ||
                    y < -adjustedSpotSize || y > transform.height + adjustedSpotSize) return;
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
                const val = values[cellIdx];
                if (val === null || val === undefined) return;
                if (!config.is_continuous) {{
                    const catIdx = Math.round(val);
                    const catName = config.categories[catIdx];
                    if (hiddenCategories.has(catName)) return;
                }}
                const x = transform.centerX + (section.x[cellIdx] - transform.dataCenterX) * transform.scale;
                const y = transform.centerY - (section.y[cellIdx] - transform.dataCenterY) * transform.scale;
                if (x < -adjustedSpotSize || x > transform.width + adjustedSpotSize ||
                    y < -adjustedSpotSize || y > transform.height + adjustedSpotSize) return;
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

    // UMAP rendering
    function renderUMAP() {{
        if (!DATA.has_umap || !umapVisible) return;

        const canvas = document.getElementById('umap-canvas');
        const ctx = canvas.getContext('2d');
        const dpr = window.devicePixelRatio || 1;
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
        const adjustedSpotSize = Math.max(1, umapSpotSize * umapZoom * 0.5);

        // Check if any categories are hidden
        const hasHidden = hiddenCategories.size > 0 && !config.is_continuous;

        // First pass: draw hidden categories as grey (if any are hidden)
        if (hasHidden) {{
            ctx.fillStyle = '#888888';
            ctx.globalAlpha = 0.2;
            DATA.sections.forEach(section => {{
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
                if (config.is_continuous) {{
                    const t = (val - config.vmin) / (config.vmax - config.vmin);
                    color = viridis(Math.max(0, Math.min(1, t)));
                }} else {{
                    const catIdx = Math.round(val);
                    color = getCategoryColor(catIdx);
                }}

                ctx.fillStyle = color;
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

        // Clear previous selection or add to it (could add shift-key support later)
        selectedCells.clear();

        // Check all cells in all sections
        DATA.sections.forEach(section => {{
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
                    selectedCells.add(`${{section.id}}:${{i}}`);
                }}
            }}
        }});

        updateSelectionInfo();
        renderUMAP();
        renderAllSections();
        if (modalSection) renderModalSection();
    }}

    // Update selection info display
    function updateSelectionInfo() {{
        const info = document.getElementById('umap-selection-info');
        if (selectedCells.size === 0) {{
            info.textContent = 'No cells selected';
        }} else {{
            info.textContent = `${{selectedCells.size.toLocaleString()}} cells selected`;
        }}
    }}

    // Clear selection
    function clearSelection() {{
        selectedCells.clear();
        updateSelectionInfo();
        renderUMAP();
        renderAllSections();
        if (modalSection) renderModalSection();
    }}

    // Toggle UMAP panel
    function toggleUMAP() {{
        umapVisible = !umapVisible;
        const panel = document.getElementById('umap-panel');
        const btn = document.getElementById('umap-toggle');
        panel.classList.toggle('visible', umapVisible);
        btn.classList.toggle('active', umapVisible);
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
            document.getElementById('umap-spot-size-label').textContent = umapSpotSize;
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
    function renderSection(section, canvas) {{
        const ctx = canvas.getContext('2d');
        const dpr = window.devicePixelRatio || 1;
        const rect = canvas.getBoundingClientRect();
        canvas.width = rect.width * dpr;
        canvas.height = rect.height * dpr;
        ctx.scale(dpr, dpr);

        const width = rect.width, height = rect.height, padding = 8;
        ctx.fillStyle = getPanelBg();
        ctx.fillRect(0, 0, width, height);

        if (section.x.length === 0) return;

        const bounds = section.bounds;
        const dataWidth = bounds.xmax - bounds.xmin;
        const dataHeight = bounds.ymax - bounds.ymin;
        const scale = Math.min((width - 2*padding) / dataWidth, (height - 2*padding) / dataHeight);
        const offsetX = padding + ((width - 2*padding) - dataWidth * scale) / 2;
        const offsetY = padding + ((height - 2*padding) - dataHeight * scale) / 2;

        const config = getColorConfig();
        const values = getSectionValues(section);

        if (showGraph && section.edges && section.edges.length) {{
            const graphColor = getGraphColor();
            ctx.strokeStyle = graphColor;
            ctx.lineWidth = Math.max(0.3, spotSize * 0.15);
            ctx.beginPath();
            for (let e = 0; e < section.edges.length; e++) {{
                const edge = section.edges[e];
                const i = edge[0];
                const j = edge[1];
                const valI = values[i];
                const valJ = values[j];
                if (valI === null || valI === undefined || valJ === null || valJ === undefined) continue;
                if (!config.is_continuous) {{
                    const catIdxI = Math.round(valI);
                    const catIdxJ = Math.round(valJ);
                    const catNameI = config.categories[catIdxI];
                    const catNameJ = config.categories[catIdxJ];
                    if (hiddenCategories.has(catNameI) || hiddenCategories.has(catNameJ)) continue;
                }}
                const x1 = offsetX + (section.x[i] - bounds.xmin) * scale;
                const y1 = height - (offsetY + (section.y[i] - bounds.ymin) * scale);
                const x2 = offsetX + (section.x[j] - bounds.xmin) * scale;
                const y2 = height - (offsetY + (section.y[j] - bounds.ymin) * scale);
                ctx.moveTo(x1, y1);
                ctx.lineTo(x2, y2);
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

                const x = offsetX + (section.x[i] - bounds.xmin) * scale;
                const y = offsetY + (section.y[i] - bounds.ymin) * scale;
                ctx.beginPath();
                ctx.arc(x, height - y, spotSize, 0, Math.PI * 2);
                ctx.fill();
            }}
            ctx.globalAlpha = 1;
        }}

        // Second pass: draw visible categories on top
        for (let i = 0; i < section.x.length; i++) {{
            const val = values[i];
            if (val === null || val === undefined) continue;

            let color;
            if (config.is_continuous) {{
                const t = (val - config.vmin) / (config.vmax - config.vmin);
                color = viridis(Math.max(0, Math.min(1, t)));
            }} else {{
                const catIdx = Math.round(val);
                const catName = config.categories[catIdx];
                if (hiddenCategories.has(catName)) continue;  // Skip hidden, already drawn as grey
                color = getCategoryColor(catIdx);
            }}

            const x = offsetX + (section.x[i] - bounds.xmin) * scale;
            const y = offsetY + (section.y[i] - bounds.ymin) * scale;
            ctx.fillStyle = color;
            ctx.beginPath();
            ctx.arc(x, height - y, spotSize, 0, Math.PI * 2);
            ctx.fill();
        }}

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

                const x = offsetX + (section.x[i] - bounds.xmin) * scale;
                const y = offsetY + (section.y[i] - bounds.ymin) * scale;
                ctx.beginPath();
                ctx.arc(x, height - y, spotSize + 1, 0, Math.PI * 2);
                ctx.stroke();
            }}
        }}
    }}

    function renderAllSections() {{
        const panels = document.querySelectorAll('.section-panel');
        let visibleCount = 0;
        let totalCells = 0;

        DATA.sections.forEach((section, idx) => {{
            const panel = panels[idx];
            if (!panel) return;

            const passes = sectionPassesFilter(section);
            panel.classList.toggle('filtered-out', !passes);

            if (passes) {{
                visibleCount++;
                totalCells += section.n_cells;
                const canvas = panel.querySelector('canvas');
                if (canvas) renderSection(section, canvas);
            }}
        }});

        // Update stats
        const colorLabel = currentGene || currentColor;
        document.getElementById('stats-text').textContent =
            `${{visibleCount}}/${{DATA.n_sections}} sections | ${{totalCells.toLocaleString()}} cells | ${{colorLabel}}`;

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

    function findNearestCell(section, mouseX, mouseY, canvasRect, transform) {{
        // transform: {{ scale, offsetX, offsetY, centerX, centerY, dataCenterX, dataCenterY, isModal }}
        const config = getColorConfig();
        const values = getSectionValues(section);
        const searchRadius = transform.isModal ? modalSpotSize * modalZoom * 2 : spotSize * 3;

        let nearestIdx = -1;
        let nearestDist = Infinity;

        for (let i = 0; i < section.x.length; i++) {{
            const val = values[i];
            if (val === null || val === undefined) continue;

            // Skip hidden categories
            if (!config.is_continuous) {{
                const catIdx = Math.round(val);
                const catName = config.categories[catIdx];
                if (hiddenCategories.has(catName)) continue;
            }}

            let screenX, screenY;
            if (transform.isModal) {{
                screenX = transform.centerX + (section.x[i] - transform.dataCenterX) * transform.scale;
                screenY = transform.centerY - (section.y[i] - transform.dataCenterY) * transform.scale;
            }} else {{
                const bounds = section.bounds;
                screenX = transform.offsetX + (section.x[i] - bounds.xmin) * transform.scale;
                screenY = transform.height - (transform.offsetY + (section.y[i] - bounds.ymin) * transform.scale);
            }}

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
            const color = viridis(Math.max(0, Math.min(1, t)));
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

    // Modal rendering
    function renderModalSection() {{
        if (!modalSection) return;

        const canvas = document.getElementById('modal-canvas');
        const ctx = canvas.getContext('2d');
        const dpr = window.devicePixelRatio || 1;
        const container = document.getElementById('modal-canvas-container');
        const rect = container.getBoundingClientRect();
        canvas.width = rect.width * dpr;
        canvas.height = rect.height * dpr;
        ctx.scale(dpr, dpr);

        const width = rect.width, height = rect.height;
        ctx.fillStyle = getPanelBg();
        ctx.fillRect(0, 0, width, height);

        if (modalSection.x.length === 0) return;

        const bounds = modalSection.bounds;
        const dataWidth = bounds.xmax - bounds.xmin;
        const dataHeight = bounds.ymax - bounds.ymin;
        const baseScale = Math.min((width - 40) / dataWidth, (height - 40) / dataHeight);
        const scale = baseScale * modalZoom;

        const centerX = width / 2 + modalPanX;
        const centerY = height / 2 + modalPanY;
        const dataCenterX = (bounds.xmin + bounds.xmax) / 2;
        const dataCenterY = (bounds.ymin + bounds.ymax) / 2;
        const adjustedSpotSize = Math.max(1, modalSpotSize * modalZoom * 0.8);

        const config = getColorConfig();
        const values = getSectionValues(modalSection);

        if (showGraph && modalSection.edges && modalSection.edges.length) {{
            const graphColor = getGraphColor();
            ctx.strokeStyle = graphColor;
            ctx.lineWidth = Math.max(0.3, modalSpotSize * modalZoom * 0.12);
            ctx.beginPath();
            for (let e = 0; e < modalSection.edges.length; e++) {{
                const edge = modalSection.edges[e];
                const i = edge[0];
                const j = edge[1];
                const valI = values[i];
                const valJ = values[j];
                if (valI === null || valI === undefined || valJ === null || valJ === undefined) continue;
                if (!config.is_continuous) {{
                    const catIdxI = Math.round(valI);
                    const catIdxJ = Math.round(valJ);
                    const catNameI = config.categories[catIdxI];
                    const catNameJ = config.categories[catIdxJ];
                    if (hiddenCategories.has(catNameI) || hiddenCategories.has(catNameJ)) continue;
                }}
                const x1 = centerX + (modalSection.x[i] - dataCenterX) * scale;
                const y1 = centerY - (modalSection.y[i] - dataCenterY) * scale;
                const x2 = centerX + (modalSection.x[j] - dataCenterX) * scale;
                const y2 = centerY - (modalSection.y[j] - dataCenterY) * scale;
                if (x1 < -adjustedSpotSize || x1 > width + adjustedSpotSize ||
                    y1 < -adjustedSpotSize || y1 > height + adjustedSpotSize) continue;
                if (x2 < -adjustedSpotSize || x2 > width + adjustedSpotSize ||
                    y2 < -adjustedSpotSize || y2 > height + adjustedSpotSize) continue;
                ctx.moveTo(x1, y1);
                ctx.lineTo(x2, y2);
            }}
            ctx.stroke();
        }}

        // First pass: draw grey background for hidden categories
        const hasHidden = hiddenCategories.size > 0 && !config.is_continuous;
        if (hasHidden) {{
            ctx.fillStyle = '#cccccc';
            ctx.globalAlpha = 0.2;
            for (let i = 0; i < modalSection.x.length; i++) {{
                const val = values[i];
                if (val === null || val === undefined) continue;
                const catIdx = Math.round(val);
                const catName = config.categories[catIdx];
                if (!hiddenCategories.has(catName)) continue;

                const x = centerX + (modalSection.x[i] - dataCenterX) * scale;
                const y = centerY - (modalSection.y[i] - dataCenterY) * scale;

                if (x < -adjustedSpotSize || x > width + adjustedSpotSize ||
                    y < -adjustedSpotSize || y > height + adjustedSpotSize) continue;

                ctx.beginPath();
                ctx.arc(x, y, adjustedSpotSize, 0, Math.PI * 2);
                ctx.fill();
            }}
            ctx.globalAlpha = 1;
        }}

        // Second pass: draw visible categories on top
        for (let i = 0; i < modalSection.x.length; i++) {{
            const val = values[i];
            if (val === null || val === undefined) continue;

            let color;
            if (config.is_continuous) {{
                const t = (val - config.vmin) / (config.vmax - config.vmin);
                color = viridis(Math.max(0, Math.min(1, t)));
            }} else {{
                const catIdx = Math.round(val);
                const catName = config.categories[catIdx];
                if (hiddenCategories.has(catName)) continue;
                color = getCategoryColor(catIdx);
            }}

            const x = centerX + (modalSection.x[i] - dataCenterX) * scale;
            const y = centerY - (modalSection.y[i] - dataCenterY) * scale;

            if (x < -adjustedSpotSize || x > width + adjustedSpotSize ||
                y < -adjustedSpotSize || y > height + adjustedSpotSize) continue;

            ctx.fillStyle = color;
            ctx.beginPath();
            ctx.arc(x, y, adjustedSpotSize, 0, Math.PI * 2);
            ctx.fill();
        }}

        // Third pass: draw selection highlights
        if (selectedCells.size > 0) {{
            ctx.strokeStyle = '#ffd700';
            ctx.lineWidth = 3;
            for (let i = 0; i < modalSection.x.length; i++) {{
                if (!isCellSelected(modalSection.id, i)) continue;

                const val = values[i];
                if (val === null || val === undefined) continue;

                // Skip hidden categories
                if (!config.is_continuous) {{
                    const catIdx = Math.round(val);
                    const catName = config.categories[catIdx];
                    if (hiddenCategories.has(catName)) continue;
                }}

                const x = centerX + (modalSection.x[i] - dataCenterX) * scale;
                const y = centerY - (modalSection.y[i] - dataCenterY) * scale;

                if (x < -adjustedSpotSize || x > width + adjustedSpotSize ||
                    y < -adjustedSpotSize || y > height + adjustedSpotSize) continue;

                ctx.beginPath();
                ctx.arc(x, y, adjustedSpotSize + 2, 0, Math.PI * 2);
                ctx.stroke();
            }}
        }}

        drawNeighborHighlights(ctx, modalSection, {{
            scale,
            centerX,
            centerY,
            dataCenterX,
            dataCenterY,
            width,
            height
        }}, adjustedSpotSize);

        document.getElementById('zoom-info').textContent = `${{Math.round(modalZoom * 100)}}%`;
    }}

    // Legend
    function renderLegend(targetId = 'legend') {{
        const legend = document.getElementById(targetId);
        const config = getColorConfig();
        const colorLabel = currentGene || currentColor;

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
            const dpr = window.devicePixelRatio || 1;
            colorbar.width = 16 * dpr;
            colorbar.height = 150 * dpr;
            ctx.scale(dpr, dpr);
            for (let i = 0; i < 150; i++) {{
                ctx.fillStyle = viridis(1 - i / 149);
                ctx.fillRect(0, i, 16, 1);
            }}
        }} else {{
            let html = `
                <div class="legend-title">${{colorLabel}}</div>
                <div class="legend-actions">
                    <button class="legend-btn" id="${{targetId}}-show-all">Show All</button>
                    <button class="legend-btn" id="${{targetId}}-hide-all">Hide All</button>
                </div>
            `;
            (config.categories || []).forEach((cat, idx) => {{
                const hiddenClass = hiddenCategories.has(cat) ? 'hidden' : '';
                html += `<div class="legend-item ${{hiddenClass}}" data-category="${{cat}}">
                    <div class="legend-color" style="background: ${{getCategoryColor(idx)}}"></div>
                    <span>${{cat}}</span>
                </div>`;
            }});
            legend.innerHTML = html;

            document.getElementById(`${{targetId}}-show-all`)?.addEventListener('click', () => {{
                hiddenCategories.clear();
                renderLegend('legend');
                renderLegend('modal-legend');
                renderAllSections();
                if (modalSection) renderModalSection();
                if (umapVisible) renderUMAP();
            }});

            document.getElementById(`${{targetId}}-hide-all`)?.addEventListener('click', () => {{
                (config.categories || []).forEach(cat => hiddenCategories.add(cat));
                renderLegend('legend');
                renderLegend('modal-legend');
                renderAllSections();
                if (modalSection) renderModalSection();
                if (umapVisible) renderUMAP();
            }});

            legend.querySelectorAll('.legend-item').forEach(item => {{
                item.addEventListener('click', () => {{
                    const cat = item.dataset.category;
                    if (hiddenCategories.has(cat)) hiddenCategories.delete(cat);
                    else hiddenCategories.add(cat);
                    renderLegend('legend');
                    renderLegend('modal-legend');
                    renderAllSections();
                    if (modalSection) renderModalSection();
                    if (umapVisible) renderUMAP();
                }});
            }});
        }}
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
                <div class="color-panel-title">Color explorer</div>
            </div>
            <div class="color-panel-section">
                <label>Search colors</label>
                <input class="color-search" id="color-search" type="text" placeholder="Type to filter...">
            </div>
            <div class="color-panel-section">
                <label>Available colors</label>
                <div class="color-list" id="color-list"></div>
            </div>
            <div class="color-panel-section">
                <label>Aggregate by</label>
                <select id="color-groupby" ${{!hasMetadata ? 'disabled' : ''}}>
                    ${{options}}
                </select>
            </div>
            <div class="color-aggregation" id="color-aggregation">
                <div class="agg-group-meta">${{hasMetadata ? 'Pick a metadata column to summarize.' : 'No metadata columns available for aggregation.'}}</div>
            </div>
        `;

        const search = document.getElementById('color-search');
        search.addEventListener('input', () => {{
            renderColorList(search.value);
        }});

        const groupBy = document.getElementById('color-groupby');
        groupBy.addEventListener('change', () => {{
            renderColorAggregation();
        }});

        renderColorList('');
        renderColorAggregation();
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
                document.getElementById('color-select').value = col;
                document.getElementById('gene-input').value = '';
                hiddenCategories.clear();
                renderLegend('legend');
                renderLegend('modal-legend');
                renderAllSections();
                if (modalSection) renderModalSection();
                if (umapVisible) renderUMAP();
                renderColorList(document.getElementById('color-search').value);
                renderColorAggregation();
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

        filterBar.innerHTML = html;

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

                renderAllSections();
            }});
        }});
    }}

    // Modal
    function openModal(sectionId) {{
        modalSection = DATA.sections.find(s => s.id === sectionId);
        if (!modalSection) return;
        modalZoom = 1; modalPanX = 0; modalPanY = 0;

        document.getElementById('modal-title').textContent = sectionId;
        const metaText = Object.entries(modalSection.metadata || {{}})
            .map(([k, v]) => `${{formatMetadataLabel(k)}}: ${{v}}`).join(' | ');
        document.getElementById('modal-meta').textContent = metaText;
        document.getElementById('modal').classList.add('active');
        renderLegend('modal-legend');
        requestAnimationFrame(renderModalSection);
    }}

    function closeModal() {{
        document.getElementById('modal').classList.remove('active');
        modalSection = null;
        hideTooltip();
    }}

    // Grid
    function initGrid() {{
        const grid = document.getElementById('grid');
        grid.innerHTML = '';

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

            panel.innerHTML = `
                <div class="section-header">
                    <div>${{section.id}}${{metaHtml}}</div>
                    <span class="expand-icon">&#x26F6;</span>
                </div>
                <canvas class="section-canvas"></canvas>
            `;
            panel.addEventListener('click', () => openModal(section.id));
            grid.appendChild(panel);
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

        colorSelect.addEventListener('change', (e) => {{
            currentColor = e.target.value;
            currentGene = null;
            document.getElementById('gene-input').value = '';
            hiddenCategories.clear();
            renderLegend('legend');
            renderLegend('modal-legend');
            renderAllSections();
            if (modalSection) renderModalSection();
            if (umapVisible) renderUMAP();
            renderColorList(document.getElementById('color-search')?.value || '');
            renderColorAggregation();
        }});

        const geneList = document.getElementById('gene-list');
        (DATA.available_genes || []).forEach(gene => {{
            const opt = document.createElement('option');
            opt.value = gene;
            geneList.appendChild(opt);
        }});

        const geneInput = document.getElementById('gene-input');
        geneInput.addEventListener('change', () => {{
            const gene = geneInput.value.trim();
            if (gene && DATA.genes_meta[gene]) {{
                currentGene = gene;
                hiddenCategories.clear();
                renderLegend('legend');
                renderLegend('modal-legend');
                renderAllSections();
                if (modalSection) renderModalSection();
                if (umapVisible) renderUMAP();
                renderColorList(document.getElementById('color-search')?.value || '');
                renderColorAggregation();
            }} else if (gene) {{
                alert(`Gene "${{gene}}" was not pre-loaded.\\nTo view it, re-export with this gene included in the genes parameter or add it to highly variable genes.`);
            }}
        }});

        document.getElementById('spot-size').addEventListener('input', (e) => {{
            spotSize = parseFloat(e.target.value);
            renderAllSections();
            if (modalSection) renderModalSection();
        }});

        document.getElementById('screenshot-btn').addEventListener('click', screenshotFullPage);

        // Legend toggle
        document.getElementById('legend-toggle').addEventListener('click', () => {{
            const legend = document.getElementById('legend');
            const btn = document.getElementById('legend-toggle');
            legend.classList.toggle('collapsed');
            btn.classList.toggle('active');
            // Re-render to adjust for new grid size
            requestAnimationFrame(renderAllSections);
        }});

        // Color explorer toggle
        buildColorPanel();
        const colorToggle = document.getElementById('color-toggle');
        const colorPanel = document.getElementById('color-panel');
        colorToggle.addEventListener('click', () => {{
            colorPanel.classList.toggle('collapsed');
            colorToggle.classList.toggle('active');
            requestAnimationFrame(renderAllSections);
        }});

        // Neighborhood graph toggle
        if (DATA.has_neighbors) {{
            const graphBtn = document.getElementById('graph-toggle');
            graphBtn.style.display = 'inline-block';
            graphBtn.addEventListener('click', () => {{
                showGraph = !showGraph;
                graphBtn.classList.toggle('active', showGraph);
                renderAllSections();
                if (modalSection) renderModalSection();
            }});

            const neighborBtn = document.getElementById('neighbor-hover-toggle');
            neighborBtn.style.display = 'inline-block';
            neighborBtn.addEventListener('click', () => {{
                neighborHoverEnabled = !neighborHoverEnabled;
                neighborBtn.classList.toggle('active', neighborHoverEnabled);
                if (!neighborHoverEnabled) {{
                    hoverNeighbors = null;
                    if (modalSection) renderModalSection();
                }}
            }});

            const hopSelect = document.getElementById('neighbor-hop-select');
            hopSelect.style.display = 'inline-block';
            hopSelect.value = neighborHopMode;
            hopSelect.addEventListener('change', () => {{
                neighborHopMode = hopSelect.value;
                if (modalSection) renderModalSection();
            }});
        }}
    }}

    function initModal() {{
        document.getElementById('modal-close').addEventListener('click', closeModal);
        document.getElementById('modal').addEventListener('click', (e) => {{
            if (e.target.id === 'modal') closeModal();
        }});
        document.addEventListener('keydown', (e) => {{ if (e.key === 'Escape') closeModal(); }});

        document.getElementById('zoom-in').addEventListener('click', () => {{
            modalZoom = Math.min(modalZoom * 1.5, 20);
            renderModalSection();
        }});
        document.getElementById('zoom-out').addEventListener('click', () => {{
            modalZoom = Math.max(modalZoom / 1.5, 0.1);
            renderModalSection();
        }});
        document.getElementById('zoom-reset').addEventListener('click', () => {{
            modalZoom = 1; modalPanX = 0; modalPanY = 0;
            renderModalSection();
        }});

        document.getElementById('modal-spot-size').addEventListener('input', (e) => {{
            modalSpotSize = parseFloat(e.target.value);
            document.getElementById('modal-spot-size-label').textContent = modalSpotSize;
            renderModalSection();
        }});

        const container = document.getElementById('modal-canvas-container');
        container.addEventListener('wheel', (e) => {{
            if (!modalSection) return;
            e.preventDefault();
            const rect = container.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;

            const bounds = modalSection.bounds;
            const dataWidth = bounds.xmax - bounds.xmin;
            const dataHeight = bounds.ymax - bounds.ymin;
            const baseScale = Math.min((rect.width - 40) / dataWidth, (rect.height - 40) / dataHeight);
            const oldScale = baseScale * modalZoom;
            const nextZoom = Math.max(0.1, Math.min(20, modalZoom * (e.deltaY > 0 ? 0.9 : 1.1)));
            const newScale = baseScale * nextZoom;

            const dataCenterX = (bounds.xmin + bounds.xmax) / 2;
            const dataCenterY = (bounds.ymin + bounds.ymax) / 2;
            const centerX = rect.width / 2 + modalPanX;
            const centerY = rect.height / 2 + modalPanY;

            const dataX = dataCenterX + (mouseX - centerX) / oldScale;
            const dataY = dataCenterY - (mouseY - centerY) / oldScale;

            const newCenterX = mouseX - (dataX - dataCenterX) * newScale;
            const newCenterY = mouseY + (dataY - dataCenterY) * newScale;
            modalPanX = newCenterX - rect.width / 2;
            modalPanY = newCenterY - rect.height / 2;
            modalZoom = nextZoom;

            renderModalSection();
        }});

        const canvas = document.getElementById('modal-canvas');
        canvas.addEventListener('mousedown', (e) => {{
            isDragging = true;
            dragStartX = e.clientX; dragStartY = e.clientY;
            lastPanX = modalPanX; lastPanY = modalPanY;
            canvas.style.cursor = 'grabbing';
            hideTooltip();
        }});
        document.addEventListener('mousemove', (e) => {{
            if (!isDragging) return;
            modalPanX = lastPanX + (e.clientX - dragStartX);
            modalPanY = lastPanY + (e.clientY - dragStartY);
            renderModalSection();
        }});
        document.addEventListener('mouseup', () => {{
            isDragging = false;
            canvas.style.cursor = 'grab';
        }});
        canvas.style.cursor = 'grab';

        // Tooltip on hover in modal
        canvas.addEventListener('mousemove', (e) => {{
            if (isDragging || !modalSection) return;

            const rect = canvas.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;

            // Calculate transform parameters (same as renderModalSection)
            const bounds = modalSection.bounds;
            const dataWidth = bounds.xmax - bounds.xmin;
            const dataHeight = bounds.ymax - bounds.ymin;
            const baseScale = Math.min((rect.width - 40) / dataWidth, (rect.height - 40) / dataHeight);
            const scale = baseScale * modalZoom;
            const centerX = rect.width / 2 + modalPanX;
            const centerY = rect.height / 2 + modalPanY;
            const dataCenterX = (bounds.xmin + bounds.xmax) / 2;
            const dataCenterY = (bounds.ymin + bounds.ymax) / 2;

            const transform = {{
                scale,
                centerX,
                centerY,
                dataCenterX,
                dataCenterY,
                isModal: true
            }};

            const cellIdx = findNearestCell(modalSection, mouseX, mouseY, rect, transform);
            if (cellIdx >= 0) {{
                const content = getCellTooltipContent(modalSection, cellIdx);
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
            modalGraphBtn.style.display = 'inline-block';
            modalGraphBtn.classList.toggle('active', showGraph);
            modalGraphBtn.addEventListener('click', () => {{
                showGraph = !showGraph;
                modalGraphBtn.classList.toggle('active', showGraph);
                const graphBtn = document.getElementById('graph-toggle');
                if (graphBtn) graphBtn.classList.toggle('active', showGraph);
                renderAllSections();
                if (modalSection) renderModalSection();
            }});

            const modalNeighborBtn = document.getElementById('modal-neighbor-hover-toggle');
            modalNeighborBtn.style.display = 'inline-block';
            modalNeighborBtn.classList.toggle('active', neighborHoverEnabled);
            modalNeighborBtn.addEventListener('click', () => {{
                neighborHoverEnabled = !neighborHoverEnabled;
                modalNeighborBtn.classList.toggle('active', neighborHoverEnabled);
                const neighborBtn = document.getElementById('neighbor-hover-toggle');
                if (neighborBtn) neighborBtn.classList.toggle('active', neighborHoverEnabled);
                if (!neighborHoverEnabled) {{
                    hoverNeighbors = null;
                    if (modalSection) renderModalSection();
                }}
            }});

            const modalHopSelect = document.getElementById('modal-neighbor-hop-select');
            modalHopSelect.style.display = 'inline-block';
            modalHopSelect.value = neighborHopMode;
            modalHopSelect.addEventListener('change', () => {{
                neighborHopMode = modalHopSelect.value;
                const hopSelect = document.getElementById('neighbor-hop-select');
                if (hopSelect) hopSelect.value = neighborHopMode;
                if (modalSection) renderModalSection();
            }});
        }}
    }}

    // Initialize
    window.addEventListener('load', () => {{
        initTheme();
        initGrid();
        initControls();
        initFilters();
        initModal();
        initUMAP();
        renderLegend('legend');
        requestAnimationFrame(renderAllSections);
        const loader = document.getElementById('loading-overlay');
        if (loader) loader.style.display = 'none';
    }});
    window.addEventListener('resize', () => {{
        renderAllSections();
        if (modalSection) renderModalSection();
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
    spot_size: float = 2,
    downsample: Optional[int] = None,
    theme: str = "light",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    outline_by: Optional[str] = "course",
    additional_colors: Optional[List[str]] = None,
    genes: Optional[List[str]] = None,
    use_hvgs: bool = True,
) -> str:
    """
    Export spatial dataset to a standalone HTML file.

    Parameters
    ----------
    dataset : SpatialDataset
        Dataset to export
    output_path : str
        Path for output HTML file
    color : str
        Initial color column or gene name
    title : str
        Page title
    min_panel_size : int
        Minimum width of each section panel in pixels (default 150).
        The grid auto-adjusts columns based on screen width.
    spot_size : float
        Default spot size
    downsample : int, optional
        Downsample cells per section (for large datasets)
    theme : str
        'light' or 'dark'
    vmin, vmax : float, optional
        Min/max for continuous color scale
    outline_by : str, optional
        Metadata column used to color panel outlines (default: "course")
    additional_colors : list, optional
        Additional obs columns to include for color switching
    genes : list, optional
        Gene names to include for expression visualization

    Returns
    -------
    str
        Path to created HTML file
    """
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

    # Prefer highly variable genes for expression if available; otherwise use provided genes
    hv_genes = None
    if use_hvgs and "highly_variable" in dataset.adata.var.columns:
        hv_mask = dataset.adata.var["highly_variable"].to_numpy(dtype=bool)
        if hv_mask.any():
            hv_genes = dataset.adata.var_names[hv_mask].tolist()[:20]
    if hv_genes is not None:
        genes = hv_genes

    if outline_by and outline_by not in dataset.metadata_columns:
        print(f"  Warning: outline_by '{outline_by}' not in metadata columns; no outlines will be shown.")

    # Get data with multiple color layers and genes
    data = dataset.to_json_data(
        color,
        downsample=downsample,
        vmin=vmin,
        vmax=vmax,
        additional_colors=additional_colors,
        genes=genes,
    )

    # Theme settings
    theme_icon = "☀️" if theme == "dark" else "🌙"
    initial_theme = theme

    # Load logo for favicon and footer
    logo_base64 = _load_logo_base64()
    if logo_base64:
        favicon_link = f'<link rel="icon" type="image/png" href="data:image/png;base64,{logo_base64}">'
        footer_logo = f'<div class="footer-logo"><img src="data:image/png;base64,{logo_base64}" alt="KaroSpace"></div>'
    else:
        favicon_link = ""
        footer_logo = ""

    # Generate HTML
    metadata_labels = {
        "last_score": "disease score",
        "last_day": "day of sacrifice",
    }

    html = HTML_TEMPLATE.format(
        title=title,
        min_panel_size=min_panel_size,
        spot_size=spot_size,
        data_json=json.dumps(data),
        palette_json=json.dumps(DEFAULT_CATEGORICAL_PALETTE),
        metadata_labels_json=json.dumps(metadata_labels),
        outline_by_json=json.dumps(outline_by),
        theme_icon=theme_icon,
        initial_theme=initial_theme,
        favicon_link=favicon_link,
        footer_logo=footer_logo,
        **colors
    )

    # Write file
    output_path = str(Path(output_path).resolve())
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)

    print(f"Exported HTML viewer to: {output_path}")
    print(f"  - {data['n_sections']} sections")
    print(f"  - {data['total_cells']:,} cells")
    print(f"  - {len(data['available_colors'])} color options")
    if genes:
        print(f"  - {len(data['genes_meta'])} genes loaded")

    return output_path
