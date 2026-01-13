"""
Export spatial data to standalone HTML viewer.

Creates self-contained HTML files with embedded data and JavaScript
for interactive visualization of spatial transcriptomics data.
"""

import json
from pathlib import Path
from typing import Optional, List

from .data_loader import SpatialDataset


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
        }}
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            background: var(--background);
            color: var(--text-color);
            min-height: 100vh;
            transition: background 0.3s, color 0.3s;
        }}
        .header {{
            padding: 12px 24px;
            background: var(--header-bg);
            border-bottom: 1px solid var(--border-color);
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
            gap: 12px;
            transition: background 0.3s, border-color 0.3s;
        }}
        .header h1 {{ font-size: 18px; font-weight: 600; }}
        .controls {{ display: flex; align-items: center; gap: 12px; flex-wrap: wrap; }}
        .control-group {{ display: flex; align-items: center; gap: 6px; }}
        .control-group label {{ font-size: 12px; color: var(--muted-color); }}
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
        select:focus, input:focus {{ outline: none; border-color: #0066cc; }}
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

        /* Filter bar */
        .filter-bar {{
            padding: 8px 24px;
            background: var(--header-bg);
            border-bottom: 1px solid var(--border-color);
            display: flex;
            align-items: center;
            gap: 16px;
            flex-wrap: wrap;
            transition: background 0.3s, border-color 0.3s;
        }}
        .filter-group {{ display: flex; align-items: center; gap: 6px; }}
        .filter-group label {{ font-size: 11px; color: var(--muted-color); text-transform: capitalize; }}
        .filter-chips {{ display: flex; gap: 4px; flex-wrap: wrap; }}
        .filter-chip {{
            padding: 3px 8px;
            font-size: 11px;
            border: 1px solid var(--border-color);
            border-radius: 12px;
            background: var(--input-bg);
            color: var(--text-color);
            cursor: pointer;
            transition: all 0.15s;
        }}
        .filter-chip:hover {{ background: var(--hover-bg); }}
        .filter-chip.active {{ background: #0066cc; color: white; border-color: #0066cc; }}
        .filter-chip.inactive {{ opacity: 0.4; }}

        .main-container {{ display: flex; height: calc(100vh - 110px); }}
        .grid-container {{
            flex: 1;
            padding: 12px;
            overflow: auto;
            display: grid;
            grid-template-columns: repeat({cols}, 1fr);
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
            padding: 6px 10px;
            background: var(--header-bg);
            border-bottom: 1px solid var(--border-color);
            font-size: 11px;
            font-weight: 500;
            display: flex;
            justify-content: space-between;
            align-items: center;
            transition: background 0.3s, border-color 0.3s;
        }}
        .section-header .expand-icon {{ font-size: 12px; opacity: 0.5; }}
        .section-meta {{ font-size: 9px; color: var(--muted-color); margin-top: 1px; }}
        .section-canvas {{ display: block; width: 100%; aspect-ratio: 1; }}

        .legend-container {{
            width: 200px;
            padding: 12px;
            background: var(--panel-bg);
            border-left: 1px solid var(--border-color);
            overflow-y: auto;
            font-size: 12px;
            transition: background 0.3s, border-color 0.3s;
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
    </style>
</head>
<body>
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
                <input type="range" id="spot-size" min="0.5" max="8" step="0.5" value="{spot_size}" style="width:80px">
            </div>
            <button class="theme-toggle" id="theme-toggle" title="Toggle dark/light mode">
                <span id="theme-icon">{theme_icon}</span>
            </button>
        </div>
        <div class="stats"><span id="stats-text"></span></div>
    </div>

    <div class="filter-bar" id="filter-bar"></div>

    <div class="main-container">
        <div class="grid-container" id="grid"></div>
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
                        <span style="margin-left: 10px; font-size: 11px; color: {muted_color};">Size:</span>
                        <input type="range" id="modal-spot-size" min="0.5" max="12" step="0.5" value="{spot_size}" style="width: 80px;">
                        <span id="modal-spot-size-label" style="font-size: 11px; min-width: 24px;">{spot_size}</span>
                    </div>
                </div>
                <div class="modal-legend" id="modal-legend"></div>
            </div>
        </div>
    </div>

    <div class="cell-tooltip" id="cell-tooltip"></div>

    <script>
    const DATA = {data_json};
    const PALETTE = {palette_json};

    // Course border color palette (rgba with 0.5 alpha for subtler borders)
    const COURSE_COLORS = {{
        'peak_I': 'rgba(228, 26, 28, 0.5)',
        'peak_II': 'rgba(55, 126, 184, 0.5)',
        'peak_III': 'rgba(77, 175, 74, 0.5)',
        'naive': 'rgba(152, 78, 163, 0.5)',
        'remission': 'rgba(255, 127, 0, 0.5)',
        'chronic': 'rgba(255, 255, 51, 0.5)',
        'acute': 'rgba(166, 86, 40, 0.5)',
        'control': 'rgba(153, 153, 153, 0.5)',
    }};

    function getCourseColor(course) {{
        if (!course) return null;
        // Check exact match first
        if (COURSE_COLORS[course]) return COURSE_COLORS[course];
        // Check case-insensitive
        const lowerCourse = course.toLowerCase();
        for (const [key, color] of Object.entries(COURSE_COLORS)) {{
            if (key.toLowerCase() === lowerCourse) return color;
        }}
        // Generate a consistent color for unknown courses
        let hash = 0;
        for (let i = 0; i < course.length; i++) {{
            hash = course.charCodeAt(i) + ((hash << 5) - hash);
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

    // Modal state
    let modalSection = null;
    let modalZoom = 1;
    let modalPanX = 0, modalPanY = 0;
    let modalSpotSize = {spot_size};
    let isDragging = false;
    let dragStartX = 0, dragStartY = 0;
    let lastPanX = 0, lastPanY = 0;

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
            }});

            document.getElementById(`${{targetId}}-hide-all`)?.addEventListener('click', () => {{
                (config.categories || []).forEach(cat => hiddenCategories.add(cat));
                renderLegend('legend');
                renderLegend('modal-legend');
                renderAllSections();
                if (modalSection) renderModalSection();
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
                }});
            }});
        }}
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
                <label>${{key}}:</label>
                <div class="filter-chips" data-filter="${{key}}">
                    ${{values.map(v => `<span class="filter-chip" data-value="${{v}}">${{v}}</span>`).join('')}}
                </div>
            </div>`;
        }}

        // Add course border legend if course metadata exists
        if (filters.course && filters.course.length > 0) {{
            html += `<div class="filter-group" style="margin-left: auto;">
                <label>Border colors:</label>
                <div style="display: flex; gap: 8px; align-items: center;">
                    ${{filters.course.map(c => `<span style="display: flex; align-items: center; gap: 3px; font-size: 10px;">
                        <span style="width: 12px; height: 12px; border: 3px solid ${{getCourseColor(c)}}; border-radius: 2px;"></span>
                        ${{c}}
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
            .map(([k, v]) => `${{k}}: ${{v}}`).join(' | ');
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

            // Apply course-based border color
            const course = section.metadata?.course;
            const borderColor = getCourseColor(course);
            if (borderColor) {{
                panel.style.borderColor = borderColor;
                panel.style.borderWidth = '3px';
            }}

            const metaParts = Object.entries(section.metadata || {{}})
                .map(([k, v]) => `${{k}}: ${{v}}`).join(' | ');
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
        }});

        const geneList = document.getElementById('gene-list');
        (DATA.all_genes || DATA.available_genes || []).forEach(gene => {{
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
            }} else if (gene && DATA.all_genes && DATA.all_genes.includes(gene)) {{
                alert(`Gene "${{gene}}" is in the dataset but was not pre-loaded.\\nTo view it, re-export with this gene included in the genes parameter.`);
            }} else if (gene) {{
                alert(`Gene "${{gene}}" not found in the dataset.`);
            }}
        }});

        document.getElementById('spot-size').addEventListener('input', (e) => {{
            spotSize = parseFloat(e.target.value);
            renderAllSections();
            if (modalSection) renderModalSection();
        }});
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
            e.preventDefault();
            modalZoom = Math.max(0.1, Math.min(20, modalZoom * (e.deltaY > 0 ? 0.9 : 1.1)));
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
            }} else {{
                hideTooltip();
            }}
        }});

        canvas.addEventListener('mouseleave', hideTooltip);
    }}

    // Initialize
    window.addEventListener('load', () => {{
        initTheme();
        initGrid();
        initControls();
        initFilters();
        initModal();
        renderLegend('legend');
        requestAnimationFrame(renderAllSections);
    }});
    window.addEventListener('resize', () => {{
        renderAllSections();
        if (modalSection) renderModalSection();
    }});
    </script>
</body>
</html>
'''


def export_to_html(
    dataset: SpatialDataset,
    output_path: str,
    color: str = "leiden",
    title: str = "Spatial Viewer",
    cols: int = 4,
    spot_size: float = 2,
    downsample: Optional[int] = None,
    theme: str = "light",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    additional_colors: Optional[List[str]] = None,
    genes: Optional[List[str]] = None,
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
    cols : int
        Number of columns in grid
    spot_size : float
        Default spot size
    downsample : int, optional
        Downsample cells per section (for large datasets)
    theme : str
        'light' or 'dark'
    vmin, vmax : float, optional
        Min/max for continuous color scale
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
        }

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

    # Generate HTML
    html = HTML_TEMPLATE.format(
        title=title,
        cols=cols,
        spot_size=spot_size,
        data_json=json.dumps(data),
        palette_json=json.dumps(DEFAULT_CATEGORICAL_PALETTE),
        theme_icon=theme_icon,
        initial_theme=initial_theme,
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
