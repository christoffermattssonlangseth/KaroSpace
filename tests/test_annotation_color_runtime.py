"""Exercise the exported legend JavaScript without scientific dependencies."""

import ast
import json
from pathlib import Path
import shutil
import subprocess
import unittest


def _viewer_source():
    tree = ast.parse((Path(__file__).parents[1] / "karospace/exporter.py").read_text())
    values = {}
    for statement in tree.body:
        if isinstance(statement, ast.Assign):
            for target in statement.targets:
                if isinstance(target, ast.Name) and target.id in {
                    "HTML_TEMPLATE", "DEFAULT_CATEGORICAL_PALETTE"
                }:
                    values[target.id] = ast.literal_eval(statement.value)
    template = values["HTML_TEMPLATE"].replace("{{", "{").replace("}}", "}")
    helpers = template[template.index("    function getCategoryColor("):
                       template.index("    function buildColorPaletteExport(")]
    legend = template[template.index("    function renderLegend("):
                      template.index("    function updateAnnotationComparisonTabVisibility(")]
    return helpers + legend, values["DEFAULT_CATEGORICAL_PALETTE"]


@unittest.skipUnless(shutil.which("node"), "Node.js is required for viewer runtime checks")
class AnnotationColorRuntimeTests(unittest.TestCase):
    def run_viewer(self, assertions):
        source, palette = _viewer_source()
        harness = r"""
const assert = require('node:assert/strict');
const DATA = {annotations_meta: {cell_type: {
    categories: ['A', 'B', 'C'], is_continuous: false
}}};
const meta = DATA.annotations_meta.cell_type;
const legend = {innerHTML: '', querySelectorAll: () => []};
// Model the browser canvas color API used to canonicalize CSS colors.
const colors = {red: '#ff0000', 'rgb(12, 34, 56)': '#0c2238',
    'hsl(120, 100%, 50%)': '#00ff00'};
const context = {
    _fillStyle: '#000000',
    set fillStyle(value) {
        if (colors[value]) this._fillStyle = colors[value];
        else if (/^#[0-9a-f]{6}$/i.test(value)) this._fillStyle = value.toLowerCase();
    },
    get fillStyle() {return this._fillStyle;},
    fillRect() {}, clearRect() {},
    getImageData() {
        return {data: [1, 3, 5].map(i => parseInt(this._fillStyle.slice(i, i + 2), 16))};
    }
};
const document = {
    getElementById: id => id === 'legend' || id === 'modal-legend' ? legend : null,
    createElement: () => ({getContext: () => context})
};
const hiddenCategories = new Set();
const currentFeature = null, currentAnnotation = 'cell_type';
const overviewBlendEnabled = false, modalSelectedCategory = null;
const linkedSpotlightEnabled = false;
const LEGEND_EYE_ICON = '', LEGEND_EYE_OFF_ICON = '', LEGEND_SPOTLIGHT_ICON = '';
const LEGEND_EXPORT_ICON = '', LEGEND_IMPORT_ICON = '';
function getColorConfig() {return {...meta, annotationCol: 'cell_type'};}
function getFeatureDisplayLabel() {return '';}
function getLinkedSpotlightCategory() {return null;}
function escapeHtml(value) {return value;}
function updateLegendSpotlightClasses() {}
function getCategoriesForColorColumn() {return meta.categories;}
function editorColors() {
    renderLegend();
    return [...legend.innerHTML.matchAll(/data-category-color="[^"]*" value="([^"]*)"/g)]
        .map(match => match[1]);
}
"""
        result = subprocess.run(
            ["node", "-e", "const PALETTE = " + json.dumps(palette) + ";\n" +
             harness + source + assertions],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_missing_palette_has_stable_random_non_grey_colors(self):
        self.run_viewer(r"""
let randomCalls = 0;
Math.random = () => {randomCalls++; return (randomCalls * 0.137) % 1;};
meta.categories = Array.from({length: 60}, (_, i) => String(i));
const initial = editorColors();
assert.ok(randomCalls > 0, 'Missing input colors should generate random colors');
assert.equal(initial.length, 60);
initial.forEach((color, i) => {
    assert.match(color, /^#[0-9a-f]{6}$/);
    assert.ok(color.slice(1, 3) !== color.slice(3, 5) ||
              color.slice(3, 5) !== color.slice(5, 7), 'Fallback must not be grey');
    assert.equal(getCategoryColor(i, 'cell_type'), color);
});
assert.deepEqual(editorColors(), initial, 'Rerendering must keep colors stable');
assert.deepEqual(ensureColorColumnPalette('cell_type'), initial);
""")

    def test_input_css_palette_matches_editor_and_survives_editing(self):
        self.run_viewer(r"""
meta.palette = ['red', 'rgb(12, 34, 56)', 'hsl(120, 100%, 50%)'];
assert.deepEqual(editorColors(), ['#ff0000', '#0c2238', '#00ff00']);
assert.equal(setCategoryColorOverride('cell_type', 'B', '#abcdef'), true);
assert.deepEqual(editorColors(), ['#ff0000', '#abcdef', '#00ff00']);
assert.equal(getCategoryColor(1, 'cell_type'), '#abcdef');
""")

    def test_valid_input_greys_are_preserved_and_missing_entries_generated(self):
        self.run_viewer(r"""
meta.palette = ['#999999', '#abc', ''];
const initial = editorColors();
assert.deepEqual(initial.slice(0, 2), ['#999999', '#aabbcc']);
assert.notEqual(initial[2], '#999999');
assert.equal(getCategoryColor(2, 'cell_type'), initial[2]);
assert.deepEqual(editorColors(), initial);
""")

    def test_invalid_css_uses_random_color_and_valid_sentinels_are_preserved(self):
        self.run_viewer(r"""
meta.palette = ['not-a-color', '#010203', '#040506'];
const initial = editorColors();
assert.deepEqual(initial.slice(1), ['#010203', '#040506']);
assert.match(initial[0], /^#[0-9a-f]{6}$/);
assert.notEqual(initial[0], '#999999');
assert.deepEqual(editorColors(), initial);
""")


if __name__ == '__main__':
    unittest.main()
