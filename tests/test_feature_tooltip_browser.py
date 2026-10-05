"""Opt-in browser regression for the exported marker heatmap tooltips.

Run with KAROSPACE_TEST_BROWSER pointing to a Chrome/Chromium executable.
No scientific dependencies or dataset are required.
"""

import ast
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


def _template():
    tree = ast.parse((Path(__file__).parents[1] / "karospace/exporter.py").read_text())
    for statement in tree.body:
        if isinstance(statement, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "HTML_TEMPLATE"
            for target in statement.targets
        ):
            return ast.literal_eval(statement.value).replace("{{", "{").replace("}}", "}")
    raise AssertionError("Viewer template missing")


@unittest.skipUnless(os.environ.get("KAROSPACE_TEST_BROWSER"), "Set KAROSPACE_TEST_BROWSER to run Chrome")
class FeatureTooltipBrowserTests(unittest.TestCase):
    def test_heatmap_hover_keeps_readable_width_at_edges_and_after_scrolling(self):
        template = _template()
        css = template[template.index("<style>") + len("<style>"):template.index("</style>")]
        source = "\n".join([
            template[template.index("    function featureGraphPanel("):
                     template.index("    function quantileSorted(")],
            template[template.index("    function heatmapZColor("):
                     template.index("    function buildSpatialMoranGraph(")],
        ])
        # Keep only the tooltip functions, without subsequent application code.
        start = template.index("    function renderPseudobulkDEDotTooltip(")
        end = template.index("\n    function ", template.index("    function bindPseudobulkDEPlotInteractions(") + 10)
        source += "\n" + template[start:end]
        harness = r"""
function escapeHtml(value) {
    return String(value).replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('"', '&quot;');
}
function getCategoryColorForValue() {return '#cc8888';}
function getCategoryColor() {return '#cc8888';}
function colorToRgbaCss(color) {return color;}
function formatCategoryLabel(col, category) {return category;}
function formatScaleNumber(value) {return String(value);}
function check(condition, message) {if (!condition) throw Error(message);}
const container = document.getElementById('fixture');
const categories = Array.from({length: 30}, (_, i) => 'Neu-Ex-' + i);
const data = {features: ['Cd3g'], categories,
    means: [categories.map(() => 1.25)], zscores: [categories.map(() => 0.1)], deStars: new Set()};
container.innerHTML = buildFeatureDEHeatmap(data, 'cell_type');
bindPseudobulkDEPlotInteractions(container, () => {});
const panel = container.querySelector('.feature-graph-panel');
const cells = [...panel.querySelectorAll('[data-tooltip-title]')];
const tip = panel.querySelector('.volcano-tooltip');
function hover(cell) {
    cell.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
    return tip.getBoundingClientRect();
}
const initial = hover(cells[0]);
check(initial.width > 100, 'Initial tooltip must be readable');
// Repeated hovers reproduce shrink-to-fit measurement using the previous left.
for (const index of [15, 20, 22, 23, 24, 25, 24, 25, 0, 25]) {
    const rect = hover(cells[index]);
    check(Math.abs(rect.width - initial.width) <= 1,
        `Tooltip collapsed at column ${index}: ${rect.width}px vs ${initial.width}px`);
    check(Math.abs(rect.height - initial.height) <= 1, 'Tooltip became a tall text column');
    const bounds = panel.getBoundingClientRect();
    check(rect.left >= bounds.left && rect.right <= bounds.right, 'Tooltip escaped visible plot');
}
panel.scrollLeft = 180;
for (const index of [8, 29, 28, 29]) {
    const rect = hover(cells[index]);
    check(Math.abs(rect.width - initial.width) <= 1,
        `Scrolled tooltip collapsed: ${rect.width}px vs ${initial.width}px`);
    const bounds = panel.getBoundingClientRect();
    check(rect.left >= bounds.left && rect.right <= bounds.right, 'Scrolled tooltip escaped plot');
}
container.style.width = '180px';
panel.scrollLeft = 500;
cells[20].setAttribute('data-tooltip-line1', 'Category: ' + 'LongCategory'.repeat(20));
const narrow = hover(cells[20]);
const bounds = panel.getBoundingClientRect();
check(narrow.width > 100 && narrow.width <= panel.clientWidth, 'Narrow plot needs readable bounded width');
check(narrow.left >= bounds.left && narrow.right <= bounds.right, 'Long text overflowed narrow plot');
"""
        html = ('<!doctype html><style>' + css +
                '</style><div id="fixture" style="width:800px"></div><pre id="result"></pre><script>' +
                source + '\ntry {\n' + harness +
                '\ndocument.getElementById("result").textContent = "PASS";' +
                '\n} catch (error) {document.getElementById("result").textContent = "FAIL: " + error.message;}</script>')
        with tempfile.TemporaryDirectory(prefix="karospace-tooltip-") as directory:
            path = Path(directory) / "fixture.html"
            path.write_text(html)
            command = [
                os.environ["KAROSPACE_TEST_BROWSER"], "--headless", "--disable-gpu",
                "--no-first-run", "--no-default-browser-check",
                "--disable-background-networking", "--disable-component-update", "--disable-extensions",
                "--timeout=5000",
                "--user-data-dir=" + str(Path(directory) / "profile"),
                "--dump-dom", path.as_uri(),
            ]
            try:
                result = subprocess.run(command, capture_output=True, text=True, timeout=15, check=False)
                self.assertEqual(result.returncode, 0, result.stderr[-2000:])
                output = result.stdout
            except subprocess.TimeoutExpired as error:
                # Some managed Chrome installations hang during shutdown after dumping the DOM.
                output = (error.stdout or b'').decode()
                self.assertIn('</html>', output, 'Chrome timed out before rendering the fixture')
        marker = '<pre id="result">'
        self.assertTrue(marker + 'PASS</pre>' in output,
                        output[output.find(marker):][:400])


if __name__ == "__main__":
    unittest.main()
