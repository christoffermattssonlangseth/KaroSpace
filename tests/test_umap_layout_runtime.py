"""Check UMAP docking while the Insights sidebar changes width."""

import ast
from pathlib import Path
import shutil
import subprocess
import unittest


def _layout_source():
    tree = ast.parse((Path(__file__).parents[1] / "karospace/exporter.py").read_text())
    for statement in tree.body:
        if isinstance(statement, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "HTML_TEMPLATE"
            for target in statement.targets
        ):
            template = ast.literal_eval(statement.value).replace("{{", "{").replace("}}", "}")
            break
    opening = template[template.index("    function openInsightsMode("):
                       template.index("    function navigateModalSection(")]
    positioning = template[template.index("    function updateUMAPPanelPosition("):
                           template.index("    function loadUMAPPanelState(")]
    return opening + positioning


@unittest.skipUnless(shutil.which("node"), "Node.js is required for viewer runtime checks")
class UMAPLayoutRuntimeTests(unittest.TestCase):
    def test_lasso_insights_opening_tracks_content_width_through_transition(self):
        harness = r"""
const assert = require('node:assert/strict');
const frames = [], timers = [], observers = [];
const properties = new Map();
let contentRight = 1200;
const content = {getBoundingClientRect: () => ({left: 0, right: contentRight})};
const umap = {style: {setProperty: (key, value) => properties.set(key, value)}};
const classes = new Set(['collapsed']);
const insights = {classList: {remove: value => classes.delete(value)}};
const toggle = {classList: {add() {}}};
const document = {
    body: {style: {setProperty() {}}},
    getElementById: id => ({'content-column': content, 'umap-panel': umap,
        'insights-panel': insights, 'insights-toggle': toggle}[id] || null),
    querySelector: selector => selector === '#content-column' ? content : null
};
class ResizeObserver {
    constructor(callback) {this.callback = callback; this.targets = []; observers.push(this);}
    observe(target) {this.targets.push(target);}
}
const window = {innerWidth: 1200, ResizeObserver, setTimeout: cb => timers.push(cb)};
function requestAnimationFrame(callback) {frames.push(callback);}
function getComputedStyle() {return {getPropertyValue: () => '100'};}
function setInsightsMode(mode) {assert.equal(mode, 'selection');}
function flushFrames() {while (frames.length) frames.shift()();}
function resizeContent(right) {
    contentRight = right;
    observers.filter(observer => observer.targets.includes(content))
        .forEach(observer => observer.callback());
    flushFrames();
}
"""
        assertions = r"""
initStickyOffsetObservers();
// performLassoSelection opens the previously collapsed sidebar through this helper.
openInsightsMode('selection');
assert.equal(classes.has('collapsed'), false);
flushFrames();
assert.equal(properties.get('--umap-fixed-right'), '8px');
for (const right of [1100, 950, 748]) {
    resizeContent(right);
    assert.equal(properties.get('--umap-fixed-right'), `${1200 - right + 8}px`,
        'UMAP must remain inside the content area throughout the Insights transition');
}
// Closing the sidebar and resizing the viewport's content must also track its boundary.
for (const right of [950, 1100, 1200]) {
    resizeContent(right);
    assert.equal(properties.get('--umap-fixed-right'), `${1200 - right + 8}px`);
}
"""
        result = subprocess.run(
            ["node", "-e", harness + _layout_source() + assertions],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
