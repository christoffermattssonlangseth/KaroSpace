# DRAFT review for PR #10 "KaroSpace for review" (Baboon61:main → main)

## === SUMMARY COMMENT (top-level review body) ===

Thanks for this — it's a big, genuinely valuable contribution. The pseudobulk DE,
pathway enrichment, tutorial overlay, gene modules, and region groups are all
features we want, and I've confirmed the image/H&E overlay layer is preserved.

Before we can merge, a few things need attention. The main point up front: even
though GitHub shows this as "mergeable / clean", that's only because it
fast-forwards — it's effectively a fork-scale rewrite (+33k/−21k, 121 files,
renamed public API, changed sidecar format, and the whole test suite removed).
So rather than fast-forward `main`, I'd like to land this via an integration
branch once the blockers below are resolved.

**Blockers (should fix before merge):**
1. **Test suite deleted with no replacement.** Commit `3651a53` removes
   `tests/test_analytics.py`, `tests/test_gene_sidecar_export.py`, and
   `tests/test_pseudobulk_interactions.py` (~3200 lines), leaving zero tests.
   Much of that covered code that still exists. Please restore the still-valid
   coverage (details inline) — at minimum the `test_analytics.py` helpers.
2. **pseudobulk log2FC / p-value come from different models** (see inline on
   `pseudobulk.py`). Reported fold-change is a mean-ratio while the p-value is
   DESeq2's — a gene can be "significant" with a mismatched log2FC.
3. **Pathway enrichment makes a network call at export time** (gseapy → Enrichr).
   Offline/CI runs will fail unless `--pathway-gmt` is supplied. Needs a cached/
   bundled fallback or clear documentation + graceful skip.

**Needs a decision (not necessarily blocking):**
4. **Public API rename** (`export_to_html`: `color→main_cell_annotation`,
   `genes→features`, `additional_colors→cell_annotations`, plus removed params)
   and the **sidecar binary format rename** (`gene→feature`). This breaks every
   existing example script and any external caller, and packages may not be
   cross-compatible between viewer versions. Let's treat it as an explicit
   breaking release (changelog + a thin back-compat shim), not a silent change.
5. **Removed features** — cluster-DE panel, dispersion insights, modal-blend
   co-expression, html2canvas screenshots. Some look intentionally superseded
   (modal-blend → OverviewSplit; html2canvas → native screenshots, fine). Please
   confirm cluster-DE and dispersion are deliberately dropped.

**Nits:**
6. `compute_pseudobulk_complex_design_de` is dead code ("no supported caller").
7. `pydeseq2` has no upper pin but relies on private internals (see inline).

Happy to pair on the integration-branch plan and on porting the tests.

---

## === INLINE COMMENTS (file : line → text) ===

### karospace/pseudobulk.py : 523  (inside _log2fc_from_pseudobulk_means / used by _fit_deseq2_pair ~599)
> This computes log2FC as a raw mean-ratio and (in the pair path) overrides
> DESeq2's model-based `log2FoldChange`, while the p-value/stat still come from
> DESeq2. Reported effect size and significance then reflect *different* models,
> which can produce genes flagged significant with a fold-change that disagrees
> with the tested contrast (confusing volcano plots). Either use DESeq2's
> (optionally shrunk) LFC consistently, or clearly label this as a display-only
> "mean-based FC" and apply it uniformly across the shared-fit path too (which
> currently doesn't re-inject it — so the two paths disagree).

### karospace/pseudobulk.py : 151  (_subset_deseq2_dataset)
> This reassigns `subset.__class__` and hand-restores ~15 private DeseqDataSet
> attributes. That's reverse-engineered PyDESeq2 internals and will likely break
> on a PyDESeq2 upgrade. Since `pydeseq2>=0.5.0` allows arbitrary newer versions,
> please add an upper bound (e.g. `>=0.5.0,<0.6`) or use a supported public API.

### karospace/pseudobulk.py : 857  (compute_pseudobulk_complex_design_de)
> Flagged "no supported caller" at line 780 and unused. ~220 lines of dead code —
> remove it, or gate behind an explicit experimental flag rather than shipping it.

### karospace/pathways.py : 89  (_load_reactome_from_gseapy → gp.get_library)
> `gp.get_library(...)` hits Enrichr over the network at export time. This makes
> exports fail offline / in CI / air-gapped unless `--pathway-gmt` is passed.
> Please cache or bundle a default Reactome GMT, and make the failure a graceful
> skip (export still succeeds, pathways omitted with a warning).

### karospace/exporter.py : 34475  (export_to_html signature — main_cell_annotation param)
> Renamed/removed params here break existing callers and all example scripts
> (`color→main_cell_annotation`, `genes→features`, `additional_colors→cell_annotations`,
> plus removed `theme`, `vmin/vmax`, `gene_encoding`, `gene_sidecar_format`,
> `marker_genes_*`, `cluster_de_*`). The sidecar/`.karospace` packaging itself is
> retained (good) but the naming/manifest format migrated gene→feature
> (`--gene-storage`→`--feature-storage`, `to_gene_sidecar_data`→`to_feature_sidecar_data`,
> `_validate_feature_sidecar_manifest_payload`), so packages written by the current
> `main` viewer may not validate/load in this one. Let's make this a deliberate
> breaking release: add a back-compat shim mapping the old kwargs (with a
> DeprecationWarning) + a manifest-version check, and note it in the changelog.

### tests/  (deleted in commit 3651a53)
> This commit removes the entire test suite (test_analytics.py,
> test_gene_sidecar_export.py, test_pseudobulk_interactions.py) with no
> replacement. Most of `test_analytics.py` covers helpers that still exist
> unchanged (`_numeric_category_perm`, `_resolve_spatial_key`,
> `_select_top_variable_genes`, `_compute_morans_i`) — please restore those
> directly. For `test_gene_sidecar_export.py`, the `cluster_de` assertions are
> obsolete but the binary sidecar / quantization / `.karospace` package contracts
> are still live — please port those against the renamed `to_feature_sidecar_data`.
