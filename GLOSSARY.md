# KaroSpace Terminology

## Feature

Use **feature** as the generic label for a displayable per-cell or per-spot value in a
modality. A feature can be an RNA gene, protein, peak, genomic range, bin, module,
or another modality-specific measurement.

Use feature in generic UI labels, CLI help, export summaries, logs, and docs when
the current modality may be anything other than RNA/gene expression.

## Gene

Use **gene** only when the context is specifically gene expression, RNA/rna data,
pathway gene sets, ORA/GSEA gene universes, or compatibility fields whose persisted
schema name is already gene-based.

Legacy public option names and payload keys that include `gene` may remain as
compatibility aliases, but new user-facing names should prefer `feature`.
