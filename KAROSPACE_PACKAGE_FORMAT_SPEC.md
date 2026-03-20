# KaroSpace Package Format Spec

Status: Draft v1  
Scope: Proposed `.karospace` shareable package format for KaroSpace sidecar exports

## 1. Purpose

This document defines a single-file packaging format for KaroSpace viewers that currently export as:

- `viewer.html`
- `viewer.genes.json`
- `viewer.genes/` shard directory

The package format is intended to make that export shareable as one file without inventing a new viewer payload schema.

The design rule for v1 is:

- `.karospace` wraps the existing sidecar export
- `.karospace` does not define a third independent viewer data model

## 2. Goals

- Preserve the current KaroSpace sidecar export contract.
- Allow one-file sharing of a sidecar-based viewer.
- Keep package contents inspectable and debuggable.
- Keep regular HTML export and regular HTML+sidecar export unchanged.
- Support browser-based loading through a package loader page.

## 3. Non-Goals

- Replace regular standalone HTML export.
- Replace regular HTML+sidecar export.
- Invent a new gene manifest format or shard schema.
- Require KaroSpace Companion / Builder as the package opener.
- Guarantee direct double-click opening from `file://` with no loader.

## 4. Export Modes

KaroSpace should conceptually support three export targets:

1. Regular HTML
   - Single embedded HTML viewer.
2. HTML + sidecar
   - HTML viewer plus `genes.json` manifest plus shard directory.
3. `.karospace`
   - Single-file package containing the same sidecar export assets as mode 2.

Mode 3 is a transport/container format for mode 2, not a new viewer format.

As a convenience, KaroSpace may also emit a sibling local opener HTML such as `viewer.loader.html` next to `viewer.karospace`. That loader is not part of the package archive itself.

## 5. File Type

### 5.1 Extension

- Custom extension: `.karospace`

### 5.2 Container

- ZIP archive
- Standard ZIP tooling should be able to inspect it
- Compression method may be `store` or `deflate`

### 5.3 Path Rules

- Archive paths must use `/`
- All paths are relative to archive root
- No absolute paths
- No `..` path traversal

## 6. Archive Layout

The v1 package root contains:

```text
/
  karospace-package.json
  index.html
  <gene manifest>.json
  <gene shard dir>/
    000.json
    001.json
    ...
```

Example:

```text
/
  karospace-package.json
  index.html
  BALO.genes.json
  BALO.genes/
    000.json
    001.json
    002.json
```

## 7. Core Invariant

If the package is unpacked as-is, the unpacked files should be equivalent to a normal KaroSpace sidecar export and should work when served over HTTP from the unpacked directory.

That means:

- `index.html` is still a normal KaroSpace HTML viewer
- the packaged gene manifest is still a normal `karospace-gene-sidecar-manifest-v2`
- the packaged gene shards are still normal `karospace-gene-sidecar-shard-v2` JSON files

The package manifest only adds packaging metadata. It must not replace the existing viewer and sidecar formats.

## 8. Required Package Manifest

Archive root must contain `karospace-package.json`.

### 8.1 Manifest Schema

Required top-level fields:

```json
{
  "format": "karospace-package-v1",
  "package_version": 1,
  "entry_html": "index.html",
  "created_at": "2026-03-20T12:34:56Z",
  "producer": {
    "name": "karospace",
    "version": "..."
  },
  "viewer": {
    "mode": "sidecar-package",
    "gene_storage": "sidecar",
    "gene_manifest_path": "BALO.genes.json",
    "gene_shard_dir": "BALO.genes"
  },
  "files": {
    "index.html": {
      "media_type": "text/html"
    },
    "BALO.genes.json": {
      "media_type": "application/json"
    },
    "BALO.genes/000.json": {
      "media_type": "application/json"
    }
  }
}
```

### 8.2 Required Semantics

- `format`
  - Must equal `karospace-package-v1`
- `package_version`
  - Integer manifest version, initially `1`
- `entry_html`
  - Relative path to the HTML entrypoint, initially expected to be `index.html`
- `created_at`
  - UTC ISO-8601 timestamp
- `producer`
  - Exporting tool identity
- `viewer.mode`
  - Must equal `sidecar-package` for v1
- `viewer.gene_storage`
  - Must equal `sidecar` for v1
- `viewer.gene_manifest_path`
  - Relative path to the packaged gene manifest
- `viewer.gene_shard_dir`
  - Relative path to the packaged shard directory
- `files`
  - Declares all packaged files needed by the loader

### 8.3 Optional Metadata

Optional but recommended:

- `title`
- `dataset_name`
- `description`
- `n_sections`
- `total_cells`
- `sha256` per file
- `size_bytes` per file

## 9. Sidecar Asset Rules

### 9.1 HTML

- The package entry HTML should be derived from the same code path as a normal sidecar export.
- It should continue to point to a relative sidecar manifest path.

### 9.2 Gene Manifest

- Must remain `karospace-gene-sidecar-manifest-v2`
- `shards` paths inside the gene manifest must be valid relative archive paths after unpacking

### 9.3 Gene Shards

- Must remain `karospace-gene-sidecar-shard-v2`

## 10. Loader Contract

The package format does not require a desktop opener.

The intended v1 opening model is:

1. User visits a KaroSpace package loader page.
2. User selects or drops a `.karospace` file.
3. Loader reads the ZIP via browser File APIs.
4. Loader validates `karospace-package.json`.
5. Loader launches `entry_html` and makes the packaged assets available to the viewer runtime.

The loader is responsible for mapping package-relative asset paths to actual bytes in the archive.

Acceptable implementation strategies include:

- rewriting package-relative asset URLs to generated object URLs
- providing a package-aware fetch shim for the viewer
- launching the viewer in a controlled iframe with package asset resolution

The package format does not require one specific browser implementation strategy, but the loader must preserve the viewer’s ability to lazy-load gene assets from the package.

The reference implementation in this repo is:

- `karospace-package-loader.html`

That loader parses the ZIP in-browser, validates `karospace-package.json`, launches `index.html` in an iframe, and resolves packaged gene-manifest and shard fetches from the archive.

## 11. Validation Rules

A `.karospace` package is valid only if:

- `karospace-package.json` exists
- `format == "karospace-package-v1"`
- `entry_html` exists in the archive
- `viewer.gene_storage == "sidecar"`
- `viewer.gene_manifest_path` exists in the archive
- `viewer.gene_shard_dir` exists in the archive
- packaged gene manifest parses and reports `format == "karospace-gene-sidecar-manifest-v2"`
- every shard path named in the gene manifest exists in the archive
- no declared path escapes archive root

## 12. Backward Compatibility

v1 does not change:

- the embedded HTML export
- the current sidecar export
- the current sidecar manifest format
- the current sidecar shard format

This is intentional. Existing users and pipelines should not be forced onto package mode.

## 13. Versioning Strategy

There are three distinct version surfaces:

1. Package manifest format
   - `karospace-package-v1`
2. HTML viewer payload/schema
   - whatever the viewer already emits
3. Gene sidecar schema
   - currently `karospace-gene-sidecar-manifest-v2`
   - currently `karospace-gene-sidecar-shard-v2`

Future package versions may wrap newer sidecar versions without changing the core design principle that package mode is a wrapper over standard sidecar export.

## 14. Implementation Phases

Suggested implementation order:

1. Package writer in `export_to_html(...)`
   - emit sidecar export as usual
   - archive it into `.karospace`
   - write `karospace-package.json`
2. Browser loader page
   - open `.karospace`
   - validate manifest
   - provide asset access to the viewer
3. CLI / GUI support
   - add package export target
4. Validation tests
   - package structure
   - manifest correctness
   - sidecar asset completeness
   - loader round-trip

## 15. Acceptance Criteria For v1

v1 should be considered complete when:

- KaroSpace can export a `.karospace` file from a sidecar export path
- the package contains a valid `karospace-package.json`
- the package can be opened from a browser loader page
- genes can still be lazy-loaded from the packaged sidecar assets
- the same viewer content works in both unpacked sidecar mode and packaged mode
