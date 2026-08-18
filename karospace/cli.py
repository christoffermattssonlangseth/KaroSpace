"""
Command-line interface for KaroSpace.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional


def _parse_section_rotations_arg(raw: str) -> Optional[Dict[str, float]]:
    text = str(raw or "").strip()
    if not text:
        return None

    rotations: Dict[str, float] = {}
    for item in text.split(","):
        token = item.strip()
        if not token:
            continue
        if ":" not in token:
            raise ValueError(
                "each section rotation must use section_id:angle format"
            )
        section_id, angle_text = token.split(":", 1)
        section_id = section_id.strip()
        angle_text = angle_text.strip()
        if not section_id:
            raise ValueError("section rotation entries must include a section_id")
        try:
            rotations[section_id] = float(angle_text)
        except ValueError as exc:
            raise ValueError(
                f"invalid angle for section {section_id!r}: {angle_text!r}"
            ) from exc

    return rotations or None


def _run_export_cli(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate HTML viewer for Xenium spatial transcriptomics data"
    )
    io_args = parser.add_argument_group("Input/output")
    dataset_args = parser.add_argument_group("Dataset loading and coordinates")
    metadata_args = parser.add_argument_group("Metadata and labels")
    viewer_args = parser.add_argument_group("Viewer layout")
    gene_args = parser.add_argument_group("Gene content and storage")
    pseudobulk_args = parser.add_argument_group("Pseudobulk DE")
    neighborhood_args = parser.add_argument_group("Neighborhoods and interactions")
    overlay_args = parser.add_argument_group("Images, deconvolution, and utilities")

    io_args.add_argument(
        "input",
        type=str,
        help="Path to input .h5ad file or SpatialData .zarr store"
    )
    io_args.add_argument(
        "-o", "--output",
        type=str,
        default="karospace.html",
        help="Output HTML file path (default: karospace.html)"
    )
    io_args.add_argument(
        "--inspect-input",
        action="store_true",
        help=(
            "Inspect input AnnData/SpatialData table metadata and exit without building a viewer or running analytics."
        ),
    )
    viewer_args.add_argument(
        "--main-cell-annotation",
        type=str,
        default="leiden",
        dest="main_cell_annotation",
        help="Main cell annotation column or gene shown first in the viewer (default: leiden)"
    )
    viewer_args.add_argument(
        "--cell-annotations",
        type=str,
        default="",
        dest="cell_annotations",
        help="Comma-separated extra cell obs annotation columns to embed as selectable annotations "
             "(e.g. a second clustering). Needed to compare annotations in the River plot."
    )
    gene_args.add_argument(
        "--features",
        type=str,
        default="",
        help="Comma-separated features to preload for expression visualization. In embedded mode, significant DE genes are embedded automatically up to the configured cap."
    )
    dataset_args.add_argument(
        "--section-key",
        type=str,
        default="sample_id",
        help="Column to identify sections (default: sample_id)"
    )
    dataset_args.add_argument(
        "--section-order",
        type=str,
        default="",
        help="Comma-separated section IDs to control section ordering."
    )
    dataset_args.add_argument(
        "--spatial-key",
        type=str,
        default="spatial",
        help="Key in obsm containing spatial coordinates (default: spatial)"
    )
    dataset_args.add_argument(
        "--spatial-x",
        type=str,
        default=None,
        help=(
            "Obs/metadata column to use as X coordinates. Must be used with "
            "--spatial-y; creates adata.obsm[spatial_key] before export."
        ),
    )
    dataset_args.add_argument(
        "--spatial-y",
        type=str,
        default=None,
        help=(
            "Obs/metadata column to use as Y coordinates. Must be used with "
            "--spatial-x; creates adata.obsm[spatial_key] before export."
        ),
    )
    dataset_args.add_argument(
        "--spatialdata-table",
        type=str,
        default=None,
        help=(
            "AnnData table key to use when the input is a SpatialData object/store. "
            "Required when multiple tables are present and no table named 'table' exists."
        ),
    )
    metadata_args.add_argument(
        "--metadata-labels",
        type=str,
        default="",
        help=(
            "JSON object mapping metadata/obs column keys to display names in the viewer "
            '(example: {"sample_id":"Sample","last_score":"Disease score"}).'
        ),
    )
    metadata_args.add_argument(
        "--section-metadata",
        type=str,
        default="",
        dest="section_metadata",
        help=(
            "Comma-separated obs columns to use as section metadata and visual filter chips "
            "(e.g. strain,region,Batch,Slide). Empty uses loader defaults."
        ),
    )
    metadata_args.add_argument(
        "--section-metadata-extra",
        type=str,
        default="",
        dest="section_metadata_extra",
        help=(
            "Comma-separated obs columns to store as section metadata without adding "
            "visual filter chips."
        ),
    )
    metadata_args.add_argument(
        "--metadata-value-order",
        type=str,
        default="",
        help='JSON object mapping metadata columns to ordered value lists (example: {"strain":["WT","KO"]}).',
    )
    metadata_args.add_argument(
        "--metadata-max-columns",
        type=int,
        default=None,
        help="Limit the number of metadata columns used, preserving order."
    )

    viewer_args.add_argument(
        "--min-panel-size",
        type=int,
        default=150,
        help="Minimum panel width in pixels (default: 150). Grid auto-adjusts columns."
    )
    viewer_args.add_argument(
        "--spot-size",
        type=str,
        default="auto",
        help="Default spot size. Use a positive number or 'auto' (default: auto)."
    )
    viewer_args.add_argument(
        "--downsample",
        type=int,
        default=None,
        help="Downsample to N cells per section (for large datasets)"
    )
    viewer_args.add_argument(
        "--title",
        type=str,
        default="KaroSpace",
        help="Page title"
    )
    viewer_args.add_argument(
        "--outlineby",
        dest="outline_by",
        type=str,
        default="course",
        help=(
            "Metadata column used to paint panel outlines. Use 'None' to disable outlines. "
            "When the column is embedded as metadata/annotation, outlines reuse that palette. (default: course)"
        )
    )
    viewer_args.add_argument(
        "--outline-by",
        dest="outline_by",
        type=str,
        help=argparse.SUPPRESS,
    )
    viewer_args.add_argument(
        "--viewer-info-html",
        type=str,
        default=None,
        help="HTML string shown in the viewer Info tab."
    )
    viewer_args.add_argument(
        "--viewer-info-html-file",
        type=str,
        default=None,
        help="Path to an HTML fragment file shown in the viewer Info tab."
    )
    viewer_args.add_argument(
        "--tutorial",
        action="store_true",
        help=(
            "(in development) Embed the interactive guided tutorial and show "
            "an in-page control to start it."
        ),
    )
    gene_args.add_argument(
        "--feature-encoding",
        choices=["auto", "dense", "sparse"],
        default="auto",
        help="Feature vector encoding. 'sparse' stores only non-zero indices/values (smaller HTML for zero-inflated data). (default: auto)"
    )
    gene_args.add_argument(
        "--feature-value-encoding",
        choices=["uint16", "uint8"],
        default="uint16",
        help="Sidecar/package feature value encoding for binary shards. (default: uint16)"
    )
    gene_args.add_argument(
        "--feature-storage",
        choices=["embedded", "sidecar"],
        default="embedded",
        help="Store requested/top DE feature vectors in the HTML (`embedded`) or write all feature vectors to a sidecar (`sidecar`). (default: embedded)"
    )
    gene_args.add_argument(
        "--feature-manifest-path",
        type=str,
        default=None,
        help="Optional output path for the feature sidecar JSON when --feature-storage sidecar."
    )
    gene_args.add_argument(
        "--feature-sidecar-shard-size",
        type=int,
        default=256,
        help="Number of features per sidecar shard. (default: 256)"
    )
    gene_args.add_argument(
        "--modalities",
        type=str,
        default=None,
        help=(
            "Comma-separated list of modalities to export (e.g. 'rna,protein'). "
            "Defaults to all detected. Non-default modalities require --feature-storage sidecar."
        ),
    )
    gene_args.add_argument(
        "--feature-sparse-zero-threshold",
        type=float,
        default=0.8,
        help="Only used with --feature-encoding auto. Use sparse encoding when zero fraction >= threshold. (default: 0.8)"
    )
    neighborhood_args.add_argument(
        "--neighbor-permutations",
        type=str,
        default="auto",
        help="Neighbor enrichment permutation count. Use 0 to disable, or 'auto' (default) which disables for very large datasets."
    )
    neighborhood_args.add_argument(
        "--neighbor-stats-annotations",
        type=str,
        default="auto",
        help="Comma-separated obs columns to compute neighbor composition stats for. Use 'auto' (default) to match --main-cell-annotation; empty disables."
    )
    pseudobulk_args.add_argument(
        "--pseudobulk",
        type=str,
        default="auto",
        help=(
            "Category pseudobulk DE mode. Use 'auto' to analyze --main-cell-annotation "
            "and --pseudobulk-additional-annotations, or 'None' to disable. (default: auto)"
        )
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-additional-annotations",
        type=str,
        default="",
        help=(
            "Comma-separated additional annotation columns to analyze when pseudobulk or "
            "interaction markers are enabled. --main-cell-annotation is included automatically."
        )
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-replicate-annotation",
        type=str,
        default=None,
        help=(
            "Obs annotation to use as the biological replicate in pseudobulk analyses. "
            "Defaults to --section-key."
        )
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-simple-constrast-categories",
        type=str,
        default="",
        help=(
            "Categories to include in Simple design category-versus-category contrasts. "
            "Use comma-separated categories only with one pseudobulk annotation. With "
            "--pseudobulk-additional-annotations, use a JSON object keyed by annotation "
            "or a nested JSON list in order [main-cell-annotation, additional...], e.g. "
            "'{\"Anno_L1\":[\"Astrocyte\",\"B cell\"],\"region\":[\"Cortex\"]}'. "
            "In zsh/bash, wrap the whole JSON value in single quotes so inner double "
            "quotes are preserved. "
            "All retained categories remain in the shared DESeq2 fit and balanced-rest "
            "contrasts still run for every category. Empty includes all categories."
        ),
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-counts-layer",
        type=str,
        default="counts",
        help="AnnData layer containing raw counts for pseudobulk DE. Use 'None' for adata.X. (default: counts)"
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-min-cell-counts",
        type=int,
        default=0,
        help=(
            "Exclude cells with fewer than this many total raw counts before pseudobulk aggregation. "
            "Use 0 to disable. (default: 0)"
        ),
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-min-gene-counts",
        type=int,
        default=0,
        help=(
            "Exclude genes with fewer than this many total raw pseudobulk counts in the shared DESeq2 fit. "
            "Use 0 to disable. (default: 0)"
        ),
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-min-cells-per-pseudobulk",
        dest="pseudobulk_min_cells_per_pseudobulk",
        type=int,
        default=20,
        help=(
            "Minimum cells required in each replicate x annotation pseudobulk sample before it can enter "
            "the shared DESeq2 fit. (default: 20)"
        ),
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-n-cpus",
        type=int,
        default=1,
        help=(
            "CPU workers used for the shared PyDESeq2 fit and as the maximum number of "
            "parallel shared-fit contrasts. (default: 1)"
        ),
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-embed-top-n-per-comparison",
        type=int,
        default=20,
        help=(
            "Maximum significant DE genes to auto-embed per category/contact comparison in embedded mode. "
            "Ignored by --feature-storage sidecar, where all gene expression vectors are written to the sidecar. (default: 20)"
        ),
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-min-replicates",
        type=int,
        default=2,
        help="Minimum paired replicates required for each reported pseudobulk contrast; at least 2 are always required. (default: 2)"
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-min-pct-expressed",
        type=float,
        default=0.0,
        help="Minimum fraction of cells expressing a gene in at least one compared group before DeseqStats. Values >1 are interpreted as percentages. (default: 0)"
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-p-adjust-method",
        choices=["fdr_bh", "bonferroni", "holm", "none"],
        default="fdr_bh",
        help="Multiple-testing correction method for pseudobulk p-values. (default: fdr_bh)"
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-padj-cutoff",
        type=float,
        default=0.05,
        help="Adjusted p-value cutoff for volcano highlighting and DE table inclusion. (default: 0.05)"
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-log2fc-cutoff",
        type=float,
        default=0.5,
        help="Absolute log2FC cutoff for volcano highlighting and DE table inclusion. (default: 0.5)"
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-deseq2-fit-type",
        choices=["parametric", "mean"],
        default="parametric",
        help="PyDESeq2 dispersion trend fit type. Use 'mean' to avoid parametric trend fallback warnings. (default: parametric)"
    )
    pseudobulk_args.add_argument(
        "--pathway-gmt",
        type=str,
        default="",
        help=(
            "Comma-separated GMT pathway files for ORA/GSEA after Simple design pseudobulk DE. "
            "When omitted, KaroSpace uses the default Reactome library via GSEApy."
        ),
    )
    pseudobulk_args.add_argument(
        "--pathway-organism",
        type=str,
        default="Human",
        help="Organism passed to GSEApy for default Reactome loading, e.g. Human or Mouse. (default: Human)",
    )
    pseudobulk_args.add_argument(
        "--pathway-top-n",
        type=int,
        default=20,
        help="Maximum ORA/GSEA pathways stored per direction and comparison. (default: 20)",
    )
    pseudobulk_args.add_argument(
        "--pathway-min-overlap",
        type=int,
        default=3,
        help="Minimum pathway/query gene overlap for ORA/GSEA reporting. (default: 3)",
    )
    pseudobulk_args.add_argument(
        "--pathway-gsea-permutations",
        type=int,
        default=100,
        help="Permutations for compact preranked GSEA p-values. Use 0 for ES/NES only. (default: 100)",
    )
    neighborhood_args.add_argument(
        "--neighbor-stats-seed",
        type=int,
        default=0,
        help="Random seed for neighbor enrichment permutations. (default: 0)"
    )
    neighborhood_args.add_argument(
        "--interaction-markers",
        type=str,
        default="auto",
        help=(
            "Contact-conditioned pseudobulk marker mode. Use 'auto' to analyze --main-cell-annotation "
            "and --pseudobulk-additional-annotations, or 'None' to disable. (default: auto)"
        )
    )
    neighborhood_args.add_argument(
        "--interaction-markers-top-targets",
        type=int,
        default=8,
        help="Number of target categories to evaluate per source for contact-conditioned markers. (default: 8)"
    )
    neighborhood_args.add_argument(
        "--interaction-markers-top-genes",
        type=int,
        default=20,
        help="Number of top DE genes to keep per source-target interaction. (default: 20)"
    )
    neighborhood_args.add_argument(
        "--interaction-markers-min-cells",
        type=int,
        default=30,
        help="Minimum cells required per replicate contact+ and contact- pseudobulk sample. (default: 30)"
    )
    neighborhood_args.add_argument(
        "--interaction-markers-min-neighbors",
        type=int,
        default=1,
        help="Minimum target neighbors required to classify source cells as contact+. (default: 1)"
    )
    overlay_args.add_argument(
        "--section-rotations",
        type=str,
        default="",
        help="Comma-separated section_id:angle pairs for initial per-section rotations with exact degree values (example: S1:37.5,S2:-90)."
    )
    gene_args.add_argument(
        "--gene-correlation-top-n",
        type=int,
        default=10,
        help="Number of top correlated genes to show per embedded gene in the discovery panel. Use 0 to disable. (default: 10)"
    )
    gene_args.add_argument(
        "--category-means-n-genes",
        type=int,
        default=500,
        help="Maximum embedded pseudobulk-DE genes to expose in category mean summaries. Use 0 to disable. (default: 500)"
    )
    gene_args.add_argument(
        "--spatial-variable-genes-n",
        type=int,
        default=200,
        help="Number of top variable genes to score with Moran's I spatial autocorrelation. Requires spatial graph in obsp. Use 0 to disable. (default: 200)"
    )
    viewer_args.add_argument(
        "--scalebar-unit",
        type=str,
        default="μm",
        help="Unit label for the scalebar (default: μm). Assumes spatial coordinates are in this unit."
    )
    overlay_args.add_argument(
        "--deconvolutions",
        type=str,
        default="",
        help='JSON object mapping deconvolution labels to obs/obsm keys.'
    )
    overlay_args.add_argument(
        "--section-images",
        type=str,
        default="",
        help="JSON object mapping section IDs to image paths or image specs."
    )
    overlay_args.add_argument(
        "--section-images-max-px",
        type=int,
        default=4096,
        help="Maximum image dimension when embedding section images. (default: 4096)"
    )

    args = parser.parse_args(argv)

    if args.pseudobulk_min_cell_counts < 0:
        parser.error("--pseudobulk-min-cell-counts must be >= 0")
    if args.pseudobulk_min_gene_counts < 0:
        parser.error("--pseudobulk-min-gene-counts must be >= 0")
    if args.pseudobulk_n_cpus < 1:
        parser.error("--pseudobulk-n-cpus must be >= 1")
    if args.pseudobulk_embed_top_n_per_comparison < 0:
        parser.error("--pseudobulk-embed-top-n-per-comparison must be >= 0")
    if args.pathway_top_n < 1:
        parser.error("--pathway-top-n must be >= 1")
    if args.pathway_min_overlap < 1:
        parser.error("--pathway-min-overlap must be >= 1")
    if args.pathway_gsea_permutations < 0:
        parser.error("--pathway-gsea-permutations must be >= 0")

    # Check input file
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if input_path.suffix.lower() not in {".h5ad", ".zarr"}:
        print(
            f"Warning: Expected .h5ad file or SpatialData .zarr store, got: {input_path.suffix}",
            file=sys.stderr,
        )

    # Import here to avoid slow startup for --help
    from .data_loader import (
        inspect_input_file,
        load_spatial_data,
        normalize_pseudobulk_simple_constrast_categories,
    )
    from .exporter import export_to_html

    if args.inspect_input:
        report = inspect_input_file(args.input, spatialdata_table=args.spatialdata_table)
        print(f"Input: {report['path']}")
        if report.get("spatialdata_table"):
            print(f"SpatialData table: {report['spatialdata_table']}")
        print(f"Cells: {report['n_cells']:,} | Genes: {report['n_genes']:,}")
        print("Available cell metadata (adata.obs):")
        for entry in report["metadata"]:
            examples = ", ".join(json.dumps(value, ensure_ascii=False) for value in entry["examples"])
            example_text = examples or "—"
            missing = f"; {entry['n_missing']:,} missing" if entry["n_missing"] else ""
            print(
                f"  - {entry['name']} [{entry['type']}; {entry['n_unique']:,} values{missing}] "
                f"examples: {example_text}"
            )
        return

    neighbor_perms: Optional[int]
    if str(args.neighbor_permutations).lower() == "auto":
        neighbor_perms = None
    else:
        try:
            neighbor_perms = int(args.neighbor_permutations)
        except ValueError:
            print("Error: --neighbor-permutations must be an integer or 'auto'", file=sys.stderr)
            sys.exit(2)

    spot_token = str(args.spot_size).strip()
    if spot_token.lower() in {"auto", "adaptive", "density"}:
        spot_size_value = "auto"
    else:
        try:
            spot_size_value = float(spot_token)
        except ValueError:
            print("Error: --spot-size must be a positive number or 'auto'", file=sys.stderr)
            sys.exit(2)
        if spot_size_value <= 0:
            print("Error: --spot-size must be a positive number or 'auto'", file=sys.stderr)
            sys.exit(2)

    def _parse_csv(value: str):
        cleaned = [v.strip() for v in str(value).split(",") if v.strip()]
        return cleaned or None

    def _parse_json_object(value: str, option_name: str):
        text = str(value or "").strip()
        if not text:
            return None
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError as exc:
            print(f"Error: {option_name} must be valid JSON: {exc}", file=sys.stderr)
            sys.exit(2)
        if not isinstance(parsed, dict):
            print(f"Error: {option_name} must be a JSON object/dictionary", file=sys.stderr)
            sys.exit(2)
        return parsed

    def _parse_optional_layer(value: str):
        text = str(value or "").strip()
        if not text or text.lower() in {"none", "null"}:
            return None
        return text

    def _parse_auto_or_none(value: str, option_name: str) -> Optional[str]:
        text = str(value or "").strip().lower()
        if text in {"", "auto"}:
            return "auto"
        if text in {"none", "null"}:
            return None
        print(f"Error: {option_name} must be 'auto' or 'None'", file=sys.stderr)
        sys.exit(2)

    if str(args.neighbor_stats_annotations).lower() == "auto":
        neighbor_stats_annotations = [args.main_cell_annotation]
    else:
        neighbor_stats_annotations = _parse_csv(args.neighbor_stats_annotations)
    pseudobulk_mode = _parse_auto_or_none(args.pseudobulk, "--pseudobulk")
    interaction_markers_mode = _parse_auto_or_none(args.interaction_markers, "--interaction-markers")
    pseudobulk_additional_annotations = _parse_csv(args.pseudobulk_additional_annotations)
    pseudobulk_annotation_columns = [
        args.main_cell_annotation,
        *(pseudobulk_additional_annotations or []),
    ]
    try:
        pseudobulk_simple_constrast_categories = normalize_pseudobulk_simple_constrast_categories(
            args.pseudobulk_simple_constrast_categories,
            pseudobulk_annotation_columns,
            option_name="--pseudobulk-simple-constrast-categories",
        )
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(2)
    pathway_gmt = _parse_csv(args.pathway_gmt) or None
    pathway_organism = str(args.pathway_organism or "Human").strip() or "Human"
    cell_annotations = _parse_csv(args.cell_annotations)
    features = _parse_csv(args.features)
    section_order = _parse_csv(args.section_order)
    section_metadata = _parse_csv(args.section_metadata)
    section_metadata_extra = _parse_csv(args.section_metadata_extra)
    metadata_labels_raw = _parse_json_object(args.metadata_labels, "--metadata-labels")
    metadata_labels = {
        str(key): str(value)
        for key, value in (metadata_labels_raw or {}).items()
        if str(key).strip() and value is not None
    } or None
    metadata_value_order_raw = _parse_json_object(args.metadata_value_order, "--metadata-value-order")
    metadata_value_order = None
    if metadata_value_order_raw:
        metadata_value_order = {}
        for key, values in metadata_value_order_raw.items():
            if not isinstance(values, list):
                print("Error: --metadata-value-order values must be lists", file=sys.stderr)
                sys.exit(2)
            metadata_value_order[str(key)] = [str(v) for v in values]
    deconvolutions_raw = _parse_json_object(args.deconvolutions, "--deconvolutions")
    deconvolutions = {
        str(key): str(value)
        for key, value in (deconvolutions_raw or {}).items()
        if str(key).strip() and value is not None
    } or None
    section_images = _parse_json_object(args.section_images, "--section-images")
    viewer_info_html = args.viewer_info_html
    if args.viewer_info_html_file:
        try:
            viewer_info_html = Path(args.viewer_info_html_file).expanduser().read_text(encoding="utf-8")
        except OSError as exc:
            print(f"Error: could not read --viewer-info-html-file: {exc}", file=sys.stderr)
            sys.exit(2)
    try:
        section_rotations = _parse_section_rotations_arg(args.section_rotations)
    except ValueError as exc:
        print(f"Error: --section-rotations {exc}", file=sys.stderr)
        sys.exit(2)

    spatial_columns = None
    if bool(args.spatial_x) != bool(args.spatial_y):
        print("Error: --spatial-x and --spatial-y must be provided together", file=sys.stderr)
        sys.exit(2)
    if args.spatial_x and args.spatial_y:
        spatial_columns = (args.spatial_x, args.spatial_y)
    outline_token = str(args.outline_by or "").strip()
    outline_by = None if outline_token.lower() in {"", "none", "null"} else outline_token

    # Load and export
    load_kwargs = {
        "section_key": args.section_key,
        "spatial_key": args.spatial_key,
    }
    if args.spatialdata_table:
        load_kwargs["spatialdata_table"] = args.spatialdata_table
    if spatial_columns is not None:
        load_kwargs["spatial_columns"] = spatial_columns
    if section_order is not None:
        load_kwargs["section_order"] = section_order
    if section_metadata is not None:
        load_kwargs["section_metadata"] = section_metadata
    if section_metadata_extra is not None:
        load_kwargs["section_metadata_extra"] = section_metadata_extra
    if metadata_value_order is not None:
        load_kwargs["metadata_value_order"] = metadata_value_order
    if args.metadata_max_columns is not None:
        load_kwargs["metadata_max_columns"] = args.metadata_max_columns
    dataset = load_spatial_data(args.input, **load_kwargs)

    modalities_arg: Optional[List[str]] = None
    if args.modalities:
        modalities_arg = [m.strip() for m in args.modalities.split(",") if m.strip()]

    output_path = export_to_html(
        dataset,
        output_path=args.output,
        main_cell_annotation=args.main_cell_annotation,
        cell_annotations=cell_annotations,
        features=features,
        title=args.title,
        modalities=modalities_arg,
        min_panel_size=args.min_panel_size,
        spot_size=spot_size_value,
        downsample=args.downsample,
        outline_by=outline_by,
        metadata_labels=metadata_labels,
        viewer_info_html=viewer_info_html,
        tutorial=args.tutorial,
        feature_encoding=args.feature_encoding,
        feature_value_encoding=args.feature_value_encoding,
        feature_storage=args.feature_storage,
        feature_manifest_path=args.feature_manifest_path,
        feature_sidecar_shard_size=args.feature_sidecar_shard_size,
        feature_sparse_zero_threshold=args.feature_sparse_zero_threshold,
        neighbor_stats_permutations=neighbor_perms,
        neighbor_stats_annotations=neighbor_stats_annotations,
        neighbor_stats_seed=args.neighbor_stats_seed,
        interaction_markers_top_targets=args.interaction_markers_top_targets,
        interaction_markers_top_genes=args.interaction_markers_top_genes,
        interaction_markers_min_cells=args.interaction_markers_min_cells,
        interaction_markers_min_neighbors=args.interaction_markers_min_neighbors,
        pseudobulk=pseudobulk_mode,
        pseudobulk_additional_annotations=pseudobulk_additional_annotations,
        pseudobulk_replicate_annotation=args.pseudobulk_replicate_annotation,
        pseudobulk_simple_constrast_categories=pseudobulk_simple_constrast_categories,
        pseudobulk_counts_layer=_parse_optional_layer(args.pseudobulk_counts_layer),
        pseudobulk_min_cell_counts=args.pseudobulk_min_cell_counts,
        pseudobulk_min_gene_counts=args.pseudobulk_min_gene_counts,
        pseudobulk_min_cells_per_pseudobulk=args.pseudobulk_min_cells_per_pseudobulk,
        pseudobulk_min_replicates=args.pseudobulk_min_replicates,
        pseudobulk_min_pct_expressed=args.pseudobulk_min_pct_expressed,
        pseudobulk_p_adjust_method=args.pseudobulk_p_adjust_method,
        pseudobulk_padj_cutoff=args.pseudobulk_padj_cutoff,
        pseudobulk_log2fc_cutoff=args.pseudobulk_log2fc_cutoff,
        pseudobulk_deseq2_fit_type=args.pseudobulk_deseq2_fit_type,
        pseudobulk_n_cpus=args.pseudobulk_n_cpus,
        pseudobulk_embed_top_n_per_comparison=args.pseudobulk_embed_top_n_per_comparison,
        pathway_gmt=pathway_gmt,
        pathway_organism=pathway_organism,
        pathway_top_n=args.pathway_top_n,
        pathway_min_overlap=args.pathway_min_overlap,
        pathway_gsea_permutations=args.pathway_gsea_permutations,
        interaction_markers=interaction_markers_mode,
        section_rotations=section_rotations,
        deconvolutions=deconvolutions,
        gene_correlation_top_n=args.gene_correlation_top_n,
        category_means_n_genes=args.category_means_n_genes,
        spatial_variable_genes_n=args.spatial_variable_genes_n,
        scalebar_unit=args.scalebar_unit,
        section_images=section_images,
        section_images_max_px=args.section_images_max_px,
    )

    if args.feature_storage == "sidecar":
        output_obj = Path(output_path).expanduser()
        print(
            "Done! Sidecar gene loading requires HTTP(S). "
            f"Serve the output directory with: python -m http.server --directory {output_obj.parent}"
        )
        print(f"Then open http://localhost:8000/{output_obj.name}")
    else:
        print(f"Done! Open {output_path} in a browser to view.")


def _run_package_sidecar_cli(argv=None):
    parser = argparse.ArgumentParser(
        description="Package an existing KaroSpace sidecar viewer into a .karospace archive"
    )
    parser.add_argument(
        "html",
        type=str,
        help="Path to an existing sidecar HTML viewer",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output .karospace file path (default: <html stem>.karospace)",
    )
    parser.add_argument(
        "--feature-manifest-path",
        type=str,
        default=None,
        help="Optional actual path to the sidecar feature manifest JSON if it differs from the path referenced in the HTML.",
    )
    parser.add_argument(
        "--feature-shard-dir",
        type=str,
        default=None,
        help="Optional actual path to the sidecar shard directory if it differs from the manifest stem directory.",
    )
    parser.add_argument(
        "--loader-output",
        type=str,
        default=None,
        help="Optional output path for the companion .loader.html file.",
    )

    args = parser.parse_args(argv)

    html_path = Path(args.html)
    if not html_path.exists():
        print(f"Error: Sidecar HTML not found: {args.html}", file=sys.stderr)
        sys.exit(1)

    from .exporter import package_sidecar_viewer

    print(f"Packaging existing sidecar viewer: {args.html}")
    package_path = package_sidecar_viewer(
        html_path,
        output_path=args.output,
        feature_manifest_path=args.feature_manifest_path,
        feature_shard_dir=args.feature_shard_dir,
        loader_output_path=args.loader_output,
    )
    print(f"Done! Share {package_path} together with its .loader.html opener.")


def _run_ome_convert_cli(argv=None):
    parser = argparse.ArgumentParser(
        prog="karospace ome-convert",
        description=(
            "Optional helper: convert stitched TIFF(s) to pyramidal tiled "
            "OME-TIFF for import into Xenium Explorer. Independent of the HTML "
            "export — use it only if you need OME-TIFFs."
        ),
    )
    parser.add_argument("inputs", nargs="+", help="Input .tif file(s)")
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: alongside each input)",
    )
    parser.add_argument(
        "--pyramid-levels",
        type=int,
        default=4,
        help="Number of pyramid levels (default: 4)",
    )
    args = parser.parse_args(argv)

    from .omeconvert import convert_to_ome

    for input_path in args.inputs:
        output_path = None
        if args.output_dir:
            p = Path(input_path)
            output_path = str(Path(args.output_dir) / (p.stem + ".ome.tif"))
        convert_to_ome(input_path, output_path, args.pyramid_levels)


def main():
    argv = list(sys.argv[1:])
    if argv and argv[0] in {"package-sidecar", "package"}:
        _run_package_sidecar_cli(argv[1:])
        return
    if argv and argv[0] in {"ome-convert", "omeconvert"}:
        _run_ome_convert_cli(argv[1:])
        return
    _run_export_cli(argv)


if __name__ == "__main__":
    main()
