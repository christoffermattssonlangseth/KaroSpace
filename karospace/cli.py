"""
Command-line interface for KaroSpace.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, Optional


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


def main():
    parser = argparse.ArgumentParser(
        description="Generate HTML viewer for Xenium spatial transcriptomics data"
    )
    parser.add_argument(
        "input",
        type=str,
        help="Path to input .h5ad file"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="karospace.html",
        help="Output HTML file path (default: karospace.html)"
    )
    parser.add_argument(
        "-c", "--color",
        type=str,
        default="leiden",
        help="Initial color column or gene (default: leiden)"
    )
    parser.add_argument(
        "-g", "--groupby",
        type=str,
        default="sample_id",
        help="Column to group sections by (default: sample_id)"
    )
    parser.add_argument(
        "--min-panel-size",
        type=int,
        default=150,
        help="Minimum panel width in pixels (default: 150). Grid auto-adjusts columns."
    )
    parser.add_argument(
        "--spot-size",
        type=str,
        default="auto",
        help="Default spot size. Use a positive number or 'auto' (default: auto)."
    )
    parser.add_argument(
        "--downsample",
        type=int,
        default=None,
        help="Downsample to N cells per section (for large datasets)"
    )
    parser.add_argument(
        "--theme",
        choices=["light", "dark"],
        default="light",
        help="Color theme (default: light)"
    )
    parser.add_argument(
        "--title",
        type=str,
        default="KaroSpace",
        help="Page title"
    )
    parser.add_argument(
        "--gene-encoding",
        choices=["auto", "dense", "sparse"],
        default="auto",
        help="Gene vector encoding. 'sparse' stores only non-zero indices/values (smaller HTML for zero-inflated data). (default: auto)"
    )
    parser.add_argument(
        "--gene-storage",
        choices=["embedded", "sidecar"],
        default="embedded",
        help="Store genes in the HTML (`embedded`) or write non-embedded genes to an auxiliary JSON sidecar (`sidecar`). (default: embedded)"
    )
    parser.add_argument(
        "--gene-aux-path",
        type=str,
        default=None,
        help="Optional output path for the gene sidecar JSON when --gene-storage sidecar."
    )
    parser.add_argument(
        "--gene-sparse-zero-threshold",
        type=float,
        default=0.8,
        help="Only used with --gene-encoding auto. Use sparse encoding when zero fraction >= threshold. (default: 0.8)"
    )
    parser.add_argument(
        "--no-pack-arrays",
        dest="pack_arrays",
        action="store_false",
        help="Disable base64 packing of large per-section arrays (coords/colors/UMAP)."
    )
    parser.add_argument(
        "--pack-arrays-min-len",
        type=int,
        default=1024,
        help="Only pack arrays when section cell count >= this value. (default: 1024)"
    )
    parser.set_defaults(pack_arrays=True)

    parser.add_argument(
        "--neighbor-permutations",
        type=str,
        default="auto",
        help="Neighbor enrichment permutation count. Use 0 to disable, or 'auto' (default) which disables for very large datasets."
    )
    parser.add_argument(
        "--neighbor-stats-groupby",
        type=str,
        default="auto",
        help="Comma-separated obs columns to compute neighbor composition stats for. Use 'auto' (default) to match the initial color; empty disables."
    )
    parser.add_argument(
        "--marker-genes-groupby",
        type=str,
        default="",
        help="Comma-separated obs columns to compute marker genes for. Empty disables (default)."
    )
    parser.add_argument(
        "--interaction-markers-groupby",
        type=str,
        default="",
        help="Comma-separated obs columns to compute contact-conditioned interaction markers for. Empty disables (default)."
    )
    parser.add_argument(
        "--cluster-de-groupby",
        type=str,
        default="",
        help="Comma-separated categorical obs columns to precompute cluster-vs-cluster differential expression for. Empty disables (default)."
    )
    parser.add_argument(
        "--cluster-de-top-n",
        type=int,
        default=20,
        help="Top N genes to keep per pairwise cluster DE result. (default: 20)"
    )
    parser.add_argument(
        "--cluster-de-method",
        type=str,
        default="wilcoxon",
        help="Method for pairwise cluster DE via scanpy.tl.rank_genes_groups (default: wilcoxon)."
    )
    parser.add_argument(
        "--cluster-de-layer",
        type=str,
        default="normalized",
        help="AnnData layer to use for pairwise cluster DE when present. (default: normalized)"
    )
    parser.add_argument(
        "--cluster-de-min-cells",
        type=int,
        default=20,
        help="Minimum cells required in both clusters to report pairwise DE. (default: 20)"
    )
    parser.add_argument(
        "--section-rotations",
        type=str,
        default="",
        help="Comma-separated section_id:angle pairs for initial per-section rotations with exact degree values (example: S1:37.5,S2:-90)."
    )
    parser.add_argument(
        "--gene-correlation-top-n",
        type=int,
        default=10,
        help="Number of top correlated genes to show per embedded gene in the discovery panel. Use 0 to disable. (default: 10)"
    )

    args = parser.parse_args()

    # Check input file
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if not input_path.suffix == ".h5ad":
        print(f"Warning: Expected .h5ad file, got: {input_path.suffix}", file=sys.stderr)

    # Import here to avoid slow startup for --help
    from .data_loader import load_spatial_data
    from .exporter import export_to_html

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

    if str(args.neighbor_stats_groupby).lower() == "auto":
        neighbor_stats_groupby = [args.color]
    else:
        neighbor_stats_groupby = _parse_csv(args.neighbor_stats_groupby)
    marker_genes_groupby = _parse_csv(args.marker_genes_groupby)
    interaction_markers_groupby = _parse_csv(args.interaction_markers_groupby)
    cluster_de_groupby = _parse_csv(args.cluster_de_groupby)
    try:
        section_rotations = _parse_section_rotations_arg(args.section_rotations)
    except ValueError as exc:
        print(f"Error: --section-rotations {exc}", file=sys.stderr)
        sys.exit(2)

    # Load and export
    print(f"Loading data from: {args.input}")
    dataset = load_spatial_data(
        args.input,
        groupby=args.groupby,
    )

    print(f"Exporting to HTML...")
    output_path = export_to_html(
        dataset,
        output_path=args.output,
        color=args.color,
        title=args.title,
        min_panel_size=args.min_panel_size,
        spot_size=spot_size_value,
        downsample=args.downsample,
        theme=args.theme,
        gene_encoding=args.gene_encoding,
        gene_storage=args.gene_storage,
        gene_aux_path=args.gene_aux_path,
        gene_sparse_zero_threshold=args.gene_sparse_zero_threshold,
        pack_arrays=args.pack_arrays,
        pack_arrays_min_len=args.pack_arrays_min_len,
        neighbor_stats_permutations=neighbor_perms,
        neighbor_stats_groupby=neighbor_stats_groupby,
        marker_genes_groupby=marker_genes_groupby,
        interaction_markers_groupby=interaction_markers_groupby,
        cluster_de_groupby=cluster_de_groupby,
        cluster_de_top_n=args.cluster_de_top_n,
        cluster_de_method=args.cluster_de_method,
        cluster_de_layer=args.cluster_de_layer,
        cluster_de_min_cells=args.cluster_de_min_cells,
        section_rotations=section_rotations,
        gene_correlation_top_n=args.gene_correlation_top_n,
    )

    print(f"Done! Open {output_path} in a browser to view.")


if __name__ == "__main__":
    main()
