"""
Command-line interface for KaroSpace.
"""

import argparse
import sys
from pathlib import Path


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
        "--cols",
        type=int,
        default=4,
        help="Number of columns in grid (default: 4)"
    )
    parser.add_argument(
        "--spot-size",
        type=float,
        default=2.0,
        help="Default spot size (default: 2.0)"
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
        cols=args.cols,
        spot_size=args.spot_size,
        downsample=args.downsample,
        theme=args.theme,
    )

    print(f"Done! Open {output_path} in a browser to view.")


if __name__ == "__main__":
    main()
