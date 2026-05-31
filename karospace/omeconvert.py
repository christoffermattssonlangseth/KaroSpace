"""
Convert stitched TIFF(s) to pyramidal tiled OME-TIFF for import into Xenium Explorer.

Usage:
    python omeconvert.py sample1.tif sample2.tif
    python omeconvert.py *.tif --output-dir ./ome_outputs/
    python omeconvert.py sample1.tif --pyramid-levels 5
"""

import tifffile
import numpy as np
import argparse
from pathlib import Path


def convert_to_ome(input_path: str, output_path: str = None, pyramid_levels: int = 4):
    input_path = Path(input_path)
    if output_path is None:
        output_path = input_path.with_suffix("").with_suffix(".ome.tif")
    else:
        output_path = Path(output_path)

    print(f"Converting: {input_path} -> {output_path}")
    img = tifffile.imread(str(input_path))
    print(f"  Shape: {img.shape}, dtype: {img.dtype}")

    # Ensure contiguous RGB array to avoid discontiguous storage warnings
    if img.ndim == 3 and img.shape[2] == 3:
        img = np.ascontiguousarray(img)
    elif img.ndim == 3 and img.shape[0] == 3:
        img = np.ascontiguousarray(np.moveaxis(img, 0, -1))

    options = {
        "tile": (1024, 1024),
        "compression": "zlib",
        "photometric": "rgb",
        "metadata": None,
    }

    with tifffile.TiffWriter(str(output_path), bigtiff=True) as tif:
        tif.write(img, subifds=pyramid_levels, **options)
        for _ in range(pyramid_levels):
            img = np.ascontiguousarray(img[::2, ::2])
            tif.write(img, subfiletype=1, **options)

    print(f"  Saved: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert TIFF(s) to pyramidal OME-TIFF for Xenium Explorer")
    parser.add_argument("inputs", nargs="+", help="Input .tif file(s)")
    parser.add_argument("--output-dir", help="Output directory (default: same as input)")
    parser.add_argument("--pyramid-levels", type=int, default=4, help="Number of pyramid levels (default: 4)")
    args = parser.parse_args()

    for input_path in args.inputs:
        output_path = None
        if args.output_dir:
            p = Path(input_path)
            output_path = Path(args.output_dir) / p.with_suffix("").with_suffix(".ome.tif").name
        convert_to_ome(input_path, output_path, args.pyramid_levels)
