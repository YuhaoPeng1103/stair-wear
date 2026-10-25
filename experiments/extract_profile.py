"""Extract a lateral wear profile from a point cloud CSV (x, y, z in mm).

Usage:
    python3 experiments/extract_profile.py \
        --points data/public_model/tread01_points.csv \
        --y-center 0 --halfwidth 20 \
        --out data/public_model/tread01_profile.csv

The output has x_mm, z_mm (raw height) and h_mm (wear depth relative to the
reference surface fitted to the outer edge fraction).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from io_field import extract_profile_from_points, reference_surface


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--points", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--y-center", type=float, default=0.0)
    p.add_argument("--halfwidth", type=float, default=20.0)
    p.add_argument("--x-min", type=float, default=-150.0)
    p.add_argument("--x-max", type=float, default=150.0)
    p.add_argument("--n", type=int, default=151)
    p.add_argument("--bin-width", type=float, default=2.0)
    p.add_argument("--edge-fraction", type=float, default=0.2)
    p.add_argument("--no-reference", action="store_true",
                   help="keep raw heights and set h_mm to the median-centred height")
    args = p.parse_args()

    points = pd.read_csv(args.points)
    missing = [c for c in ("x", "y", "z") if c not in points.columns]
    if missing:
        raise SystemExit(f"points CSV missing columns: {missing}")

    grid = np.linspace(args.x_min, args.x_max, args.n)
    z = extract_profile_from_points(points, args.y_center, args.halfwidth,
                                    grid, args.bin_width)
    if args.no_reference:
        h = z - np.nanmedian(z)
        poly = None
    else:
        h, poly = reference_surface(grid, z, args.edge_fraction)

    out = pd.DataFrame({"x_mm": grid, "z_mm": z, "h_mm": h})
    out.to_csv(args.out, index=False)
    print(f"wrote {args.out}: {len(out)} points")
    print(f"h_mm range: {np.nanmin(h):.3f} to {np.nanmax(h):.3f} mm")
    if poly is not None:
        print(f"reference polynomial: {poly}")


if __name__ == "__main__":
    main()
