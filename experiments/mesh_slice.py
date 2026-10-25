"""Slice a staircase mesh and extract lateral wear profiles.

This is a semi-automatic tool. It assumes the staircase is axis-aligned enough
that one axis is 'up', one is 'width' (across the tread), and one is 'depth'
(along the walking direction). Inspect the printed bounding box first.

Example:
    python3 experiments/mesh_slice.py \
        --mesh data/public_model/umm_qais.glb \
        --out data/public_model/profiles \
        --up z --width x --depth y \
        --section-depth 0 --depth-halfwidth 30 \
        --n 151 --scale 1000

If the mesh is a Scene, all geometries are concatenated.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from io_field import reference_surface


def load_mesh(path: Path):
    import trimesh
    obj = trimesh.load(path, force="scene")
    if isinstance(obj, trimesh.Scene):
        # to_mesh() applies the scene-graph transforms and concatenates the
        # 65k-vertex chunks that Sketchfab uses.
        mesh = obj.to_mesh()
    elif isinstance(obj, trimesh.Trimesh):
        mesh = obj
    else:
        raise SystemExit(f"unsupported object: {type(obj)}")
    if not isinstance(mesh, trimesh.Trimesh) or len(mesh.faces) == 0:
        raise SystemExit("no triangle geometry in the file")
    return mesh


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--up", default="z", choices=["x", "y", "z"])
    ap.add_argument("--width", default="x", choices=["x", "y", "z"])
    ap.add_argument("--depth", default="y", choices=["x", "y", "z"])
    ap.add_argument("--section-depth", type=float, default=0.0,
                    help="coordinate of the depth slice (mesh units)")
    ap.add_argument("--depth-halfwidth", type=float, default=25.0)
    ap.add_argument("--height-min", type=float, default=None)
    ap.add_argument("--height-max", type=float, default=None)
    ap.add_argument("--n", type=int, default=151)
    ap.add_argument("--bin-width", type=float, default=2.0)
    ap.add_argument("--edge-fraction", type=float, default=0.2)
    ap.add_argument("--scale", type=float, default=1.0,
                    help="multiply mesh coordinates by this factor (e.g. 1000 for m->mm)")
    ap.add_argument("--n-samples", type=int, default=2_000_000)
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    mesh = load_mesh(Path(args.mesh))
    print(f"vertices={len(mesh.vertices)} faces={len(mesh.faces)}")
    print(f"bounds:\n{mesh.bounds}")
    print(f"extents: {mesh.extents}")

    points, _ = mesh.sample(args.n_samples, return_index=True)
    coords = {"x": 0, "y": 1, "z": 2}
    up, width, depth = coords[args.up], coords[args.width], coords[args.depth]
    pts = points * args.scale

    mask = np.abs(pts[:, depth] - args.section_depth * args.scale) <= (args.depth_halfwidth * args.scale)
    if args.height_min is not None:
        mask &= pts[:, up] >= args.height_min * args.scale
    if args.height_max is not None:
        mask &= pts[:, up] <= args.height_max * args.scale
    sub = pts[mask]
    if len(sub) < 100:
        raise SystemExit(f"only {len(sub)} points in the slice; adjust the section/height bounds")

    x = sub[:, width]
    z = sub[:, up]
    xmin, xmax = np.percentile(x, [0.5, 99.5])
    grid = np.linspace(xmin, xmax, args.n)
    bins = np.arange(xmin, xmax + args.bin_width, args.bin_width)
    idx = np.digitize(x, bins)
    zz = np.full(len(grid), np.nan)
    for j, gx in enumerate(grid):
        b = int(np.argmin(np.abs(bins - gx)))
        vals = z[idx == b]
        if len(vals):
            zz[j] = np.median(vals)
    good = np.isfinite(zz)
    if good.sum() < 10:
        raise SystemExit("too few populated bins; increase --depth-halfwidth or --n-samples")
    zz = np.interp(grid, grid[good], zz[good])

    h, poly = reference_surface(grid, zz, args.edge_fraction)
    outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame({"x_mm": grid, "z_mm": zz, "h_mm": h})
    out.to_csv(outdir / "profile.csv", index=False)
    print(f"wrote {outdir/'profile.csv'}  n={len(out)}")
    print(f"wear range: {h.min():.3f} to {h.max():.3f} (mesh units x scale)")
    print(f"reference polynomial: {poly}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot(grid, zz, color="#0072B2", lw=1, label="surface")
        ax.plot(grid, zz - h, color="#D55E00", lw=1, ls="--", label="reference")
        ax.set_xlabel("width (mesh units x scale)")
        ax.set_ylabel("height")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(outdir / "profile.png", dpi=200)
        print(f"wrote {outdir/'profile.png'}")


if __name__ == "__main__":
    main()
