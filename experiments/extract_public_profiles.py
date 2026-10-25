"""Extract lateral wear profiles from a segmented public staircase model.

The script segments horizontal tread faces by height, clusters them in plan
view, and extracts a median-height profile across each tread. It is tuned for
the Paestum and Roman marble models, whose scene coordinates are Y-up after
trimesh applies the glTF transforms.

Outputs:
  data/public_model/profiles/<model>_tread<k>.csv
  data/public_model/profiles/<model>_profiles.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from io_field import reference_surface

MODELS = {
    "paestum": {
        "path": ROOT / "data/public_model/extracted/paestum/scene.gltf",
        "width_axis": 2,   # Z
        "depth_axis": 0,   # X
        "y_values": [792.5, 1062.5, 1292.5, 1512.5, 1742.5],
    },
    "roman": {
        "path": ROOT / "data/public_model/extracted/roman/scene.gltf",
        "width_axis": 0,   # X
        "depth_axis": 2,   # Z
        "y_values": [147.5, 347.5, 562.5, 787.5, 927.5],
    },
}


def load_scaled(path: Path, scale: float):
    import trimesh
    obj = trimesh.load(path, force="scene")
    mesh = obj.to_mesh() if isinstance(obj, trimesh.Scene) else obj
    mesh.apply_scale(scale)
    return mesh


def segment_tread(mesh, y0, width_axis, depth_axis, height_halfwidth=12,
                  eps=25, min_samples=15):
    c = mesh.triangles_center
    n = mesh.face_normals
    sel = (np.abs(n[:, 1]) > 0.9) & (np.abs(c[:, 1] - y0) < height_halfwidth)
    if sel.sum() < 50:
        return None
    plan = c[sel][:, [depth_axis, width_axis]]
    labels = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(plan)
    clusters = [lab for lab in set(labels) if lab != -1]
    if not clusters:
        return None
    best = max(clusters, key=lambda lab: (labels == lab).sum())
    sub = plan[labels == best]
    return {
        "depth_min": float(sub[:, 0].min()), "depth_max": float(sub[:, 0].max()),
        "width_min": float(sub[:, 1].min()), "width_max": float(sub[:, 1].max()),
        "n_faces": int((labels == best).sum()),
    }


def extract(mesh, y0, box, width_axis, depth_axis, depth_halfwidth=20,
            height_halfwidth=10, bin_width=5.0, n_samples=1_500_000,
            seed=1):
    rng = np.random.default_rng(seed)
    pts, _ = mesh.sample(n_samples, return_index=True)
    depth = pts[:, depth_axis]
    width = pts[:, width_axis]
    height = pts[:, 1]
    depth_center = 0.5 * (box["depth_min"] + box["depth_max"])
    mask = ((np.abs(depth - depth_center) <= depth_halfwidth) &
            (np.abs(height - y0) <= height_halfwidth + 5) &
            (width >= box["width_min"] - 5) & (width <= box["width_max"] + 5))
    sub = pts[mask]
    if len(sub) < 200:
        return None
    x = sub[:, width_axis]
    z = sub[:, 1]
    xmin, xmax = np.percentile(x, [0.5, 99.5])
    grid = np.linspace(xmin, xmax, max(30, int((xmax - xmin) / bin_width) + 1))
    bins = np.arange(xmin, xmax + bin_width, bin_width)
    idx = np.digitize(x, bins)
    zz = np.full(len(grid), np.nan)
    for j, gx in enumerate(grid):
        b = int(np.argmin(np.abs(bins - gx)))
        vals = z[idx == b]
        if len(vals):
            zz[j] = np.median(vals)
    good = np.isfinite(zz)
    if good.sum() < 20:
        return None
    zz = np.interp(grid, grid[good], zz[good])
    h, poly = reference_surface(grid, zz, edge_fraction=0.2, degree=2)
    return pd.DataFrame({"x_mm": grid, "z_mm": zz, "h_mm": h}), poly


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", choices=sorted(MODELS), default="roman")
    ap.add_argument("--scale", type=float, default=1000.0)
    ap.add_argument("--bin-width", type=float, default=5.0)
    args = ap.parse_args()
    cfg = MODELS[args.model]
    mesh = load_scaled(cfg["path"], args.scale)
    outdir = ROOT / "data" / "public_model" / "profiles"
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for k, y0 in enumerate(cfg["y_values"], start=1):
        box = segment_tread(mesh, y0, cfg["width_axis"], cfg["depth_axis"])
        if box is None:
            print(f"tread {k} (y={y0}): no cluster, skip")
            continue
        res = extract(mesh, y0, box, cfg["width_axis"], cfg["depth_axis"],
                      bin_width=args.bin_width)
        if res is None:
            print(f"tread {k} (y={y0}): extraction failed, skip")
            continue
        df, poly = res
        df.to_csv(outdir / f"{args.model}_tread{k}.csv", index=False)
        rows.append({"tread": k, "y0_mm": y0, "n": len(df),
                     "width_min_mm": box["width_min"], "width_max_mm": box["width_max"],
                     "wear_min_mm": float(df.h_mm.min()), "wear_max_mm": float(df.h_mm.max()),
                     "wear_range_mm": float(df.h_mm.max() - df.h_mm.min())})
        print(f"tread {k} (y={y0}): width {box['width_min']:.0f}..{box['width_max']:.0f} mm, "
              f"wear {df.h_mm.min():.2f}..{df.h_mm.max():.2f} mm")

    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / f"{args.model}_summary.csv", index=False)

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.0, 3.0))
    for k in summary["tread"]:
        df = pd.read_csv(outdir / f"{args.model}_tread{k}.csv")
        ax.plot(df.x_mm - df.x_mm.mean(), df.h_mm, lw=1.2, label=f"tread {k}")
    ax.axhline(0, color="0.5", lw=0.6, ls=":")
    ax.set_xlabel("lateral position (mm, centred)")
    ax.set_ylabel("wear depth (mm)")
    ax.legend(frameon=False, ncol=2, fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir / f"{args.model}_profiles.png", dpi=300)
    fig.savefig(outdir / f"{args.model}_profiles.pdf")
    plt.close(fig)
    print("wrote", outdir / f"{args.model}_profiles.png")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
