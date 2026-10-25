"""Real-geometry case study on the Roman marble staircase model.

The point is not to estimate traffic (the model has no known traffic). The point
is to show how much the inferred wear depth depends on the reference-surface
model when the same scan is used. That sensitivity is the real-data counterpart
of the simulation's identifiability result.

Outputs:
  data/public_model/case/roman_reference_sensitivity.csv
  data/public_model/case/roman_case.png / .pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "data/public_model" / "extracted" / "roman" / "scene.gltf"
OUT = ROOT / "data" / "public_model" / "case"
OUT.mkdir(parents=True, exist_ok=True)

Y_VALUES = [147.5, 347.5, 562.5, 787.5, 927.5]
WIDTH_AXIS, DEPTH_AXIS = 0, 2   # X width, Z depth


def load_scaled(path, scale=1000.0):
    import trimesh
    obj = trimesh.load(path, force="scene")
    mesh = obj.to_mesh() if isinstance(obj, trimesh.Scene) else obj
    mesh.apply_scale(scale)
    return mesh


def segment(mesh, y0):
    c = mesh.triangles_center
    n = mesh.face_normals
    sel = (np.abs(n[:, 1]) > 0.9) & (np.abs(c[:, 1] - y0) < 12)
    if sel.sum() < 50:
        return None
    plan = c[sel][:, [DEPTH_AXIS, WIDTH_AXIS]]
    lab = DBSCAN(eps=25, min_samples=15).fit_predict(plan)
    clusters = [l for l in set(lab) if l != -1]
    if not clusters:
        return None
    best = max(clusters, key=lambda l: (lab == l).sum())
    sub = plan[lab == best]
    return {
        "depth_min": float(sub[:, 0].min()), "depth_max": float(sub[:, 0].max()),
        "width_min": float(sub[:, 1].min()), "width_max": float(sub[:, 1].max()),
    }


def fit_reference(P, method):
    x = P[:, WIDTH_AXIS]
    d = P[:, DEPTH_AXIS]
    y = P[:, 1]
    wx = (x - x.min()) / (x.max() - x.min())
    wd = (d - d.min()) / (d.max() - d.min())
    if method == "plane":
        ring = (wx < 0.15) | (wx > 0.85) | (wd < 0.15) | (wd > 0.85)
        A = np.c_[np.ones(ring.sum()), x[ring], d[ring]]
        coef, *_ = np.linalg.lstsq(A, y[ring], rcond=None)
        return coef[0] + coef[1] * x + coef[2] * d
    if method == "quad":
        ring = (wx < 0.2) | (wx > 0.8)
        A = np.c_[np.ones(ring.sum()), x[ring], x[ring] ** 2]
        coef, *_ = np.linalg.lstsq(A, y[ring], rcond=None)
        return coef[0] + coef[1] * x + coef[2] * x ** 2
    if method == "edge":
        ring = (wx < 0.1) | (wx > 0.9)
        return np.full(len(x), np.median(y[ring]))
    raise ValueError(method)


def profile(P, method, bin_width=10.0):
    x = P[:, WIDTH_AXIS]
    d = P[:, DEPTH_AXIS]
    y = P[:, 1]
    dc = 0.5 * (d.min() + d.max())
    band = np.abs(d - dc) <= 30
    xb, yb, db = x[band], y[band], d[band]
    # reference is fitted on the full tread, not only the band
    ref_all = fit_reference(P, method)
    ref_band = np.interp(np.arange(band.sum()), np.arange(len(P)), ref_all[band]) if False else None
    # fit reference on the band's x/d with the same ring criterion is unstable;
    # use the full-tread coefficients evaluated on the band
    if method == "plane":
        x_all = x; d_all = d; y_all = y
        wx = (x_all - x_all.min()) / (x_all.max() - x_all.min())
        wd = (d_all - d_all.min()) / (d_all.max() - d_all.min())
        ring = (wx < 0.15) | (wx > 0.85) | (wd < 0.15) | (wd > 0.85)
        A = np.c_[np.ones(ring.sum()), x_all[ring], d_all[ring]]
        coef, *_ = np.linalg.lstsq(A, y_all[ring], rcond=None)
        ref = coef[0] + coef[1] * xb + coef[2] * db
    elif method == "quad":
        x_all, y_all = x, y
        wx = (x_all - x_all.min()) / (x_all.max() - x_all.min())
        ring = (wx < 0.2) | (wx > 0.8)
        A = np.c_[np.ones(ring.sum()), x_all[ring], x_all[ring] ** 2]
        coef, *_ = np.linalg.lstsq(A, y_all[ring], rcond=None)
        ref = coef[0] + coef[1] * xb + coef[2] * xb ** 2
    else:
        wx = (x - x.min()) / (x.max() - x.min())
        ring = (wx < 0.1) | (wx > 0.9)
        ref = np.full(len(xb), np.median(y[ring]))
    h = ref - yb  # positive where the tread is below the reference
    bins = np.arange(xb.min(), xb.max() + bin_width, bin_width)
    idx = np.digitize(xb, bins)
    xs, hs = [], []
    for b in range(len(bins)):
        v = h[idx == b]
        if len(v) >= 15:
            xs.append(bins[b - 1] + bin_width / 2 if b > 0 else bins[0])
            hs.append(np.median(v))
    xs, hs = np.array(xs), np.array(hs)
    # light smoothing to suppress photogrammetric spikes
    if len(hs) >= 5:
        hs = np.array([np.median(hs[max(0, i - 1):min(len(hs), i + 2)])
                       for i in range(len(hs))])
    return xs, hs


def main():
    mesh = load_scaled(MODEL)
    rows = []
    curves = {}
    for k, y0 in enumerate(Y_VALUES, start=1):
        box = segment(mesh, y0)
        if box is None:
            print("skip tread", k)
            continue
        pts, _ = mesh.sample(2_000_000, return_index=True)
        mask = ((pts[:, DEPTH_AXIS] >= box["depth_min"] - 5) &
                (pts[:, DEPTH_AXIS] <= box["depth_max"] + 5) &
                (pts[:, WIDTH_AXIS] >= box["width_min"] - 5) &
                (pts[:, WIDTH_AXIS] <= box["width_max"] + 5) &
                (np.abs(pts[:, 1] - y0) <= 10))
        P = pts[mask]
        if len(P) < 500:
            print("skip tread", k, "few points")
            continue
        curves[k] = {}
        for method in ("plane", "quad", "edge"):
            xs, hs = profile(P, method)
            if len(xs) < 10:
                continue
            curves[k][method] = (xs, hs)
            # central wear metric: median of the central 50% of the profile
            xc = 0.5 * (xs.min() + xs.max())
            central = np.abs(xs - xc) <= 0.25 * (xs.max() - xs.min())
            h_center = float(np.median(hs[central]))
            rows.append({"tread": k, "y0_mm": y0, "reference": method,
                         "h_center_mm": h_center,
                         "h_min_mm": float(np.min(hs)), "h_max_mm": float(np.max(hs)),
                         "width_min_mm": box["width_min"], "width_max_mm": box["width_max"]})

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "roman_reference_sensitivity.csv", index=False)
    # ratio of inferred N across reference models, per tread (N proportional to h)
    ratios = []
    for k, g in df.groupby("tread"):
        pos = g[g.h_center_mm > 0]
        if len(pos) >= 2:
            ratios.append({"tread": k, "N_ratio_max_min": float(pos.h_center_mm.max() / pos.h_center_mm.min())})
    pd.DataFrame(ratios).to_csv(OUT / "roman_N_ratio.csv", index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9),
                             gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    for k, methods in curves.items():
        if "plane" in methods:
            xs, hs = methods["plane"]
            ax.plot(xs - np.mean(xs), hs, lw=1.0, label=f"tread {k}")
    ax.axhline(0, color="0.5", lw=0.6, ls=":")
    ax.set_xlabel("lateral position (mm, centred)")
    ax.set_ylabel("wear depth (mm)")
    ax.legend(frameon=False, ncol=2, fontsize=6.5)
    ax.set_title("a  Profiles, plane reference", loc="left")

    ax = axes[1]
    widths = 0.25
    for i, method in enumerate(("plane", "quad", "edge")):
        sub = df[df.reference == method]
        ax.bar(sub.tread + (i - 1) * widths, sub.h_center_mm, width=widths,
               label=method)
    ax.axhline(0, color="0.5", lw=0.6, ls=":")
    ax.set_xlabel("tread")
    ax.set_ylabel("central wear depth (mm)")
    ax.legend(frameon=False, fontsize=6.5)
    ax.set_title("b  Reference-model sensitivity", loc="left")
    for a in axes:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "roman_case.png", dpi=300)
    fig.savefig(OUT / "roman_case.pdf")
    print(df.to_string(index=False))
    print(pd.DataFrame(ratios).to_string(index=False))


if __name__ == "__main__":
    main()
