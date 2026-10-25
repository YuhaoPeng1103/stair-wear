"""Inspect the three public staircase models: bounds, normal alignment, projections."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import trimesh

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "public_model" / "extracted"
OUT = ROOT / "data" / "public_model" / "inspection"
OUT.mkdir(exist_ok=True)


def load(path):
    obj = trimesh.load(path, force="scene")
    return obj.to_mesh() if isinstance(obj, trimesh.Scene) else obj


for name in ["paestum", "roman", "umm_qais"]:
    mesh = load(BASE / name / "scene.gltf")
    print(f"\n=== {name} ===")
    print("verts", len(mesh.vertices), "faces", len(mesh.faces))
    print("bounds", np.round(mesh.bounds, 3).tolist())
    print("extents", np.round(mesh.extents, 3).tolist())
    # fraction of area whose normal is nearly aligned with each axis
    normals = mesh.face_normals
    areas = mesh.area_faces
    for i, axis in enumerate("xyz"):
        frac = areas[np.abs(normals[:, i]) > 0.9].sum() / areas.sum()
        print(f"  horizontal-face area fraction (normal // {axis}): {frac:.3f}")
    pts, _ = mesh.sample(250_000, return_index=True)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, (a, b, lab) in zip(axes, [(0, 2, "X-Z"), (1, 2, "Y-Z"), (0, 1, "X-Y")]):
        ax.scatter(pts[:, a], pts[:, b], s=0.05, c=pts[:, 1], cmap="viridis", alpha=0.5)
        ax.set_xlabel("XYZ"[a]); ax.set_ylabel("XYZ"[b])
        ax.set_title(f"{name} {lab}")
        ax.set_aspect("equal", adjustable="datalim")
    fig.tight_layout()
    fig.savefig(OUT / f"{name}_projections.png", dpi=150)
    plt.close(fig)
    print("wrote", OUT / f"{name}_projections.png")
