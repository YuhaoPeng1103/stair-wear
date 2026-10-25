"""Graphical abstract for the JCH submission (13 x 5 cm, 300 dpi)."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "manuscript" / "figures"
C = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
     "red": "#D55E00", "grey": "#666666"}

plt.rcParams.update({
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7.5,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.bbox": "tight",
})

sl = np.load(ROOT / "results" / "slices_identifiability.npz")
pr = np.load(ROOT / "results" / "profiles_identifiability.npz")

fig, axes = plt.subplots(1, 3, figsize=(13 / 2.54, 5 / 2.54),
                         gridspec_kw={"width_ratios": [1, 1.15, 0.9]})

# (a) ridge
ax = axes[0]
logN, logk = sl["logN_grid"], sl["logk_grid"]
Z = sl["slice_free"]; Z = Z - np.nanmax(Z)
ax.contourf(logN, logk, Z.T, levels=[-30, -5, -2, -1, 0], cmap="Blues_r")
ax.contour(logN, logk, Z.T, levels=[-1.92], colors="black", linewidths=0.6)
ax.plot([6], [-6], marker="*", ms=6, color=C["red"])
ax.set_xlabel(r"$\log_{10}N$"); ax.set_ylabel(r"$\log_{10}\kappa$")
ax.set_title("free $\\kappa$: 6 dex", fontsize=7.5)
ax.set_xticks([3, 5, 7, 9]); ax.set_yticks([-9, -7, -5, -3])
ax.tick_params(length=2, pad=1)

# (b) profiles
ax = axes[1]
for name, col, lab in [("free", C["red"], "free"),
                       ("weak", C["orange"], "weak"),
                       ("strong", C["blue"], "calibrated")]:
    g = pr[f"grid_{name}"]; p = pr[f"profile_{name}"]
    ax.plot(g, p - np.nanmax(p), color=col, lw=1.0, label=lab)
ax.axhline(-1.92, color="black", lw=0.5, ls="--")
ax.set_xlim(3, 9); ax.set_ylim(-12, 0.3)
ax.set_xlabel(r"$\log_{10}N$"); ax.set_ylabel(r"$\Delta$ log-posterior")
ax.set_title("calibrated: 0.1 dex", fontsize=7.5)
ax.legend(frameon=False, fontsize=5.2, loc="lower left")
ax.set_xticks([3, 5, 7, 9]); ax.tick_params(length=2, pad=1)

# (c) three tiers
ax = axes[2]
labels = ["none", "literature", "direct"]
vals = [6.0, 1.1, 0.1]
notes = ["", "+172%", "-0.5%"]
cols = [C["grey"], C["orange"], C["blue"]]
bars = ax.bar(labels, vals, color=cols, width=0.6)
for b, v, note in zip(bars, vals, notes):
    x = b.get_x() + b.get_width() / 2
    ax.text(x, v + 0.12, f"{v:g}", ha="center", va="bottom", fontsize=6)
    if note:
        ax.text(x, v + 0.65, note, ha="center", va="bottom", fontsize=5.0,
                color="0.25")
ax.set_ylim(0, 7.2)
ax.set_ylabel("95% CI width (dex)")
ax.set_title("calibration tiers", fontsize=7.5)
ax.tick_params(axis="x", length=0, pad=1)
ax.tick_params(axis="y", length=2, pad=1)

for ax in axes:
    ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout(w_pad=1.2)
fig.savefig(FIG / "graphical_abstract.png", dpi=300)
fig.savefig(FIG / "graphical_abstract.pdf")
print("wrote", FIG / "graphical_abstract.png")
