"""Figure 9: sensitivity of the identifiability result to model assumptions."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "manuscript" / "figures"
C = {"blue": "#0072B2", "orange": "#E69F00", "red": "#D55E00", "grey": "#666666"}

plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 9, "axes.titlesize": 9,
    "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.bbox": "tight",
})

df = pd.read_csv(ROOT / "results" / "sensitivity_summary.csv")
panels = [
    ("sigma_eff", r"contact-kernel width $\sigma_{\mathrm{eff}}$ (mm)", "a"),
    ("lane_spacing", "lane spacing (mm)", "b"),
    ("direction_offset", r"direction offset $\Delta$ (mm)", "c"),
]
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))
for ax, (param, xlabel, letter) in zip(axes, panels):
    sub = df[df.parameter == param]
    for name, col, label in (("free", C["red"], "free $\\kappa$"),
                             ("calibrated", C["blue"], "calibrated")):
        g = sub[sub.scenario == name].sort_values("value")
        ax.plot(g.value, g.width_median, marker="o", ms=3.5, lw=1.2,
                color=col, label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("95% CI width of $\\log_{10}N$ (dex)")
    ax.set_yscale("log")
    ax.set_title(f"{letter}  {param.replace('_', ' ')}", loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2, pad=1)
axes[0].legend(frameon=False)
fig.tight_layout(w_pad=1.2)
fig.savefig(FIG / "fig9_sensitivity.png", dpi=300)
fig.savefig(FIG / "fig9_sensitivity.pdf")
print("wrote", FIG / "fig9_sensitivity.png")
