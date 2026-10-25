"""Generate publication figures from the simulation outputs.

Run after experiments/run_simulation.py:
    python3 experiments/make_figures.py

Figures are written to manuscript/figures/ as PDF and PNG.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import model as M

RESULTS = ROOT / "results"
FIG = ROOT / "manuscript" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# Okabe-Ito colourblind-safe palette
C = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#666666",
    "sky": "#56B4E9",
}

plt.rcParams.update({
    "font.size": 8.5,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "legend.fontsize": 7.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png", dpi=300)
    plt.close(fig)
    print("wrote", name)


def load_csv(name):
    path = RESULTS / name
    return pd.read_csv(path) if path.exists() else None


# ---------------------------------------------------------------------------
# Figure 1: model schematic and example profiles
# ---------------------------------------------------------------------------
def fig1():
    geom = M.Geometry(direction_offset_mm=36.0)
    x = geom.x
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.4),
                             gridspec_kw={"width_ratios": [1.15, 1, 1]})

    ax = axes[0]
    ax.set_xlim(-150, 150)
    ax.set_ylim(-1.1, 1.1)
    ax.axis("off")
    ax.plot([-150, 150], [0, 0], color="black", lw=1)
    ax.add_patch(Rectangle((-150, -0.55), 300, 0.55, color="#e8e2d6",
                               ec="black", lw=0.6))
    for c in M.lane_centers(3, geom):
        ax.annotate("", xy=(c, 0.75), xytext=(c, 0.05),
                    arrowprops=dict(arrowstyle="-|>", color=C["blue"], lw=1))
        ax.text(c, 0.82, "foot", ha="center", fontsize=7, color=C["blue"])
    ax.annotate("ascending", xy=(60, 0.35), xytext=(78, 0.50),
                fontsize=7, color=C["red"],
                arrowprops=dict(arrowstyle="->", color=C["red"], lw=0.8))
    ax.annotate("descending", xy=(-60, 0.35), xytext=(-148, 0.50),
                fontsize=7, color=C["green"],
                arrowprops=dict(arrowstyle="->", color=C["green"], lw=0.8))
    ax.text(0, -0.85, "tread width  W", ha="center", fontsize=7)
    ax.annotate("", xy=(-150, -0.75), xytext=(150, -0.75),
                arrowprops=dict(arrowstyle="<->", color="black", lw=0.6))
    ax.set_title("a  Foot placement and contact kernels", loc="left")

    ax = axes[1]
    for L, col in zip([1, 2, 3], [C["blue"], C["orange"], C["green"]]):
        ax.plot(x, M.shape_profile(0.0, L, geom), color=col, lw=1.2,
                label=f"{L} lane{'s' if L > 1 else ''}")
    ax.set_xlabel("lateral position x (mm)")
    ax.set_ylabel("unit wear depth")
    ax.legend(frameon=False)
    ax.set_title("b  Side-by-side lanes", loc="left")

    ax = axes[2]
    for d, col in zip([-0.6, 0.0, 0.6], [C["green"], C["grey"], C["red"]]):
        ax.plot(x, M.shape_profile(d, 2, geom), color=col, lw=1.2,
                label=f"δ = {d:+.1f}")
    ax.set_xlabel("lateral position x (mm)")
    ax.set_ylabel("unit wear depth")
    ax.legend(frameon=False)
    ax.set_title("c  Direction bias", loc="left")

    for ax in axes[1:]:
        ax.set_xlim(-150, 150)
        ax.set_ylim(0, 1.05)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "fig1_model")


# ---------------------------------------------------------------------------
# Figure 2: identifiability ridge and profiles
# ---------------------------------------------------------------------------
def fig2():
    slice_file = RESULTS / "slices_identifiability.npz"
    prof_file = RESULTS / "profiles_identifiability.npz"
    if not slice_file.exists() or not prof_file.exists():
        print("skip fig2: missing files")
        return
    sl = np.load(slice_file)
    pr = np.load(prof_file)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9),
                             gridspec_kw={"width_ratios": [1.05, 1]})

    ax = axes[0]
    logN = sl["logN_grid"]
    logk = sl["logk_grid"]
    Z = sl["slice_free"]
    Z = Z - np.nanmax(Z)
    cs = ax.contourf(logN, logk, Z.T, levels=[-20, -10, -5, -2, -1, 0],
                     cmap="viridis", extend="both")
    ax.contour(logN, logk, Z.T, levels=[-1.92], colors="white", linewidths=1.2)
    ax.plot([6.0], [-6.0], marker="*", ms=9, color=C["red"], mec="white",
            label="truth")
    ax.set_xlabel(r"$\log_{10} N$")
    ax.set_ylabel(r"$\log_{10}\kappa$")
    ax.set_title("a  Free $\\kappa$: a multiplicative ridge", loc="left")
    ax.legend(frameon=False, loc="upper right")
    cb = fig.colorbar(cs, ax=ax, pad=0.02)
    cb.set_label(r"$\Delta$ log-posterior", fontsize=7)
    cb.ax.tick_params(labelsize=7)

    ax = axes[1]
    for name, col, label in [("free", C["red"], "free $\\kappa$"),
                             ("weak", C["orange"], "weak calibration"),
                             ("strong", C["blue"], "strong calibration")]:
        g = pr[f"grid_{name}"]
        p = pr[f"profile_{name}"]
        ax.plot(g, p - np.nanmax(p), color=col, lw=1.3, label=label)
    ax.axhline(-1.92, color="black", lw=0.7, ls="--")
    ax.axvline(6.0, color=C["grey"], lw=0.7, ls=":")
    ax.set_xlabel(r"$\log_{10} N$")
    ax.set_ylabel(r"$\Delta$ log-posterior")
    ax.set_ylim(-12, 0.4)
    ax.set_xlim(3, 9)
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("b  Profile likelihood", loc="left")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "fig2_identifiability")


# ---------------------------------------------------------------------------
# Figure 3: noise sweep
# ---------------------------------------------------------------------------
def fig3():
    df = load_csv("noise.csv")
    if df is None:
        print("skip fig3")
        return
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7))
    ax = axes[0]
    for name, col, label in [("free", C["red"], "free $\\kappa$"),
                             ("strong", C["blue"], "calibrated $\\kappa$")]:
        g = df[df.scenario == name].groupby("sigma_obs")["width"]
        med = g.median()
        lo = g.quantile(0.25)
        hi = g.quantile(0.75)
        ax.plot(med.index, med.values, marker="o", ms=3.5, color=col, lw=1.2,
                label=label)
        ax.fill_between(med.index, lo.values, hi.values, color=col, alpha=0.18,
                        lw=0)
    ax.set_xscale("log")
    ax.set_xlabel("profile noise $\\sigma_{obs}$ (mm)")
    ax.set_ylabel(r"95% CI width of $\log_{10} N$ (dex)")
    ax.legend(frameon=False)
    ax.set_title("a  Identified range", loc="left")

    ax = axes[1]
    sub = df[df.scenario == "strong"]
    g = sub.groupby("sigma_obs")["N_map_rel"]
    med = g.median() * 100
    lo = g.quantile(0.25) * 100
    hi = g.quantile(0.75) * 100
    ax.plot(med.index, med.values, marker="o", ms=3.5, color=C["blue"], lw=1.2)
    ax.fill_between(med.index, lo.values, hi.values, color=C["blue"], alpha=0.18,
                    lw=0)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("profile noise $\\sigma_{obs}$ (mm)")
    ax.set_ylabel("relative error in $N$ (%)")
    ax.set_title("b  Calibrated bias", loc="left")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "fig3_noise")


# ---------------------------------------------------------------------------
# Figure 4: lanes and direction
# ---------------------------------------------------------------------------
def fig4():
    lanes = load_csv("lanes.csv")
    direc = load_csv("direction.csv")
    if lanes is None or direc is None:
        print("skip fig4")
        return
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7))
    ax = axes[0]
    acc = lanes.groupby("true_lanes")["correct"].mean() * 100
    ax.bar(acc.index.astype(str), acc.values, color=C["blue"], width=0.6)
    for i, v in enumerate(acc.values):
        ax.text(i, v + 2, f"{v:.0f}%", ha="center", fontsize=7.5)
    ax.set_ylim(0, 110)
    ax.set_xlabel("true number of lanes")
    ax.set_ylabel("correct BIC selection (%)")
    ax.set_title("a  Lane count", loc="left")

    ax = axes[1]
    for d, col in zip(sorted(direc.delta_true.unique()),
                      [C["red"], C["orange"], C["grey"], C["sky"], C["blue"]]):
        sub = direc[direc.delta_true == d]
        ax.errorbar(np.full(len(sub), d), sub.delta_est,
                    yerr=[sub.delta_est - sub.ci_low, sub.ci_high - sub.delta_est],
                    fmt="o", ms=3, color=col, ecolor=col, elinewidth=0.8,
                    capsize=2, alpha=0.8)
    ax.plot([-1, 1], [-1, 1], color="black", lw=0.7, ls="--")
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_xlabel("true direction bias $\\delta$")
    ax.set_ylabel("estimated $\\delta$ (95% CI)")
    ax.set_title("b  Direction bias", loc="left")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "fig4_lanes_direction")


# ---------------------------------------------------------------------------
# Figure 5: calibration transfer
# ---------------------------------------------------------------------------
def fig5():
    df = load_csv("calibration.csv")
    if df is None:
        print("skip fig5")
        return
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7))
    ax = axes[0]
    ax.errorbar(df.index, df.logk_cal, yerr=1.96 * df.logk_sd_used,
                fmt="o", ms=3, color=C["blue"], ecolor=C["blue"], elinewidth=0.8,
                capsize=2)
    ax.axhline(-6.0, color=C["red"], lw=0.8, ls="--", label="truth")
    ax.set_xlabel("calibration replicate")
    ax.set_ylabel(r"estimated $\log_{10}\kappa$")
    ax.legend(frameon=False)
    ax.set_title("a  Calibrated wear coefficient", loc="left")

    ax = axes[1]
    ax.hist(df.N_map_rel * 100, bins=8, color=C["blue"], alpha=0.8)
    ax.axvline(0, color=C["red"], lw=0.8, ls="--")
    ax.set_xlabel("relative error in $N$ (%)")
    ax.set_ylabel("count")
    ax.set_title("b  Transfer to the target stair", loc="left")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "fig5_calibration")


# ---------------------------------------------------------------------------
# Figure 6: reference-surface bias
# ---------------------------------------------------------------------------
def fig6():
    df = load_csv("reference.csv")
    if df is None:
        print("skip fig6")
        return
    fig, ax = plt.subplots(figsize=(3.8, 2.8))
    quads = sorted(df.quad_true.unique())
    width = 0.34
    for i, use_quad in enumerate((0, 1)):
        med, iqr = [], []
        for q in quads:
            sub = df[(df.quad_true == q) & (df.use_quad == use_quad)]
            med.append(sub.N_rel_error.median() * 100)
            iqr.append((sub.N_rel_error.quantile(0.75) - sub.N_rel_error.quantile(0.25)) * 100)
        xpos = np.arange(len(quads)) + (i - 0.5) * width
        ax.bar(xpos, med, width=width, color=C["orange"] if use_quad == 0 else C["blue"],
               alpha=0.85, label="without quadratic term" if use_quad == 0 else "with quadratic term")
        ax.errorbar(xpos, med, yerr=iqr, fmt="none", ecolor="black", elinewidth=0.7,
                    capsize=2)
        for x, m in zip(xpos, med):
            ax.plot(x, m, marker="o", ms=2.5, color="black")
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(np.arange(len(quads)))
    ax.set_xticklabels([f"{q:g}" for q in quads])
    ax.set_xlabel("quadratic reference-surface bias (mm)")
    ax.set_ylabel("relative error in $N$ (%)")
    ax.legend(frameon=False, fontsize=7, loc="lower left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "fig6_reference")


# ---------------------------------------------------------------------------
# Figure 7: measurement-precision nomogram
# ---------------------------------------------------------------------------
def fig7():
    df = load_csv("noise.csv")
    if df is None:
        print("skip fig7")
        return
    fig, ax = plt.subplots(figsize=(3.6, 2.7))
    sub = df[df.scenario == "strong"]
    g = sub.groupby("sigma_obs")["width"]
    med = g.median()
    ax.loglog(med.index, med.values, marker="o", ms=4, color=C["blue"], lw=1.2)
    for target, col in [(0.1, C["green"]), (0.2, C["orange"]), (0.3, C["red"])]:
        ax.axhline(target, color=col, lw=0.8, ls="--")
        ax.text(med.index.min() * 1.05, target * 1.05, f"{target:.1f} dex",
                fontsize=7, color=col)
    ax.set_xlabel("profile noise $\\sigma_{obs}$ (mm)")
    ax.set_ylabel(r"95% CI width of $\log_{10} N$")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "fig7_precision")


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    fig7()
