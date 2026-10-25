"""Evaluate a literature anchor for kappa on the synthetic benchmark.

Anchor: Borobudur Conservation Office projection of 1.8 cm wear per 50 million
visitors -> kappa_eff = 3.6e-7 mm per visitor (log10 = -6.444). This is an
average over the staircase, not a single-point calibration, and the material
(andesite) differs from the simulation. We therefore evaluate it as a prior,
not as a calibration.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))
import model as M
from run_simulation import GEOM, BASE, profile_case

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

LOGK_ANCHOR = float(np.log10(18.0 / 5.0e7))  # -6.4437

rows = []
for sd in (0.3, 0.5):
    prior = M.Prior(logk_mean=LOGK_ANCHOR, logk_sd=sd)
    for rep in range(20):
        x, h = M.simulate_profile(BASE, GEOM, 0.2, seed=12000 + rep)
        s, _, _ = profile_case(x, h, 0.2, BASE.lanes, prior)
        rows.append({"logk_prior_mean": LOGK_ANCHOR, "logk_prior_sd": sd,
                     "rep": rep, "N_map_rel": s["N_map_rel"],
                     "N_rel_error": s["N_rel_error"], "ci_low": s["ci_low"],
                     "ci_high": s["ci_high"], "width": s["width"],
                     "covered": s["covered"], "identified": s["identified"]})
        print(f"sd={sd} rep={rep+1}/20 width={s['width']:.2f} "
              f"Nerr={s['N_map_rel']*100:.1f}% covered={s['covered']}")

df = pd.DataFrame(rows)
df.to_csv(OUT / "literature_anchor.csv", index=False)
summary = df.groupby("logk_prior_sd").agg(
    width_median=("width", "median"),
    Nerr_median=("N_map_rel", "median"),
    Nerr_q25=("N_map_rel", lambda x: x.quantile(.25)),
    Nerr_q75=("N_map_rel", lambda x: x.quantile(.75)),
    coverage=("covered", "mean"),
)
print(summary.to_string())
summary.to_csv(OUT / "literature_anchor_summary.csv")
