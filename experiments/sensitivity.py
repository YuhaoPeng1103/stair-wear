"""Sensitivity of the identifiability result to model assumptions.

Varies the contact-kernel width, lane spacing, and direction offset, and
reports the 95% interval for N under (i) a free wear coefficient and (ii) a
calibrated prior. The goal is to show that the product confound is structural,
not an artefact of the baseline geometry.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import model as M
from inference import profile_parameter, summarize_profile

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
GRID = np.linspace(3.0, 9.0, 61)


def fit(geom, prior, seed):
    truth = M.Truth(log10_N=6.0, delta=0.3, log10_kappa=-6.0, lanes=2)
    x, h = M.simulate_profile(truth, geom, 0.2, seed=seed)
    theta, _ = M.map_estimate(x, h, 0.2, 2, geom, prior, False)
    res = profile_parameter(x, h, 0.2, 2, geom, prior, "logN", GRID,
                            use_quad=False, theta_start=theta)
    s = summarize_profile(res, 6.0)
    return {
        "width": s["width"], "ci_low": s["ci_low"], "ci_high": s["ci_high"],
        "covered": s["covered"], "N_map_rel": float(10 ** (theta[0] - 6.0) - 1.0),
        "identified": bool(s["width"] < 2.0),
    }


rows = []
# 1) contact-kernel width: foot_sigma varied, pressure sigma fixed at 10 mm
for foot_sigma in (8.0, 15.0, 22.0, 28.0):
    geom = M.Geometry(foot_sigma_mm=foot_sigma)
    sigma_eff = float(np.hypot(foot_sigma, geom.pressure_sigma_mm))
    for rep in range(4):
        for name, prior in (("free", M.Prior()), ("calibrated", M.Prior(logk_mean=-6.0, logk_sd=0.05))):
            s = fit(geom, prior, 20000 + rep)
            rows.append({"parameter": "sigma_eff", "value": sigma_eff,
                         "scenario": name, "rep": rep, **s})
    print("sigma", sigma_eff)

# 2) lane spacing
for spacing in (45.0, 65.0, 85.0, 105.0):
    geom = M.Geometry(lane_spacing_mm=spacing)
    for rep in range(4):
        for name, prior in (("free", M.Prior()), ("calibrated", M.Prior(logk_mean=-6.0, logk_sd=0.05))):
            s = fit(geom, prior, 21000 + rep)
            rows.append({"parameter": "lane_spacing", "value": spacing,
                         "scenario": name, "rep": rep, **s})
    print("spacing", spacing)

# 3) direction offset
for offset in (18.0, 24.0, 36.0, 48.0):
    geom = M.Geometry(direction_offset_mm=offset)
    for rep in range(4):
        for name, prior in (("free", M.Prior()), ("calibrated", M.Prior(logk_mean=-6.0, logk_sd=0.05))):
            s = fit(geom, prior, 22000 + rep)
            rows.append({"parameter": "direction_offset", "value": offset,
                         "scenario": name, "rep": rep, **s})
    print("offset", offset)

df = pd.DataFrame(rows)
df.to_csv(OUT / "sensitivity.csv", index=False)
summary = df.groupby(["parameter", "value", "scenario"]).agg(
    width_median=("width", "median"),
    width_min=("width", "min"), width_max=("width", "max"),
    N_err_median=("N_map_rel", "median"),
    N_err_absmax=("N_map_rel", lambda x: np.max(np.abs(x))),
    coverage=("covered", "mean"),
).reset_index()
summary.to_csv(OUT / "sensitivity_summary.csv", index=False)
print(summary.to_string(index=False))
