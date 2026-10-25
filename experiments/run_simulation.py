"""Profile-likelihood simulation study for calibrating stair-use inference.

Run from the project root:
    python3 experiments/run_simulation.py

The study has six blocks:
  A  identifiability: free vs weak vs strong calibration of kappa
  B  measurement-noise sweep
  C  lane-count selection by BIC
  D  direction-bias sweep
  E  calibration transfer from known-traffic stairs
  F  reference-surface (quadratic) bias

Outputs are written to results/.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import model as M
from calibration import estimate_kappa
from inference import (likelihood_grid_2d, profile_parameter, summarize_profile)

RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

GEOM = M.Geometry(direction_offset_mm=36.0)
BASE = M.Truth(log10_N=6.0, delta=0.3, log10_kappa=-6.0, lanes=2)
PRIORS = {
    "free": M.Prior(),
    "weak": M.Prior(logk_mean=-6.0, logk_sd=0.20),
    "strong": M.Prior(logk_mean=-6.0, logk_sd=0.05),
}
GRID_LOGN = np.linspace(3.0, 9.0, 61)
GRID_DELTA = np.linspace(-1.0, 1.0, 41)


def profile_case(x, h, sigma_obs, lanes, prior, use_quad=False,
                 theta_start=None):
    theta_map, lp = M.map_estimate(x, h, sigma_obs, lanes, GEOM, prior, use_quad,
                                   starts=None if theta_start is None else [theta_start])
    res = profile_parameter(x, h, sigma_obs, lanes, GEOM, prior, "logN",
                            GRID_LOGN, use_quad=use_quad, theta_start=theta_map)
    s = summarize_profile(res, np.log10(BASE.N))
    s.update({
        "N_est": float(10.0 ** res.estimate),
        "N_rel_error": float(10.0 ** res.estimate / BASE.N - 1.0),
        "N_map": float(10.0 ** theta_map[0]),
        "N_map_rel": float(10.0 ** (theta_map[0] - np.log10(BASE.N)) - 1.0),
        "delta_map": float(theta_map[1]),
        "logk_map": float(theta_map[2]),
        "env_map": float(theta_map[3]),
        "identified": bool(res.width < 2.0),
        "log_post_map": float(lp),
        "theta_map": theta_map.tolist(),
    })
    if use_quad:
        s["quad_map"] = float(theta_map[4])
    return s, res, theta_map


def main():
    t0 = time.time()
    rows_ident = []
    profile_store = {}
    grid_store = {}

    # ------------------------------------------------------------------
    # A. Identifiability
    # ------------------------------------------------------------------
    n_rep = 12
    for rep in range(n_rep):
        x, h = M.simulate_profile(BASE, GEOM, 0.2, seed=1000 + rep)
        for name, prior in PRIORS.items():
            s, res, theta_map = profile_case(x, h, 0.2, BASE.lanes, prior)
            rows_ident.append({"scenario": name, "rep": rep, **s})
            if rep == 0:
                profile_store[f"profile_{name}"] = res.profile
                profile_store[f"grid_{name}"] = res.grid
                theta_free = theta_map.copy()
                logN_grid = np.linspace(3.0, 9.0, 41)
                logk_grid = np.linspace(-9.0, -3.0, 41)
                grid_store[f"slice_{name}"] = likelihood_grid_2d(
                    x, h, 0.2, BASE.lanes, GEOM, prior, logN_grid, logk_grid,
                    theta_map)
                grid_store[f"logN_grid"] = logN_grid
                grid_store[f"logk_grid"] = logk_grid
        print(f"[A identifiability] rep {rep + 1}/{n_rep}  ({time.time() - t0:.0f}s)")
    pd.DataFrame(rows_ident).to_csv(RESULTS / "identifiability.csv", index=False)
    np.savez_compressed(RESULTS / "profiles_identifiability.npz", **profile_store)
    np.savez_compressed(RESULTS / "slices_identifiability.npz", **grid_store)

    # ------------------------------------------------------------------
    # B. Noise sweep
    # ------------------------------------------------------------------
    rows_noise = []
    for sigma in [0.05, 0.1, 0.2, 0.5, 1.0]:
        for rep in range(6):
            x, h = M.simulate_profile(BASE, GEOM, sigma, seed=3000 + rep)
            for name in ("free", "strong"):
                s, _, _ = profile_case(x, h, sigma, BASE.lanes, PRIORS[name])
                rows_noise.append({"sigma_obs": sigma, "scenario": name,
                                   "rep": rep, **s})
        print(f"[B noise] sigma={sigma}  ({time.time() - t0:.0f}s)")
    pd.DataFrame(rows_noise).to_csv(RESULTS / "noise.csv", index=False)

    # ------------------------------------------------------------------
    # C. Lane-count selection by BIC
    # ------------------------------------------------------------------
    rows_lanes = []
    for true_lanes in [1, 2, 3, 4]:
        truth = M.Truth(log10_N=6.0, delta=0.3, log10_kappa=-6.0, lanes=true_lanes)
        for rep in range(10):
            x, h = M.simulate_profile(truth, GEOM, 0.2, seed=5000 + rep)
            bic_rows = M.bic_for_lanes(x, h, 0.2, GEOM, PRIORS["strong"],
                                       use_quad=False, max_lanes=4)
            best = min(bic_rows, key=lambda r: r["BIC"])
            rows_lanes.append({"true_lanes": true_lanes, "rep": rep,
                               "selected_lanes": best["lanes"],
                               "correct": int(best["lanes"] == true_lanes),
                               "bic": json.dumps({str(r["lanes"]): r["BIC"] for r in bic_rows})})
        print(f"[C lanes] true={true_lanes}  ({time.time() - t0:.0f}s)")
    pd.DataFrame(rows_lanes).to_csv(RESULTS / "lanes.csv", index=False)

    # ------------------------------------------------------------------
    # D. Direction-bias sweep
    # ------------------------------------------------------------------
    rows_dir = []
    for delta_true in [-0.8, -0.4, 0.0, 0.4, 0.8]:
        truth = M.Truth(log10_N=6.0, delta=delta_true, log10_kappa=-6.0, lanes=2)
        for rep in range(6):
            x, h = M.simulate_profile(truth, GEOM, 0.2, seed=6000 + rep)
            theta_map, _ = M.map_estimate(x, h, 0.2, 2, GEOM, PRIORS["strong"], False)
            res = profile_parameter(x, h, 0.2, 2, GEOM, PRIORS["strong"], "delta",
                                    GRID_DELTA, theta_start=theta_map)
            rows_dir.append({"delta_true": delta_true, "rep": rep,
                             "delta_est": res.estimate,
                             "ci_low": res.ci_low, "ci_high": res.ci_high,
                             "covered": bool(res.ci_low <= delta_true <= res.ci_high),
                             "width": res.width})
        print(f"[D direction] delta={delta_true}  ({time.time() - t0:.0f}s)")
    pd.DataFrame(rows_dir).to_csv(RESULTS / "direction.csv", index=False)

    # ------------------------------------------------------------------
    # E. Calibration transfer
    # ------------------------------------------------------------------
    rows_cal = []
    N_cal = [2e5, 5e5, 1e6]
    for rep in range(12):
        profiles = []
        for j, N in enumerate(N_cal):
            truth_cal = M.Truth(log10_N=np.log10(N), delta=0.2,
                                log10_kappa=-6.0, lanes=2)
            _, h = M.simulate_profile(truth_cal, GEOM, 0.1, seed=8000 + 10 * rep + j)
            profiles.append(h)
        cal = estimate_kappa(profiles, N_cal, 0.1, [2, 2, 2], GEOM, PRIORS["free"])
        # include site-to-site material variability, not only the standard error
        sd_cal = float(np.hypot(cal.logk_sd, 0.05))
        prior_cal = M.Prior(logk_mean=cal.logk_mean, logk_sd=max(sd_cal, 0.03))
        x, h = M.simulate_profile(BASE, GEOM, 0.2, seed=9000 + rep)
        s, _, _ = profile_case(x, h, 0.2, 2, prior_cal)
        rows_cal.append({"rep": rep, "logk_true": -6.0,
                         "logk_cal": cal.logk_mean, "logk_sd_used": prior_cal.logk_sd,
                         "logk_error": cal.logk_mean + 6.0,
                         "N_rel_error": s["N_rel_error"],
                         "N_map_rel": s["N_map_rel"],
                         "N_ci_low": s["ci_low"], "N_ci_high": s["ci_high"],
                         "N_covered": s["covered"], "delta_map": s["delta_map"]})
        print(f"[E calibration] rep {rep + 1}/12  ({time.time() - t0:.0f}s)")
    pd.DataFrame(rows_cal).to_csv(RESULTS / "calibration.csv", index=False)

    # ------------------------------------------------------------------
    # F. Reference-surface bias
    # ------------------------------------------------------------------
    rows_ref = []
    for quad in [0.3, 0.6, 1.2]:
        for rep in range(6):
            truth = M.Truth(log10_N=6.0, delta=0.3, log10_kappa=-6.0, lanes=2,
                            quad_mm=quad)
            x, h = M.simulate_profile(truth, GEOM, 0.2, seed=10000 + rep,
                                       reference_bias=True)
            for use_quad in (False, True):
                s, _, _ = profile_case(x, h, 0.2, 2, PRIORS["strong"],
                                       use_quad=use_quad)
                rows_ref.append({"quad_true": quad, "use_quad": int(use_quad),
                                 "rep": rep,
                                 "N_rel_error": s["N_rel_error"],
                                 "N_covered": s["covered"],
                                 "quad_map": s.get("quad_map", 0.0)})
        print(f"[F reference] quad={quad}  ({time.time() - t0:.0f}s)")
    pd.DataFrame(rows_ref).to_csv(RESULTS / "reference.csv", index=False)

    summary = {
        "runtime_s": round(time.time() - t0, 1),
        "geometry": {k: v for k, v in GEOM.__dict__.items()},
        "truth": {k: v for k, v in BASE.__dict__.items()},
        "priors": {k: {"logk_mean": v.logk_mean, "logk_sd": v.logk_sd}
                   for k, v in PRIORS.items()},
        "n_rep_identifiability": n_rep,
        "files": sorted(p.name for p in RESULTS.glob("*.*")),
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
