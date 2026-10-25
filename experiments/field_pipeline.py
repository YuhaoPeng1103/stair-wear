"""End-to-end field pipeline: metadata + profiles -> calibration -> inference.

Usage:
    python3 experiments/field_pipeline.py \
        --metadata data/field/metadata.csv \
        --profiles data/field/profiles.csv \
        --out results/field

The pipeline expects the schema in docs/data_collection_SOP.md. It is designed
to be run after the field campaign; until then, `data/field/` contains a
template with two example rows (clearly marked as examples, not measurements).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import model as M
from calibration import estimate_kappa, fit_kappa_known_N
from inference import profile_parameter, summarize_profile
from io_field import (load_metadata, load_profiles, mean_profile_for_stair,
                      validate_dataset)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--metadata", default=str(ROOT / "data" / "field" / "metadata.csv"))
    p.add_argument("--profiles", default=str(ROOT / "data" / "field" / "profiles.csv"))
    p.add_argument("--out", default=str(ROOT / "results" / "field"))
    p.add_argument("--max-lanes", type=int, default=4)
    p.add_argument("--sigma-default", type=float, default=0.2)
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    geom = M.Geometry(direction_offset_mm=36.0)
    grid = geom.x

    meta = load_metadata(args.metadata)
    profiles = load_profiles(args.profiles)
    problems = validate_dataset(meta, profiles)
    for msg in problems:
        print("[warning]", msg)

    # ------------------------------------------------------------------
    # Calibration staircases
    # ------------------------------------------------------------------
    cal_rows = []
    cal_meta = meta[meta["role"] == "calibration"]
    cal_profiles, cal_N, cal_lanes = [], [], []
    for _, row in cal_meta.iterrows():
        sid = row["stair_id"]
        try:
            h, noise = mean_profile_for_stair(profiles, sid, grid)
        except Exception as exc:  # noqa: BLE001
            print(f"[skip] {sid}: {exc}")
            continue
        sigma = float(row.get("scan_accuracy_mm", args.sigma_default) or args.sigma_default)
        sigma = max(sigma, noise if np.isfinite(noise) else 0.0, 1e-3)
        lanes = int(row["lanes"]) if pd.notna(row["lanes"]) else 2
        theta, lp = fit_kappa_known_N(grid, h, sigma, float(row["N_known"]),
                                      lanes, geom, M.Prior())
        logk = float(theta[1])
        # propagate the uncertainty in N_known into kappa
        n_rel_sd = float(row.get("footfall_rel_sd", 0.0) or 0.0)
        logk_sd = float(max(np.hypot(0.02, n_rel_sd / np.log(10.0)), 0.01))
        cal_rows.append({
            "stair_id": sid, "N_known": float(row["N_known"]),
            "lanes": lanes, "delta": float(theta[0]),
            "logk": logk, "kappa": float(10.0 ** logk),
            "logk_sd": logk_sd, "sigma_obs": sigma,
            "material": row["material"], "log_post": lp,
        })
        cal_profiles.append(h); cal_N.append(float(row["N_known"])); cal_lanes.append(lanes)

    if not cal_rows:
        raise SystemExit("No usable calibration staircases. See data/field/README.md.")

    cal = estimate_kappa(cal_profiles, cal_N, args.sigma_default, cal_lanes, geom,
                         M.Prior())
    # combine the per-stair standard error with the mean stated N uncertainty
    n_sd_mean = float(np.mean([r["logk_sd"] for r in cal_rows]))
    logk_sd = float(max(np.hypot(cal.logk_sd, n_sd_mean), 0.03))
    prior_cal = M.Prior(logk_mean=cal.logk_mean, logk_sd=logk_sd)
    print(f"[calibration] log10 kappa = {cal.logk_mean:.3f} +/- {logk_sd:.3f} "
          f"({len(cal_rows)} stairs)")

    # ------------------------------------------------------------------
    # Application staircases
    # ------------------------------------------------------------------
    app_rows = []
    grid_logN = np.linspace(3.0, 9.0, 81)
    for _, row in meta[meta["role"] == "application"].iterrows():
        sid = row["stair_id"]
        try:
            h, noise = mean_profile_for_stair(profiles, sid, grid)
        except Exception as exc:  # noqa: BLE001
            print(f"[skip] {sid}: {exc}")
            continue
        sigma = float(row.get("scan_accuracy_mm", args.sigma_default) or args.sigma_default)
        sigma = max(sigma, noise if np.isfinite(noise) else 0.0, 1e-3)

        # lane selection by BIC under the calibrated prior
        bic_rows = M.bic_for_lanes(grid, h, sigma, geom, prior_cal,
                                   use_quad=True, max_lanes=args.max_lanes)
        best = min(bic_rows, key=lambda r: r["BIC"])
        lanes = int(best["lanes"])
        theta_map, lp = M.map_estimate(grid, h, sigma, lanes, geom, prior_cal, True)
        res = profile_parameter(grid, h, sigma, lanes, geom, prior_cal, "logN",
                                grid_logN, use_quad=True, theta_start=theta_map)
        s = summarize_profile(res, np.nan)
        app_rows.append({
            "stair_id": sid, "lanes_selected": lanes,
            "N_estimate": float(10.0 ** res.estimate),
            "N_ci_low": float(10.0 ** res.ci_low),
            "N_ci_high": float(10.0 ** res.ci_high),
            "logN_estimate": float(res.estimate),
            "logN_width": float(res.width),
            "delta_estimate": float(theta_map[1]),
            "logk_estimate": float(theta_map[2]),
            "sigma_obs": sigma, "material": row["material"],
        })

    pd.DataFrame(cal_rows).to_csv(outdir / "calibration_stairs.csv", index=False)
    pd.DataFrame(app_rows).to_csv(outdir / "application_stairs.csv", index=False)
    summary = {
        "n_calibration": len(cal_rows),
        "n_application": len(app_rows),
        "logk_mean": cal.logk_mean,
        "logk_sd_used": logk_sd,
        "calibration_per_stair": cal.per_stair,
        "problems": problems,
    }
    (outdir / "field_summary.json").write_text(json.dumps(summary, indent=2),
                                               encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
