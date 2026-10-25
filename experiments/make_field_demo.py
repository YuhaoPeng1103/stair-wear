"""Create a synthetic 'field demo' dataset to test the end-to-end pipeline.

This is NOT field data. It is a known-truth integration test for
experiments/field_pipeline.py. The generated files live in data/field_demo/ and
are labelled as synthetic in the metadata.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import model as M

OUT = ROOT / "data" / "field_demo"
OUT.mkdir(parents=True, exist_ok=True)

GEOM = M.Geometry(direction_offset_mm=36.0)
rng = np.random.default_rng(42)

STAIRS = [
    # stair_id, role, N, delta, lanes, material, install_year, service_days
    ("C1", "calibration", 2.0e5, 0.20, 2, "marble", 2010, 5000),
    ("C2", "calibration", 5.0e5, 0.20, 2, "granite", 2012, 4500),
    ("C3", "calibration", 1.0e6, 0.20, 2, "terrazzo", 2008, 6000),
    ("A1", "application", 5.0e5, -0.30, 1, "limestone", 1900, 40000),
    ("A2", "application", 1.0e6, 0.30, 2, "marble", 1850, 60000),
    ("A3", "application", 2.0e6, 0.00, 2, "granite", 1920, 35000),
]

meta_rows = []
profile_rows = []
for i, (sid, role, N, delta, lanes, material, year, days) in enumerate(STAIRS):
    sigma = 0.1 if role == "calibration" else 0.2
    truth = M.Truth(log10_N=np.log10(N), delta=delta, log10_kappa=-6.0,
                    lanes=lanes, env_mm=0.0, quad_mm=0.0)
    clean = M.forward(GEOM.x, truth.log10_N, truth.delta, truth.log10_kappa,
                      truth.lanes, env_mm=0.0, quad_mm=0.0, geom=GEOM)
    for tread in range(1, 6):
        for prof in ("L1", "L2", "L3"):
            h = clean + rng.normal(0.0, sigma, size=len(GEOM.x))
            for x, val in zip(GEOM.x, h):
                profile_rows.append({
                    "stair_id": sid, "tread_index": tread, "profile_id": prof,
                    "x_mm": round(float(x), 3), "h_mm": round(float(val), 4),
                    "valid": 1,
                })
    meta_rows.append({
        "stair_id": sid, "role": role, "install_year": year,
        "service_days": days, "mean_daily_footfall": N / days,
        "footfall_rel_sd": 0.10, "material": material,
        "hardness_leeb": 600 if material == "granite" else 300,
        "lanes": lanes, "scan_accuracy_mm": sigma,
        "notes": "SYNTHETIC DEMO — not field data",
    })

pd.DataFrame(meta_rows).to_csv(OUT / "metadata.csv", index=False)
pd.DataFrame(profile_rows).to_csv(OUT / "profiles.csv", index=False)
print(f"wrote {OUT/'metadata.csv'} and {OUT/'profiles.csv'}")
print(pd.DataFrame(meta_rows)[["stair_id", "role", "lanes", "mean_daily_footfall"]])
