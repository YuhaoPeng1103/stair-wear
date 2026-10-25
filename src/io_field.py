"""Field-data ingestion and profile extraction.

The field schema is defined in docs/data_collection_SOP.md:

  metadata.csv:
    stair_id, role, install_year, service_days, mean_daily_footfall,
    footfall_rel_sd, material, hardness_leeb, lanes, scan_accuracy_mm

  profiles.csv:
    stair_id, tread_index, profile_id, x_mm, h_mm, valid

`h_mm` is wear depth (positive downward) relative to the reconstructed
reference surface. If only raw heights `z_mm` are available, use
`extract_profile_from_points` or fit the reference with `reference_surface`.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REQUIRED_METADATA = [
    "stair_id", "role", "install_year", "service_days",
    "mean_daily_footfall", "footfall_rel_sd", "material",
    "hardness_leeb", "lanes", "scan_accuracy_mm",
]
REQUIRED_PROFILES = ["stair_id", "tread_index", "profile_id", "x_mm", "h_mm", "valid"]


def load_metadata(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = [c for c in REQUIRED_METADATA if c not in df.columns]
    if missing:
        raise ValueError(f"metadata.csv missing columns: {missing}")
    df["role"] = df["role"].astype(str).str.lower().str.strip()
    if not set(df["role"]).issubset({"calibration", "application"}):
        raise ValueError("role must be 'calibration' or 'application'")
    df["N_known"] = df["mean_daily_footfall"] * df["service_days"]
    return df


def load_profiles(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = [c for c in REQUIRED_PROFILES if c not in df.columns]
    if missing:
        raise ValueError(f"profiles.csv missing columns: {missing}")
    df = df[df["valid"].astype(int) == 1].copy()
    return df


def reference_surface(x: np.ndarray, h: np.ndarray, edge_fraction: float = 0.2,
                      degree: int = 2) -> tuple[np.ndarray, np.poly1d]:
    """Fit a smooth reference to the unworn outer edges and return wear depth."""
    order = np.argsort(x)
    x, h = x[order], h[order]
    n = len(x)
    k = max(3, int(round(edge_fraction * n)))
    mask = np.zeros(n, dtype=bool)
    mask[:k] = True
    mask[-k:] = True
    coef = np.polyfit(x[mask], h[mask], degree)
    ref = np.polyval(coef, x)
    # positive wear depth: reference minus measured height
    return ref - h, np.poly1d(coef)


def resample(x: np.ndarray, y: np.ndarray, grid: np.ndarray) -> np.ndarray:
    order = np.argsort(x)
    return np.interp(grid, x[order], y[order], left=np.nan, right=np.nan)


def mean_profile_for_stair(profiles: pd.DataFrame, stair_id: str,
                           grid: np.ndarray) -> tuple[np.ndarray, float]:
    """Average all valid longitudinal profiles of a staircase onto a common grid."""
    sub = profiles[profiles["stair_id"] == stair_id]
    if sub.empty:
        raise ValueError(f"no profiles for {stair_id}")
    stack = []
    for _, g in sub.groupby(["tread_index", "profile_id"]):
        g = g.sort_values("x_mm")
        stack.append(resample(g["x_mm"].to_numpy(float), g["h_mm"].to_numpy(float), grid))
    arr = np.vstack(stack)
    mean = np.nanmean(arr, axis=0)
    # fill gaps at the edges by nearest valid value
    idx = np.arange(len(mean))
    good = np.isfinite(mean)
    if not good.all():
        mean = np.interp(idx, idx[good], mean[good])
    noise = float(np.nanstd(arr - mean[None, :], ddof=1))
    return mean, noise


def extract_profile_from_points(points: pd.DataFrame, y_center: float,
                                y_halfwidth: float, grid: np.ndarray,
                                bin_width: float = 2.0) -> np.ndarray:
    """Extract a lateral profile from a point cloud with columns x,y,z (mm).

    Points within a longitudinal band around y_center are binned in x; the
    median z per bin is used as a robust surface height. The returned array is
    the height on `grid`, without reference-surface subtraction.
    """
    sub = points[(points["y"] >= y_center - y_halfwidth) &
                 (points["y"] <= y_center + y_halfwidth)]
    if sub.empty:
        raise ValueError("no points in the requested y band")
    bins = np.arange(grid.min(), grid.max() + bin_width, bin_width)
    idx = np.digitize(sub["x"].to_numpy(float), bins)
    z = np.full(len(grid), np.nan)
    for j, xg in enumerate(grid):
        b = np.argmin(np.abs(bins - xg))
        vals = sub.loc[idx == b, "z"].to_numpy(float)
        if len(vals):
            z[j] = np.median(vals)
    good = np.isfinite(z)
    if good.sum() < 5:
        raise ValueError("too few bins with data")
    return np.interp(grid, grid[good], z[good])


def validate_dataset(meta: pd.DataFrame, profiles: pd.DataFrame) -> list[str]:
    problems = []
    cal = meta[meta["role"] == "calibration"]
    app = meta[meta["role"] == "application"]
    if len(cal) < 3:
        problems.append(f"only {len(cal)} calibration staircases; protocol requires >=3")
    if len(app) < 3:
        problems.append(f"only {len(app)} application staircases; protocol requires >=3")
    for sid in meta["stair_id"]:
        if sid not in set(profiles["stair_id"]):
            problems.append(f"no profiles for {sid}")
    if (meta["N_known"] <= 0).any():
        problems.append("non-positive cumulative footfall in metadata")
    return problems
