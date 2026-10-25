"""Calibration of the wear coefficient on staircases with known traffic.

The calibration staircase must have:
  * a known cumulative footfall count N_cal (door counter, turnstile, or a
    documented visitor record), and
  * a measured tread profile h_cal(x).

Because the forward model is linear in kappa, a profile with known N gives a
direct estimate of kappa once the shape parameters (delta, lane count) are
fitted. Combining several calibration staircases yields a log-normal posterior
for kappa that can be plugged into the inference of an unlabelled historic
staircase.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from model import Geometry, Prior, forward, log_prior


@dataclass
class CalibrationResult:
    logk_mean: float
    logk_sd: float
    logk_values: np.ndarray
    n_stairs: int
    per_stair: list[dict]


def fit_kappa_known_N(x: np.ndarray, h_obs: np.ndarray, sigma_obs: float,
                      N_known: float, lanes: int, geom: Geometry,
                      prior: Prior, use_quad: bool = False) -> tuple[np.ndarray, float]:
    """MAP fit of [delta, logk, env, (quad)] with N fixed at N_known."""
    logN = np.log10(N_known)

    def neg(theta: np.ndarray) -> float:
        delta, logk, env = theta[:3]
        quad = theta[3] if use_quad else 0.0
        if not (prior.delta_bounds[0] <= delta <= prior.delta_bounds[1]):
            return 1e12
        if not (prior.logk_bounds[0] <= logk <= prior.logk_bounds[1]):
            return 1e12
        if env < 0 or (use_quad and quad < 0):
            return 1e12
        pred = forward(x, logN, delta, logk, lanes, env_mm=env,
                       quad_mm=quad, geom=geom)
        resid = (h_obs - pred) / sigma_obs
        lp = -0.5 * np.sum(resid ** 2) - 0.5 * (env / prior.env_scale) ** 2
        if use_quad:
            lp += -0.5 * (quad / prior.quad_scale) ** 2
        return -lp

    starts = [np.array([0.0, -6.0, 0.0, 0.0])[: (4 if use_quad else 3)],
              np.array([0.4, -5.5, 0.1, 0.2])[: (4 if use_quad else 3)]]
    bounds = [prior.delta_bounds, prior.logk_bounds, (0.0, 3.0)]
    if use_quad:
        bounds.append((0.0, 2.0))
    best = None
    for s in starts:
        res = minimize(neg, s, method="L-BFGS-B", bounds=bounds,
                       options={"maxiter": 3000, "ftol": 1e-12})
        if best is None or res.fun < best.fun:
            best = res
    assert best is not None
    return np.asarray(best.x), float(-best.fun)


def estimate_kappa(profiles: list[np.ndarray], N_known: list[float],
                   sigma_obs: float, lanes: list[int], geom: Geometry,
                   prior: Prior | None = None, use_quad: bool = False) -> CalibrationResult:
    prior = prior or Prior()
    logks = []
    per_stair = []
    for h, N, L in zip(profiles, N_known, lanes):
        theta, lp = fit_kappa_known_N(geom.x, h, sigma_obs, N, L, geom, prior, use_quad)
        logk = float(theta[1])
        logks.append(logk)
        per_stair.append({"N_known": float(N), "lanes": int(L),
                          "delta": float(theta[0]), "logk": logk,
                          "kappa": float(10.0 ** logk), "log_post": lp})
    logks_arr = np.asarray(logks)
    mean = float(logks_arr.mean())
    # Standard error of the mean in log space; floor prevents overconfidence.
    sd = float(max(logks_arr.std(ddof=1) / np.sqrt(len(logks_arr)), 0.01))
    return CalibrationResult(mean, sd, logks_arr, len(logks), per_stair)


def propagate_to_N(h_peak: float, logk_mean: float, logk_sd: float,
                   shape_peak: float = 1.0) -> dict:
    """First-order propagation from kappa uncertainty to N for a peak measurement."""
    logN = np.log10(h_peak) - logk_mean - np.log10(shape_peak)
    # dN/N = ln(10) * sigma_logk
    rel_sd = np.log(10.0) * logk_sd
    return {
        "N_median": float(10.0 ** logN),
        "N_rel_sd": float(rel_sd),
        "N_ci95": [float(10.0 ** (logN - 1.96 * logk_sd)),
                    float(10.0 ** (logN + 1.96 * logk_sd))],
    }
