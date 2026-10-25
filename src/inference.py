"""Profile-likelihood inference for the stair-wear inverse problem.

Why profile likelihood rather than MCMC?
----------------------------------------
When kappa is free, the posterior has an exact multiplicative ridge:
h(x) depends on N and kappa only through the product kappa*N. A random-walk
MCMC mixes poorly along that ridge, and posterior summaries can look
artificially precise. Profiling one parameter at a time shows the ridge
directly and yields honest confidence intervals: if the profile is flat, the
parameter is not identified.

The same machinery handles a calibration prior on kappa: the prior enters the
objective as a penalty, and the profile is a penalised profile likelihood.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

from model import Geometry, Prior, log_posterior, map_estimate


@dataclass
class ProfileResult:
    param: str
    grid: np.ndarray
    profile: np.ndarray
    map_theta: np.ndarray
    best_value: float
    ci_low: float
    ci_high: float
    bounded_low: bool
    bounded_high: bool
    covered: bool | None = None

    @property
    def estimate(self) -> float:
        return float(self.grid[int(np.argmax(self.profile))])

    @property
    def width(self) -> float:
        return float(self.ci_high - self.ci_low)


PARAM_INDEX = {"logN": 0, "delta": 1, "logk": 2, "env": 3, "quad": 4}


def _bounds_for(prior: Prior, use_quad: bool) -> list[tuple[float, float]]:
    b = [prior.logN_bounds, prior.delta_bounds, prior.logk_bounds, (0.0, 3.0)]
    if use_quad:
        b.append((0.0, 2.0))
    return b


def _neg_logpost(theta, x, h, sigma, lanes, geom, prior, use_quad):
    lp = log_posterior(theta, x, h, sigma, lanes, geom, prior, use_quad)
    return 1e12 if not np.isfinite(lp) else -lp


def profile_parameter(x: np.ndarray, h: np.ndarray, sigma_obs: float, lanes: int,
                      geom: Geometry, prior: Prior, param: str,
                      grid: np.ndarray, use_quad: bool = False,
                      theta_start: np.ndarray | None = None) -> ProfileResult:
    """Profile the negative log posterior over one parameter."""
    idx = PARAM_INDEX[param]
    bounds = _bounds_for(prior, use_quad)
    if theta_start is None:
        theta_start, _ = map_estimate(x, h, sigma_obs, lanes, geom, prior, use_quad)
    theta_start = np.asarray(theta_start, float)

    free_idx = [i for i in range(len(theta_start)) if i != idx]
    free_bounds = [bounds[i] for i in free_idx]
    lower = np.array([b[0] for b in free_bounds])
    upper = np.array([b[1] for b in free_bounds])
    prof = np.empty(len(grid))
    sol = np.empty((len(grid), len(free_idx)))

    # Walk outward from the MAP grid point and warm-start each fit from the
    # neighbouring solution. This is the standard trick that makes profile
    # likelihood cheap; without it the flat ridge costs thousands of
    # L-BFGS-B iterations per grid point.
    order = sorted(range(len(grid)), key=lambda k: abs(grid[k] - theta_start[idx]))
    last = np.clip(theta_start[free_idx], lower, upper)
    for rank, k in enumerate(order):
        value = grid[k]
        if not (bounds[idx][0] <= value <= bounds[idx][1]):
            prof[k] = -np.inf
            sol[k] = last
            continue

        def objective(v):
            theta = theta_start.copy()
            theta[idx] = value
            theta[free_idx] = v
            return _neg_logpost(theta, x, h, sigma_obs, lanes, geom, prior, use_quad)

        starts = [last]
        if rank < 2:  # at the MAP, also try the unconstrained MAP start
            starts.append(np.clip(theta_start[free_idx], lower, upper))
        best = None
        for s in starts:
            res = minimize(objective, s, method="L-BFGS-B", bounds=free_bounds,
                           options={"maxiter": 400, "ftol": 1e-10})
            if best is None or res.fun < best.fun:
                best = res
        prof[k] = -best.fun
        last = np.clip(best.x, lower, upper)
        sol[k] = last

    prof = np.asarray(prof)
    i_best = int(np.argmax(prof))
    best_value = float(prof[i_best])
    threshold = best_value - 0.5 * chi2.ppf(0.95, 1)

    # contiguous region around the maximum where prof >= threshold
    mask = prof >= threshold
    lo = i_best
    while lo > 0 and mask[lo - 1]:
        lo -= 1
    hi = i_best
    while hi < len(grid) - 1 and mask[hi + 1]:
        hi += 1

    return ProfileResult(
        param=param, grid=np.asarray(grid), profile=prof,
        map_theta=theta_start, best_value=best_value,
        ci_low=float(grid[lo]), ci_high=float(grid[hi]),
        bounded_low=bool(lo > 0), bounded_high=bool(hi < len(grid) - 1),
    )


def summarize_profile(res: ProfileResult, truth: float) -> dict:
    return {
        "estimate": res.estimate,
        "ci_low": res.ci_low,
        "ci_high": res.ci_high,
        "width": res.width,
        "bounded_low": res.bounded_low,
        "bounded_high": res.bounded_high,
        "covered": bool(res.ci_low <= truth <= res.ci_high),
    }


def likelihood_grid_2d(x: np.ndarray, h: np.ndarray, sigma_obs: float, lanes: int,
                       geom: Geometry, prior: Prior,
                       logN_grid: np.ndarray, logk_grid: np.ndarray,
                       theta_map: np.ndarray, use_quad: bool = False) -> np.ndarray:
    """2-D penalised log-posterior slice with the other parameters at the MAP."""
    out = np.empty((len(logN_grid), len(logk_grid)))
    theta = np.asarray(theta_map, float).copy()
    for i, logN in enumerate(logN_grid):
        for j, logk in enumerate(logk_grid):
            theta[0] = logN
            theta[2] = logk
            out[i, j] = -_neg_logpost(theta, x, h, sigma_obs, lanes, geom, prior, use_quad)
    return out
