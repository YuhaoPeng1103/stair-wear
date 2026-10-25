"""Forward wear model, likelihood, MAP estimation, and MCMC for stair-use inference.

Physical model
--------------
For a tread of width W, let x be the lateral coordinate. A footfall at lateral
position c removes material with a contact-pressure kernel p(x - c). If foot
placements follow a density g(c), the expected wear depth after N footfalls is

    h(x) = kappa * N * (g * p)(x) + e(x),

where kappa = k/H absorbs the wear coefficient and surface hardness, and e(x)
is a smooth environmental/reference-surface term. Directional preference enters
through the relative weight of the ascending and descending contact kernels;
side-by-side walking enters through the number of lateral lanes L.

The observable is the measured tread profile h_obs(x). The inverse problem
estimates N, the direction bias delta, L, and (optionally) kappa.

Key identifiability fact: when kappa is free, only the product kappa*N is
identified from the profile. Calibrating kappa on a staircase with known
traffic breaks the degeneracy.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np


# ---------------------------------------------------------------------------
# Geometry and default parameters
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Geometry:
    width_mm: float = 300.0
    n_points: int = 151          # 2 mm sampling across the tread
    lane_spacing_mm: float = 65.0
    direction_offset_mm: float = 36.0
    foot_sigma_mm: float = 22.0
    pressure_sigma_mm: float = 10.0

    @property
    def x(self) -> np.ndarray:
        return np.linspace(-self.width_mm / 2, self.width_mm / 2, self.n_points)

    @property
    def sigma_eff(self) -> float:
        return float(np.hypot(self.foot_sigma_mm, self.pressure_sigma_mm))


@dataclass(frozen=True)
class Truth:
    log10_N: float = 6.0        # N = 1e6 footfalls
    delta: float = 0.30         # +1 = all ascending, -1 = all descending
    log10_kappa: float = -6.0   # kappa in mm per footfall
    lanes: int = 2
    env_mm: float = 0.0
    quad_mm: float = 0.0        # quadratic reference-surface bias at the edges

    @property
    def N(self) -> float:
        return 10.0 ** self.log10_N

    @property
    def kappa(self) -> float:
        return 10.0 ** self.log10_kappa


# ---------------------------------------------------------------------------
# Forward model
# ---------------------------------------------------------------------------
def lane_centers(lanes: int, geom: Geometry) -> np.ndarray:
    if lanes <= 1:
        return np.zeros(1)
    return (np.arange(lanes) - (lanes - 1) / 2.0) * geom.lane_spacing_mm


def shape_profile(delta: float, lanes: int, geom: Geometry) -> np.ndarray:
    """Unit-peak lateral wear shape for a given direction bias and lane count."""
    x = geom.x
    centres = lane_centers(lanes, geom)
    f_a = 0.5 * (1.0 + delta)
    f_d = 1.0 - f_a
    mu_a = +geom.direction_offset_mm / 2.0
    mu_d = -geom.direction_offset_mm / 2.0
    s = geom.sigma_eff

    h = np.zeros_like(x)
    for c in centres:
        h += f_a * np.exp(-0.5 * ((x - (c + mu_a)) / s) ** 2)
        h += f_d * np.exp(-0.5 * ((x - (c + mu_d)) / s) ** 2)
    h /= max(lanes, 1)
    peak = h.max()
    if peak <= 0:
        raise ValueError("Degenerate wear shape")
    return h / peak


def forward(
    x: np.ndarray,
    log10_N: float,
    delta: float,
    log10_kappa: float,
    lanes: int,
    env_mm: float = 0.0,
    quad_mm: float = 0.0,
    geom: Geometry | None = None,
) -> np.ndarray:
    geom = geom or Geometry()
    shape = shape_profile(delta, lanes, geom)
    # shape is defined on geom.x; interpolate if x differs
    if x.shape != geom.x.shape or not np.allclose(x, geom.x):
        shape = np.interp(x, geom.x, shape)
    h = (10.0 ** log10_kappa) * (10.0 ** log10_N) * shape
    if env_mm:
        h = h + env_mm
    if quad_mm:
        half = geom.width_mm / 2.0
        h = h + quad_mm * (x / half) ** 2
    return h


# ---------------------------------------------------------------------------
# Priors and likelihood
# ---------------------------------------------------------------------------
@dataclass
class Prior:
    logN_bounds: tuple[float, float] = (3.0, 9.0)
    delta_bounds: tuple[float, float] = (-1.0, 1.0)
    logk_bounds: tuple[float, float] = (-9.0, -3.0)
    logk_mean: float | None = None      # calibration prior; None -> uniform
    logk_sd: float = 0.15
    env_scale: float = 0.75             # half-normal scale, mm
    quad_scale: float = 0.5             # half-normal scale, mm at edge


def log_prior(theta: np.ndarray, prior: Prior, use_quad: bool) -> float:
    logN, delta, logk, env = theta[:4]
    if not (prior.logN_bounds[0] <= logN <= prior.logN_bounds[1]):
        return -np.inf
    if not (prior.delta_bounds[0] <= delta <= prior.delta_bounds[1]):
        return -np.inf
    if prior.logk_mean is None:
        if not (prior.logk_bounds[0] <= logk <= prior.logk_bounds[1]):
            return -np.inf
    else:
        if not (prior.logk_bounds[0] <= logk <= prior.logk_bounds[1]):
            return -np.inf
    if env < 0:
        return -np.inf
    lp = -0.5 * (env / prior.env_scale) ** 2
    if prior.logk_mean is not None:
        lp += -0.5 * ((logk - prior.logk_mean) / prior.logk_sd) ** 2
    if use_quad:
        quad = theta[4]
        if quad < 0:
            return -np.inf
        lp += -0.5 * (quad / prior.quad_scale) ** 2
    return float(lp)


def log_likelihood(theta: np.ndarray, x: np.ndarray, h_obs: np.ndarray,
                   sigma_obs: float, lanes: int, geom: Geometry,
                   use_quad: bool) -> float:
    pred = forward(
        x, theta[0], theta[1], theta[2], lanes,
        env_mm=theta[3], quad_mm=(theta[4] if use_quad else 0.0), geom=geom,
    )
    resid = (h_obs - pred) / sigma_obs
    return float(-0.5 * np.sum(resid ** 2) - len(x) * np.log(sigma_obs))


def log_posterior(theta: np.ndarray, x: np.ndarray, h_obs: np.ndarray,
                  sigma_obs: float, lanes: int, geom: Geometry,
                  prior: Prior, use_quad: bool) -> float:
    lp = log_prior(theta, prior, use_quad)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, h_obs, sigma_obs, lanes, geom, use_quad)


# ---------------------------------------------------------------------------
# MAP and MCMC
# ---------------------------------------------------------------------------
def map_estimate(x: np.ndarray, h_obs: np.ndarray, sigma_obs: float, lanes: int,
                 geom: Geometry, prior: Prior, use_quad: bool,
                 starts: list[np.ndarray] | None = None) -> tuple[np.ndarray, float]:
    from scipy.optimize import minimize

    ndim = 5 if use_quad else 4
    if starts is None:
        base = np.array([6.0, 0.0, -6.0, 0.0, 0.0])[:ndim]
        starts = [base + np.array([d, 0, 0, 0, 0])[:ndim]
                  for d in (-1.0, 0.0, 1.0)]

    def neg(theta: np.ndarray) -> float:
        val = log_posterior(theta, x, h_obs, sigma_obs, lanes, geom, prior, use_quad)
        return 1e12 if not np.isfinite(val) else -val

    bounds = [prior.logN_bounds, prior.delta_bounds, prior.logk_bounds,
              (0.0, 3.0)]
    if use_quad:
        bounds.append((0.0, 2.0))

    best = None
    for s in starts:
        res = minimize(neg, np.asarray(s, float), method="L-BFGS-B", bounds=bounds,
                       options={"maxiter": 4000, "ftol": 1e-12})
        if best is None or res.fun < best.fun:
            best = res
    assert best is not None
    return np.asarray(best.x), float(-best.fun)


def run_mcmc(logpost_fn: Callable[[np.ndarray], float], theta0: np.ndarray,
             proposal_sd: np.ndarray, nsteps: int = 12000, burn: int = 3000,
             thin: int = 5, seed: int = 1) -> np.ndarray:
    rng = np.random.default_rng(seed)
    ndim = len(theta0)
    chain = np.empty((nsteps, ndim))
    current = np.asarray(theta0, float).copy()
    current_lp = logpost_fn(current)
    accepted = 0

    for i in range(nsteps):
        proposal = current + rng.normal(0.0, proposal_sd, size=ndim)
        proposal_lp = logpost_fn(proposal)
        if np.log(rng.random()) < proposal_lp - current_lp:
            current, current_lp = proposal, proposal_lp
            accepted += 1
        chain[i] = current

    # simple adaptation during the first half of burn-in
    if burn < nsteps:
        acc_rate = accepted / nsteps
    else:
        acc_rate = accepted / nsteps
    samples = chain[burn::thin]
    return samples


def run_chains(logpost_fn: Callable[[np.ndarray], float], theta0: np.ndarray,
               proposal_sd: np.ndarray, nchains: int = 4, nsteps: int = 12000,
               burn: int = 3000, thin: int = 5, seed: int = 1) -> np.ndarray:
    rng = np.random.default_rng(seed)
    chains = []
    for c in range(nchains):
        init = np.asarray(theta0, float) + rng.normal(0, proposal_sd * 0.5)
        chains.append(run_mcmc(logpost_fn, init, proposal_sd, nsteps, burn, thin,
                               seed=seed + 100 + c))
    return np.stack(chains, axis=0)  # (nchains, nsamples, ndim)


def numerical_hessian(f: Callable[[np.ndarray], float], theta: np.ndarray,
                      eps: float = 1e-4) -> np.ndarray:
    """Central finite-difference Hessian of a scalar function."""
    theta = np.asarray(theta, float)
    d = len(theta)
    H = np.zeros((d, d))
    f0 = f(theta)
    for i in range(d):
        ei = np.zeros(d)
        ei[i] = eps
        for j in range(i, d):
            ej = np.zeros(d)
            ej[j] = eps
            fpp = f(theta + ei + ej)
            fpm = f(theta + ei - ej)
            fmp = f(theta - ei + ej)
            fmm = f(theta - ei - ej)
            H[i, j] = H[j, i] = (fpp - fpm - fmp + fmm) / (4.0 * eps * eps)
    # regularise tiny numerical noise on the diagonal
    H[np.diag_indices(d)] += 1e-8 * max(1.0, abs(f0))
    return H


def laplace_covariance(logpost_fn: Callable[[np.ndarray], float],
                       theta_map: np.ndarray, scale: float = 2.38,
                       eigen_floor: float = 1e-6) -> np.ndarray:
    """Proposal covariance from the inverse Hessian at the MAP.

    The Hessian of the negative log posterior is the observed information.
    Eigenvalues are floored so a flat (non-identified) direction still gets a
    finite proposal width.
    """
    H = numerical_hessian(lambda t: -logpost_fn(t), theta_map)
    H = 0.5 * (H + H.T)
    vals, vecs = np.linalg.eigh(H)
    vals = np.clip(vals, eigen_floor, None)
    cov = vecs @ np.diag(1.0 / vals) @ vecs.T
    cov = 0.5 * (cov + cov.T)
    return cov * (scale ** 2) / len(theta_map)


def run_mcmc_cov(logpost_fn: Callable[[np.ndarray], float], theta0: np.ndarray,
                 cov: np.ndarray, nsteps: int = 12000, burn: int = 3000,
                 thin: int = 5, seed: int = 1) -> np.ndarray:
    """Random-walk Metropolis with a full-covariance Gaussian proposal."""
    rng = np.random.default_rng(seed)
    ndim = len(theta0)
    chain = np.empty((nsteps, ndim))
    current = np.asarray(theta0, float).copy()
    current_lp = logpost_fn(current)
    L = np.linalg.cholesky(cov + 1e-12 * np.eye(ndim))
    for i in range(nsteps):
        proposal = current + L @ rng.normal(size=ndim)
        proposal_lp = logpost_fn(proposal)
        if np.log(rng.random()) < proposal_lp - current_lp:
            current, current_lp = proposal, proposal_lp
        chain[i] = current
    return chain[burn::thin]


def run_chains_cov(logpost_fn: Callable[[np.ndarray], float], theta0: np.ndarray,
                   cov: np.ndarray, nchains: int = 4, nsteps: int = 12000,
                   burn: int = 3000, thin: int = 5, seed: int = 1) -> np.ndarray:
    rng = np.random.default_rng(seed)
    chains = []
    for c in range(nchains):
        init = np.asarray(theta0, float) + rng.multivariate_normal(
            np.zeros(len(theta0)), 0.25 * cov)
        chains.append(run_mcmc_cov(logpost_fn, init, cov, nsteps, burn, thin,
                                   seed=seed + 100 + c))
    return np.stack(chains, axis=0)


def rhat(chains: np.ndarray) -> np.ndarray:
    """Gelman-Rubin R-hat for a (nchains, nsamples, ndim) array."""
    m, n, d = chains.shape
    means = chains.mean(axis=1)
    variances = chains.var(axis=1, ddof=1)
    W = variances.mean(axis=0)
    B = n * means.var(axis=0, ddof=1)
    var_hat = (n - 1) / n * W + B / n
    return np.sqrt(var_hat / W)


# ---------------------------------------------------------------------------
# Model selection over lane counts
# ---------------------------------------------------------------------------
def bic_for_lanes(x: np.ndarray, h_obs: np.ndarray, sigma_obs: float,
                  geom: Geometry, prior: Prior, use_quad: bool = False,
                  max_lanes: int = 4) -> list[dict]:
    rows = []
    n = len(x)
    k = 5 if use_quad else 4
    for lanes in range(1, max_lanes + 1):
        theta, lp = map_estimate(x, h_obs, sigma_obs, lanes, geom, prior, use_quad)
        bic = -2.0 * lp + k * np.log(n)
        rows.append({"lanes": lanes, "log_post": lp, "BIC": bic, "theta": theta.tolist()})
    return rows


def simulate_profile(truth: Truth, geom: Geometry, sigma_obs: float,
                     seed: int = 0, reference_bias: bool = False) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    x = geom.x
    h = forward(x, truth.log10_N, truth.delta, truth.log10_kappa, truth.lanes,
                env_mm=truth.env_mm,
                quad_mm=(truth.quad_mm if reference_bias else 0.0), geom=geom)
    noisy = h + rng.normal(0.0, sigma_obs, size=h.shape)
    return x, noisy


# ---------------------------------------------------------------------------
# Summary helpers
# ---------------------------------------------------------------------------
def summarize_samples(samples: np.ndarray, truth: Truth, use_quad: bool = False) -> dict:
    """Summaries for parameter vector [logN, delta, logk, env, (quad)]."""
    logN = samples[:, 0]
    delta = samples[:, 1]
    logk = samples[:, 2]
    env = samples[:, 3]
    N = 10.0 ** logN
    out = {
        "logN_mean": float(logN.mean()),
        "logN_sd": float(logN.std(ddof=1)),
        "N_mean": float(N.mean()),
        "N_median": float(np.median(N)),
        "N_rel_error_median": float(np.median(N) / truth.N - 1.0),
        "N_ci95": [float(np.percentile(N, 2.5)), float(np.percentile(N, 97.5))],
        "N_covered": bool(np.percentile(N, 2.5) <= truth.N <= np.percentile(N, 97.5)),
        "delta_mean": float(delta.mean()),
        "delta_sd": float(delta.std(ddof=1)),
        "delta_ci95": [float(np.percentile(delta, 2.5)), float(np.percentile(delta, 97.5))],
        "delta_covered": bool(np.percentile(delta, 2.5) <= truth.delta <= np.percentile(delta, 97.5)),
        "logk_mean": float(logk.mean()),
        "logk_sd": float(logk.std(ddof=1)),
        "env_mean": float(env.mean()),
        "corr_logN_logk": float(np.corrcoef(logN, logk)[0, 1]),
    }
    if use_quad:
        out["quad_mean"] = float(samples[:, 4].mean())
    return out
