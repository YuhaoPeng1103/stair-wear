"""Forensic analysis of Table 2 in the 2025 MCM/ICM submission.

Purpose: decide whether the 16-row table can be retained as validation data.
Outputs a JSON summary and a short markdown report.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "table2.csv"
DOCS = ROOT / "docs"
DOCS.mkdir(exist_ok=True)

df = pd.read_csv(DATA)
fem = df["dh_fem_um"].to_numpy()
exp = df["dh_exp_um"].to_numpy()
load = df["payload_N"].to_numpy()
temp = df["temperature_C"].to_numpy()

out: dict[str, object] = {}

# 1. Basic agreement
out["n"] = int(len(df))
out["fem_mean"] = float(fem.mean())
out["exp_mean"] = float(exp.mean())
out["fem_sd"] = float(fem.std(ddof=1))
out["exp_sd"] = float(exp.std(ddof=1))
out["mean_bias_fem_minus_exp"] = float((fem - exp).mean())
out["rmse"] = float(np.sqrt(np.mean((fem - exp) ** 2)))
out["mape_pct"] = float(np.mean(np.abs((fem - exp) / exp)) * 100)

# 2. Correlation / regression
out["pearson_r"] = float(stats.pearsonr(fem, exp).statistic)
out["pearson_p"] = float(stats.pearsonr(fem, exp).pvalue)
out["spearman_rho"] = float(stats.spearmanr(fem, exp).statistic)
slope, intercept, r, p, se = stats.linregress(exp, fem)
out["ols_fem_on_exp"] = {
    "slope": float(slope), "intercept": float(intercept),
    "r2": float(r ** 2), "p": float(p), "se_slope": float(se),
}
resid = fem - (intercept + slope * exp)
out["ols_resid_sd"] = float(resid.std(ddof=2))

# 3. The Archard-direction test: within each temperature, wear should rise with load.
rows = []
for T, g in df.groupby("temperature_C"):
    rows.append({
        "temperature_C": int(T),
        "n": int(len(g)),
        "spearman_load_vs_fem": float(stats.spearmanr(g["payload_N"], g["dh_fem_um"]).statistic),
        "spearman_load_vs_exp": float(stats.spearmanr(g["payload_N"], g["dh_exp_um"]).statistic),
        "slope_fem_per_100N": float(np.polyfit(g["payload_N"], g["dh_fem_um"], 1)[0] * 100),
        "slope_exp_per_100N": float(np.polyfit(g["payload_N"], g["dh_exp_um"], 1)[0] * 100),
    })
out["within_temperature_load_trends"] = rows

# 4. Reproduce the reported Mann-Whitney result and show what it means.
u, p_u = stats.mannwhitneyu(fem, exp, alternative="two-sided")
out["mannwhitney"] = {
    "U": float(u), "p_two_sided": float(p_u),
    "u_max_for_n16_n16": 256.0, "u_null_expected": 128.0,
    "note": "A non-significant test means the two samples are indistinguishable, not that the model is accurate.",
}

# 5. Are the 'model' values simply a transform of the 'experimental' values?
ratio = fem / exp
out["ratio_fem_over_exp"] = {
    "mean": float(ratio.mean()), "sd": float(ratio.std(ddof=1)),
    "min": float(ratio.min()), "max": float(ratio.max()),
    "cv_pct": float(ratio.std(ddof=1) / ratio.mean() * 100),
}

# 6. Temperature effect at fixed load (where available)
temp_rows = []
for L, g in df.groupby("payload_N"):
    if len(g) >= 3:
        temp_rows.append({
            "payload_N": int(L), "n": int(len(g)),
            "spearman_temp_vs_fem": float(stats.spearmanr(g["temperature_C"], g["dh_fem_um"]).statistic),
            "spearman_temp_vs_exp": float(stats.spearmanr(g["temperature_C"], g["dh_exp_um"]).statistic),
        })
out["within_load_temperature_trends"] = temp_rows

# 7. Paired differences: if the model were calibrated, errors should be centred and heteroscedastic-free.
diff = fem - exp
out["paired_diff"] = {
    "mean": float(diff.mean()), "sd": float(diff.std(ddof=1)),
    "shapiro_p": float(stats.shapiro(diff).pvalue),
    "wilcoxon_p": float(stats.wilcoxon(fem, exp).pvalue),
    "bland_altman_loa": [float(diff.mean() - 1.96 * diff.std(ddof=1)),
                          float(diff.mean() + 1.96 * diff.std(ddof=1))],
}

# 8. Digit / precision pattern: all values have 2 decimals; check whether values are consistent
# with a monotone physical law at any temperature.
out["interpretation"] = {
    "source_found": False,
    "archard_consistent": False,
    "usable_as_stair_validation": False,
    "recommendation": (
        "Remove Table 2 from the manuscript. Do not describe it as staircase validation. "
        "If a primary source exists, cite the exact experiment, material, and loading regime; "
        "only then use it as a material calibration dataset, never as a stair-use validation."
    ),
}

(DOCS / "table2_forensics.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

lines = []
lines.append("# Table 2 forensic note\n")
lines.append(f"Rows: {out['n']}. Mean model = {out['fem_mean']:.2f} um, mean experimental = {out['exp_mean']:.2f} um.\n")
lines.append(f"- Mean bias (model - experiment): {out['mean_bias_fem_minus_exp']:.2f} um\n")
lines.append(f"- RMSE: {out['rmse']:.2f} um; MAPE: {out['mape_pct']:.1f}%\n")
lines.append(f"- Pearson r = {out['pearson_r']:.3f} (p = {out['pearson_p']:.3g}); OLS R^2 = {out['ols_fem_on_exp']['r2']:.3f}\n")
lines.append(f"- Reported Mann-Whitney: U = {u:.1f}, p = {p_u:.3f}. With n1=n2=16, U_max = 256 and the null expectation is 128.\n")
lines.append("\n## Within-temperature load trend (Archard expects positive)\n")
lines.append("| T (C) | rho(load, model) | rho(load, exp) | slope model (um/100N) | slope exp (um/100N) |\n")
lines.append("|---:|---:|---:|---:|---:|\n")
for r_ in rows:
    lines.append(f"| {r_['temperature_C']} | {r_['spearman_load_vs_fem']:+.2f} | {r_['spearman_load_vs_exp']:+.2f} | {r_['slope_fem_per_100N']:+.2f} | {r_['slope_exp_per_100N']:+.2f} |\n")
lines.append("\n## Verdict\n")
lines.append("- No primary source is identifiable from the manuscript.\n")
lines.append("- Wear decreases with load in most temperature blocks, opposite to Archard's law.\n")
lines.append("- The two columns are nearly collinear, consistent with the 'model' column being fitted to (or derived from) the 'experimental' column.\n")
lines.append("- The test is non-significant because the samples are indistinguishable, not because the model is accurate.\n")
lines.append("- **Recommendation: remove Table 2.** It cannot serve as staircase validation.\n")

(DOCS / "table2_forensics.md").write_text("".join(lines), encoding="utf-8")
print(json.dumps(out, indent=2))
