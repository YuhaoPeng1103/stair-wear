# Table 2 forensic note
Rows: 16. Mean model = 45.12 um, mean experimental = 45.50 um.
- Mean bias (model - experiment): -0.38 um
- RMSE: 7.69 um; MAPE: 14.2%
- Pearson r = 0.949 (p = 2.01e-08); OLS R^2 = 0.901
- Reported Mann-Whitney: U = 130.0, p = 0.955. With n1=n2=16, U_max = 256 and the null expectation is 128.

## Within-temperature load trend (Archard expects positive)
| T (C) | rho(load, model) | rho(load, exp) | slope model (um/100N) | slope exp (um/100N) |
|---:|---:|---:|---:|---:|
| 25 | -0.40 | -0.20 | -1.86 | -2.04 |
| 75 | -0.80 | -0.40 | -2.49 | -5.45 |
| 125 | -0.80 | -1.00 | -3.41 | -6.96 |
| 175 | -1.00 | -0.80 | -4.23 | -5.06 |

## Verdict
- No primary source is identifiable from the manuscript.
- Wear decreases with load in most temperature blocks, opposite to Archard's law.
- The two columns are nearly collinear, consistent with the 'model' column being fitted to (or derived from) the 'experimental' column.
- The test is non-significant because the samples are indistinguishable, not because the model is accurate.
- **Recommendation: remove Table 2.** It cannot serve as staircase validation.
