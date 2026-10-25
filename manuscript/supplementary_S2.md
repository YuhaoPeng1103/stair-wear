# Supplementary Material S2 — Reproducible code and synthetic data

## S2.1 Environment

- Python 3.10
- numpy 2.2, scipy 1.15, pandas 2.3, matplotlib 3.10
- No proprietary solver is required. The finite-element claim of the original competition paper is not used anywhere in this manuscript.

Run the scripts from the project root. Set `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`; the profile-likelihood fits are small and threading only adds overhead.

## S2.2 Repository layout

```
src/model.py          forward model, priors, MAP, MCMC utilities
src/inference.py      profile likelihood, 2-D slices
src/calibration.py    known-traffic calibration and propagation
src/io_field.py       field schema, reference surface, profile extraction
experiments/run_simulation.py    main simulation study
experiments/make_figures.py      Figures 1-9
experiments/table2_forensics.py  audit of the original Table 2
experiments/field_pipeline.py    calibration + application on field data
experiments/make_field_demo.py   synthetic known-truth integration test
```

## S2.3 Reproduce the simulation

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 experiments/run_simulation.py
```

This writes `results/identifiability.csv`, `noise.csv`, `lanes.csv`,
`direction.csv`, `calibration.csv`, `reference.csv`, and
`results/summary_stats.json`. Runtime is about 4–6 minutes on a laptop.
All random seeds are fixed in the script.

Key parameters: \(N=10^6\), \(\delta=0.3\), \(\kappa=10^{-6}\) mm per footfall,
\(L=2\), \(W=300\) mm, lane spacing 65 mm, direction offset 36 mm,
\(\sigma=24\) mm, 2 mm sampling, default \(\sigma_{\mathrm{obs}}=0.2\) mm.

## S2.4 Reproduce the figures

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 experiments/make_figures.py
```

Figures are written to `manuscript/figures/` as PDF and 300 dpi PNG.

## S2.5 Field pipeline

```bash
python3 experiments/field_pipeline.py \
  --metadata data/field/metadata.csv \
  --profiles data/field/profiles.csv \
  --out results/field
```

The pipeline (i) averages valid profiles per staircase onto a common grid,
(ii) fits \(\kappa\) on calibration staircases with known \(N\),
(iii) selects the lane count by BIC, and (iv) profiles \(N\) for application
staircases with the calibrated prior. The output includes confidence intervals
and the reference-surface fit.

## S2.6 Known-truth integration test

`experiments/make_field_demo.py` generates a synthetic dataset with known
truth and the same schema as field data. This is a test of the pipeline, **not**
field data. Three calibration staircases (\(N=2\times10^5,5\times10^5,10^6\))
and three application staircases were generated with 0.1–0.2 mm noise. The
pipeline recovered:

| stair | true \(N\) | estimated \(N\) | 95% interval | true lanes | selected lanes | true \(\delta\) | estimated \(\delta\) |
|---|---:|---:|---|---:|---:|---:|---:|
| A1 | \(5.0\times10^5\) | \(5.01\times10^5\) | \(4.22\times10^5\)–\(5.96\times10^5\) | 1 | 1 | -0.30 | -0.27 |
| A2 | \(1.0\times10^6\) | \(1.00\times10^6\) | \(8.41\times10^5\)–\(1.19\times10^6\) | 2 | 2 | +0.30 | +0.30 |
| A3 | \(2.0\times10^6\) | \(2.00\times10^6\) | \(1.68\times10^6\)–\(2.37\times10^6\) | 2 | 2 | 0.00 | 0.00 |

Calibrated \(\log_{10}\kappa=-6.007\) against a true value of \(-6.000\)
(site-to-site SD used: 0.049 dex). The test confirms that the schema, the
profile averaging, the calibration, and the profile-likelihood inversion are
consistent. It does not validate the model against real stone.

## S2.7 Audit of the original Table 2

```bash
python3 experiments/table2_forensics.py
```

The audit reproduces the reported Mann–Whitney result (\(U=130.0\),
\(p=0.955\)) and shows that, within every temperature block, wear depth
decreases as load increases. Under Archard's law the opposite is expected.
The table has no provenance in the original manuscript and is not used in this
paper. See `docs/table2_forensics.md`.

## S2.8 Public-model case study

The Roman marble staircase model is distributed by Sketchfab under CC BY 4.0
(author: chabychab). Download the model while logged in, unzip it under
`data/public_model/extracted/roman/`, and run:

```bash
python3 experiments/extract_public_profiles.py --model roman
python3 experiments/public_model_case.py
```

The first script segments five treads and writes individual profiles to
`data/public_model/profiles/`. The second computes the apparent central wear
under three reference models (plane, quadratic in width, edge median) and
writes `data/public_model/case/roman_reference_sensitivity.csv` and
`roman_case.png`. The central estimate ranges from -5.3 to +8.8 mm and changes
sign; the implied N therefore changes by at least a factor of 6.3 on the one
tread where all three estimates are positive.

The Paestum and Umm Qais models are provided for comparison; the Paestum stair
is a curved stair with narrow treads, and the Umm Qais model includes the
surrounding structure, so both require more segmentation than the Roman model.

Attribution (required by CC BY 4.0):

> This work is based on "Roman marble stairs" (https://sketchfab.com/3d-models/roman-marble-stairs-3a26ffc0044a4ca4a41d244f2f528a5a) by chabychab (https://sketchfab.com/chabychab), licensed under CC BY 4.0. It is also based on "Paestum old stairs" (https://sketchfab.com/3d-models/paestum-old-stairs-0c2cf3f7048343368d14e2c798721abd) by chabychab, and "Stair over an arch - Umm Qais Museum" (https://sketchfab.com/3d-models/stair-over-an-arch-umm-qais-museum-d088a03863ad40cabbb17326e4efc0d4) by SJoudeh, both licensed under CC BY 4.0.

## S2.9 Literature anchor

```bash
python3 experiments/literature_anchor.py
```

The script applies the Borobudur projection (18 mm per 50 million visitors,
log10 kappa = -6.444) as a prior with 0.30 and 0.50 dex standard deviations,
and writes `results/literature_anchor.csv` and
`results/literature_anchor_summary.csv`. The anchor narrows the N interval but
biases N by +172% because it is 0.44 dex below the simulation's true kappa. The
input table is `data/literature/borobudur_anchor.csv`.

## S2.10 Licence and citation

The code and synthetic data are released under the MIT licence (to be confirmed
by the corresponding author). If you use the code, cite the manuscript and the
repository DOI (to be assigned at acceptance).
