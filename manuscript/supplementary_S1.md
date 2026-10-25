# Supplementary Material S1 — Field protocol for calibrated stair-wear inference

This protocol accompanies the manuscript “What worn stair treads can and cannot tell: identifiability, calibration, and measurement design for wear-based inspection of historic stairs.” It is written so that a team of two or three people with a camera, a scale bar, and a door counter can execute it.

## S1.1 Design overview

Two kinds of staircases are needed.

- **Calibration staircases** have a documented installation year and a known cumulative footfall \(N_{\mathrm{cal}}\). They are used to estimate the wear coefficient \(\kappa=k/H\). They need not be historic; modern library, office, or institutional stairs are preferable because their traffic is countable and their material is homogeneous.
- **Application staircases** are the historic stairs whose use history is unknown. They are measured in the same way, and the calibrated \(\kappa\) is used to infer \(N\) with a confidence interval.

The minimum publishable dataset is three calibration staircases and three application staircases, with five treads per staircase and three longitudinal profiles per tread.

## S1.2 Calibration staircase criteria

A calibration staircase should satisfy all of the following:

1. Installation or construction year documented by a contract, plaque, property record, or dated photograph.
2. At least five years of service since installation, with no resurfacing, cladding, or replacement of the tread surface.
3. A single homogeneous material (marble, granite, limestone, terrazzo, concrete, or wood) with a known surface finish.
4. Cumulative footfall countable through a door counter, turnstile, ticket system, or a 14-day manual count.
5. Permission for non-destructive measurement.

Exclude staircases with carpets, anti-slip strips, recent coatings, or a documented change of use that makes the cumulative count unreliable.

## S1.3 Traffic counting

**Automatic counts.** Export daily counts from the access-control system. If the counter covers a building entrance rather than the staircase, estimate the staircase share \(p\) by manual counting in three periods (morning peak, midday, evening peak), 30 minutes each, and record the share.

**Manual counts.** Count for 14 consecutive days, three 30-minute periods per day (e.g. 08:00–08:30, 12:00–12:30, 17:30–18:00). Record ascending and descending separately. Extrapolate to a daily total with period weights; keep weekdays and weekends separate. Note closures, holidays, and construction.

**Cumulative footfall.**

\[
N=\bar n_{\mathrm{day}}\,D_{\mathrm{service}}\,p,
\qquad
\sigma_N/N=\sqrt{(\sigma_{\bar n}/\bar n)^2+(\sigma_p/p)^2},
\]

where \(\bar n_{\mathrm{day}}\) is the mean daily count, \(D_{\mathrm{service}}\) the number of service days, and \(p\) the staircase share. Record the relative SD in `metadata.csv`.

## S1.4 Geometric measurement

**Photogrammetry (recommended, low cost).** Take 40–80 overlapping photographs per tread, covering the tread, riser, and nosing; avoid direct specular reflection. Place at least two scale bars at right angles and four to six coded targets. Reconstruct in Agisoft Metashape, RealityCapture, or COLMAP. Export a scaled dense point cloud or mesh. Typical precision is 0.5–2 mm; it must be verified by repeat measurement, not assumed.

**Structured light or laser scanning (high precision).** Use a handheld or tripod scanner with a point spacing of 1 mm or less. Scan each tread from at least two stations with 30% overlap. Typical precision is 0.1–0.5 mm.

**Profile gauge (fallback).** Use a profile gauge or depth gauge along three longitudinal lines at 10 mm intervals. Precision is about 1 mm; use only as a supplement or for a pilot survey.

## S1.5 Sampling layout

Select five consecutive treads per staircase, avoiding the top and bottom steps to reduce boundary effects. For each tread:

- three longitudinal profiles at 50 mm from the left edge, at the centre line, and 50 mm from the right edge, from the nosing to the rear edge, at 5–10 mm intervals;
- one transverse profile at mid-depth, across the full width, at 5–10 mm intervals;
- record step height, tread depth, nosing radius, material, and surface treatment.

## S1.6 Reference surface

The reference surface is the largest practical uncertainty. Use, in order of preference:

1. the unworn nosing or the unworn outer 20% of the tread;
2. an unused or blocked neighbouring tread;
3. a plane fitted to the outer 20% of the tread and extrapolated to the centre;
4. a contemporary staircase of the same material and construction as a control.

Report the method, the fit residual, and the sensitivity of \(N\) to the reference-surface choice. The model includes a quadratic term \(q(x/(W/2))^2\); fit it rather than assuming it is zero.

## S1.7 Material characterisation

Measure Leeb hardness (or Schmidt rebound) at five points per tread and report the median and SD. If a sample of the same stone is available from a previous intervention, measure density and porosity. Do not assume a literature value of \(k\); calibrate \(\kappa\) on the calibration staircases.

## S1.8 Data schema

Two CSV files are required. The complete column definitions are in `data/field/README.md`. In summary:

- `metadata.csv`: `stair_id`, `role`, `install_year`, `service_days`, `mean_daily_footfall`, `footfall_rel_sd`, `material`, `hardness_leeb`, `lanes`, `scan_accuracy_mm`, `notes`.
- `profiles.csv`: `stair_id`, `tread_index`, `profile_id`, `x_mm`, `h_mm`, `valid`.

## S1.9 Quality control

- [ ] At least five treads per staircase and three longitudinal profiles per tread.
- [ ] Scale check against the step height: relative error < 0.5%.
- [ ] Repeat measurement of one tread three times: SD < 0.3 mm.
- [ ] Reference-surface fit residual < 0.5 mm.
- [ ] Traffic count of at least 14 days, including weekdays and weekends.
- [ ] Installation year supported by a document.
- [ ] Site permission and data-use permission on file.
- [ ] Raw point clouds and photographs archived with a DOI.

## S1.10 Ethics and permissions

Obtain written permission from the building owner or heritage authority before entering, photographing, or scanning. Do not take destructive samples. Choose inconspicuous locations for hardness measurements. Do not record personally identifying information during traffic counts; avoid identifiable faces in photographs. Agree anonymisation with the owner before publishing site data.

## S1.11 Minimum dataset and expected uncertainty

With three calibration staircases and a site-to-site SD of 0.05 dex in \(\log_{10}\kappa\), Eq. (8) of the manuscript predicts a relative SD in \(N\) of about 12%, before adding profile noise. At 0.2 mm profile noise, the simulation gives a median bias of \(-3.0\%\) and a 95% interval width of 0.1 dex. These are planning values, not guarantees; the field data should be used to estimate the actual calibration variance.


## S1.12 Literature anchor (not a substitute for calibration)

If no calibration staircase is available, a published wear rate can be used as a
prior. The Borobudur Conservation Office reports a projection of 1.8 cm of wear
per 50 million visitors, i.e. an effective kappa of about 3.6 x 10^-7 mm per
visitor. This is an average over the staircase surface, not a single-point
measurement, and the stone is andesite. A prior of log10 kappa = -6.44 with a
standard deviation of 0.30 dex narrows the 95% interval for N from six orders of
magnitude to about 1.1 dex, but biases N by +172% when the target material has a
different kappa. Report the anchor, its source, and the material mismatch; do not
call it a calibration.
