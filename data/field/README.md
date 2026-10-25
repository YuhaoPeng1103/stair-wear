# Field data directory

Put the real field data here. Do **not** copy `data/field_demo/` into this
folder: the demo is synthetic and is labelled as such.

## Files

- `metadata.csv` — one row per staircase.
- `profiles.csv` — long-format profiles, one row per point.

The full protocol is in `docs/data_collection_SOP.md`.

## metadata.csv

| column | type | required | notes |
|---|---|---|---|
| `stair_id` | string | yes | unique, e.g. `C1`, `A1` |
| `role` | string | yes | `calibration` or `application` |
| `install_year` | int | yes | documented construction/installation year |
| `service_days` | int | yes | days since installation, excluding closures |
| `mean_daily_footfall` | float | yes (calibration) | mean daily users of this staircase |
| `footfall_rel_sd` | float | yes | relative SD of the count (0.1 = 10%) |
| `material` | string | yes | marble, granite, limestone, terrazzo, wood, ... |
| `hardness_leeb` | float | recommended | median Leeb hardness |
| `lanes` | int | recommended | number of side-by-side lanes; use 1 if unknown |
| `scan_accuracy_mm` | float | yes | repeat-measurement SD |
| `notes` | string | optional | permissions, anomalies, missing data |

## profiles.csv

| column | type | required | notes |
|---|---|---|---|
| `stair_id` | string | yes | must match `metadata.csv` |
| `tread_index` | int | yes | 1–5 (or more) |
| `profile_id` | string | yes | `L1`–`L3` for longitudinal, `T1` for transverse |
| `x_mm` | float | yes | lateral coordinate, 0 at tread centre |
| `h_mm` | float | yes | wear depth, positive downward, relative to reference |
| `valid` | int | yes | 1 = use, 0 = exclude |

## Minimum dataset

- 3 calibration staircases with known cumulative traffic
- 3 application staircases
- 5 treads per staircase
- 3 longitudinal profiles per tread
- 14 days of traffic counting per calibration staircase

## Run

```bash
python3 experiments/field_pipeline.py \
  --metadata data/field/metadata.csv \
  --profiles data/field/profiles.csv \
  --out results/field
```

The pipeline writes `calibration_stairs.csv`, `application_stairs.csv`, and
`field_summary.json`.
