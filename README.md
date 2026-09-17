# HSCI-H exposure before MCI and between MCI and AD

Uses existing HeatIndex daily NetCDF outputs to supply exposure measures for
patients diagnosed with mild cognitive impairment (MCI) and Alzheimer's disease
(AD). This repository computes exposures; it does not estimate disease risk.

## Methodology and provenance

The unmodified files under `vendor/` come from
https://github.com/Mostafa-Rezaali/HeatIndex at commit
`74114851b81a77db4fda6c4e6daca370cbaf001d`.
`vendor/upstream_hospital.py` is the original hospital appender and
`vendor/heatindex/` contains its original helpers. They supply NetCDF date decoding,
ZIP centroid lookup, grid masks, and the daily neighborhood average.

- Reuse HI-based HSCI-H products with monthly P90 and P95 climatological thresholds.
- Require detector metadata declaring three valid heatwave days, one post-rule
  grace day, and source variable `HI_EXCDMAG`. Spatial and persistence detection
  have already been performed upstream; this workflow does not rerun the detector.
- HSCI-H is the upstream **domain-wide** daily `HSCI` value. It is assigned by
  date, not multiplied by ZIP exposure. People with identical date windows have
  identical HSCI-H totals regardless of ZIP. Summation units are degree-Celsius days.
- Local heatwave days follow the current source repository: a date must appear
  in the domain HW-days file AND have positive finite mean HI exceedance in the
  ZIP centroid's 3x3 PRISM neighborhood. There is no population or polygon weighting.
- The 1-, 3-, and 5-year windows use calendar-year offsets and exclude MCI day.
  February 29 anniversaries use February 28 in a non-leap year.
- The interval between diagnoses includes MCI day and excludes AD day. Only
  May through September (MJJAS) dates contribute. Annual rows are clipped to the
  interval and marked as partial when they do not contain all of May-September.
- Missing AD does not imply follow-up through the present: interval outputs are
  missing. AD before MCI is flagged. Equal dates yield an empty, zero-day interval.
- The supplied ZIP is assumed to apply throughout the exposure history; residential
  moves cannot be inferred. Input rows are retained, including duplicate IDs.

## Inputs

Install dependencies with `python -m pip install -r requirements.txt` in your
chosen environment. The existing HiPerGator `torch_b200` environment is the default.

Place the patient CSV in `/blue/nessie/mostafarezaali/400M_PRISM`.
Default headers: `ID,zip5,MCI_DATE,AD_DATE`; default dates: `YYYY-MM-DD`.
IDs and ZIPs are read as text. Short numeric ZIPs are padded to five digits.
Other headers and date formats can be supplied through the options below.

Required existing files in that directory:

* `USZipsWithLatLon_20231227.csv`
* `HI_EXCDMAG_daily_1981_2025_90.nc` and `_95.nc`
* `HI_EXCD_MJJAS_HWdays_90.nc` and `_95.nc`

These must be **completed** upstream products. Their paths can be changed with
`--mag-template`, `--hw-template`, and `--zip-csv`. The date axis defines coverage;
this workflow cannot detect an upstream file whose missing slices were incorrectly
encoded as valid zeros. Absent dates in the complete sparse HW-day product are
non-heatwave days, provided the daily MAG time axis covers them.

## Run on HiPerGator

Clone or pull:

```bash
cd /blue/nessie/mostafarezaali/400M_PRISM && module load git && if [ -d HSCI_Alzheimer/.git ]; then git -C HSCI_Alzheimer pull --ff-only; else git clone https://github.com/Mostafa-Rezaali/HSCI_Alzheimer.git; fi
```

Submit after replacing `patients.csv` with the actual input filename:

```bash
cd /blue/nessie/mostafarezaali/400M_PRISM && PATIENT_CSV=patients.csv sbatch HSCI_Alzheimer/submit_mci_ad.slurm
```

One CPU job, 64 CPUs, 500 GB, maximum 24 hours, P90 and P95 sequentially.
ZIP groups are processed in parallel. No download or changes to existing climate
files are performed. Re-running recomputes exposure tables; it does not resume a
partial table. Final CSVs are replaced only after a complete percentile run.

Custom columns/date format example:

```bash
PATIENT_CSV=patients.csv ID_COL=patient_id ZIP_COL=ZIP MCI_COL=mci_date AD_COL=ad_date DATE_FORMAT='%m/%d/%Y' sbatch HSCI_Alzheimer/submit_mci_ad.slurm
```

## Outputs and dictionary

Under `400M_PRISM/HSCI_Alzheimer_outputs`, for each percentile:

* `MCI_AD_exposures_90.csv` / `_95.csv`: one row per input record.
* `MCI_AD_summers_90.csv` / `_95.csv`: one row per input record and overlapping
  MJJAS summer between MCI and AD. No invented rows for missing/invalid intervals.
* `run_manifest.json`: source revision, arguments, and boundary conventions.

Common keys: `source_row` is the original CSV line number (header=1); `ID`,
`zip5`, `MCI_DATE`, `AD_DATE` preserve patient linkage; `percentile` is 90 or 95.
`mci_valid`, `ad_valid`, `zip_linked` are 0/1 flags. `interval_status` is `ok`,
`invalid_or_missing_mci`, `invalid_or_missing_ad`, or `ad_before_mci`.

The following measures use suffixes `_1y_prior_mci`, `_3y_prior_mci`,
`_5y_prior_mci`, and `_mci_to_ad` in the patient table. The summer table uses
the same names without suffixes, plus `summer_year` and `partial_summer` (0/1).

| Measure | Definition |
|---|---|
| `HSCIH_accumulated` | Sum of daily domain HSCI-H over the interval's MJJAS dates |
| `heatwave_days` | Number of domain heatwave days with positive local HI exceedance |
| `start`, `end_exclusive` | Calendar interval boundaries |
| `expected_summer_days` | Number of MJJAS dates requested |
| `hsci_observed_days` | Requested dates with available finite daily HSCI-H, including known non-HW zeros |
| `hw_observed_days` | Requested dates with an evaluable local heatwave indicator |
| `hsci_complete`, `hw_complete` | 1 if observed count equals expected count, otherwise 0 |

Totals are missing when the corresponding coverage is incomplete, not sums of
partial histories. Missing ZIP linkage affects local counts but not domain HSCI-H.
An all-winter or zero-length valid interval has zero expected summer days and
zero accumulated exposure. For complete intervals, annual totals reconcile with
the MCI-to-AD total (within floating point summation tolerance).

Run tests with `python -m unittest discover -s tests -v`.
# Repairing legacy empty heatwave counts

`repair_hw_columns.py` repairs CSVs produced before the sparse-MAG fix.
The original MAG writer stores only positive exceedances; other values are
NaN. The old exposure code incorrectly treated those non-exceedances as
missing observations, then suppressed window totals. The repair subtracts
the number of non-domain-heatwave dates from saved `hw_observed_days` to
recover the positive local heatwave-day count. This reproduces the upstream
finite-and-positive rule without rereading spatial grids. It does not turn
unlinked ZIPs or dates outside the climate time axis into known exposures.
The sparse product cannot distinguish genuinely unavailable local source
data from non-exceedance; this limitation also applies to upstream counts.

Submit `submit_repair_hw.slurm` from HiPerGator after pulling. Defaults read
`400M_PRISM/HSCI_Alzheimer_outputs` and write separate files to
`400M_PRISM/HSCI_Alzheimer_outputs_HW_fixed`. Override `INPUT_DIR` if needed.
Only `heatwave_days*`, `hw_observed_days*`, and `hw_complete*` change; every
other field is checked for equality as strings. Originals are preserved.
Use the same climate products used in the original run. This lightweight
repair requests one CPU, 8 GB, and one hour; it reads time axes and HSCI only.
