"""MCI/AD exposure windows using pinned HeatIndex daily products."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent / 'vendor'))
import numpy as np
import pandas as pd
import netCDF4
import upstream_hospital as source

SOURCE_COMMIT = '74114851b81a77db4fda6c4e6daca370cbaf001d'


def summer_dates(start, end):
    """Calendar dates in the half-open interval [start, end), restricted to MJJAS."""
    if end <= start:
        return pd.DatetimeIndex([])
    dates = pd.date_range(start, end, inclusive='left')
    return dates[(dates.month >= 5) & (dates.month <= 9)]


def summarize(start, end, daily):
    expected = summer_dates(start, end)
    subset = daily.reindex(expected)
    n = len(expected)
    nh = int(subset.hsci.notna().sum())
    nw = int(subset.hw.notna().sum())
    return dict(start=start.date().isoformat(), end_exclusive=end.date().isoformat(),
                expected_summer_days=n, hsci_observed_days=nh, hw_observed_days=nw,
                hsci_complete=int(nh == n), hw_complete=int(nw == n),
                HSCIH_accumulated=float(subset.hsci.sum()) if nh == n else np.nan,
                heatwave_days=int(subset.hw.sum()) if nw == n else np.nan)


def normalize_zip(value):
    s = str(value).strip()
    if s.endswith('.0'):
        s = s[:-2]
    return s.zfill(5) if s.isdigit() and len(s) <= 5 else None


def zip_daily(task):
    zc, mask, mag, dates, hsci, domain = task
    hw = np.full(len(dates), np.nan)
    try:
        if mask is not None:
            index = source.make_date_index(dates)
            for i, day in enumerate(dates):
                if not domain[i]:
                    hw[i] = 0
                else:
                    v = source.read_zip_avg(mag, 'HI_EXCDMAG', index, mask, day, False)
                    # Upstream MAG is sparse: non-exceedances are NaN.
                    hw[i] = int(np.isfinite(v) and v > 0)
        return zc, pd.DataFrame({'hsci': hsci, 'hw': hw}, index=dates)
    finally:
        for ds in source._NC_DATASET_CACHE.values():
            ds.close()
        source._NC_DATASET_CACHE.clear()
        source._ZIP_AVG_CACHE.clear()


def load_products(mag, hw, pct):
    ctx = source.load_hi_context(pct, mag, hw)
    dates = ctx['dates_daily']
    hd = pd.DatetimeIndex(ctx['hw_dates'])
    if dates.has_duplicates or hd.has_duplicates or not dates.is_monotonic_increasing:
        raise ValueError('Climate dates must be unique and daily MAG dates sorted.')
    if not hd.isin(dates).all():
        raise ValueError('HW dates are absent from the daily magnitude time axis.')
    for path in (mag, hw):
        with netCDF4.Dataset(path) as ds:
            if hasattr(ds, 'percentile') and int(ds.percentile) != pct:
                raise ValueError(f'Percentile mismatch: {path}')
    with netCDF4.Dataset(hw) as ds:
        if int(getattr(ds, 'min_duration_days', -1)) != 3 or int(getattr(ds, 'grace_days', -1)) != 1:
            raise ValueError('HSCI-H product must declare min_duration_days=3 and grace_days=1.')
        if getattr(ds, 'source_exceedance_variable', '') != 'HI_EXCDMAG':
            raise ValueError('Expected an HI-based HSCI-H product.')
    values = np.array([ctx['hsci_by_date'].get(d, 0.0) for d in dates], dtype=float)
    return dates, values, dates.isin(hd)


def atomic_csv(frame, path):
    tmp = path.with_suffix('.csv.tmp')
    frame.to_csv(tmp, index=False)
    os.replace(tmp, path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir', type=Path, default=Path('/blue/nessie/mostafarezaali/400M_PRISM'))
    p.add_argument('--patient-csv', type=Path, required=True)
    p.add_argument('--id-col', default='ID')
    p.add_argument('--zip-col', default='zip5')
    p.add_argument('--mci-col', default='MCI_DATE')
    p.add_argument('--ad-col', default='AD_DATE')
    p.add_argument('--date-format', default='%Y-%m-%d')
    p.add_argument('--pcts', default='90,95')
    p.add_argument('--workers', type=int, default=64)
    p.add_argument('--zip-csv', default='USZipsWithLatLon_20231227.csv')
    p.add_argument('--mag-template', default='HI_EXCDMAG_daily_1981_2025_{pct}.nc')
    p.add_argument('--hw-template', default='HI_EXCD_MJJAS_HWdays_{pct}.nc')
    p.add_argument('--output-dir', type=Path)
    args = p.parse_args()
    pct_list = [int(x) for x in args.pcts.split(',')]
    if not pct_list or len(set(pct_list)) != len(pct_list) or any(x not in (90,95) for x in pct_list):
        p.error('--pcts must select distinct values from 90,95')
    patient_path = args.data_dir / args.patient_csv
    patients = pd.read_csv(patient_path, dtype=str, keep_default_na=False)
    required = [args.id_col,args.zip_col,args.mci_col,args.ad_col]
    if any(c not in patients for c in required):
        p.error(f'Required input columns: {required}. Available: {list(patients.columns)}')
    if len(patients) == 0:
        p.error('Patient CSV is empty.')
    out = args.output_dir or args.data_dir / 'HSCI_Alzheimer_outputs'
    out.mkdir(parents=True, exist_ok=True)
    mci = pd.to_datetime(patients[args.mci_col], format=args.date_format, errors='coerce').dt.normalize()
    ad = pd.to_datetime(patients[args.ad_col], format=args.date_format, errors='coerce').dt.normalize()
    # Avoid pandas coercing missing normalized ZIPs from None to float NaN.
    zips = [normalize_zip(value) for value in patients[args.zip_col]]
    unique_zips = sorted(set(z for z in zips if z is not None))
    for pct in pct_list:
        mag = args.data_dir / args.mag_template.format(pct=pct)
        hw = args.data_dir / args.hw_template.format(pct=pct)
        dates, hsci, domain = load_products(mag, hw, pct)
        lat, lon = source.read_lat_lon(mag)
        masks = source.build_zip_masks(argparse.Namespace(mask_cache='', zip_buffer_cells=1,
                                       zip_csv=str(args.data_dir / args.zip_csv)), unique_zips, lat, lon)
        groups = sorted(set(zips), key=lambda z: '' if z is None else z)
        tasks = [(z, masks.get('z'+z) if z else None, str(mag), dates, hsci, domain) for z in groups]
        rows, summers = [], []
        workers = min(max(1,args.workers),len(tasks))
        pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
        try:
            results = pool.map(zip_daily, tasks) if pool else map(zip_daily,tasks)
            for number, (zc,daily) in enumerate(results,1):
                indices = [i for i,z in enumerate(zips) if z == zc]
                for i in indices:
                    identity = dict(source_row=i+2, ID=patients.iloc[i][args.id_col], zip5=zc,
                                    MCI_DATE=patients.iloc[i][args.mci_col], AD_DATE=patients.iloc[i][args.ad_col],
                                    percentile=pct)
                    row = dict(identity, mci_valid=int(pd.notna(mci[i])), ad_valid=int(pd.notna(ad[i])),
                               zip_linked=int(zc is not None and 'z'+zc in masks))
                    status = 'invalid_or_missing_mci' if pd.isna(mci[i]) else ('invalid_or_missing_ad' if pd.isna(ad[i]) else ('ad_before_mci' if ad[i]<mci[i] else 'ok'))
                    row['interval_status'] = status
                    for years in (1,3,5):
                        if pd.notna(mci[i]):
                            stats = summarize(mci[i]-pd.DateOffset(years=years),mci[i],daily)
                            row.update({f'{key}_{years}y_prior_mci':v for key,v in stats.items()})
                    if status == 'ok':
                        row.update({f'{key}_mci_to_ad':v for key,v in summarize(mci[i],ad[i],daily).items()})
                        for year in range(mci[i].year,ad[i].year+1):
                            start = max(mci[i],pd.Timestamp(year,5,1))
                            end = min(ad[i],pd.Timestamp(year,10,1))
                            if start < end:
                                summers.append(dict(identity, summer_year=year,
                                    partial_summer=int(start!=pd.Timestamp(year,5,1) or end!=pd.Timestamp(year,10,1)),
                                    **summarize(start,end,daily)))
                    rows.append(row)
                print(f'P{pct}: ZIP groups {number}/{len(tasks)} complete',flush=True)
        finally:
            if pool:
                pool.shutdown()
        # Fix column schema even when every patient has an invalid date.
        stats_keys = list(summarize(pd.Timestamp('2000-01-01'),pd.Timestamp('2000-01-01'),
                                   pd.DataFrame(columns=['hsci','hw'],index=pd.DatetimeIndex([]))))
        identity_keys = ['source_row','ID','zip5','MCI_DATE','AD_DATE','percentile']
        columns = identity_keys+['mci_valid','ad_valid','zip_linked','interval_status']
        for suffix in ['1y_prior_mci','3y_prior_mci','5y_prior_mci','mci_to_ad']:
            columns += [f'{k}_{suffix}' for k in stats_keys]
        atomic_csv(pd.DataFrame(rows).reindex(columns=columns).sort_values('source_row'),out/f'MCI_AD_exposures_{pct}.csv')
        atomic_csv(pd.DataFrame(summers).reindex(columns=identity_keys+['summer_year','partial_summer']+stats_keys),
                   out/f'MCI_AD_summers_{pct}.csv')
    manifest = dict(source_repo='Mostafa-Rezaali/HeatIndex', source_commit=SOURCE_COMMIT,
                    arguments={k:str(v) for k,v in vars(args).items()}, patient_rows=len(patients),
                    boundaries='[MCI minus calendar years, MCI); [MCI, AD); MJJAS only')
    (out/'run_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(f'Finished. Outputs: {out}',flush=True)


if __name__ == '__main__':
    main()
