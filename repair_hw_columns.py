"""Repair legacy HW counts using saved observation counts, without grid reads.

Only for outputs from the original sparse-MAG observation-count bug. The
upstream MAG writer retains strictly positive exceedances and writes NaN
otherwise. Thus legacy observations = positive local HW days + non-domain
HW dates. Subtracting the latter recovers the exact upstream heat-day count.
"""
import argparse
from pathlib import Path

import pandas as pd

from exposures import atomic_csv, load_products, summer_dates


def repair_frame(frame, dates, domain, linked, annual=False):
    result = frame.copy()
    suffixes = [''] if annual else ['_1y_prior_mci', '_3y_prior_mci',
                                   '_5y_prior_mci', '_mci_to_ad']
    non_hw = dates[~domain]
    for i, row in frame.iterrows():
        for suffix in suffixes:
            start, end = row['start'+suffix], row['end_exclusive'+suffix]
            if not start or not end:
                continue
            expected = summer_dates(pd.Timestamp(start), pd.Timestamp(end))
            n = len(expected)
            is_linked = linked.get(str(row['source_row']), False)
            # Unlinked ZIPs are genuinely unknown, except empty intervals.
            if not is_linked:
                continue
            if row['hw_complete'+suffix] == '1':
                continue
            legacy = int(float(row['hw_observed_days'+suffix]))
            count = legacy - int(expected.isin(non_hw).sum())
            eligible = int(expected.isin(dates[domain]).sum())
            if not 0 <= count <= eligible:
                raise ValueError(f'Not a legacy sparse-MAG result: row {i}, {suffix}')
            covered = int(expected.isin(dates).sum())
            result.at[i, 'hw_observed_days'+suffix] = str(covered)
            result.at[i, 'hw_complete'+suffix] = str(int(covered == n))
            result.at[i, 'heatwave_days'+suffix] = str(count) if covered == n else ''
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir', type=Path, required=True)
    p.add_argument('--input-dir', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--pcts', default='90,95')
    args = p.parse_args()
    if args.input_dir.resolve() == args.output_dir.resolve():
        p.error('Use a separate output directory to preserve originals.')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for pct in [int(x) for x in args.pcts.split(',')]:
        dates, _, domain = load_products(
            args.data_dir / f'HI_EXCDMAG_daily_1981_2025_{pct}.nc',
            args.data_dir / f'HI_EXCD_MJJAS_HWdays_{pct}.nc', pct)
        name = f'MCI_AD_exposures_{pct}.csv'
        patients = pd.read_csv(args.input_dir/name, dtype=str, keep_default_na=False)
        if patients.source_row.duplicated().any():
            raise ValueError('Duplicate source_row in exposure file')
        linked = dict(zip(patients.source_row, patients.zip_linked.eq('1')))
        for kind, annual in [('exposures', False), ('summers', True)]:
            name = f'MCI_AD_{kind}_{pct}.csv'
            frame = pd.read_csv(args.input_dir/name, dtype=str, keep_default_na=False)
            if not frame.source_row.isin(linked).all():
                raise ValueError('Summer rows do not match exposure rows')
            fixed = repair_frame(frame, dates, domain, linked, annual)
            untouched = [c for c in frame if not c.startswith(('hw_', 'heatwave_days'))]
            pd.testing.assert_frame_equal(frame[untouched], fixed[untouched])
            atomic_csv(fixed, args.output_dir/name)
            print(f'Repaired {name}; all non-HW fields unchanged.', flush=True)


if __name__ == '__main__':
    main()
