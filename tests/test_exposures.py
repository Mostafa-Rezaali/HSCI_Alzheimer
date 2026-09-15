import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import exposures as e


class Windows(unittest.TestCase):
    def test_calendar_and_missing(self):
        self.assertEqual(len(e.summer_dates(pd.Timestamp('2020-06-01'),pd.Timestamp('2020-06-01'))),0)
        end = pd.Timestamp('2020-02-29')
        self.assertEqual((end-pd.DateOffset(years=1)).date().isoformat(),'2019-02-28')
        dates = e.summer_dates(pd.Timestamp('2019-01-01'),pd.Timestamp('2020-01-01'))
        self.assertEqual(len(dates),153)
        daily = pd.DataFrame({'hsci':1.,'hw':1.},index=dates)
        stats = e.summarize(pd.Timestamp('2019-05-01'),pd.Timestamp('2019-05-04'),daily)
        self.assertEqual(stats['HSCIH_accumulated'],3)
        self.assertEqual(stats['heatwave_days'],3)
        stats = e.summarize(pd.Timestamp('2018-05-01'),pd.Timestamp('2019-05-04'),daily)
        self.assertTrue(np.isnan(stats['HSCIH_accumulated']))

    def test_end_to_end(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            dates=e.summer_dates(pd.Timestamp('2014-01-01'),pd.Timestamp('2021-01-01'))
            hot=dates.day%2 == 0
            for pct,value in [(90,2.),(95,1.)]:
                for hw in (False,True):
                    path=root/(f'HI_EXCD_MJJAS_HWdays_{pct}.nc' if hw else f'HI_EXCDMAG_daily_1981_2025_{pct}.nc')
                    selected=dates[hot] if hw else dates
                    with netCDF4.Dataset(path,'w') as ds:
                        ds.percentile=pct
                        ds.createDimension('time',len(selected))
                        ds.createVariable('time','f8',('time',))[:]=[int(d.strftime('%Y%m%d')) for d in selected]
                        if hw:
                            ds.min_duration_days=3
                            ds.grace_days=1
                            ds.source_exceedance_variable='HI_EXCDMAG'
                            ds.createVariable('HSCI','f4',('time',))[:]=value
                        else:
                            for axis,coords in [('lat',[29.,30.,31.]),('lon',[-83.,-82.,-81.])]:
                                ds.createDimension(axis,3)
                                ds.createVariable(axis,'f8',(axis,))[:]=coords
                            ds.createVariable('HI_EXCDMAG','f4',('lat','lon','time'))[:]=value
            pd.DataFrame({'postal code':['01234'],'latitude':[30.],'longitude':[-82.]}).to_csv(root/'USZipsWithLatLon_20231227.csv',index=False)
            pd.DataFrame({'ID':['001','002','003','004'], 'zip5':['01234','99999','01234','01234'],
                          'MCI_DATE':['2019-06-01']*4,'AD_DATE':['2020-07-01','2020-07-01','','2018-01-01']}).to_csv(root/'patients.csv',index=False)
            subprocess.run([sys.executable,str(Path(e.__file__)), '--data-dir',str(root),
                            '--patient-csv','patients.csv','--workers','2'],check=True,capture_output=True,text=True)
            for pct,value in [(90,2.),(95,1.)]:
                table=pd.read_csv(root/'HSCI_Alzheimer_outputs'/f'MCI_AD_exposures_{pct}.csv',dtype={'ID':str,'zip5':str})
                annual=pd.read_csv(root/'HSCI_Alzheimer_outputs'/f'MCI_AD_summers_{pct}.csv')
                self.assertEqual(table.iloc[0].ID,'001')
                self.assertEqual(table.iloc[0].zip5,'01234')
                self.assertTrue(np.isnan(table.iloc[1].heatwave_days_mci_to_ad))
                self.assertEqual(table.iloc[0].HSCIH_accumulated_mci_to_ad,table.iloc[1].HSCIH_accumulated_mci_to_ad)
                self.assertEqual(table.iloc[2].interval_status,'invalid_or_missing_ad')
                self.assertEqual(table.iloc[3].interval_status,'ad_before_mci')
                for years in (1,3,5):
                    expected=e.summer_dates(pd.Timestamp('2019-06-01')-pd.DateOffset(years=years),pd.Timestamp('2019-06-01'))
                    n=int((expected.day%2==0).sum())
                    self.assertEqual(table.iloc[0][f'HSCIH_accumulated_{years}y_prior_mci'],n*value)
                    self.assertEqual(table.iloc[0][f'heatwave_days_{years}y_prior_mci'],n)
                first=annual[annual.source_row==2]
                self.assertEqual(first.HSCIH_accumulated.sum(),table.iloc[0].HSCIH_accumulated_mci_to_ad)
                self.assertEqual(first.heatwave_days.sum(),table.iloc[0].heatwave_days_mci_to_ad)
                self.assertTrue(first.partial_summer.eq(1).all())


if __name__=='__main__':
    unittest.main()
