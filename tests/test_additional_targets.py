"""Check requested target geometry independently of the payload generator."""
import json
import re
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body
from astropy.time import Time
from astropy.utils import iers

ROOT=Path(__file__).resolve().parents[1]
REVIEW=ROOT/'data/additional_targets/H5MUSES-327'


@unittest.skipUnless((REVIEW/'target_payload.json').exists(), 'Local archive review unavailable')
class AdditionalTargetTests(unittest.TestCase):
    def test_visible_catalogue_entry_keeps_identity_uncertainty(self):
        page=(ROOT/'docs/index.html').read_text()
        data=json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>',page,re.S).group(1))
        target=next(t for t in data['targets'] if t['name']=='H5MUSES-327')
        self.assertIsNone(target['z'])
        self.assertEqual(target['pool_role'],'requested')
        self.assertIn('stellar',target['science']['state_classification'])
        self.assertTrue(np.isfinite([target['ux'],target['uy']]).all())
        self.assertGreater(len(target['ztf']['g']),100)
        self.assertEqual(len(target['wise']['W1']),27)
        self.assertEqual(len(target['neo']),22)
        self.assertIsNotNone(target['cut40'])
        visits=[r for r in data['october']['rows'] if r['name']==target['name']]
        self.assertEqual(len(visits),1)
        self.assertEqual(visits[0]['night'],'oct26')

    def test_first_and_latest_starts_fit_full_visit_in_independent_frame(self):
        iers.conf.auto_download=False
        iers.conf.auto_max_age=None
        location=EarthLocation.from_geodetic(-116.865*u.deg,33.3563*u.deg,1712*u.m)
        coord=SkyCoord(261.137*u.deg,59.44636111111111*u.deg)
        standards=pd.read_csv(ROOT/'observing/oct26_27_2026/standards.csv')
        windows=json.loads((REVIEW/'visibility.json').read_text())
        for night, tiers in windows.items():
            for w in tiers['1.8']:
                for key in ['first_start_pdt','last_start_pdt']:
                    start=pd.Timestamp(w[key]);end=start+pd.Timedelta(minutes=14)
                    times=Time(start.to_pydatetime())+np.r_[np.arange(0,840,7),840]*u.s
                    moon=get_body('moon',times,location)
                    separation=coord.transform_to(moon.frame).separation(moon).deg
                    aa=coord.transform_to(AltAz(obstime=times,location=location,pressure=0*u.hPa))
                    self.assertTrue((aa.alt.deg>0).all())
                    self.assertLessEqual(aa.secz.value.max(),1.8)
                    self.assertGreaterEqual(separation.min(),40)
                    for r in standards[standards.NIGHT==night].itertuples():
                        self.assertFalse(start<pd.Timestamp(r.END_PDT) and end>pd.Timestamp(r.START_PDT))


if __name__=='__main__':unittest.main()
