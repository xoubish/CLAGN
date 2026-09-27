"""Validate the October packet against serialized telescope coordinates.

Local observing products are intentionally untracked; skip in a fresh checkout.
"""
from pathlib import Path
import csv
import json
import unittest
import warnings

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body
from astropy.time import Time
from astropy.utils import iers

ROOT=Path(__file__).resolve().parents[1]
PACKET=ROOT/'observing/oct26_27_2026'
iers.conf.auto_download=False
iers.conf.auto_max_age=None


def packet_csv(name):
    return pd.read_csv(PACKET/name).rename(columns=str.lower)


@unittest.skipUnless((PACKET/'validation.json').exists(),'Local October packet not built')
class OctoberPacketTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.queue=packet_csv('all_80_targets.csv').fillna('')
        cls.validation=json.loads((PACKET/'validation.json').read_text())

    def test_cardinality_settings_and_jwst_coverage(self):
        q=self.queue
        self.assertEqual(q.groupby('night').size().to_dict(),{'oct26':40,'oct27':40})
        self.assertEqual(q.name.nunique(),80)
        jwst=pd.read_csv(ROOT/'jwst_proposal/inputs/jwst_sample_cycle6.csv')
        audit=packet_csv('jwst_coverage.csv').fillna('')
        self.assertEqual(set(audit.name),set(jwst.internal_id))
        for r in audit.itertuples():
            self.assertTrue(r.name in set(q.name) or r.observed_sep23)
        observed=set(pd.read_csv(ROOT/'sep23_data/reduction_20260924/products_p330e/target_summary.csv').NAME)
        self.assertTrue((set(q.name)&observed)<=set(jwst.internal_id))
        windows=packet_csv('visibility_windows.csv')
        for name in set(audit.name)-set(q.name):
            self.assertFalse((windows.name==name).any())
        for col,v in [('seconds_each',300),('exposures',2),('slit_arcsec',1.5),('binspat',2),('binspect',3),('visit_minutes',12)]:
            self.assertTrue(q[col].eq(v).all(),col)

    def test_no_overlap_or_twilight_overrun_and_time_budget(self):
        standards=packet_csv('standards.csv')
        gaps=packet_csv('free_time.csv')
        for night,m in self.validation['nights'].items():
            q=self.queue[self.queue.night==night]
            self.assertEqual(q.order.tolist(),list(range(1,41)))
            visits=pd.concat([q[['start_pdt','end_pdt']],standards[standards.night==night][['start_pdt','end_pdt']]])
            visits=visits.sort_values('start_pdt')
            previous=pd.Timestamp(m['dusk_pdt'])
            for r in visits.itertuples():
                a,b=pd.Timestamp(r.start_pdt),pd.Timestamp(r.end_pdt)
                self.assertGreaterEqual(a,previous)
                self.assertLessEqual(b,pd.Timestamp(m['dawn_pdt']))
                previous=b
            self.assertAlmostEqual(gaps[gaps.night==night].minutes.sum()+500,m['minutes'],places=5)

    def test_upload_coordinates_and_independent_moon_frame(self):
        location=EarthLocation.from_geodetic(-116.865*u.deg,33.3563*u.deg,1712*u.m)
        for night in ['oct26','oct27']:
            with (PACKET/f'{night}_ngps.csv').open() as f:
                rows=list(csv.DictReader(f))
            q=self.queue[self.queue.night==night]
            self.assertEqual([r['NAME'] for r in rows],q.name.tolist())
            self.assertTrue(all(r['EXPTIME']=='SET 300' and r['NEXP']=='2' and r['SLITWIDTH']=='SET 1.5' and r['BINSPAT']=='2' and r['BINSPECT']=='3' for r in rows))
            for upload,r in zip(rows,q.itertuples()):
                coord=SkyCoord(upload['RA'],upload['DECL'],unit=(u.hourangle,u.deg))
                original=SkyCoord(r.ra*u.deg,r.dec*u.deg)
                self.assertLess(coord.separation(original).arcsec,.02)
                times=Time(pd.date_range(r.start_pdt,r.end_pdt,freq='30s').to_pydatetime())
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore',message='Tried to get polar motions')
                    moon=get_body('moon',times,location)
                    # Compare in the topocentric GCRS frame, independently of
                    # the planner's shared AltAz target/Moon separation.
                    separation=coord.transform_to(moon.frame).separation(moon).deg
                    aa=coord.transform_to(AltAz(obstime=times,location=location,pressure=0*u.hPa))
                self.assertGreaterEqual(separation.min(),40,r.name)
                self.assertTrue(np.all(aa.alt.deg>0),r.name)
                self.assertLessEqual(aa.secz.max().value,2 if r.jwst_id else 1.8,r.name)
                self.assertAlmostEqual(separation.min(),r.moon_min,places=2)

    def test_backups_are_distinct_from_primaries_and_fit_slots(self):
        b=packet_csv('slot_backups.csv')
        self.assertFalse(set(b.name)&set(self.queue.name))
        self.assertEqual(b.groupby(['night','replaces_order']).size().min(),2)
        self.assertTrue(b.r_adopted.lt(19).all())
        self.assertTrue(b.airmass_max_actual.le(1.8).all())
        self.assertTrue(b.moon_min.ge(40).all())
        lookup=self.queue.set_index(['night','order'])
        for r in b.itertuples():
            primary=lookup.loc[(r.night,r.replaces_order)]
            self.assertEqual(r.start_pdt,primary.start_pdt)
            self.assertEqual(r.end_pdt,primary.end_pdt)


if __name__=='__main__':
    unittest.main()
