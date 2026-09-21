"""Check the local telescope packet against geometry, slot limits and page exports."""
from pathlib import Path
import csv
import json
import re
import unittest
import warnings

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body
from astropy.time import Time
from astropy.utils import iers
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT/'observing/sep23'
iers.conf.auto_download = False
iers.conf.auto_max_age = None


@unittest.skipUnless((DEST/'snr5_plan.json').exists(), 'Requires the local observing packet')
class SeptemberPacketTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet = json.loads((DEST/'packet.json').read_text())
        cls.plan = json.loads((DEST/'snr5_plan.json').read_text())
        cls.rows = list(csv.DictReader((DEST/'sep23_primaries_ngps.csv').read_text().splitlines()))
        cls.backup_rows = list(csv.DictReader((DEST/'sep23_backups_ngps.csv').read_text().splitlines()))

    def test_sequence_preserves_originals_and_only_promotes_former_backups(self):
        p = self.packet
        names = [v['name'] for v in p['primaries']]
        self.assertEqual(len(names), len(set(names)))
        self.assertGreater(len(names), len(self.plan['original_primaries']))
        self.assertTrue(set(self.plan['original_primaries']) <= set(names))
        self.assertEqual(set(names)-set(self.plan['original_primaries']), set(p['promoted']))
        self.assertTrue(set(p['promoted']) <= set(self.plan['former_backups']))
        self.assertEqual(p['settings']['goal_snr'], 5)
        for a,b in zip(p['sequence'],p['sequence'][1:]):
            self.assertLessEqual(pd.Timestamp(a['end_utc']), pd.Timestamp(b['start_utc']))
        science_end = pd.Timestamp(p['primaries'][-1]['end_utc'])
        reserve_start = pd.Timestamp(p['reserved']['start_pdt'],tz='America/Los_Angeles')
        reserve_end = pd.Timestamp(p['reserved']['end_pdt'],tz='America/Los_Angeles')
        self.assertGreaterEqual((reserve_end-reserve_start).total_seconds(),1200)
        for v in p['sequence']:
            self.assertTrue(pd.Timestamp(v['end_utc'])<=reserve_start or pd.Timestamp(v['start_utc'])>=reserve_end)
        late = p['primaries'][-1]
        self.assertEqual(late['name'],'P12457')
        self.assertEqual(late['start_pdt'],'2026-09-24 00:08')
        self.assertLessEqual(late['airmass_max_actual'],1.5)
        self.assertLessEqual(science_end,pd.Timestamp(p['sequence'][-1]['start_utc']))
        for v in p['primaries']:
            if v['name'] in {'P8548','P12457'}:
                self.assertTrue(v['host_contaminated'])
                self.assertIn('AGN-only S/N is lower',v['caution'])
        self.assertEqual([v['plan']['seconds_each'] for v in p['sequence'] if v['role']=='standard'],[30,5])
        for v in p['primaries']+p['backups']:
            self.assertEqual(v['plan']['exposures'],2)
            self.assertEqual(v['plan']['seconds_each'],300)
            self.assertEqual(v['plan']['visit_minutes'],20)
            self.assertGreaterEqual(v['plan']['predicted_snr'],4.99)
            self.assertGreaterEqual(v['plan']['visit_minutes'],2*v['plan']['seconds_each']/60+10)

    def test_backup_rows_fit_their_primary_slots_and_are_not_primaries(self):
        primaries = {f"P{v['rank']:02}":v for v in self.packet['primaries']}
        names = {v['name'] for v in primaries.values()}
        self.assertEqual(len(self.backup_rows), len(self.packet['backups']))
        for v,row in zip(self.packet['backups'],self.backup_rows):
            p = primaries[v['replaces']]
            self.assertNotIn(v['name'],names)
            self.assertEqual(v['start_utc'],p['start_utc'])
            self.assertLessEqual(v['end_utc'],p['end_utc'])
            self.assertEqual(row['name'],v['name'])
            self.assertEqual(row['exptime'],f"SET {v['plan']['seconds_each']}")
            self.assertIn(v['replaces'],row['Note'])
        for p in primaries.values():
            self.assertEqual(len(p['backups']),3)
            self.assertEqual(len(set(p['backups'])),3)

    def test_geometry_from_csv_coordinates_over_entire_visits(self):
        site = EarthLocation.from_geodetic(-116.865*u.deg,33.3563*u.deg,1712*u.m)
        pairs = list(zip(self.packet['sequence'],self.rows))+list(zip(self.packet['backups'],self.backup_rows))
        for v,row in pairs:
            self.assertLessEqual(len(row['Note']),24)
            self.assertLessEqual(len(row['Comment']),1024)
            self.assertEqual(len(row),12)
            times = Time(pd.date_range(v['start_utc'],v['end_utc'],freq='10s').to_pydatetime())
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore',message='Tried to get polar motions')
                frame = AltAz(obstime=times,location=site,pressure=0*u.hPa)
                coord = SkyCoord(row['RA'],row['DECL'],unit=(u.hourangle,u.deg)).transform_to(frame)
                moon = get_body('moon',times,site).transform_to(frame)
            self.assertTrue(np.all(coord.alt.deg>0),v['name'])
            self.assertLessEqual(float(coord.secz.max()),float(row['airmass_max'])+1e-5,v['name'])
            self.assertGreaterEqual(float(coord.separation(moon).deg.min()),40-1e-5,v['name'])

    def test_page_downloads_and_public_spectra_match_authoritative_files(self):
        def payload(path,kind):
            return json.loads(re.search(r'<script id="'+kind+r'" type="application/json">(.*?)</script>',path.read_text(),re.S).group(1))
        public_source = payload(ROOT/'docs/index.html','candidate-data')
        public_targets = {t['name']:t for t in public_source['targets']}
        for path,private in [(ROOT/'docs/sep23_primaries.html',False),(DEST/'sep23_primaries_local.html',True)]:
            page = payload(path,'primary-data')
            self.assertEqual([t['name'] for t in page['targets']],[v['name'] for v in self.packet['primaries']])
            expected = ['sep23_primaries_ngps.csv']+(['sep23_backups_ngps.csv'] if private else [])
            self.assertEqual(set(page['packet']['files']),set(expected))
            for name in expected:
                self.assertEqual(page['packet']['files'][name],(DEST/name).read_text())
            if not private:
                self.assertNotIn('backups',page['packet'])
                for target in page['targets']:
                    self.assertNotIn('science',target)
                    self.assertNotIn('exposure_plans',target)
                    self.assertEqual(target['spec'],public_targets[target['name']]['spec'])
        for name in ['sep23_primaries_ngps.csv','sep23_backups_ngps.csv']:
            self.assertEqual((ROOT/'observing'/name).resolve(),(DEST/name).resolve())


if __name__=='__main__':
    unittest.main()
