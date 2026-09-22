"""Check the local telescope packet against geometry, slot limits, one setting and page exports."""
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

    def test_sequence_keeps_protected_primaries_and_one_setting(self):
        p = self.packet
        S = p['settings']
        names = [v['name'] for v in p['primaries']]
        self.assertEqual(len(names), len(set(names)))
        self.assertGreaterEqual(len(names), len(self.plan['protected']))
        self.assertTrue(set(self.plan['protected']) <= set(names))
        if not self.plan.get('user_selection'):
            self.assertTrue(set(self.plan['original_primaries']) <= set(self.plan['protected']))
        self.assertEqual(set(names)-set(self.plan['protected']), set(p['promoted']))
        self.assertEqual(S['goal_snr'], 5)
        self.assertEqual((S['slit_arcsec'], S['binspat'], S['binspect']), (1.5, 2, 3))
        for a, b in zip(p['sequence'], p['sequence'][1:]):
            self.assertLessEqual(pd.Timestamp(a['end_utc']), pd.Timestamp(b['start_utc']))
        science_end = pd.Timestamp(p['primaries'][-1]['end_utc'])
        self.assertLessEqual(science_end, pd.Timestamp(p['sequence'][-1]['start_utc']))
        # Unscheduled time inside the science block is the delay reserve.
        block = pd.Timestamp(S['science_end_pdt'], tz='America/Los_Angeles')-pd.Timestamp(S['science_start_pdt'], tz='America/Los_Angeles')
        used = sum(v['plan']['visit_minutes'] for v in p['primaries'])
        self.assertGreaterEqual(block.total_seconds()/60-used, S['reserve_minutes'])
        self.assertEqual(p['reserved']['minutes'], block.total_seconds()/60-used)
        for v in p['primaries']:
            if v['name'] in {'P8548', 'P12457'}:
                self.assertTrue(v['host_contaminated'])
                self.assertIn('AGN-only S/N is lower', v['caution'])
        self.assertEqual([v['plan']['seconds_each'] for v in p['sequence'] if v['role'] == 'standard'], [30, 5])
        for v in p['primaries']+p['backups']:
            n = v['plan']['exposures']
            self.assertIn(n, (2, 3, 4))
            self.assertEqual(v['plan']['seconds_each'], S['seconds_each'])
            self.assertEqual(v['plan']['visit_minutes'], int(n*S['seconds_each']/60+(n-1)*S['readout_minutes']+S['overhead_minutes']))
            if v['role'] == 'backup' or v['name'] not in self.plan['protected']:
                self.assertGreaterEqual(v['snr_per_angstrom'], S['goal_snr']-0.01)
            else:
                self.assertIn('S/N', v['caution']) if v['snr_per_angstrom'] < S['goal_snr'] else None
            self.assertLessEqual(v['airmass_max_actual'], v['plan']['airmass']+1e-6)
            self.assertGreaterEqual(v['moon_min'], S['moon_min_deg']-1e-6)
            self.assertGreaterEqual(v['plan']['visit_minutes'], n*S['seconds_each']/60+S['overhead_minutes'])

    def test_order_matches_the_plan_and_slot_table(self):
        table = pd.read_csv(DEST/'snr5_slot_table.csv')
        for v, choice in zip(self.packet['primaries'], self.plan['sequence']):
            self.assertEqual(v['name'], choice['name'])
            self.assertEqual(v['start_pdt'], choice['start_pdt'])
            row = table[(table.name == v['name']) & (table.start_pdt == v['start_pdt'])]
            self.assertEqual(len(row), 1)
            self.assertAlmostEqual(float(row.snr_per_angstrom.iloc[0]), v['snr_per_angstrom'], places=6)
            if v['name'] not in self.plan['protected']:
                self.assertTrue(bool(row.meets_goal.iloc[0]))
        starts = [pd.Timestamp(v['start_pdt']) for v in self.packet['primaries']]
        self.assertEqual(starts, sorted(starts))

    def test_backup_rows_fit_their_primary_slots_and_are_not_primaries(self):
        primaries = {f"P{v['rank']:02}": v for v in self.packet['primaries']}
        names = {v['name'] for v in primaries.values()}
        self.assertEqual(len(self.backup_rows), len(self.packet['backups']))
        for v, row in zip(self.packet['backups'], self.backup_rows):
            p = primaries[v['replaces']]
            self.assertNotIn(v['name'], names)
            self.assertEqual(v['start_utc'], p['start_utc'])
            self.assertLessEqual(v['end_utc'], p['end_utc'])
            self.assertEqual(row['name'], v['name'])
            self.assertEqual(row['exptime'], f"SET {v['plan']['seconds_each']}")
            self.assertIn(v['replaces'], row['Note'])
        for p in primaries.values():
            self.assertGreaterEqual(len(p['backups']), 1)
            self.assertLessEqual(len(p['backups']), 3)
            self.assertEqual(len(set(p['backups'])), len(p['backups']))

    def test_geometry_and_setting_from_csv_rows_over_entire_visits(self):
        site = EarthLocation.from_geodetic(-116.865*u.deg, 33.3563*u.deg, 1712*u.m)
        S = self.packet['settings']
        pairs = list(zip(self.packet['sequence'], self.rows))+list(zip(self.packet['backups'], self.backup_rows))
        for v, row in pairs:
            self.assertLessEqual(len(row['Note']), 24)
            self.assertLessEqual(len(row['Comment']), 1024)
            self.assertEqual(len(row), 12)
            self.assertEqual(row['slitwidth'], f"SET {S['slit_arcsec']}")
            self.assertEqual((row['binspat'], row['binspect'], row['slitangle']), (str(S['binspat']), str(S['binspect']), 'PA'))
            times = Time(pd.date_range(v['start_utc'], v['end_utc'], freq='10s').to_pydatetime())
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', message='Tried to get polar motions')
                frame = AltAz(obstime=times, location=site, pressure=0*u.hPa)
                coord = SkyCoord(row['RA'], row['DECL'], unit=(u.hourangle, u.deg)).transform_to(frame)
                moon = get_body('moon', times, site).transform_to(frame)
            self.assertTrue(np.all(coord.alt.deg > 0), v['name'])
            self.assertLessEqual(float(coord.secz.max()), float(row['airmass_max'])+1e-5, v['name'])
            self.assertGreaterEqual(float(coord.separation(moon).deg.min()), 40-1e-5, v['name'])

    def test_page_downloads_and_public_spectra_match_authoritative_files(self):
        def payload(path, kind):
            return json.loads(re.search(r'<script id="'+kind+r'" type="application/json">(.*?)</script>', path.read_text(), re.S).group(1))
        public_source = payload(ROOT/'docs/index.html', 'candidate-data')
        public_targets = {t['name']: t for t in public_source['targets']}
        for path, private in [(ROOT/'docs/sep23_primaries.html', False), (DEST/'sep23_primaries_local.html', True)]:
            page = payload(path, 'primary-data')
            expected_names = [v['name'] for v in self.packet['primaries'] if private or (not v['plan'].get('reference_private') and v['name'] in public_targets)]
            self.assertEqual([t['name'] for t in page['targets']], expected_names)
            self.assertEqual(page['packet']['settings']['slit_arcsec'], self.packet['settings']['slit_arcsec'])
            expected = ['sep23_primaries_ngps.csv']+(['sep23_backups_ngps.csv'] if private else [])
            self.assertEqual(set(page['packet']['files']), set(expected))
            for name in expected:
                self.assertEqual(page['packet']['files'][name], (DEST/name).read_text())
            if not private:
                self.assertNotIn('backups', page['packet'])
                for target in page['targets']:
                    self.assertNotIn('science', target)
                    self.assertNotIn('exposure_plans', target)
                    self.assertEqual(target['spec'], public_targets[target['name']]['spec'])
        for name in ['sep23_primaries_ngps.csv', 'sep23_backups_ngps.csv']:
            self.assertEqual((ROOT/'observing'/name).resolve(), (DEST/name).resolve())


if __name__ == '__main__':
    unittest.main()
