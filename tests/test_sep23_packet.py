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
        self.assertEqual([v['plan']['seconds_each'] for v in p['sequence'] if v['role'] == 'standard'], [60, 10])
        for v in p['primaries']+p['backups']:
            n = v['plan']['exposures']
            self.assertIn(n, (2, 3, 4))
            self.assertEqual(v['plan']['seconds_each'], S['seconds_each'])
            self.assertEqual(v['plan']['visit_minutes'], {2: 16, 3: 22, 4: 28}[n])
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
        locked = self.plan.get('user_selection', {}).get('preserve_sequence', [])
        if locked:
            self.assertEqual([(v['name'], v['start_pdt']) for v in self.packet['primaries']], [(v['name'], v['start_pdt']) for v in locked])
        starts = [pd.Timestamp(v['start_pdt']) for v in self.packet['primaries']]
        self.assertEqual(starts, sorted(starts))

    def test_displayed_windows_cover_the_full_primary_visits(self):
        options = json.loads((ROOT/'data/reselection_2026-09-20/airmass_options.json').read_text())
        for visit in self.packet['primaries']:
            night = next(n for n in options['windows'][visit['name']] if n['night'] == 'sep23')
            self.assertTrue(any(pd.Timestamp(r['start_utc'], tz='UTC') <= pd.Timestamp(visit['start_utc'])
                                and pd.Timestamp(r['end_utc'], tz='UTC') >= pd.Timestamp(visit['end_utc'])
                                for r in night['tier_ranges']['extended']), visit['name'])

    def test_backup_rows_fit_their_primary_slots_and_are_not_primaries(self):
        primaries = {f"P{v['rank']:02}": v for v in self.packet['primaries']}
        names = {v['name'] for v in primaries.values()}
        self.assertEqual(len(self.backup_rows), len(self.packet['backups']))
        for v, row in zip(self.packet['backups'], self.backup_rows):
            p = primaries[v['replaces']]
            self.assertNotIn(v['name'], names)
            self.assertEqual(v['start_utc'], p['start_utc'])
            self.assertEqual(v['end_utc'], p['end_utc'])
            self.assertEqual(v['plan']['exposures'], p['plan']['exposures'])
            self.assertLess(v['airmass_max_actual'], 1.5)
            self.assertGreater(v['moon_min'], 40)
            self.assertLess(v['eligibility']['r_mag'], 19)
            self.assertTrue(set(v['eligibility']['quasar_basis'].split(';')) & {'expanded_DR16_QSO', 'full_DR16Q_catalog'})
            self.assertFalse(v['plan']['reference_private'])
            self.assertEqual(v['eligibility']['reference_mjd'], v['plan']['reference_mjd'])
            self.assertEqual(row['name'], v['name'])
            self.assertEqual(row['exptime'], f"SET {v['plan']['seconds_each']}")
            self.assertIn(v['replaces'], row['Note'])
        for p in primaries.values():
            self.assertEqual(len(p['backups']), 2)
            self.assertEqual(len(set(p['backups'])), len(p['backups']))

    def test_backups_have_real_public_spectral_references_and_magnitude_provenance(self):
        import sys
        sys.path.insert(0, str(ROOT/'scripts'))
        from spectral_utils import archival_records, accepted_reference
        targets = pd.read_csv(ROOT/'data/reselection_2026-09-20/compact_review_objects.csv').set_index('name', drop=False)
        science = pd.read_csv(ROOT/'data/reselection_2026-09-20/three_night_review/science_and_sensitivity.csv').set_index('name')
        for visit in self.packet['backups']:
            name = visit['name']
            ref = accepted_reference(targets.loc[name], [r for r in archival_records(name) if not r.get('proprietary')])
            self.assertIsNotNone(ref, name)
            self.assertEqual(ref[0]['mjd'], visit['eligibility']['reference_mjd'])
            latest = science.loc[name, 'ztf_r_latest180_mag']
            expected = latest if pd.notna(latest) else targets.loc[name, 'r_planning']
            self.assertEqual(expected, visit['eligibility']['r_mag'])

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
            if v['role'] == 'backup':
                self.assertLess(float(coord.secz.max()), 1.5, v['name'])
                self.assertGreater(float(coord.separation(moon).deg.min()), 40, v['name'])

    def test_single_page_embeds_the_sequence_and_csvs(self):
        def payload(path):
            return json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>', path.read_text(), re.S).group(1))
        primaries = self.packet['primaries']
        for path, private in [(ROOT/'docs/index.html', False), (ROOT/'data/reselection_2026-09-20/observer_page_local.html', True)]:
            page = payload(path)
            public_names = {t['name'] for t in payload(ROOT/'docs/index.html')['targets']}
            expected = {v['name'] for v in primaries if private or v['name'] in public_names}
            self.assertEqual(set(page['sep23_sequence']), expected)
            for name, s in page['sep23_sequence'].items():
                v = next(x for x in primaries if x['name'] == name)
                self.assertEqual((s['rank'], s['start_pdt'], s['exposures']), (v['rank'], v['start_pdt'], v['plan']['exposures']))
            self.assertEqual(page['files']['sep23_primaries_ngps.csv'], self.packet['files' if private else 'public_files']['sep23_primaries_ngps.csv'])
            if private:
                self.assertEqual(page['files']['sep23_backups_ngps.csv'], (DEST/'sep23_backups_ngps.csv').read_text())
                self.assertEqual(len(page['backups']), len(self.packet['backups']))
            else:
                self.assertEqual(page['files']['sep23_backups_ngps.csv'], self.packet['public_files']['sep23_backups_ngps.csv'])
                self.assertEqual(len(page['backups']), 2*len(primaries))
                self.assertTrue(all(b['name'] in public_names for b in page['backups']))
                self.assertNotIn('see the local page', path.read_text())
                self.assertNotIn('proprietary', path.read_text())
                for t in page['targets']:
                    if any(v.get('reference_private') for v in self.plan['plans'].get(t['name'], {}).values()):
                        self.assertEqual(t['sep23_slots'], [])
                    if t.get('science', {}).get('reference_private'):
                        self.assertNotIn('continuum_AB', t['science'])
                        self.assertNotIn('post_spectrum_optical_trigger', t['science'])
                for v in self.packet['primaries']:
                    if v['plan'].get('reference_private'):
                        self.assertIsNone(page['sep23_sequence'][v['name']]['snr_per_angstrom'])
                        row = next(r for r in csv.DictReader(page['files']['sep23_primaries_ngps.csv'].splitlines()) if r['name']==v['name'])
                        self.assertNotIn('model continuum S/N', row['Comment'])
            self.assertEqual(len(page['sequence_rows']), len(self.packet['sequence']))
            self.assertEqual(page['decisions']['setting']['slit_arcsec'], self.packet['settings']['slit_arcsec'])
            self.assertTrue(page['run']['nights'] and page['run']['calibrations'] and page['run']['procedure'])
            self.assertEqual([t['code'] for t in page['targets']], [f'{i:03d}' for i in range(1, len(page['targets'])+1)])
        for name in ['sep23_primaries_ngps.csv', 'sep23_backups_ngps.csv']:
            self.assertEqual((ROOT/'observing'/name).resolve(), (DEST/name).resolve())

if __name__ == '__main__':
    unittest.main()
