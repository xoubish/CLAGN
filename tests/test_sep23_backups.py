"""Boundary and provenance checks for backup selection, without the local data cache."""
import importlib
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import sep23_backups as backups


class BackupSelectionTests(unittest.TestCase):
    def setUp(self):
        self.target = dict(name='Q', origin='expanded_DR16_QSO', r_planning=18.,
                           r_source='historical SDSS catalog', review_region='Zeltyn',
                           manifold_cl_neighbor_fraction=.2, balmer_pair_in_range=True)

    def test_strict_brightness_and_no_fallback_over_a_faint_recent_measurement(self):
        for mag in [19., 19.1]:
            self.assertIsNone(backups.catalog_evidence(self.target, {'ztf_r_latest180_mag': mag}))
        self.assertIsNotNone(backups.catalog_evidence(self.target, {'ztf_r_latest180_mag': 18.99}))
        evidence = backups.catalog_evidence(self.target, {'ztf_r_latest180_mag': float('nan')})
        self.assertEqual(evidence['r_mag'], 18.)
        self.assertEqual(evidence['r_source'], 'historical SDSS catalog')
        self.assertIsNone(backups.catalog_evidence(dict(self.target, r_planning=float('nan')), {}))

    def test_galaxy_or_generic_agn_identity_does_not_establish_a_quasar(self):
        galaxy = dict(self.target, origin='expanded_DR16_GALAXY', agn_identity_eligible=True)
        self.assertIsNone(backups.catalog_evidence(galaxy, {}))
        self.assertIsNotNone(backups.catalog_evidence(dict(galaxy, origins='full_DR16Q_catalog'), {}))

    def test_manifold_priority_public_spectrum_and_primary_exclusion(self):
        targets = pd.DataFrame([
            dict(self.target, name='region', field_status='clear', r_planning=18.9),
            dict(self.target, name='bright_control', field_status='clear', r_planning=14., review_region='Other manifold region'),
            dict(self.target, name='primary', field_status='clear'),
            dict(self.target, name='private_identity', field_status='clear'),
            dict(self.target, name='no_spectrum', field_status='clear'),
        ])
        science = pd.DataFrame({'name': targets.name, 'ztf_r_latest180_mag': float('nan')}).set_index('name')
        def reference(target, records):
            self.assertTrue(all(not r.get('proprietary') for r in records))
            return None if target['name'] == 'no_spectrum' else ({'mjd': 55000}, 18., 1.)
        with patch.object(backups, 'archival_records', return_value=[{'proprietary': True}, {'mjd': 55000}]), patch.object(backups, 'accepted_reference', side_effect=reference):
            selected = backups.candidates(targets, science, set(targets.name)-{'private_identity'}, {'primary'})
        self.assertEqual([c['target']['name'] for c in selected], ['region', 'bright_control'])

    def test_full_primary_duration_and_strict_geometry_before_etc(self):
        planner = importlib.import_module('41_september_snr5_plan')
        candidate = dict(target=self.target, reference=(55000, 18., False, 1.))
        for exposures, duration in [(2, 16), (3, 22), (4, 28)]:
            primary = dict(start_pdt='2026-09-23 21:20', plan=dict(exposures=exposures, visit_minutes=duration))
            for x, moon in [(1.5, 41), (1.4, 40), (1.51, 50)]:
                slot = dict(start_pdt=primary['start_pdt'], airmass_max=x, moon_deg=moon)
                with patch.object(planner, 'candidate_slots', return_value=([slot], {})) as geometry, patch.object(planner, 'evaluate') as etc:
                    self.assertEqual(backups.select([candidate], primary, {}, {'Q': 0}, None, '/tmp'), [])
                    geometry.assert_called_once_with({}, 0, duration)
                    etc.assert_not_called()


if __name__ == '__main__':
    unittest.main()
