"""Regression checks for the September calculation and ingestion audit."""
import importlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from spectral_utils import bin_indices, centered_record, selected_records, accepted_reference


class SpectralTests(unittest.TestCase):
    def test_center_samples_stay_in_their_bins(self):
        np.testing.assert_array_equal(bin_indices([3000, 3006, 3012], [2996, 3000, 3006, 3012, 3016]), [-1, 0, 1, 2, 3])
        legacy = dict(wave=[3600., 3606., 3612.], flux=[1., 2., 3.])
        corrected = centered_record(legacy)
        self.assertEqual(corrected['wave'], [3603., 3609., 3615.])
        self.assertEqual(centered_record(corrected), corrected)
        self.assertEqual(legacy['wave'][0], 3600.)

    def test_display_and_reference_use_same_reduction(self):
        wave = np.arange(4600., 5200., 6.)
        old = dict(source='SDSS', mjd=59000., wave=wave.tolist(), flux=np.ones(len(wave)).tolist(), grid_version=2)
        current = dict(old, run2d='v6_2_1', flux=(np.ones(len(wave))*2).tolist())
        self.assertEqual(selected_records([old, current])[0]['flux'][0], 2.)
        reference, _, continuum = accepted_reference(dict(name='test', z=0), [old, current])
        self.assertEqual(reference['run2d'], 'v6_2_1')
        self.assertAlmostEqual(continuum, 2.)

    def test_missing_measurements_are_unclassified_and_both_sides_required(self):
        ingest = importlib.import_module('12_ngps_ingest')
        self.assertIn('unclassified', ingest.verdict('T', {}, {}))
        self.assertNotIn('no change', ingest.verdict('T', {'Hb': 20}, None))
        w = np.arange(4700., 5000., 6.)
        self.assertNotIn('Hb', ingest.ew_indices(w, np.ones(len(w)), 0))
        verdict = ingest.verdict('T', {'Hb': 2., 'Ha': 30.}, {'EW_Hb_rest': 10., 'EW_Ha_rest': 10.})
        self.assertIn('Hb decreased', verdict)
        self.assertIn('Ha increased', verdict)
        self.assertNotIn('turn-on', verdict)

    def test_current_catalog_precedes_legacy(self):
        ingest = importlib.import_module('12_ngps_ingest')
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)/'reselection_2026-09-20'
            p.mkdir()
            pd.DataFrame([dict(name='P2190', z=.3)]).to_csv(p/'compact_review_objects.csv', index=False)
            pd.DataFrame([dict(name='P2190', z=.9), dict(name='Old', z=.2)]).to_csv(Path(directory)/'targets_old.csv', index=False)
            with patch.object(ingest, 'DATA', directory):
                self.assertEqual(ingest.target_catalog().loc['P2190', 'z'], .3)
                self.assertIn('Old', ingest.target_catalog().index)

    def test_ingest_writes_centered_spectrum_and_actual_exposure_date(self):
        ingest = importlib.import_module('12_ngps_ingest')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            catalog = root/'reselection_2026-09-20'
            catalog.mkdir()
            pd.DataFrame([dict(name='Synthetic', z=.3)]).to_csv(catalog/'compact_review_objects.csv', index=False)
            source = root/'Synthetic_R.csv'
            pd.DataFrame(dict(wave=[6000., 6006., 6012.], flux=[1., 2., 3.])).to_csv(source, index=False)
            output = root/'ngps_spectra'
            with patch.object(ingest, 'DATA', directory), patch.object(ingest, 'OUT', str(output)), patch.object(sys, 'argv', ['12', str(source), '--mjd', '61210.25']):
                ingest.main()
            frame = pd.read_csv(output/'Synthetic.csv')
            np.testing.assert_array_equal(frame.wave_A, [6000., 6006., 6012.])
            np.testing.assert_array_equal(frame.flux, [1., 2., 3.])
            self.assertEqual(json.loads((output/'Synthetic.json').read_text())['mjd'], 61210.25)
            self.assertIn('unclassified', pd.read_csv(root/'ngps_lines.csv').verdict.iloc[0])

    def test_hbeta_fit_masks_oxygen_and_has_red_continuum(self):
        fit = importlib.import_module('46_line_history')
        w = np.arange(4700., 5131., 2.)
        flux = 5 + 10*np.exp(-.5*((w-4862.68)/25)**2)
        flux += 100*np.exp(-.5*((w-5008.24)/2)**2)
        result = fit.fit_line(w, flux, 0, 'Hbeta')
        self.assertIsNotNone(result)
        self.assertTrue((result['x'] >= 5080).any())
        self.assertFalse(((result['x'] >= 4995) & (result['x'] <= 5022)).any())
        self.assertAlmostEqual(result['flux'], 10*25*np.sqrt(2*np.pi), delta=1)


class PlanningTests(unittest.TestCase):
    def test_long_visits_do_not_truncate_added_readout(self):
        planner = importlib.import_module('41_september_snr5_plan')
        self.assertEqual(planner.visit_minutes(2), 16)
        self.assertGreaterEqual(planner.visit_minutes(3), 21.6)
        self.assertGreaterEqual(planner.visit_minutes(4), 27.2)

    @unittest.skipUnless((ROOT/'data/reselection_2026-09-20/ngps_etc').exists(), 'Requires cached official ETC')
    def test_official_etc_applies_zenith_seeing_once(self):
        model = importlib.import_module('22_september_etc')
        args = model.ETC.parser.parse_args(['G', '500', '508', 'EXPTIME', '300', '-slit', 'SET', '1.5', '-seeing', '1.3', '500', '-airmass', '1.8', '-skymag', '18.5', '-mag', '18', '-magsystem', 'AB', '-magfilter', 'match', '-noslicer'])
        model.ETC.check_inputs_add_units(args)
        self.assertAlmostEqual(args.seeing[0].value, 1.3*1.8**.6)


if __name__ == '__main__':
    unittest.main()
