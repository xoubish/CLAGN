"""Scientific invariants: flux units, gaps, calibration, and propagated errors."""
import unittest
import numpy as np
from analysis import integration_weights, line_index, common_support_pair


def spectrum(z=0., scale=1.):
    rest = np.arange(4600., 5300., .8)
    flux = 10 + .002*(rest-4862.68) + 7*np.exp(-.5*((rest-4862.68)/24)**2)
    return dict(wave=rest*(1+z), flux=flux*scale, error=np.full(len(rest), .1*scale), z=z,
                quality=np.zeros(len(rest), int), channel=np.full(len(rest), 'R'), source='NGPS')


class MeasurementTests(unittest.TestCase):
    def test_exact_linear_integral_on_irregular_grid(self):
        w = np.array([0., .6, 1.5, 2.1, 3.])
        h, cov = integration_weights(w, np.ones(len(w), bool), [(.25, 2.75)])
        self.assertAlmostEqual(h @ (2*w+3), (2.75**2+3*2.75)-(.25**2+3*.25))
        self.assertAlmostEqual(cov, 1.)

    def test_masked_gap_is_not_interpolated(self):
        h, cov = integration_weights(np.arange(6.), [1, 1, 0, 0, 1, 1], [(0, 5)])
        self.assertAlmostEqual(cov, .4)
        self.assertAlmostEqual(h.sum(), 2.)

    def test_missing_wavelength_gap_is_not_interpolated(self):
        w = np.array([0, 1, 2, 20, 21, 22.])
        h, cov = integration_weights(w, np.ones(6, bool), [(0, 22)])
        self.assertAlmostEqual(h.sum(), 4.)

    def test_grey_calibration_preserves_ew(self):
        a, b = line_index(spectrum(), 'Hb'), line_index(spectrum(scale=2.7), 'Hb')
        self.assertAlmostEqual(a['wing_ew'], b['wing_ew'])
        self.assertAlmostEqual(b['wing_flux']/a['wing_flux'], 2.7)
        self.assertAlmostEqual(a['wing_ew_err_diagonal'], b['wing_ew_err_diagonal'])

    def test_rest_ew_and_observed_flux_differential(self):
        a, b = line_index(spectrum(), 'Hb'), line_index(spectrum(z=.5), 'Hb')
        self.assertAlmostEqual(a['wing_ew'], b['wing_ew'])
        self.assertAlmostEqual(b['wing_flux']/a['wing_flux'], 1.5)

    def test_missing_errors_never_invented(self):
        s = spectrum(); s['error'][:] = np.nan
        self.assertTrue(np.isnan(line_index(s, 'Hb')['wing_ew_err_diagonal']))

    def test_sideband_uncertainty_propagation(self):
        rng = np.random.default_rng(439)
        s = spectrum(); expected = line_index(s, 'Hb')['wing_ew_err_diagonal']
        draws = [line_index(dict(s, flux=s['flux']+rng.normal(0, s['error'])), 'Hb')['wing_ew'] for _ in range(400)]
        self.assertLess(abs(np.std(draws)/expected-1), .15)

    def test_missing_line_pixels_fail_coverage(self):
        s = spectrum(); rest = s['wave']
        s['flux'][(rest > 4800) & (rest < 4820)] = np.nan
        self.assertFalse(line_index(s, 'Hb')['valid'])

    def test_repeat_pair_uses_identical_support_without_invented_errors(self):
        a, b = spectrum(), spectrum()
        a['flux'][(a['wave'] > 4800) & (a['wave'] < 4820)] = np.nan
        x, y = common_support_pair(a, b)
        np.testing.assert_array_equal(np.isfinite(x['flux']), np.isfinite(y['flux']))
        gap = (x['wave'] > 4800) & (x['wave'] < 4820)
        self.assertTrue(np.isnan(x['flux'][gap]).all())
        self.assertTrue(np.isnan(x['error']).all())


if __name__ == '__main__':
    unittest.main()
