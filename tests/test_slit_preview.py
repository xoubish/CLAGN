"""Physical direction and scheduling regressions for the image slit overlay."""
import importlib
from pathlib import Path
import sys
import unittest

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, AltAz
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
slit = importlib.import_module('slit_preview')


class SlitPreviewTests(unittest.TestCase):
    def test_axis_points_toward_zenith_in_image_frame(self):
        times = Time(['2026-09-24T04:12:00', '2026-09-24T04:56:00'])
        for ra, dec in [(245.00754, 29.590792), (330, 20), (10, 60)]:
            pa, x = slit.slit_angles(ra, dec, times)
            target = SkyCoord(ra*u.deg, dec*u.deg)
            frame = AltAz(obstime=times, location=slit.PALOMAR.location)
            zenith = SkyCoord(az=np.zeros(2)*u.deg, alt=np.full(2,90)*u.deg, frame=frame).icrs
            # Independent great-circle direction to the zenith. Aberration
            # makes this differ slightly from the local differential direction.
            expected = target.position_angle(zenith).deg % 180
            delta = (pa-expected+90) % 180-90
            np.testing.assert_allclose(delta, 0, atol=0.03)
            self.assertTrue(np.isfinite(x).all())

    def test_midpoint_and_endpoints_survive_utc_date_rollover(self):
        start = Time('2026-09-24T06:48:00').mjd
        target = dict(name='example', ra=330, dec=20, nights=[dict(
            night='sep23', curve=[[start,1.2,60],[start+5/1440,1.1,60]])])
        slot = dict(start_utc='2026-09-24T06:48:00Z', end_utc='2026-09-24T07:04:00Z')
        preview = slit.previews(target, {'example':slot})['sep23']
        sample = preview['samples'][preview['default_index']]
        self.assertTrue(preview['scheduled'])
        self.assertAlmostEqual(sample[0], Time('2026-09-24T06:56:00').mjd)
        self.assertEqual(sample[3], 'visit midpoint')
        self.assertEqual({s[3] for s in preview['samples'] if s[3]},
                         {'visit start','visit midpoint','visit end'})

    def test_unscheduled_preview_uses_observable_samples(self):
        m = Time('2026-10-27T04:00:00').mjd
        target = dict(name='example',ra=330,dec=20,nights=[dict(night='oct26',curve=[
            [m,1.0,20], [m+.01,None,60], [m+.02,2.1,60],
            [m+.03,1.5,60], [m+.04,1.3,60]])])
        p = slit.previews(target, {'example':{}})['oct26']
        self.assertFalse(p['scheduled'])
        self.assertEqual(len(p['samples']),2)
        self.assertEqual(p['samples'][p['default_index']][0],m+.04)


if __name__ == '__main__':
    unittest.main()
