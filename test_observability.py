"""Regression checks for the apparent-Moon frame error and window integration."""
import importlib
import unittest
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, get_body
from astropy.time import Time
from astropy.utils import iers

iers.conf.auto_download=False
iers.conf.auto_max_age=None
obs=importlib.import_module('05_observability')


class ApparentGeometryTests(unittest.TestCase):
    def test_september_p5125_regression_and_independent_frame(self):
        p=pd.read_csv('data/targets_sep23.csv').set_index('name').loc['P5125']
        t=Time('2026-09-24T06:05:00')
        c=SkyCoord(p.ra*u.deg,p.dec*u.deg)
        actual=obs.apparent_moon_separation(c,t).deg
        moon=get_body('moon',t,obs.PALOMAR.location)
        independent=c.transform_to(moon.frame).separation(moon).deg
        self.assertAlmostEqual(actual,18.14,places=1)
        self.assertAlmostEqual(actual,independent,places=3)

    def test_broadcast_preserves_target_order(self):
        c=SkyCoord([0,90,180]*u.deg,[0,30,50]*u.deg)
        t=Time(['2026-10-27T04:00:00','2026-10-27T08:00:00'])
        x=obs.apparent_moon_separation(c[:,None],t).deg
        self.assertEqual(x.shape,(3,2))
        for i in range(3):
            np.testing.assert_allclose(x[i],obs.apparent_moon_separation(c[i],t).deg,atol=1e-8)

    def test_always_visible_target_does_not_gain_extra_grid_interval(self):
        # North pole altitude at Palomar is >30 degrees throughout every window.
        r=obs.score(pd.DataFrame({'ra':[0.0],'dec':[90.0]}),step_min=10)
        for name,(date,part) in obs.NIGHTS.items():
            a,b,_,_=obs.night_window(date,part)
            self.assertAlmostEqual(r.loc[0,f'hrs_{name}'],(b-a).to_value(u.hour),places=8)


if __name__=='__main__':unittest.main()
