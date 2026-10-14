"""Check complete-visit geometry and queue order, including window boundaries."""
import json
import re
import sys
import unittest
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, SkyCoord, get_body
from astropy.time import Time
from astropy.utils import iers

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from observer_visibility import LOCATION


class VisibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>',
                                    (ROOT/'docs/index.html').read_text(),re.S)[1])

    def test_queue_order_and_every_scheduled_visit_is_available(self):
        for night,V in self.data['visibility'].items():
            if night=='sep23':
                expected=[n for n,r in sorted(self.data['sep23_sequence'].items(),key=lambda p:p[1]['rank'])]
            else:expected=[r['name'] for r in self.data['october']['rows'] if r['night']==night]
            self.assertEqual([r['name'] for r in V['rows']],expected)
            for r in V['rows']:
                self.assertEqual(r['planned_end']-r['planned_start'],r['duration']*60000)
                index=round((r['planned_start']-V['sample_start'])/V['step'])
                self.assertGreater(r['samples'][index][2],0,(night,r['name']))
                for j,(x,moon,tier) in enumerate(r['samples']):
                    if not tier:continue
                    self.assertLessEqual(x,r['airmass_limit'])
                    self.assertGreaterEqual(moon,40)
                    self.assertLessEqual(x,{1:1.5,2:1.8,3:2}[tier])
                    start=V['sample_start']+j*V['step'];end=start+r['duration']*60000
                    self.assertGreaterEqual(start,V['start'])
                    self.assertLessEqual(end,V['end'])
                    self.assertFalse(any(start<s['end'] and end>s['start'] for s in V['standards']))

    def test_window_edges_independently_at_twenty_second_cadence(self):
        iers.conf.auto_download=False;iers.conf.auto_max_age=None
        coordinates={t['name']:(t['ra'],t['dec']) for t in self.data['targets']}
        coordinates.update({t['name']:(float(t['ra']),float(t['dec'])) for t in self.data['october']['rows']})
        for night,V in self.data['visibility'].items():
            for duration in sorted({r['duration'] for r in V['rows']}):
                starts=[];ras=[];decs=[];limits=[]
                for r in V['rows']:
                    if r['duration']!=duration:continue
                    ra,dec=coordinates[r['name']]
                    for w in r['windows']:
                        for start in [w['start'],w['end']]:
                            starts.append(start/1000);ras.append(ra);decs.append(dec)
                            limits.append(min(r['airmass_limit'],{1:1.5,2:1.8,3:2}[w['tier']]))
                times=Time(np.array(starts)[:,None]+np.arange(0,duration*60+1,20),format='unix')
                frame=AltAz(obstime=times,location=LOCATION,pressure=0*u.hPa)
                aa=SkyCoord(ras*u.deg,decs*u.deg)[:,None].transform_to(frame)
                self.assertTrue(np.all(aa.alt.deg>0))
                self.assertTrue(np.all(aa.secz.value.max(axis=1)<=limits),night)
                separation=aa.separation(get_body('moon',times,LOCATION).transform_to(frame)).deg
                self.assertGreaterEqual(float(separation.min()),40,night)


if __name__=='__main__':unittest.main()
