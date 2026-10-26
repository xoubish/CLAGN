"""Full-visit start windows for the three observing queues at Palomar."""
import csv
import io
from datetime import datetime

import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body
from astropy.time import Time
from astropy.utils import iers

LOCATION = EarthLocation.from_geodetic(-116.865*u.deg, 33.3563*u.deg, 1712*u.m)
STEP = 30000


def millis(value):
    return round(datetime.fromisoformat(value.replace('Z', '+00:00')).timestamp()*1000)


def attach_visibility(payload):
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    targets = {t['name']: t for t in payload['targets']}
    result = {}
    for night in ['sep23', 'oct26', 'oct27']:
        if night == 'sep23':
            meta = payload['nights'][night]
            start, end = [(meta[k]-40587)*86400000 for k in ['start_mjd', 'end_mjd']]
            queue = [dict(r, name=name, order=r['rank'], ra=targets[name]['ra'], dec=targets[name]['dec'])
                     for name, r in sorted(payload['sep23_sequence'].items(), key=lambda p:p[1]['rank'])]
            standards = [dict(name=r['name'], start=millis(r['start_utc']), end=millis(r['end_utc']))
                         for r in payload['sequence_rows'] if r['role']=='standard']
        else:
            meta = payload['october']['nights'][night]
            start, end = millis(meta['dusk_pdt']), millis(meta['dawn_pdt'])
            limits = {r['NAME']:float(r['AIRMASS_MAX']) for r in csv.DictReader(io.StringIO(payload['files'][night+'_ngps.csv']))}
            queue = [dict(r, airmass_limit=limits[r['name']]) for r in payload['october']['rows'] if r['night']==night]
            standards = [dict(name=r['name'], start=millis(r['start_pdt']), end=millis(r['end_pdt']))
                         for r in payload['october']['standards'] if r['night']==night]
            standards += [dict(name=r['name'], start=millis(r['start_pdt']), end=millis(r['end_pdt']), kind='manual visit')
                          for r in payload['october'].get('manual_visits',[]) if r['night']==night]
            queue += [dict(r, manual=True, airmass_limit=1.8)
                      for r in payload['october'].get('manual_visits',[]) if r['night']==night]
            queue.sort(key=lambda r:r['start_pdt'])
            for order,r in enumerate(queue,1):r['order']=order
        grid = np.arange(np.ceil(start/STEP)*STEP, np.floor(end/STEP)*STEP+1, STEP)
        times = Time(grid/1000, format='unix')
        frame = AltAz(obstime=times, location=LOCATION, pressure=0*u.hPa)
        coords = SkyCoord([float(r['ra']) for r in queue]*u.deg, [float(r['dec']) for r in queue]*u.deg)
        aa = coords[:,None].transform_to(frame)
        airmass = np.where(aa.alt.deg>0, aa.secz.value, np.inf)
        moon = aa.separation(get_body('moon', times, LOCATION).transform_to(frame)).deg
        rows = []
        for i, r in enumerate(queue):
            duration = float(r['visit_minutes'])
            count = int(round(duration*60000/STEP))+1
            xmax = np.lib.stride_tricks.sliding_window_view(airmass[i], count).max(axis=1)
            mmin = np.lib.stride_tricks.sliding_window_view(moon[i], count).min(axis=1)
            starts = grid[:len(xmax)]
            free = np.ones(len(starts), dtype=bool)
            for s in standards:
                if r.get('manual') and s['name']==r['name']:continue
                free &= ~((starts < s['end']) & (starts+duration*60000 > s['start']))
            limit = float(r['airmass_limit'])
            # Small margins protect boundaries between samples; evaluate the entire
            # acquisition + exposure block, not just the target at its start.
            allowed = free & (mmin >= 40.02) & (xmax <= limit-.001)
            tiers = np.where(allowed, np.where(xmax<=1.499, 1, np.where(xmax<=1.799, 2, 3)), 0)
            windows = []
            boundaries = np.r_[0, np.where(np.diff(tiers)!=0)[0]+1, len(tiers)]
            for a,b in zip(boundaries[:-1], boundaries[1:]):
                if tiers[a]:
                    windows.append(dict(start=int(starts[a]), end=int(starts[b-1]), tier=int(tiers[a])))
            samples = [[round(float(x),3) if np.isfinite(x) and x<100 else None,
                        round(float(m),2), int(t)] for x,m,t in zip(xmax,mmin,tiers)]
            rows.append(dict(name=r['name'], order=int(r['order']), duration=duration, airmass_limit=limit,
                manual=bool(r.get('manual')),
                planned_start=millis(r['start_utc']), planned_end=millis(r.get('end_utc') or r['end_pdt']),
                windows=windows, samples=samples))
        result[night] = dict(start=round(start), end=round(end), sample_start=int(grid[0]), step=STEP,
                             count=len(grid), standards=standards, rows=rows)
    payload['visibility'] = result
