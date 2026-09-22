"""Slit previews in the ICRS north-up/east-left SDSS cutout frame.

These are planning illustrations, not slit-angle commands. NGPS PA requests
are evaluated for the actual exposure by the observing software.
"""
import importlib

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, SkyCoord
from astropy.time import Time
from astropy.utils import iers

iers.conf.auto_download = False
iers.conf.auto_max_age = None
PALOMAR = importlib.import_module('05_observability').PALOMAR


def slit_angles(ra, dec, times):
    """Position angle east of ICRS north, modulo 180 degrees; sec(z).

    Transform a short step toward the zenith back into the image frame rather
    than mixing apparent sidereal time with catalogue J2000 RA. Refraction is
    off, as in the observability calculations. PA is undefined at the zenith.
    """
    target = SkyCoord(ra*u.deg, dec*u.deg)
    frame = AltAz(obstime=times, location=PALOMAR.location, pressure=0*u.hPa)
    horizontal = target.transform_to(frame)
    altitude = horizontal.alt.to_value(u.deg)
    step = np.minimum(1/3600, (90-altitude)/2)
    toward_zenith = SkyCoord(az=horizontal.az, alt=(altitude+step)*u.deg,
                            frame=frame).icrs
    angles = target.position_angle(toward_zenith).to_value(u.deg) % 180
    return np.where(altitude > 89.9999, np.nan, angles), horizontal.secz.value


def previews(target, sequence):
    """Five-minute observable samples plus exact visit start/midpoint/end."""
    result = {}
    for night in target.get('nights', []):
        key = night['night']
        # Reuse the page's Moon/airmass screen; do not assign unscheduled visits.
        candidates = [(m, x) for m, x, moon in night.get('curve', [])
                      if x is not None and 1 <= x <= 2 and moon >= 40]
        if not candidates:
            continue
        preferred = [(m, x) for m, x in candidates if x <= 1.8] or candidates
        default = min(preferred, key=lambda row: row[1])[0]
        labels = {}
        slot = sequence.get(target['name']) if key == 'sep23' else None
        mjds = [m for m, _ in candidates]
        if slot:
            start, end = Time([slot['start_utc'], slot['end_utc']]).mjd
            default = (start+end)/2
            labels = {float(start): 'visit start', float(default): 'visit midpoint',
                      float(end): 'visit end'}
            mjds.extend(labels)
        mjds = sorted(set(mjds))
        times = Time(mjds, format='mjd', scale='utc')
        angles, airmasses = slit_angles(target['ra'], target['dec'], times)
        samples = [[float(m), round(float(pa), 3) if np.isfinite(pa) else None,
                    round(float(x), 3), labels.get(m, '')]
                   for m, pa, x in zip(mjds, angles, airmasses)]
        result[key] = dict(samples=samples, default_index=mjds.index(default),
                           scheduled=bool(slot))
    return result
