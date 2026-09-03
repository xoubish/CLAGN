"""
05_observability.py  --  Stage 4 of FOLLOWUP_PLAN.md: Palomar/NGPS observability for the 2026B runs.

Runs (from S. Hemmati, 2026-09-03):  2026-09-23 first half of the night;  2026-10-26 and 2026-10-27 full nights.
All bright time -> prefer bright targets and large moon separation.

Usage:  /opt/anaconda3/bin/python 05_observability.py <in.csv> <out.csv> [ra_col dec_col]
Adds per-night columns:  hrs_<night> (hours inside the window with airmass < AIRMASS_MAX),
                          minX_<night> (best airmass in window), moonsep_<night> (deg, at window midpoint),
                          plus hrs_any (max over nights) and moonsep_min.
Also prints the night summary (twilights, window, moon illumination/position, RA range that is observable).
"""
import sys, os
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body
from astroplan import Observer, FixedTarget, moon_illumination

PALOMAR = Observer(location=EarthLocation.from_geodetic(-116.8650 * u.deg, 33.3563 * u.deg, 1712 * u.m),
                   timezone='US/Pacific', name='Palomar')
AIRMASS_MAX = 2.0
# night label -> (local calendar date of the evening, fraction of the night to use: 'first' or 'full')
NIGHTS = {'sep23': ('2026-09-23', 'first'), 'oct26': ('2026-10-26', 'full'), 'oct27': ('2026-10-27', 'full')}


def night_window(date_str, part):
    """Return (t_start, t_end) as astropy Times in UTC for the observing window."""
    # local noon of that date, then the following evening twilight
    noon_local = pd.Timestamp(f'{date_str} 12:00', tz='US/Pacific')
    t_noon = Time(noon_local.tz_convert('UTC').to_pydatetime())
    dusk = PALOMAR.twilight_evening_astronomical(t_noon, which='next')
    dawn = PALOMAR.twilight_morning_astronomical(dusk, which='next')
    mid = dusk + (dawn - dusk) / 2
    if part == 'first':
        return dusk, mid, dusk, dawn
    return dusk, dawn, dusk, dawn


def describe_night(label, t0, t1, dusk, dawn):
    times = t0 + np.linspace(0, (t1 - t0).to(u.hour).value, 25) * u.hour
    moon_mid = get_body('moon', t0 + (t1 - t0) / 2, PALOMAR.location)
    illum = moon_illumination(t0 + (t1 - t0) / 2)
    moon_alt = PALOMAR.altaz(times, moon_mid).alt.deg
    lst0 = PALOMAR.local_sidereal_time(t0).hour; lst1 = PALOMAR.local_sidereal_time(t1).hour
    to_local = lambda t: pd.Timestamp(t.utc.datetime, tz='UTC').tz_convert('US/Pacific').strftime('%H:%M')
    print(f'\n=== {label}: {t0.utc.datetime:%Y-%m-%d} UTC ===')
    print(f'  astronomical twilight: {to_local(dusk)} -> {to_local(dawn)} local; window used {to_local(t0)} -> {to_local(t1)} '
          f'({(t1 - t0).to(u.hour).value:.1f} h)')
    print(f'  LST at window start/end: {lst0:.1f}h / {lst1:.1f}h  -> meridian RA range {lst0:.1f}h..{lst1:.1f}h; '
          f'airmass<2 roughly RA {lst0-3.5:.1f}h .. {lst1+3.5:.1f}h (Dec ~ +30)')
    print(f'  Moon: illumination {100*illum:.0f}%, RA {moon_mid.ra.hour:.1f}h Dec {moon_mid.dec.deg:+.0f} deg, '
          f'altitude in window {moon_alt.min():.0f}..{moon_alt.max():.0f} deg (above horizon {100*(moon_alt>0).mean():.0f}% of window)')
    return moon_mid, illum


def score(df, ra_col='ra', dec_col='dec', step_min=10):
    coords = SkyCoord(df[ra_col].values * u.deg, df[dec_col].values * u.deg)
    out = df.copy()
    for label, (date_str, part) in NIGHTS.items():
        t0, t1, dusk, dawn = night_window(date_str, part)
        moon_mid, illum = describe_night(label, t0, t1, dusk, dawn)
        n = int((t1 - t0).to(u.min).value // step_min) + 1
        times = t0 + np.arange(n) * step_min * u.min
        # airmass grid: (ntargets, ntimes)
        altaz = PALOMAR.altaz(times[np.newaxis, :], coords[:, np.newaxis])
        secz = altaz.secz.value
        good = (altaz.alt.deg > 0) & (secz < AIRMASS_MAX) & (secz > 0)
        out[f'hrs_{label}'] = good.sum(axis=1) * step_min / 60.0
        best = np.where(good, secz, np.inf).min(axis=1)
        out[f'minX_{label}'] = np.where(np.isfinite(best), best, np.nan).round(2)
        out[f'moonsep_{label}'] = coords.separation(moon_mid).deg.round(0)
    hrs_cols = [c for c in out.columns if c.startswith('hrs_')]
    out['hrs_any'] = out[hrs_cols].max(axis=1)
    out['moonsep_min'] = out[[c for c in out.columns if c.startswith('moonsep_')]].min(axis=1)
    return out


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    inp, outp = sys.argv[1], sys.argv[2]
    ra_col, dec_col = (sys.argv[3:5] if len(sys.argv) >= 5 else ('ra', 'dec'))
    df = pd.read_csv(inp)
    res = score(df, ra_col, dec_col)
    res.to_csv(outp, index=False)
    print(f'\nwrote {outp}: {len(res)} rows; observable >=1h on some night: {(res.hrs_any >= 1).sum()}')
