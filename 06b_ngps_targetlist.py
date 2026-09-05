"""
06b_ngps_targetlist.py  --  write the per-night target lists in the NGPS observing-software CSV format.

Format per https://caltechopticalobservatories.github.io/NGPS/users-manual/target-lists.html (read 2026-09-05):
  name, RA (HH:MM:SS.S, J2000), DECL (+DD:MM:SS), slitwidth ("SET X" arcsec | "PSF X" | "SNR X"), exptime ("SET s" | "SNR X"),
  nexp, binspect, binspat, slitangle (deg | "PA"), airmass_max, mag, magsystem, magfilter, channel, wrange (Å, "a:b"),
  Note (<= 24 chars), Comment (<= 1024 chars). Unknown headers are ignored; a header row is required.
Two files per night:
  finders/<night>/ngps_<night>_fixed.csv  exposure times from our model (SET seconds per sub-exposure, nexp >= 2, <= 900 s)
  finders/<night>/ngps_<night>_snr.csv    exptime "SNR 7" so the sequencer/ETC solves the time from the target's magnitude,
                                          the channel and a wavelength window on the diagnostic line (Hα if z <= 0.55, else Hβ)
Slit 1.3" fixed (the software default; single slice is what the ETC currently models), binning 2x3 (BINSPAT x BINSPEC) per the
observing page's guidance for a ~1.3" slit, slit angle PA (parallactic, no ADC), airmass_max 2.0.
Backups are included after the primaries with Note 'backup'.
"""
import os, sys, glob
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, 'data')
CHANNELS = [('U', 3050, 4430), ('G', 4250, 5960), ('R', 5620, 7950), ('I', 7530, 10400)]
SLIT, BINSPAT, BINSPEC, AIRMASS_MAX, SNR_TARGET = 1.3, 2, 3, 2.0, 7


def channel_for(wave):
    """channel whose bandpass contains `wave`, preferring the one where it sits furthest from an edge."""
    best, margin = None, -1
    for c, a, b in CHANNELS:
        if a <= wave <= b:
            m = min(wave - a, b - wave)
            if m > margin:
                best, margin = c, m
    return best


def diagnostic(z):
    """observed wavelength of the line the exposure is planned on: Hα if it is in range, else Hβ, else Mg II."""
    for name, w0 in [('Ha', 6563.0), ('Hb', 4861.0), ('MgII', 2798.0)]:
        w = w0 * (1 + z)
        if 3100 <= w <= 10300:
            return name, w
    return 'Hb', 4861.0 * (1 + z)


def rows_for(night):
    t = pd.read_csv(os.path.join(DATA, f'targets_{night}.csv'))
    m = pd.read_csv(os.path.join(DATA, 'master_list_scored.csv'), low_memory=False).set_index('name')
    out_fixed, out_snr = [], []
    for r in t.itertuples():
        c = SkyCoord(r.ra * u.deg, r.dec * u.deg)
        ra = c.ra.to_string(u.hour, sep=':', precision=1, pad=True); dec = c.dec.to_string(u.deg, sep=':', precision=0, alwayssign=True, pad=True)
        z = float(r.z) if pd.notna(r.z) else 0.3
        line, w = diagnostic(z); ch = channel_for(w) or 'R'
        wrange = f'{int(w - 40)}:{int(w + 40)}'
        rmag = float(r.r_mag) if pd.notna(r.r_mag) else 19.0
        texp = float(r.t_exp_min) if pd.notna(getattr(r, 't_exp_min', np.nan)) else 10.0
        nexp = int(max(2, np.ceil(texp / 10.0))); per = min(900, int(round(texp * 60 / nexp / 10) * 10))
        status = 'primary' if r.rank > 0 else 'backup'
        note = (f'{r.tier} #{int(r.rank)} p{r.priority_night:.1f}' if r.rank > 0 else f'{r.tier} backup')[:24]
        mm = m.loc[r.name] if r.name in m.index else None
        comment = (f'{status}; {r.tier}; z={z:.3f}; r={rmag:.1f}; line {line} at {w:.0f} A ({ch}); trend {getattr(r, "trend", "")}; '
                   f'last spectrum {getattr(r, "years_since_last_spec", np.nan):.1f} yr ago; exposure model {texp:.0f} min; '
                   f'{str(getattr(r, "notes", ""))[:400]}')[:1024]
        base = dict(name=r.name, RA=ra, DECL=dec, slitwidth=f'SET {SLIT}', nexp=nexp, binspect=BINSPEC, binspat=BINSPAT, slitangle='PA',
                    airmass_max=AIRMASS_MAX, mag=round(rmag, 2), magsystem='AB', magfilter='r', channel=ch, wrange=wrange, Note=note, Comment=comment)
        out_fixed.append({**base, 'exptime': f'SET {per}'})
        out_snr.append({**base, 'exptime': f'SNR {SNR_TARGET}'})
    cols = ['name', 'RA', 'DECL', 'slitwidth', 'exptime', 'nexp', 'binspect', 'binspat', 'slitangle', 'airmass_max', 'mag', 'magsystem', 'magfilter', 'channel', 'wrange', 'Note', 'Comment']
    return pd.DataFrame(out_fixed)[cols], pd.DataFrame(out_snr)[cols]


if __name__ == '__main__':
    nights = sys.argv[1:] or ['sep23', 'oct26', 'oct27']
    for n in nights:
        fx, sn = rows_for(n)
        d = os.path.join(HERE, 'finders', n); os.makedirs(d, exist_ok=True)
        fx.to_csv(os.path.join(d, f'ngps_{n}_fixed.csv'), index=False); sn.to_csv(os.path.join(d, f'ngps_{n}_snr.csv'), index=False)
        print(f'{n}: {len(fx)} rows -> finders/{n}/ngps_{n}_fixed.csv and ngps_{n}_snr.csv ({int((fx.Note.str.contains("backup")).sum())} backups)')
    print(fx.head(3).to_string(index=False))
