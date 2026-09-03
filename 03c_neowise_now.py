"""
03c_neowise_now.py  --  Stage 3d extension: NEOWISE-R single-exposure W1/W2 photometry (IRSA, 2013-12 .. 2024-02)
so the mid-IR record continues past the unWISE time-domain catalog (which stops in 2020-12).

Usage: /opt/anaconda3/bin/python 03c_neowise_now.py <in.csv> <tag> [ra_col dec_col name_col]
Bulk IRSA Gator upload (spatial=Upload, 3" radius), 100 positions per request, cached per batch in
data/neowise_cache/<tag>/.  Frames with qual_frame > 0, cc_flags starting '00' and finite w1mpro are kept and
grouped into ~6-month visits.

Outputs
  data/neowise_visits_<tag>.csv   name, mjd (visit median), w1, w1err, w2, n   (one row per visit)
  data/neowise_now_<tag>.csv      per object: n_frames, n_visits, mjd_first_neo, mjd_last_neo, w1_first, w1_last,
                                  dw1_neowise (last - first, mag, + = fading), w1_slope_2yr (mag/yr), w1_amp_visits,
                                  w1_flux_last_mjy (Vega ZP 309.54 Jy), w2_last
"""
import io, os, sys, time
import numpy as np
import pandas as pd
import requests
from astropy.table import Table
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, 'data')
GATOR = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query'
BATCH = 100
W1_ZP_JY = 309.54


def gator(pos, tries=3):
    up = Table({'ra': pos.ra.values, 'dec': pos.dec.values})
    buf = io.StringIO(); ascii.write(up, buf, format='ipac')
    for attempt in range(tries):
        try:
            r = requests.post(GATOR, data={'catalog': 'neowiser_p1bs_psd', 'spatial': 'Upload', 'uplist': 'pos.tbl', 'radius': 3,
                                           'radunits': 'arcsec', 'outfmt': 1,
                                           'selcols': 'ra,dec,mjd,w1mpro,w1sigmpro,w2mpro,w2sigmpro,qual_frame,cc_flags'},
                              files={'filename': ('pos.tbl', buf.getvalue())}, timeout=1800)
            r.raise_for_status()
            if 'stat="ERROR"' in r.text[:300]:
                raise RuntimeError(r.text[:200])
            t = ascii.read(r.text, format='ipac').to_pandas()
            return t
        except Exception as e:
            print(f'   gator attempt {attempt+1}: {str(e)[:120]}', flush=True); time.sleep(30 * (attempt + 1))
    raise RuntimeError('Gator failed')


def main():
    inp, tag = sys.argv[1], sys.argv[2]
    ra_col, dec_col, name_col = (sys.argv[3:6] if len(sys.argv) >= 6 else ('ra', 'dec', 'name'))
    t = pd.read_csv(inp).drop_duplicates(name_col)
    pos = pd.DataFrame({'name': t[name_col].astype(str), 'ra': t[ra_col].astype(float), 'dec': t[dec_col].astype(float)}).reset_index(drop=True)
    cache = os.path.join(DATA, 'neowise_cache', tag); os.makedirs(cache, exist_ok=True)
    t0 = time.time(); frames = []
    for i in range(0, len(pos), BATCH):
        f = os.path.join(cache, f'batch_{i:05d}.csv')
        if os.path.exists(f):
            frames.append(pd.read_csv(f)); continue
        sub = pos.iloc[i:i + BATCH]
        g = gator(sub)
        # assign every returned frame to the nearest uploaded position (within 3")
        if len(g):
            cs = SkyCoord(g.ra.values * u.deg, g.dec.values * u.deg); cp = SkyCoord(sub.ra.values * u.deg, sub.dec.values * u.deg)
            idx, sep, _ = cs.match_to_catalog_sky(cp)
            g = g[sep.arcsec < 3.0].copy(); g['name'] = sub.name.values[idx[sep.arcsec < 3.0]]
        g.to_csv(f, index=False); frames.append(g)
        print(f'[{time.time()-t0:5.0f}s] batch {i//BATCH+1}/{(len(pos)-1)//BATCH+1}: {len(g)} frames', flush=True)
    fr = pd.concat(frames, ignore_index=True)
    fr = fr[(fr.qual_frame > 0) & fr.cc_flags.astype(str).str[:2].eq('00') & np.isfinite(fr.w1mpro) & np.isfinite(fr.w1sigmpro)]
    fr['visit'] = np.round(fr.mjd / 180.0)
    v = fr.groupby(['name', 'visit']).agg(mjd=('mjd', 'median'), w1=('w1mpro', 'median'), w1err=('w1mpro', lambda x: 1.2533 * x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan),
                                          w2=('w2mpro', 'median'), n=('w1mpro', 'size')).reset_index()
    v = v[v.n >= 3]
    v.to_csv(os.path.join(DATA, f'neowise_visits_{tag}.csv'), index=False)
    rows = []
    for name, g in v.groupby('name'):
        g = g.sort_values('mjd')
        first = g[g.mjd < g.mjd.min() + 400].w1.median(); last = g[g.mjd > g.mjd.max() - 400].w1.median()
        rr = g[g.mjd > g.mjd.max() - 730]
        slope = np.polyfit((rr.mjd - rr.mjd.min()) / 365.25, rr.w1, 1)[0] if len(rr) >= 3 and rr.mjd.max() - rr.mjd.min() > 300 else np.nan
        rows.append(dict(name=name, n_frames=int(fr[fr.name == name].shape[0]), n_visits=len(g), mjd_first_neo=g.mjd.min(), mjd_last_neo=g.mjd.max(),
                         w1_first=first, w1_last=last, dw1_neowise=last - first, w1_slope_2yr=slope, w1_amp_visits=g.w1.max() - g.w1.min(),
                         w1_flux_last_mjy=W1_ZP_JY * 1e3 * 10 ** (-0.4 * last), w2_last=g[g.mjd > g.mjd.max() - 400].w2.median()))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(DATA, f'neowise_now_{tag}.csv'), index=False)
    print(f'wrote data/neowise_now_{tag}.csv ({len(out)} of {len(pos)} objects) and data/neowise_visits_{tag}.csv ({len(v)} visits); '
          f'last visit median MJD {out.mjd_last_neo.median():.0f}; faded >0.3 mag since 2014: {(out.dw1_neowise > 0.3).sum()}, brightened >0.3: {(out.dw1_neowise < -0.3).sum()}')


if __name__ == '__main__':
    main()
