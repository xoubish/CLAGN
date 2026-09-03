"""
03b_ztf_now.py  --  Stage 3b/3d of FOLLOWUP_PLAN.md: current ZTF photometry and change since a reference epoch.

Usage: /opt/anaconda3/bin/python 03b_ztf_now.py <in.csv> <tag> [ra_col dec_col id_col mjd_ref_col]
  Queries the IRSA ZTF light-curve service (PSF photometry, all DRs to date) within 3 arcsec of each position,
  g and r bands, catflags < 32768.  Per-object cache in data/ztf_cache/<tag>/<id>.csv (safe to re-run).
Output data/ztf_now_<tag>.csv with, per object:
  n_ztf_g, n_ztf_r, mjd_first_ztf, mjd_last_ztf, r_first (median r over first 365 d), r_last (median r over last 365 d),
  g_last, dr_last_minus_first, r_at_ref (median r within +-180 d of mjd_ref, if given), dr_since_ref (= r_last - r_at_ref),
  slope_r_2yr (mag/yr over the last 2 yr; positive = fading), r_amp (95th-5th percentile of r).
"""
import io, os, sys, time
import numpy as np
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
URL = 'https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves'


def fetch(ra, dec, path, tries=3):
    if os.path.exists(path):
        return pd.read_csv(path) if os.path.getsize(path) > 5 else pd.DataFrame()
    for attempt in range(tries):
        try:
            r = requests.get(URL, params={'POS': f'CIRCLE {ra:.6f} {dec:.6f} 0.00083', 'BANDNAME': 'g,r', 'FORMAT': 'csv',
                                          'BAD_CATFLAGS_MASK': 32768}, timeout=600)
            r.raise_for_status()
            df = pd.read_csv(io.StringIO(r.text)) if r.text.strip() and not r.text.lstrip().startswith('<') else pd.DataFrame()
            (df if len(df) else pd.DataFrame(columns=['mjd'])).to_csv(path, index=False)
            return df
        except Exception as e:
            print(f'   {os.path.basename(path)}: {str(e)[:80]} (attempt {attempt+1})', flush=True); time.sleep(10 * (attempt + 1))
    return pd.DataFrame()


def summarise(df, mjd_ref=np.nan):
    out = dict(n_ztf_g=0, n_ztf_r=0, mjd_first_ztf=np.nan, mjd_last_ztf=np.nan, r_first=np.nan, r_last=np.nan, g_last=np.nan,
               dr_last_minus_first=np.nan, r_at_ref=np.nan, dr_since_ref=np.nan, slope_r_2yr=np.nan, r_amp=np.nan)
    if not len(df) or 'mjd' not in df or df.mjd.isna().all():
        return out
    df = df[np.isfinite(df.mag)]
    g = df[df.filtercode == 'zg']; r = df[df.filtercode == 'zr']
    out['n_ztf_g'], out['n_ztf_r'] = len(g), len(r)
    out['mjd_first_ztf'], out['mjd_last_ztf'] = df.mjd.min(), df.mjd.max()
    if len(r) >= 5:
        out['r_first'] = r[r.mjd < r.mjd.min() + 365].mag.median()
        out['r_last'] = r[r.mjd > r.mjd.max() - 365].mag.median()
        out['dr_last_minus_first'] = out['r_last'] - out['r_first']
        out['r_amp'] = np.percentile(r.mag, 95) - np.percentile(r.mag, 5)
        rr = r[r.mjd > r.mjd.max() - 730]
        if len(rr) >= 5 and rr.mjd.max() - rr.mjd.min() > 200:
            out['slope_r_2yr'] = np.polyfit((rr.mjd - rr.mjd.min()) / 365.25, rr.mag, 1)[0]
        if np.isfinite(mjd_ref):
            near = r[np.abs(r.mjd - mjd_ref) < 180]
            if len(near) >= 3:
                out['r_at_ref'] = near.mag.median(); out['dr_since_ref'] = out['r_last'] - out['r_at_ref']
    if len(g) >= 5:
        out['g_last'] = g[g.mjd > g.mjd.max() - 365].mag.median()
    return out


def main():
    inp, tag = sys.argv[1], sys.argv[2]
    ra_col, dec_col, id_col = (sys.argv[3:6] if len(sys.argv) >= 6 else ('ra', 'dec', 'name'))
    mjd_ref_col = sys.argv[6] if len(sys.argv) >= 7 else None
    t = pd.read_csv(inp)
    cache = os.path.join(DATA, 'ztf_cache', tag); os.makedirs(cache, exist_ok=True)
    t0 = time.time()
    def one(row):
        path = os.path.join(cache, f'{str(row[id_col]).replace("/", "_")}.csv')
        df = fetch(float(row[ra_col]), float(row[dec_col]), path)
        ref = float(row[mjd_ref_col]) if mjd_ref_col and mjd_ref_col in row and pd.notna(row[mjd_ref_col]) else np.nan
        return {id_col: row[id_col], **summarise(df, ref)}
    rows = []
    with ThreadPoolExecutor(max_workers=4) as ex:
        for k, res in enumerate(ex.map(one, [r for _, r in t.iterrows()])):
            rows.append(res)
            if (k + 1) % 25 == 0:
                print(f'[{time.time()-t0:5.0f}s] {k+1}/{len(t)}', flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(DATA, f'ztf_now_{tag}.csv'), index=False)
    print(f'wrote data/ztf_now_{tag}.csv ({len(out)} rows); with r light curve: {(out.n_ztf_r >= 5).sum()}; '
          f'median mjd_last {out.mjd_last_ztf.median():.0f}; faded >0.5 mag since first year: {(out.dr_last_minus_first > 0.5).sum()}; '
          f'brightened >0.5: {(out.dr_last_minus_first < -0.5).sum()}')


if __name__ == '__main__':
    main()
