"""
14b_turnon_arm.py  --  turn-on discovery arm from SDSS narrow-line AGN galaxies (Yang+2025 / Zhu+2026 recipe on our data).

Parent: data/turnon_parent_galaxies.csv (SDSS GALAXY spectra, subclass AGN, z 0.02-0.35, r < 19.5, in the run windows; 11,608).
Triggers
  MIR   NEOWISE W1 AND W2 brightened by >= 0.2 mag between the first and last 500 d of 2014-2024 (visit medians), each at >= 3 sigma
  OPT   ZTF DR24 median g brighter than the synthetic g of the archival spectrum by > 0.3 mag after removing the population offset
        (ZTF PSF photometry of an extended galaxy captures less host light than the 3" fibre; the offset is taken out per petroRad bin)
  COL   W1-W2 reddened by > 0.1 (the nucleus emerging)
Score = calibrated success rate of the matching literature criterion: MIR+OPT 0.70 (Yang+2025: 82/115), MIR only 0.45, OPT only 0.30,
x1.2 with COL (cap 0.85).  Known CLAGNs (Camus & Panda 2026) are flagged and left out of the discovery arm.
Observability from 05_observability.py for the triggered galaxies only.  Output: data/turnon_candidates.csv (read by 14_select_targets.py).
"""
import os, sys, glob, subprocess, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from astropy.coordinates import SkyCoord; import astropy.units as u
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, 'data'); PY = sys.executable


def neowise_summary():
    """per-galaxy NEOWISE W1/W2 change 2014 -> 2024 from 03c's output or, while it is still running, from its cached batches."""
    f = os.path.join(DATA, 'neowise_visits_turnon.csv')
    if os.path.exists(f):
        v = pd.read_csv(f)
    else:
        fr = pd.concat([pd.read_csv(p) for p in sorted(glob.glob(os.path.join(DATA, 'neowise_cache', 'turnon', 'batch_*.csv'))) if os.path.getsize(p) > 100], ignore_index=True)
        fr = fr[(fr.qual_frame > 0) & fr.cc_flags.astype(str).str[:2].eq('00') & np.isfinite(fr.w1mpro) & np.isfinite(fr.w1sigmpro)]
        fr['visit'] = np.round(fr.mjd / 180.0)
        v = fr.groupby(['name', 'visit']).agg(mjd=('mjd', 'median'), w1=('w1mpro', 'median'), w1err=('w1mpro', lambda x: 1.2533 * x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan),
                                              w2=('w2mpro', 'median'), w2err=('w2mpro', lambda x: 1.2533 * x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan), n=('w1mpro', 'size')).reset_index()
        v = v[v.n >= 3]
    if 'w2err' not in v: v['w2err'] = v.w1err
    v = v.sort_values(['name', 'mjd']); g = v.groupby('name'); mn, mx = g.mjd.transform('min'), g.mjd.transform('max')
    first = v[v.mjd < mn + 500].groupby('name').agg(w1_first=('w1', 'median'), w2_first=('w2', 'median'), e1_first=('w1err', 'median'), e2_first=('w2err', 'median'), n_first=('w1', 'size'))
    last = v[v.mjd > mx - 500].groupby('name').agg(w1_last=('w1', 'median'), w2_last=('w2', 'median'), e1_last=('w1err', 'median'), e2_last=('w2err', 'median'), n_last=('w1', 'size'), mjd_last_neo=('mjd', 'max'))
    s = first.join(last); s['dw1'] = s.w1_last - s.w1_first; s['dw2'] = s.w2_last - s.w2_first
    s['sig1'] = -s.dw1 / np.sqrt(s.e1_first ** 2 + s.e1_last ** 2).replace(0, np.nan); s['sig2'] = -s.dw2 / np.sqrt(s.e2_first ** 2 + s.e2_last ** 2).replace(0, np.nan)
    s['dcolor'] = (s.w1_last - s.w2_last) - (s.w1_first - s.w2_first); s['w1w2_now'] = s.w1_last - s.w2_last; s['n_visits'] = g.size()
    return s


def main():
    gal = pd.read_csv(os.path.join(DATA, 'turnon_parent_galaxies.csv')).set_index('name')
    zo = pd.read_csv(os.path.join(DATA, 'ztf_objects_turnon.csv')).set_index('name'); syn = pd.read_csv(os.path.join(DATA, 'synthetic_ztf_turnon.csv')).set_index('name')
    neo = neowise_summary()
    x = gal.join(zo.drop(columns=['ra', 'dec'], errors='ignore')).join(syn[['syn_g', 'syn_r', 'snr_g', 'snr_r']]).join(neo)
    print(f'parent {len(x)}; with ZTF g {x.ztf_g_medianmag.notna().sum()}, synthetic g {x.syn_g.notna().sum()}, NEOWISE {x.dw1.notna().sum()}')
    for b in ('g', 'r'):
        bad = (x[f'ztf_{b}_medianmag'] < 10) | (x[f'ztf_{b}_medianmag'] > 23) | (x[f'ztf_{b}_ngoodobs'] < 20); x.loc[bad, [f'ztf_{b}_medianmag', f'ztf_{b}_amp']] = np.nan
        x[f'd_{b}'] = x[f'ztf_{b}_medianmag'] - x[f'syn_{b}']
    x.loc[(x.snr_g < 3) | (x.snr_r < 3), ['d_g', 'd_r']] = np.nan
    # population offset of ZTF (PSF) minus synthetic (3" fibre) as a function of galaxy size, removed so that d_*_corr = 0 means "no change"
    rb = pd.qcut(x.petroRad_r.clip(0.5, 20), 8, duplicates='drop')
    for b in ('g', 'r'):
        x[f'd_{b}_corr'] = x[f'd_{b}'] - x.groupby(rb)[f'd_{b}'].transform('median')
    x['d_opt'] = np.where(x.d_g_corr.notna() & x.d_r_corr.notna() & ((x.d_g_corr - x.d_r_corr).abs() <= 0.8), 0.5 * (x.d_g_corr + x.d_r_corr), np.where(x.d_g_corr.notna(), x.d_g_corr, x.d_r_corr))
    x['mir'] = (x.dw1 <= -0.2) & (x.dw2 <= -0.2) & (x.sig1 >= 3) & (x.sig2 >= 3)
    x['opt'] = (x.d_opt < -0.3) & (x.ztf_g_amp > 0.3)
    x['col'] = x.dcolor > 0.1
    # the MIR trigger is mandatory (Yang+2025 is MIR-first; ZTF PSF photometry of extended hosts makes an optical-only trigger unreliable):
    # MIR + optical 0.70 (82/115 confirmed), MIR + some ZTF variability 0.45, MIR alone 0.30; colour reddening x1.2 (cap 0.85)
    x['score'] = np.select([x.mir & x.opt, x.mir & (x.ztf_g_amp > 0.3), x.mir], [0.70, 0.45, 0.30], 0.0); x.loc[x.col & (x.score > 0), 'score'] = (x.score[x.col & (x.score > 0)] * 1.2).clip(upper=0.85)
    cat = pd.read_csv(os.path.join(DATA, 'external', 'clagn_catalog_camus_panda2026.csv'), low_memory=False); cat = cat[np.isfinite(cat.ra_deg) & np.isfinite(cat.dec_deg)]
    i, s, _ = SkyCoord(x.ra.values * u.deg, x.dec.values * u.deg).match_to_catalog_sky(SkyCoord(cat.ra_deg.values * u.deg, cat.dec_deg.values * u.deg)); x['known_clagn'] = s.arcsec < 2
    print(f'triggers: MIR {int(x.mir.sum())}, OPT {int(x.opt.sum())}, both {int((x.mir & x.opt).sum())}, with colour reddening {int((x.col & (x.score > 0)).sum())}; known CLAGNs among triggered: {int((x.known_clagn & (x.score > 0)).sum())}')
    cand = x[(x.score > 0) & ~x.known_clagn].copy()
    cand['r_now'] = cand.ztf_r_medianmag.where(cand.ztf_r_medianmag.notna(), cand.modelMag_r)
    cand['why'] = ['TURN-ON candidate (galaxy parent, ' + str(r.subclass) + f', z={r.z:.3f}): ' + ('; '.join(p for p in [
        f'W1 brightened {-r.dw1:.2f} mag and W2 {-r.dw2:.2f} mag 2014-24 ({r.sig1:.0f}/{r.sig2:.0f} sigma)' if r.mir else (f'W1 change {-r.dw1:+.2f} mag' if pd.notna(r.dw1) else 'no NEOWISE yet'),
        f'ZTF g brighter than the archival spectrum by {-r.d_opt:.2f} mag (host-offset corrected), ZTF amplitude {r.ztf_g_amp:.2f}' if r.opt else (f'optical change {r.d_opt:+.2f} mag' if pd.notna(r.d_opt) else 'no optical measure'),
        f'W1-W2 reddened by {r.dcolor:+.2f} (nucleus emerging)' if r.col else ''] if p)) + f'; score {r.score:.2f} = literature success rate of this trigger set' for r in cand.itertuples()]
    pos = os.path.join(DATA, 'turnon_candidates_pos.csv'); cand.reset_index()[['name', 'ra', 'dec']].to_csv(pos, index=False)
    obs = os.path.join(DATA, 'turnon_candidates_obs.csv')
    subprocess.run([PY, os.path.join(HERE, '05_observability.py'), pos, obs], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    o = pd.read_csv(obs).set_index('name'); ocols = [c for c in o.columns if c.startswith(('hrs_', 'minX_', 'moonsep_'))]
    cand = cand.join(o[ocols])
    out = cand.reset_index()[['name', 'ra', 'dec', 'z', 'r_now', 'modelMag_r', 'petroRad_r', 'subclass', 'score', 'mir', 'opt', 'col', 'dw1', 'dw2', 'sig1', 'sig2', 'dcolor', 'w1w2_now', 'd_g_corr', 'd_r_corr', 'd_opt', 'ztf_g_amp', 'known_clagn', 'why'] + ocols]
    out.to_csv(os.path.join(DATA, 'turnon_candidates.csv'), index=False)
    print(f'{len(out)} turn-on candidates written; observable >= 1.5 h on some night: {int((out[[c for c in ocols if c.startswith("hrs_")]].max(axis=1) >= 1.5).sum())}; score >= 0.7: {int((out.score >= 0.7).sum())}')
    print(out.sort_values('score', ascending=False).head(15)[['name', 'z', 'r_now', 'score', 'dw1', 'dw2', 'dcolor', 'd_opt', 'ztf_g_amp', 'subclass']].round(2).to_string(index=False))


if __name__ == '__main__':
    main()
