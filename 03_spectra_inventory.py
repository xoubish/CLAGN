"""
03_spectra_inventory.py  --  Stage 3a of FOLLOWUP_PLAN.md: archival spectra per target.

Usage:  /opt/anaconda3/bin/python 03_spectra_inventory.py <input.csv> <tag> [ra_col dec_col name_col]
        input.csv needs one row per target with ra, dec (deg) and a name column.

Sources
  * SDSS DR19 `allspec` (every SDSS-I..V optical spectrum, incl. SDSS-V BOSS daily epochs to late 2022),
    with class/subclass/z from `SpecObjAll` (SDSS-I-IV) and `mos_sdssv_boss_spall` (SDSS-V; also gives
    spectroflux_g/r/i = synthetic flux at the spectral epoch, in nanomaggies).
    Queried through the SkyServer REST SQL endpoint in batches (UNION ALL of cone searches).
  * DESI DR1 `desi_dr1.zpix` via the NOIRLab Data Lab TAP service (healpix coadds; mean/min/max MJD).
  * LAMOST: not yet (TODO).

Outputs
  data/spectra_epochs_<tag>.csv    one row per (target, spectrum epoch)
  data/spectra_summary_<tag>.csv   one row per target: n_spec, mjd_first, mjd_last, years_since_last_spec,
                                   surveys, class_list, subclass_list, class_change_flag, broadline_change_flag,
                                   last_class, last_subclass, r_nmgy_last_sdssv
"""
import io, os, sys, time
import numpy as np
import pandas as pd
import requests
import pyvo
from astropy.time import Time

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
SKY = 'https://skyserver.sdss.org/dr19/SkyServerWS/SearchTools/SqlSearch'
TAP = pyvo.dal.TAPService('https://datalab.noirlab.edu/tap')
RADIUS_ARCMIN = 2.0 / 60.0          # 2 arcsec
BATCH = 20
TODAY_MJD = Time.now().mjd


def skyserver(sql, tries=4):
    for attempt in range(tries):
        try:
            r = requests.get(SKY, params={'cmd': sql, 'format': 'csv'}, timeout=300)
            r.raise_for_status()
            txt = r.text
            if txt.startswith('#Table1'):
                txt = txt.split('\n', 1)[1]
            if not txt.strip() or txt.lstrip().startswith('<'):
                raise RuntimeError(txt[:200])
            return pd.read_csv(io.StringIO(txt))
        except Exception as e:
            print(f'   skyserver attempt {attempt+1}: {str(e)[:120]}', flush=True)
            time.sleep(10 * (attempt + 1))
    raise RuntimeError('SkyServer query failed')


def sdss_allspec(targets):
    """All allspec rows within 2 arcsec of each target. targets: DataFrame with idx, ra, dec."""
    cols = ('allspec_id, sdss_phase, instrument, sdss_id, catalogid, fiberid, plate_or_fps_field, mjd, run2d, coadd, '
            'programname, survey, ra, dec, specobjid, sas_url')
    # SkyServer's REST endpoint rejects UNION ALL (multi-SELECT), so run one indexed cone query per target,
    # a few in parallel.
    from concurrent.futures import ThreadPoolExecutor
    def one(r):
        q = (f"SELECT {r.idx} AS idx, {cols} FROM allspec WHERE dbo.fDistanceArcMinEq({r.ra:.6f}, {r.dec:.6f}, ra, dec) < {RADIUS_ARCMIN:.5f}")
        try:
            return skyserver(q)
        except Exception as e:
            print(f'   target {r.idx} failed: {str(e)[:80]}', flush=True); return pd.DataFrame()
    out = []
    with ThreadPoolExecutor(max_workers=3) as ex:
        for k, df in enumerate(ex.map(one, list(targets.itertuples()))):
            out.append(df)
            if (k + 1) % 25 == 0:
                print(f'   allspec: {k+1}/{len(targets)} targets done', flush=True)
    out = [d for d in out if len(d)]
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def sdss_classes(ep):
    """Attach class/subclass/z/zwarning to allspec epochs."""
    ep = ep.copy()
    for c in ['class', 'subclass', 'firstcarton']:
        ep[c] = pd.Series([None] * len(ep), index=ep.index, dtype=object)
    for c in ['z', 'zwarning', 'sn_median_all', 'spectroflux_g', 'spectroflux_r', 'spectroflux_i']:
        ep[c] = np.nan
    # SDSS-I..IV: SpecObjAll by specobjid
    legacy = ep[(ep.sdss_phase < 5) & ep.specobjid.notna()]
    ids = [str(int(x)) for x in legacy.specobjid.unique()]
    for i in range(0, len(ids), 200):
        q = ("SELECT specObjID AS specobjid, class, subclass, z, zWarning AS zwarning, snMedian AS sn_median_all, "
             "spectroFlux_g AS spectroflux_g, spectroFlux_r AS spectroflux_r, spectroFlux_i AS spectroflux_i "
             f"FROM SpecObjAll WHERE specObjID IN ({','.join(ids[i:i+200])})")
        so = skyserver(q)
        so['specobjid'] = so['specobjid'].astype('int64')
        m = ep.specobjid.isin(so.specobjid)
        ep.loc[m, so.columns[1:]] = ep.loc[m, ['specobjid']].merge(so, on='specobjid', how='left')[so.columns[1:]].values
    # SDSS-V daily epochs: mos_sdssv_boss_spall by catalogid + mjd
    v = ep[(ep.sdss_phase == 5) & (ep.coadd != 'allepoch') & ep.catalogid.notna() & (ep.catalogid > 0)]
    if len(v):
        cats = [str(int(x)) for x in v.catalogid.unique()]
        rows = []
        for i in range(0, len(cats), 100):
            q = ("SELECT catalogid, mjd, field, class, subclass, z, zwarning, sn_median_all, spectroflux_g, spectroflux_r, spectroflux_i, firstcarton "
                 f"FROM mos_sdssv_boss_spall WHERE catalogid IN ({','.join(cats[i:i+100])})")
            rows.append(skyserver(q))
        sp = pd.concat(rows, ignore_index=True).drop_duplicates(['catalogid', 'mjd'])
        key = ['catalogid', 'mjd']
        ep['catalogid'] = ep['catalogid'].astype('Int64'); sp['catalogid'] = sp['catalogid'].astype('Int64')
        mrg = ep.loc[v.index, key].merge(sp, on=key, how='left')
        for c in ['class', 'subclass', 'z', 'zwarning', 'sn_median_all', 'spectroflux_g', 'spectroflux_r', 'spectroflux_i']:
            ep.loc[v.index, c] = mrg[c].values
        ep.loc[v.index, 'firstcarton'] = mrg['firstcarton'].values
    return ep


def desi_zpix(targets):
    out = []
    for i in range(0, len(targets), BATCH):
        sub = targets.iloc[i:i + BATCH]
        # Data Lab's ADQL layer does not translate CONTAINS/POINT for this table; use small boxes instead
        rad = RADIUS_ARCMIN / 60.0
        cond = ' OR '.join(f"(mean_fiber_ra BETWEEN {r.ra - rad/np.cos(np.radians(r.dec)):.6f} AND {r.ra + rad/np.cos(np.radians(r.dec)):.6f}"
                           f" AND mean_fiber_dec BETWEEN {r.dec - rad:.6f} AND {r.dec + rad:.6f})" for r in sub.itertuples())
        q = ("SELECT targetid, mean_fiber_ra AS ra, mean_fiber_dec AS dec, z, zwarn AS zwarning, spectype AS class, survey, program, "
             f"mean_mjd AS mjd, min_mjd, max_mjd, coadd_numexp, coadd_numnight FROM desi_dr1.zpix WHERE {cond}")
        t = None
        for attempt in range(3):
            try:
                t = TAP.search(q).to_table().to_pandas(); break
            except Exception as e:
                print(f'   DESI attempt {attempt+1}: {str(e)[:120]}', flush=True); time.sleep(10)
        if t is None:
            print('   DESI batch failed; continuing without it', flush=True); continue
        # assign each row to the nearest target in this batch
        if len(t):
            d = np.array([[np.hypot((t.ra - r.ra) * np.cos(np.radians(r.dec)), t.dec - r.dec) for r in sub.itertuples()]]).squeeze(0)
            t['idx'] = sub.idx.values[np.argmin(d, axis=0)]
        out.append(t)
        print(f'   DESI batch {i//BATCH+1}/{(len(targets)-1)//BATCH+1}: {len(t)} rows', flush=True)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def summarise(ep, targets):
    ep = ep[ep['is_coadd'] != True] if 'is_coadd' in ep else ep
    rows = []
    for r in targets.itertuples():
        e = ep[ep.idx == r.idx].sort_values('mjd')
        cls = [str(c) for c in e['class'].fillna('?')]
        sub = [str(s) for s in e['subclass'].fillna('')]
        known = [c for c in cls if c not in ('?', 'nan')]
        bl = [('BROADLINE' in s.upper()) for s, c in zip(sub, cls) if c not in ('?', 'nan')]
        last_sdssv = e[(e.source == 'SDSS') & (e.sdss_phase == 5) & e.spectroflux_r.notna()]
        rows.append(dict(idx=r.idx, name=r.name, ra=r.ra, dec=r.dec, n_spec=len(e),
                         n_sdss=(e.source == 'SDSS').sum(), n_desi=(e.source == 'DESI').sum(),
                         mjd_first=e.mjd.min() if len(e) else np.nan, mjd_last=e.mjd.max() if len(e) else np.nan,
                         years_since_last_spec=(TODAY_MJD - e.mjd.max()) / 365.25 if len(e) else np.nan,
                         baseline_yr=(e.mjd.max() - e.mjd.min()) / 365.25 if len(e) else np.nan,
                         surveys='|'.join(sorted(set(e.survey.astype(str)))),
                         mjd_list='|'.join(f'{m:.0f}' for m in e.mjd), class_list='|'.join(cls), subclass_list='|'.join(sub),
                         class_change_flag=len(set(known)) > 1, broadline_change_flag=len(set(bl)) > 1,
                         last_class=cls[-1] if cls else '', last_subclass=sub[-1] if sub else '',
                         z_spec_median=e.z.median() if len(e) else np.nan,
                         r_nmgy_last_sdssv=last_sdssv.spectroflux_r.iloc[-1] if len(last_sdssv) else np.nan,
                         mjd_last_sdssv=last_sdssv.mjd.iloc[-1] if len(last_sdssv) else np.nan))
    return pd.DataFrame(rows)


def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    inp, tag = sys.argv[1], sys.argv[2]
    ra_col, dec_col, name_col = (sys.argv[3:6] if len(sys.argv) >= 6 else ('ra', 'dec', 'name'))
    tin = pd.read_csv(inp)
    targets = pd.DataFrame({'idx': np.arange(len(tin)), 'ra': tin[ra_col].astype(float), 'dec': tin[dec_col].astype(float),
                            'name': tin[name_col].astype(str)})
    t0 = time.time()
    print(f'{len(targets)} targets. SDSS DR19 allspec ...', flush=True)
    ep = sdss_allspec(targets)
    ep = sdss_classes(ep) if len(ep) else ep
    if len(ep):
        ep['source'] = 'SDSS'
        ep['is_coadd'] = ep['coadd'].astype(str).str.contains('allepoch')
    print(f'[{time.time()-t0:.0f}s] SDSS epochs: {len(ep)}  (coadd rows: {int(ep.is_coadd.sum()) if len(ep) else 0}). DESI DR1 ...', flush=True)
    de = desi_zpix(targets)
    if len(de):
        de['source'] = 'DESI'; de['sdss_phase'] = np.nan; de['is_coadd'] = False
        de['spectroflux_r'] = np.nan; de['subclass'] = ''
    allep = pd.concat([ep, de], ignore_index=True, sort=False)
    allep = allep.merge(targets[['idx', 'name']], on='idx', how='left')
    allep.to_csv(os.path.join(DATA, f'spectra_epochs_{tag}.csv'), index=False)
    summ = summarise(allep, targets)
    summ.to_csv(os.path.join(DATA, f'spectra_summary_{tag}.csv'), index=False)
    print(f'[{time.time()-t0:.0f}s] wrote data/spectra_epochs_{tag}.csv ({len(allep)} rows) and data/spectra_summary_{tag}.csv')
    print(summ[['n_spec', 'n_sdss', 'n_desi', 'years_since_last_spec', 'class_change_flag', 'broadline_change_flag']].describe(include='all').T.to_string())
    print('targets with 0 spectra:', (summ.n_spec == 0).sum(), '| with DESI:', (summ.n_desi > 0).sum(),
          '| last spectrum < 2 yr ago:', (summ.years_since_last_spec < 2).sum(),
          '| class change:', summ.class_change_flag.sum(), '| broad-line flag change:', summ.broadline_change_flag.sum())


if __name__ == '__main__':
    main()
