"""
03g_sdssv_discovery.py  --  spectroscopic discovery pass over the proprietary SDSS-V spAll (objects NOT limited to our master list).

Our master list is built from DR16 quasars, so comparing it with new epochs can only reveal turn-offs and repeat changes.
This pass looks at every SDSS-V spectrum inside the observable RA/Dec windows and finds objects whose *own* SDSS-V epochs
disagree: class GALAXY -> QSO (turn-on), QSO -> GALAXY (turn-off), or a spectrophotometric r change larger than --dr mag.
These are selected by their spectra, not by the WISE manifold: keep them as a separate tier when evaluating the selection.

Two streaming passes over the cached spAll-lite (03f_sdssv_internal.py downloads it; ~22 s per pass from local disk):
  1. compact pass: all science rows in the windows -> per-object epoch history -> candidate flags
  2. full pass for the candidates only -> data/sdssv_discovery_epochs.csv (03d schema, with sas_url) and
     data/sdssv_discovery.csv (one row per object: flags, z, r at first/last epoch, in_master, J-name)
Both outputs are PROPRIETARY (git-ignored).  Credentials: same login file as 03f (see its docstring).

Usage: /opt/anaconda3/bin/python 03g_sdssv_discovery.py [--version v6_2_1] [--coadd daily] [--zmax 0.9] [--dr 0.8] [--snmin 2] [--baseline 180]
"""
import os, sys, time, argparse, importlib.util
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('f03', os.path.join(HERE, '03f_sdssv_internal.py'))
f03 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f03)
DATA = f03.DATA
RA_WINDOWS_H = [(15.5, 24.0), (0.0, 4.5), (6.5, 11.0)]      # same as 02_parent_pool.py (airmass < 2 on the run nights)
DEC_MIN = -15.0
CLS = {b'QSO': 1, b'GALAXY': 2, b'STAR': 3}
CLS_NAME = {0: '?', 1: 'QSO', 2: 'GALAXY', 3: 'STAR'}
COMPACT = np.dtype([('sdss_id', 'i8'), ('catalogid', 'i8'), ('mjd', 'i4'), ('ra', 'f8'), ('dec', 'f8'), ('cls', 'i1'),
                    ('z', 'f4'), ('zwarning', 'i4'), ('sn', 'f4'), ('fr', 'f4'), ('nspecobs', 'i2'), ('bhm', 'i1'), ('fq', 'i1')])
# bhm: 1 = Black Hole Mapper carton (bhm_*) or an open-fibre/extragalactic target type; 0 = Milky Way Mapper star.  MWM stars
# (OB, CV, YSO, WD, A-star cartons) routinely flip between the GALAXY and QSO templates at z ~ 0 and would swamp the list.


def in_windows(ra, dec):
    ok = np.zeros(len(ra), bool)
    for h0, h1 in RA_WINDOWS_H:
        ok |= (ra >= h0 * 15) & (ra < h1 * 15)
    return ok & (dec > DEC_MIN)


def select_science(arr):
    ra, dec = f03.radec(arr)
    ok = np.isfinite(ra) & np.isfinite(dec) & in_windows(ra, dec) & (arr['CATALOGID'] > 0)
    if 'OBJTYPE' in arr.dtype.names:
        ot = np.char.strip(arr['OBJTYPE'])
        ok &= ~(np.char.startswith(ot, b'SKY') | np.char.startswith(ot, b'SPECTROPHOTO'))
    return ok


def compact(arr):
    out = np.zeros(len(arr), COMPACT)
    out['sdss_id'] = arr['SDSS_ID'] if 'SDSS_ID' in arr.dtype.names else -1
    out['catalogid'] = arr['CATALOGID']; out['mjd'] = arr['MJD']
    out['ra'], out['dec'] = f03.radec(arr)
    cls = np.char.strip(arr['CLASS']); code = np.zeros(len(arr), 'i1')
    for k, v in CLS.items():
        code[cls == k] = v
    out['cls'] = code; out['z'] = arr['Z']; out['zwarning'] = arr['ZWARNING']; out['sn'] = arr['SN_MEDIAN_ALL']
    sf = np.asarray(arr['SPECTROFLUX'], dtype=float); out['fr'] = sf[:, 2] if sf.ndim == 2 else np.nan
    out['nspecobs'] = arr['NSPECOBS'] if 'NSPECOBS' in arr.dtype.names else -1
    carton = np.char.strip(arr['FIRSTCARTON']); ot = np.char.strip(arr['OBJTYPE']) if 'OBJTYPE' in arr.dtype.names else np.full(len(arr), b'')
    out['bhm'] = (np.char.startswith(carton, b'bhm') | np.char.startswith(carton, b'open') | (ot == b'QSO') | (ot == b'GALAXY')) & ~np.char.startswith(carton, b'mwm')
    out['fq'] = (np.char.strip(arr['FIELDQUALITY']) == b'good') if 'FIELDQUALITY' in arr.dtype.names else 1
    return out


def jname(ra, dec):
    c = SkyCoord(ra * u.deg, dec * u.deg)
    return ['J' + c[i].ra.to_string(unit=u.hour, sep='', precision=2, pad=True) + c[i].dec.to_string(sep='', precision=1, alwayssign=True, pad=True) for i in range(len(c))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default='v6_2_1'); ap.add_argument('--coadd', default='daily', choices=['daily', 'epoch'])
    ap.add_argument('--zmax', type=float, default=0.9, help='keep objects with median z below this (Hb/Ha reachable by NGPS)')
    ap.add_argument('--dr', type=float, default=0.8, help='spectrophotometric r change (mag) between first and last good epoch that flags a candidate')
    ap.add_argument('--snmin', type=float, default=2.0, help='SN_MEDIAN_ALL floor for an epoch to count')
    ap.add_argument('--baseline', type=float, default=180, help='minimum days between the first and last good epoch for a change to count')
    ap.add_argument('--sn-ends', dest='sn_ends', type=float, default=3.0, help='SN_MEDIAN_ALL floor for the two epochs that define the change')
    a = ap.parse_args()
    path = os.path.join(f03.CACHE, a.coadd, os.path.basename(f03.spall_url(a.version, a.coadd)))
    if not os.path.exists(path + '.ok'):
        sys.exit(f'{path} not downloaded yet: run 03f_sdssv_internal.py first')

    # ---- pass 1: compact history of every science spectrum in the windows
    t0 = time.time(); c = pd.DataFrame(f03.stream_spall(path, select_science, convert=compact))
    c['key'] = np.where(c.sdss_id > 0, c.sdss_id, -c.catalogid)
    print(f'{len(c):,} science spectra in the windows (RA {RA_WINDOWS_H} h, Dec > {DEC_MIN}) for {c.key.nunique():,} objects  [{time.time()-t0:.0f}s]')
    # a mis-positioned fibre records its own position and a faint, low-S/N spectrum (P9434, MJD 61213: 2.4" off, S/N 3,
    # r 1.5 mag too faint, "GALAXY"): drop epochs > 1" from the object's median position before anything else
    med = c.groupby('key')[['ra', 'dec']].transform('median')
    off = np.hypot((c.ra - med.ra) * np.cos(np.radians(c.dec)), c.dec - med.dec) * 3600
    n_off = int((off > 1.0).sum()); c = c[off <= 1.0]
    print(f'   {n_off:,} epochs dropped for sitting > 1\" from their object (mis-positioned fibres)')
    print(f'   of which {int(c.bhm.sum()):,} BHM/extragalactic-targeted spectra ({c[c.bhm == 1].key.nunique():,} objects); MWM stellar cartons are dropped here')
    good = c[(c.sn >= a.snmin) & (c.cls > 0) & (c.bhm == 1)].sort_values(['key', 'mjd'])
    g = good.groupby('key')
    h = pd.DataFrame(dict(n_epochs=g.size(), mjd_first=g.mjd.first(), mjd_last=g.mjd.last(), cls_first=g.cls.first(), cls_last=g.cls.last(),
                          n_qso=g.cls.apply(lambda s: int((s == 1).sum())), n_gal=g.cls.apply(lambda s: int((s == 2).sum())),
                          ra=g.ra.median(), dec=g.dec.median(), sn_min=g.sn.min(), catalogid=g.catalogid.last(), sdss_id=g.sdss_id.last(),
                          zw_first=g.zwarning.first(), zw_last=g.zwarning.last(), fq_first=g.fq.first(), fq_last=g.fq.last(),
                          z_first=g.z.first(), z_last=g.z.last(), sn_first=g.sn.first(), sn_last=g.sn.last()))
    # the change is judged between the FIRST and LAST good epoch, so only those two must be clean (intermediate visits with a
    # redshift warning or a bad field are common and would otherwise disqualify most objects)
    h['n_zwarn'] = (h.zw_first != 0).astype(int) + (h.zw_last != 0).astype(int); h['n_badfield'] = (h.fq_first == 0).astype(int) + (h.fq_last == 0).astype(int)
    h['z'] = h[['z_first', 'z_last']].median(axis=1); h['z_spread'] = (h.z_last - h.z_first).abs()
    fr = good[good.fr > 0].groupby('key').fr; h['r_first'] = 22.5 - 2.5 * np.log10(fr.first().reindex(h.index)); h['r_last'] = 22.5 - 2.5 * np.log10(fr.last().reindex(h.index))
    h['dr'] = h.r_last - h.r_first
    multi = h[h.n_epochs >= 2].copy()
    print(f'{len(h):,} objects with a good epoch (S/N >= {a.snmin}); {len(multi):,} with >= 2 good epochs; '
          f'{int(((h.n_epochs == 1) & (h.cls_last == 1) & (h.z < 0.8)).sum()):,} single-epoch QSOs at z<0.8 (pool for a DR17 turn-on match)')
    # class flips between consecutive nights are template noise, not physics: a change needs a baseline of >= --baseline days
    span_ok = (multi.mjd_last - multi.mjd_first) >= a.baseline
    multi['turn_on'] = (multi.cls_first == 2) & (multi.cls_last == 1) & span_ok
    multi['turn_off'] = (multi.cls_first == 1) & (multi.cls_last == 2) & span_ok
    multi['brightened'] = (multi.dr < -a.dr) & span_ok
    multi['dimmed'] = (multi.dr > a.dr) & span_ok
    # a real change keeps its redshift: consistent z (spread < 1 % of 1+z) between the good epochs, z > 0.02 (not a star),
    # no redshift warning and a good field on the epochs used; spectrophotometric jumps need S/N >= 5 on both ends
    ok = (multi.z > 0.02) & (multi.z < a.zmax) & (multi.z_spread < 0.01 * (1 + multi.z)) & (multi.n_zwarn == 0) & (multi.n_badfield == 0) \
         & (multi.cls_last != 3) & (multi.cls_first != 3) & (multi[['sn_first', 'sn_last']].min(axis=1) >= a.sn_ends)
    sn_ends = multi[['sn_first', 'sn_last']].min(axis=1); multi['brightened'] &= sn_ends >= 5; multi['dimmed'] &= sn_ends >= 5
    cand = multi[ok & (multi.turn_on | multi.turn_off | multi.brightened | multi.dimmed)].copy()
    print(f'   {int(ok.sum()):,} multi-epoch BHM objects pass the quality gate (consistent z in 0.02-{a.zmax}, zwarning 0, good fields)')

    # ---- known or new?  (2" match to the master list)
    m = pd.read_csv(os.path.join(DATA, 'master_list_scored.csv'), low_memory=False, usecols=['name', 'ra', 'dec', 'tier'])
    cm = SkyCoord(m.ra.values * u.deg, m.dec.values * u.deg); cc = SkyCoord(cand.ra.values * u.deg, cand.dec.values * u.deg)
    idx, sep, _ = cc.match_to_catalog_sky(cm); hit = sep.arcsec < 2
    cand['in_master'] = hit; cand['master_name'] = np.where(hit, m.name.values[idx], ''); cand['master_tier'] = np.where(hit, m.tier.values[idx], '')
    cand['name'] = jname(cand.ra.values, cand.dec.values)
    cand['cls_first'] = cand.cls_first.map(CLS_NAME); cand['cls_last'] = cand.cls_last.map(CLS_NAME)
    cand['baseline_yr'] = (cand.mjd_last - cand.mjd_first) / 365.25

    def counts(df):
        return dict(turn_on=int(df.turn_on.sum()), turn_off=int(df.turn_off.sum()), brightened=int(df.brightened.sum()), dimmed=int(df.dimmed.sum()), total=len(df))
    print(f'candidates (z < {a.zmax}, |dr| > {a.dr} or class change, S/N >= {a.snmin}):')
    print(f'   in master list : {counts(cand[cand.in_master])}')
    print(f'   NEW            : {counts(cand[~cand.in_master])}')
    bright = cand[~cand.in_master & (cand.r_last < 19.5)]
    print(f'   NEW and r_last < 19.5 : {counts(bright)}')

    # ---- pass 2: full rows for the candidates -> 03d-style epoch table
    ids = set(cand.catalogid.astype(np.int64))
    ep = f03.spall_table(f03.stream_spall(path, lambda arr: np.isin(arr['CATALOGID'], list(ids))))
    ep['key'] = np.where(ep.sdss_id > 0, ep.sdss_id, -ep.catalogid)
    ep = ep.merge(cand[['name', 'in_master', 'master_name']].rename_axis('key').reset_index(), on='key', how='inner')
    ep['run2d'] = a.version; ep['coadd'] = a.coadd; ep['sdss_phase'] = 5; ep['source'] = 'SDSS'; ep['proprietary'] = True; ep['is_coadd'] = False
    ep['sas_url'] = [f03.spec_url(a.version, a.coadd, f, j, c, sf) for f, j, c, sf in zip(ep.field, ep.mjd, ep.catalogid, ep.spec_file)]
    ep = ep.rename(columns={'cls': 'class'}).sort_values(['name', 'mjd'])
    cols03 = ['name', 'ra', 'dec', 'mjd', 'sdss_phase', 'run2d', 'coadd', 'programname', 'firstcarton', 'class', 'subclass', 'z', 'zwarning',
              'sn_median_all', 'spectroflux_g', 'spectroflux_r', 'spectroflux_i', 'catalogid', 'field', 'sas_url', 'is_coadd', 'source', 'proprietary', 'in_master', 'master_name']
    ep[cols03].rename(columns={'field': 'plate_or_fps_field'}).to_csv(os.path.join(DATA, 'sdssv_discovery_epochs.csv'), index=False)
    carton = ep.groupby('name').firstcarton.first()
    cand = cand.set_index('name'); cand['carton'] = carton.reindex(cand.index); cand = cand.reset_index()
    out = ['name', 'ra', 'dec', 'z', 'z_spread', 'n_epochs', 'mjd_first', 'mjd_last', 'baseline_yr', 'cls_first', 'cls_last', 'n_qso', 'n_gal', 'r_first', 'r_last', 'dr',
           'sn_min', 'turn_on', 'turn_off', 'brightened', 'dimmed', 'in_master', 'master_name', 'master_tier', 'carton', 'sdss_id', 'catalogid']
    cand[out].sort_values(['in_master', 'turn_on', 'turn_off', 'dr'], ascending=[True, False, False, True]).to_csv(os.path.join(DATA, 'sdssv_discovery.csv'), index=False)
    print(f'wrote data/sdssv_discovery.csv ({len(cand)} objects) and data/sdssv_discovery_epochs.csv ({len(ep)} epochs)')
    show = cand[~cand.in_master & (cand.turn_on | cand.turn_off)].sort_values(['turn_on', 'r_last'], ascending=[False, True])
    if len(show):
        print('NEW class-change candidates (turn-ons first):')
        print(show[['name', 'z', 'n_epochs', 'mjd_first', 'mjd_last', 'cls_first', 'cls_last', 'r_first', 'r_last', 'sn_min', 'carton']].head(40).to_string(index=False, float_format=lambda x: f'{x:.2f}'))


if __name__ == '__main__':
    main()
