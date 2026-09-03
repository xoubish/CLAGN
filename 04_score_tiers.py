"""
04_score_tiers.py  --  Stages 3d/4 of FOLLOWUP_PLAN.md: merge, score, tier, allocate to nights.

Inputs (all optional except the Zeltyn ones; missing files are skipped with a message)
  data/zeltyn_candidates_obs.csv       Zeltyn CL-AGNs/EVQs: embedding, spectra summary, observability (05_observability.py)
  data/parent_pool_scored.csv          DR16 quasar pool with manifold scores (02_parent_pool.py project)
  data/parent_pool_obs.csv             same pool with observability columns (05_observability.py)
  data/spectra_summary_pool.csv        spectra inventory for the pool subset (03_spectra_inventory.py)
  data/ztf_now_<tag>.csv               current ZTF photometry (03b_ztf_now.py), tags 'zeltyn' and 'pool'
  data/wise_cache/pool/*.parquet       unWISE W1 light curves -> W1 change since the archival SDSS spectrum

Priority (recorded in FOLLOWUP_PLAN.md Section 4):
  M = clip(max(zeltyn_density_ratio/3, clagn_score/0.15), 0, 2)           manifold CLAGN-likeness (1 = threshold)
  P = clip(|log10(flux_now/flux_at_spec)| / log10(1.5), 0, 2)            photometric change since last spectrum
      (W1 from unWISE to 2020-12, ZTF r to 2025 when available; the larger of the two)
  S = 1 if years_since_last_spec >= 3 else 0.3
  B = 1.0 (r<18.5) / 0.8 (18.5-19) / 0.5 (19-19.5) / 0.25 (fainter)
  priority = B * S * (M + P) + 0.5 * class_change_flag
Tiers: T3 confirmed Zeltyn CL-AGN; T2 Zeltyn EVQ with M >= 1 or in the Zeltyn region; T1 pool object with M >= 1
       (not a Zeltyn object); T4 control: pool object with clagn_score == 0 and zeltyn_density_ratio <= 0.3, matched in
       z and r to the T1 picks; T0 = everything else (kept in the master table, never allocated).
Hard cuts for allocation: hrs_<night> >= 1.5, moonsep_<night> >= 40, z <= 0.8 (T3 exempt from z), r <= 19.5.
"""
import os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
NIGHTS = {'sep23': 14, 'oct26': 34, 'oct27': 34}          # target slots per night (~20 min each incl. overhead)
QUOTA = {'T3': 3, 'T4': 4}                                # per-night caps for revisit / control tiers
MIN_HRS, MIN_MOONSEP, RMAX, ZMAX = 1.5, 40.0, 19.5, 0.8


def load(name, required=False):
    p = os.path.join(DATA, name)
    if os.path.exists(p):
        return pd.read_csv(p)
    if required:
        raise FileNotFoundError(p)
    print(f'   (no {name}; skipping)'); return None


def w1_change_from_cache(cache_dir, mjd_ref_by_id):
    """flux ratio median(last 3 W1 epochs) / W1 near mjd_ref, per objectid, from cached unWISE chunks."""
    if not os.path.isdir(cache_dir):
        return {}
    lc = pd.concat([pd.read_parquet(os.path.join(cache_dir, f)) for f in sorted(os.listdir(cache_dir)) if f.endswith('.parquet')])
    lc = lc[lc.index.get_level_values('band') == 'WISE_W1'].reset_index()[['objectid', 'time', 'flux']]
    out = {}
    for oid, g in lc.groupby('objectid'):
        g = g.sort_values('time')
        if len(g) < 6:
            continue
        now = g.flux.tail(3).median()
        ref = mjd_ref_by_id.get(oid, np.nan)
        near = g[np.abs(g.time - ref) < 400] if np.isfinite(ref) else g.head(0)
        at_ref = near.flux.median() if len(near) >= 2 else g.flux.head(3).median()   # fallback: earliest epochs
        out[oid] = dict(w1_now_mjy=now, w1_at_spec_mjy=at_ref, w1_ratio_now_over_spec=now / at_ref if at_ref > 0 else np.nan,
                        w1_mjd_last=g.time.max())
    return out


def brightness_weight(r):
    return np.select([r < 18.5, r < 19.0, r < 19.5], [1.0, 0.8, 0.5], 0.25)


def build_master():
    # ---------------- Zeltyn tiers
    z = load('zeltyn_candidates_obs.csv', required=True)
    zt = pd.DataFrame({'name': z['name'], 'source_catalog': 'Zeltyn24', 'zeltyn_class': z['class_zeltyn'], 'ra': z.ra, 'dec': z.dec,
                       'z': z.z, 'r_mag': z.r_mag_last_sdssv, 'umap_x': z.umap_x, 'umap_y': z.umap_y, 'clagn_score': z.clagn_score,
                       'in_region_zeltyn': z.in_region_zeltyn, 'in_region_clagn': z.in_region_clagn, 'zeltyn_density_ratio': np.nan,
                       'n_spec': z.n_spec, 'mjd_last_spec': z.mjd_last, 'years_since_last_spec': z.years_since_last_spec,
                       'class_change_flag': z.class_change_flag.fillna(False).astype(bool), 'last_class': z.last_class})
    for c in [c for c in z.columns if c.startswith(('hrs_', 'minX_', 'moonsep_'))]:
        zt[c] = z[c]
    zn = load('ztf_now_zeltyn.csv')
    if zn is not None:
        zt = zt.merge(zn, on='name', how='left')
        # Zeltyn r_mag from SDSS-V synthetic flux; prefer current ZTF r when available
        zt['r_mag'] = zt['r_last'].fillna(zt['r_mag'])
    zt['tier'] = np.where(z.is_clagn, 'T3', np.where((z.clagn_score >= 0.15) | z.in_region_zeltyn, 'T2', 'T0'))
    zt['notes'] = 'Zeltyn+2024 ' + z['class_zeltyn'].astype(str) + ', SDSS-V last epoch MJD ' + z.mjd_last.round(0).astype('Int64').astype(str)
    frames = [zt]

    # ---------------- pool tiers
    ps = load('parent_pool_scored.csv')
    if ps is not None:
        po = load('parent_pool_obs.csv')
        if po is not None:
            ps = ps.merge(po[['poolid'] + [c for c in po.columns if c.startswith(('hrs_', 'minX_', 'moonsep_'))]], on='poolid', how='left')
        ps = ps[ps.projected.fillna(False)]
        pt = pd.DataFrame({'name': ['P' + str(i) for i in ps.poolid], 'poolid': ps.poolid, 'source_catalog': 'DR16_QSO', 'ra': ps.ra, 'dec': ps.dec,
                           'z': ps.z, 'r_mag': ps.psfmag_r, 'umap_x': ps.umap_x, 'umap_y': ps.umap_y, 'clagn_score': ps.clagn_score,
                           'in_region_zeltyn': ps.in_region_zeltyn, 'in_region_clagn': ps.in_region_clagn,
                           'zeltyn_density_ratio': ps.zeltyn_density_ratio, 'sdss_plate_mjd_fiber': ps.plate.astype(str) + '-' + ps.mjd.astype(str) + '-' + ps.fiberid.astype(str),
                           'n_spec': 1, 'mjd_last_spec': ps.mjd, 'class_change_flag': False, 'last_class': 'QSO'})
        for c in [c for c in ps.columns if c.startswith(('hrs_', 'minX_', 'moonsep_'))]:
            pt[c] = ps[c].values
        sp = load('spectra_summary_pool.csv')
        if sp is not None:
            sp = sp.rename(columns={'name': 'name'})
            pt = pt.merge(sp[['name', 'n_spec', 'mjd_last', 'class_change_flag', 'last_class', 'n_desi']].rename(
                columns={'n_spec': 'n_spec_inv', 'mjd_last': 'mjd_last_inv', 'class_change_flag': 'ccf_inv', 'last_class': 'last_class_inv'}),
                on='name', how='left')
            has = pt.n_spec_inv.notna()
            pt.loc[has, 'n_spec'] = pt.loc[has, 'n_spec_inv']; pt.loc[has, 'mjd_last_spec'] = pt.loc[has, 'mjd_last_inv']
            pt.loc[has, 'class_change_flag'] = pt.loc[has, 'ccf_inv'].astype(bool); pt.loc[has, 'last_class'] = pt.loc[has, 'last_class_inv']
            pt = pt.drop(columns=['n_spec_inv', 'mjd_last_inv', 'ccf_inv', 'last_class_inv'])
        pt['years_since_last_spec'] = (61286.0 - pt.mjd_last_spec) / 365.25          # 61286 = 2026-09-03
        w1 = w1_change_from_cache(os.path.join(DATA, 'wise_cache', 'pool'), dict(zip(ps.poolid, ps.mjd)))
        if w1:
            w1df = pd.DataFrame.from_dict(w1, orient='index'); w1df.index.name = 'poolid'
            pt = pt.merge(w1df.reset_index(), on='poolid', how='left')
        zn = load('ztf_now_pool.csv')
        if zn is not None:
            pt = pt.merge(zn, on='name', how='left')
            pt['r_mag'] = pt['r_last'].fillna(pt['r_mag'])          # current ZTF r when available (bright time!)
        # exclude pool objects that are Zeltyn objects (already tiered there) or published CLAGNs
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        cz = SkyCoord(zt.ra.values * u.deg, zt.dec.values * u.deg); cp = SkyCoord(pt.ra.values * u.deg, pt.dec.values * u.deg)
        idx, sep, _ = cp.match_to_catalog_sky(cz)
        pt['is_zeltyn_dup'] = sep.arcsec < 2
        pt['notes'] = 'DR16 QSO, SDSS spectrum ' + pt.sdss_plate_mjd_fiber
        lit = load('literature_clagn.csv')
        pt['is_literature_clagn'] = False
        if lit is not None and len(lit):
            cl = SkyCoord(lit.ra.values * u.deg, lit.dec.values * u.deg)
            li, ls, _ = cp.match_to_catalog_sky(cl)
            pt['is_literature_clagn'] = ls.arcsec < 2
            pt.loc[pt.is_literature_clagn, 'notes'] = 'published CLAGN (' + lit.ref.values[li[pt.is_literature_clagn.values]] + ')'
        M = np.clip(np.fmax(pt.zeltyn_density_ratio / 3.0, pt.clagn_score / 0.15), 0, 2)
        pt['tier'] = np.where(pt.is_zeltyn_dup | pt.is_literature_clagn, 'T0', np.where(M >= 1, 'T1',
                     np.where((pt.clagn_score == 0) & (pt.zeltyn_density_ratio <= 0.3) & ~pt.in_region_zeltyn & ~pt.in_region_clagn, 'T4', 'T0')))
        frames.append(pt)

    m = pd.concat(frames, ignore_index=True, sort=False)

    # ---------------- proprietary SDSS-V epochs, if the user exported them (data/sdssv_internal_epochs.csv)
    # expected columns (case-insensitive): ra, dec (or plug_ra/plug_dec), mjd, class, subclass, z, zwarning; optional spectroflux_r
    sv = load('sdssv_internal_epochs.csv')
    if sv is not None and len(sv):
        sv.columns = [c.lower() for c in sv.columns]
        sv = sv.rename(columns={'plug_ra': 'ra', 'plug_dec': 'dec', 'fiber_ra': 'ra', 'fiber_dec': 'dec'})
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        cm = SkyCoord(m.ra.values * u.deg, m.dec.values * u.deg); cs = SkyCoord(sv.ra.values * u.deg, sv.dec.values * u.deg)
        idx, sep, _ = cs.match_to_catalog_sky(cm)
        sv = sv[sep.arcsec < 2].assign(mi=idx[sep.arcsec < 2])
        n_new = 0
        for mi, g in sv.groupby('mi'):
            newer = g[g.mjd > m.at[mi, 'mjd_last_spec']] if pd.notna(m.at[mi, 'mjd_last_spec']) else g
            if not len(newer):
                continue
            n_new += 1
            last = newer.sort_values('mjd').iloc[-1]
            m.at[mi, 'n_spec'] = (m.at[mi, 'n_spec'] if pd.notna(m.at[mi, 'n_spec']) else 0) + len(newer)
            m.at[mi, 'mjd_last_spec'] = last.mjd
            m.at[mi, 'years_since_last_spec'] = (61286.0 - last.mjd) / 365.25
            if 'class' in newer and pd.notna(last.get('class')):
                prev = str(m.at[mi, 'last_class'])
                m.at[mi, 'class_change_flag'] = bool(m.at[mi, 'class_change_flag']) or (prev not in ('', '?', 'nan') and prev != str(last['class']))
                m.at[mi, 'last_class'] = str(last['class'])
            if 'spectroflux_r' in newer and pd.notna(last.get('spectroflux_r')) and last.spectroflux_r > 0:
                m.at[mi, 'r_at_last_sdssv_internal'] = 22.5 - 2.5 * np.log10(last.spectroflux_r)
            m.at[mi, 'notes'] = str(m.at[mi, 'notes']) + f'; SDSS-V internal epoch MJD {last.mjd:.0f}'
        print(f'   SDSS-V internal epochs: {len(sv)} rows matched, {n_new} objects gained a newer epoch')
        if 'r_at_last_sdssv_internal' in m and 'r_last' in m:
            newer_ref = m.r_at_last_sdssv_internal.notna() & m.r_last.notna()
            m.loc[newer_ref, 'dr_since_ref'] = m.loc[newer_ref, 'r_last'] - m.loc[newer_ref, 'r_at_last_sdssv_internal']

    # ---------------- scores
    m['M'] = np.clip(np.fmax(m.zeltyn_density_ratio.fillna(0) / 3.0, m.clagn_score.fillna(0) / 0.15), 0, 2)
    p_w1 = np.abs(np.log10(m.get('w1_ratio_now_over_spec', pd.Series(np.nan, index=m.index)))) / np.log10(1.5)
    p_ztf = np.abs(m.get('dr_since_ref', pd.Series(np.nan, index=m.index))) / (2.5 * np.log10(1.5))
    m['P'] = np.clip(np.fmax(p_w1.fillna(0), p_ztf.fillna(0)), 0, 2)
    m['S'] = np.where(m.years_since_last_spec >= 3, 1.0, 0.3)
    m['B'] = brightness_weight(m.r_mag.fillna(20.5))
    m['priority'] = (m.B * m.S * (m.M + m.P) + 0.5 * m.class_change_flag.fillna(False).astype(float)).round(3)
    # expected transition direction from the photometry (fading type-1 -> turn-off candidate, brightening -> turn-on)
    ratio = m.get('w1_ratio_now_over_spec', pd.Series(np.nan, index=m.index))
    dr = m.get('dr_since_ref', pd.Series(np.nan, index=m.index))
    m['trend'] = np.select([(ratio < 1 / 1.3) | (dr > 0.3), (ratio > 1.3) | (dr < -0.3)], ['fading', 'brightening'], 'flat/unknown')
    return m


def allocate(m):
    picked = set(); lists = {}
    for night, nslots in NIGHTS.items():
        hrs, sep = m.get(f'hrs_{night}'), m.get(f'moonsep_{night}')
        if hrs is None:
            print(f'   no observability columns for {night}'); continue
        elig = (hrs >= MIN_HRS) & (sep >= MIN_MOONSEP) & (m.r_mag.fillna(99) <= RMAX) & ((m.z <= ZMAX) | (m.tier == 'T3')) \
               & m.tier.isin(['T1', 'T2', 'T3', 'T4']) & ~m.name.isin(picked)
        cand = m[elig].sort_values(['priority', 'r_mag'], ascending=[False, True])
        chosen = []
        counts = {'T1': 0, 'T2': 0, 'T3': 0, 'T4': 0}
        for _, row in cand.iterrows():
            if len(chosen) >= nslots:
                break
            if row.tier in QUOTA and counts[row.tier] >= QUOTA[row.tier]:
                continue
            chosen.append(row.name); counts[row.tier] += 1
        sel = m.loc[chosen].copy(); sel['night'] = night
        sel['rank'] = np.arange(1, len(sel) + 1)
        backups = cand[~cand.index.isin(chosen)].head(nslots).copy(); backups['night'] = night; backups['rank'] = 0
        lists[night] = pd.concat([sel, backups])
        picked |= set(sel.name)
        print(f'   {night}: {len(sel)} targets ({counts}), {len(backups)} backups; eligible pool {elig.sum()}')
    return lists


if __name__ == '__main__':
    m = build_master()
    m.to_csv(os.path.join(DATA, 'master_list_scored.csv'), index=False)
    print(f'master list: {len(m)} rows; tiers: {m.tier.value_counts().to_dict()}')
    lists = allocate(m)
    cols = ['rank', 'night', 'tier', 'name', 'ra', 'dec', 'z', 'r_mag', 'priority', 'M', 'P', 'trend', 'years_since_last_spec', 'n_spec',
            'last_class', 'clagn_score', 'zeltyn_density_ratio', 'in_region_zeltyn', 'notes']
    for night, df in lists.items():
        cols_n = cols + [f'hrs_{night}', f'minX_{night}', f'moonsep_{night}']
        df[[c for c in cols_n if c in df.columns]].to_csv(os.path.join(DATA, f'targets_{night}.csv'), index=False)
        print(f'\n=== {night}: top targets ===')
        print(df[df['rank'] > 0][[c for c in cols_n if c in df.columns]].head(40).to_string(index=False))
