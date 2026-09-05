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
NIGHTS = {'sep23': 20, 'oct26': 45, 'oct27': 45}          # loose upper bound on list length; the time budget below governs
HOURS = {'sep23': 4.3, 'oct26': 9.5, 'oct27': 9.5}        # usable science hours per window (twilight-to-mid / full night, minus ~0.5 h standards)
QUOTA = {'T3': 3, 'T4': 4}                                # per-night caps for revisit / control tiers
OVERHEAD_MIN = 5.0                                        # slew + acquisition + readout per target


def exposure_minutes(r, z, moonsep):
    """On-source minutes for S/N ~ 7 per bin on the diagnostic broad line, scaled from the proposal's NGPS ETC anchor
    (r = 18.5, dark sky: ~9 min for S/N ~ 10). Sky-limited scaling: t ~ (S/N)^2 * 10^(0.8 dr) * 10^(0.4 dsky).
    Bright moon adds ~2 mag of sky in the blue-green (Hb region) and ~1 mag in the red (Ha region) at > 60 deg,
    +0.5 mag more at 40-60 deg and +1 mag more inside 40 deg. Ha is usable to z ~ 0.55 in NGPS; beyond that the
    diagnostic is Hb, which costs the blue-green sky penalty. Floor 8 min (practical minimum), cap 60 min."""
    r = np.where(np.isfinite(r), r, 20.0); z = np.where(np.isfinite(z), z, 0.5)
    dsky = np.where(z <= 0.55, 1.0, 2.0) + np.select([moonsep >= 60, moonsep >= 40], [0.0, 0.5], 1.0)
    t = 9.0 * (7.0 / 10.0) ** 2 * 10 ** (0.8 * (r - 18.5)) * 10 ** (0.4 * dsky)
    return np.clip(t, 8.0, 60.0)
# per-night floors so the science tiers of the proposal are all represented (filled first, from each tier's own
# priority order); the remaining slots go by overall priority
FLOOR = {'sep23': {'T3': 1, 'T2': 1, 'T4': 2}, 'oct26': {'T3': 3, 'T2': 3, 'T4': 4}, 'oct27': {'T3': 3, 'T2': 3, 'T4': 4}}
# No magnitude cut: brightness only weights the priority (B). Moon separation is a weight too (W_moon), with a
# hard floor at MIN_MOONSEP degrees.
MIN_HRS, MIN_MOONSEP, ZMAX = 1.5, 30.0, 0.8


def moon_weight(sep):
    return np.select([sep >= 60, sep >= 40, sep >= MIN_MOONSEP], [1.0, 0.7, 0.4], 0.0)


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

    # ---------------- NEOWISE-R W1 to 2024 (03c_neowise_now.py), both tags
    neo = [x for x in (load('neowise_now_pool.csv'), load('neowise_now_zeltyn.csv')) if x is not None]
    if neo:
        neo = pd.concat(neo).drop_duplicates('name')
        m = m.merge(neo[['name', 'n_visits', 'mjd_last_neo', 'w1_first', 'w1_last', 'dw1_neowise', 'w1_slope_2yr', 'w1_flux_last_mjy']], on='name', how='left')
        if 'w1_at_spec_mjy' in m:
            m['w1_ratio_2024_over_spec'] = m.w1_flux_last_mjy / m.w1_at_spec_mjy

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

    # ---------------- archival spectral direction (03d_fetch_spectra.py) -> state-reversal candidates
    # spec_dir: what the archival Hβ EW did between first and last SDSS epoch; a target whose *current* photometry moves the
    # other way (turned off, now brightening; turned on, now fading) is a state-reversal candidate. The recurrent-CLAGN
    # literature (dim-state plateaus of ~4-7 yr; turn-on flares that fade within months to years) makes these the best
    # bets for catching a change in 2026, so they get the same +0.5 bonus as an archival class change.
    m['spec_dir'] = ''; m['reversal_candidate'] = False
    sl = load('spectra_lines.csv')
    if sl is not None and len(sl):
        S = sl[sl.source == 'SDSS'].sort_values('mjd'); g = S.groupby('name')
        f, l = g.first(), g.last()
        span = (l.mjd - f.mjd) > 365
        off = span & (l.EW_Hb_rest < 0.5 * f.EW_Hb_rest) & (f.EW_Hb_rest > 8)
        on = span & (l.EW_Hb_rest > 2 * f.EW_Hb_rest.clip(lower=1)) & (l.EW_Hb_rest > 8)
        sd = pd.Series('', index=f.index); sd[off] = 'turned off'; sd[on] = 'turned on'
        m['spec_dir'] = m.name.map(sd).fillna('')
    # ---------------- scores
    m['M'] = np.clip(np.fmax(m.zeltyn_density_ratio.fillna(0) / 3.0, m.clagn_score.fillna(0) / 0.15), 0, 2)
    nanS = pd.Series(np.nan, index=m.index)
    p_w1 = np.abs(np.log10(m.get('w1_ratio_now_over_spec', nanS))) / np.log10(1.5)          # unWISE to 2020 vs at-spectrum
    p_w1b = np.abs(np.log10(m.get('w1_ratio_2024_over_spec', nanS))) / np.log10(1.5)        # NEOWISE 2024 vs at-spectrum
    p_neo = np.abs(m.get('dw1_neowise', nanS)) / (2.5 * np.log10(1.5))                      # NEOWISE 2024 vs 2014 (all objects)
    p_ztf = np.abs(m.get('dr_since_ref', nanS)) / (2.5 * np.log10(1.5))
    m['P'] = np.clip(np.fmax.reduce([p_w1.fillna(0), p_w1b.fillna(0), p_neo.fillna(0), p_ztf.fillna(0)]), 0, 2)
    m['S'] = np.where(m.years_since_last_spec >= 3, 1.0, 0.3)
    m['B'] = brightness_weight(m.r_mag.fillna(20.5))
    m['priority'] = (m.B * m.S * (m.M + m.P) + 0.5 * m.class_change_flag.fillna(False).astype(float)).round(3)
    # controls (T4) are the opposite case: we want bright, photometrically QUIET objects off the CLAGN regions
    isT4 = m.tier == 'T4'
    m.loc[isT4, 'priority'] = (m.B[isT4] * m.S[isT4] * (1.0 - np.clip(m.P[isT4], 0, 1))).round(3)
    # diagnostic features inside the NGPS range (3200-10400 A) at each redshift
    NGPS = (3050.0, 10400.0)          # U 3050-4430, G 4250-5960, R 5620-7950, I 7530-10400 Å (technical specifications page)
    LINES = [('MgII', 2798.0), ('Hb+[OIII]', 5007.0), ('Ha', 6563.0), ('CaII_trip', 8600.0), ('[SIII]9531', 9531.0)]
    def lines_in_range(z):
        if not np.isfinite(z):
            return ''
        return '+'.join(n for n, w in LINES if NGPS[0] <= w * (1 + z) <= NGPS[1])
    m['lines_in_ngps'] = m.z.map(lines_in_range)
    # expected transition direction from the photometry (fading type-1 -> turn-off candidate, brightening -> turn-on)
    ratio = m.get('w1_ratio_now_over_spec', nanS)
    dr = m.get('dr_since_ref', nanS)
    dw1 = m.get('dw1_neowise', nanS)
    slope = m.get('w1_slope_2yr', nanS)
    fading = (ratio < 1 / 1.3) | (dr > 0.3) | (dw1 > 0.3) | (slope > 0.1)
    brightening = (ratio > 1.3) | (dr < -0.3) | (dw1 < -0.3) | (slope < -0.1)
    m['trend'] = np.select([fading & brightening, fading, brightening], ['mixed (IR vs optical)', 'fading', 'brightening'], 'flat/unknown')
    # second, near-independent vote from the ZTF g,r + W1 manifold (10_combined_rescore.py; Spearman 0.18 with the W1 score,
    # AUC 0.70 on held-out Zeltyn CL-AGNs): +0.5 when both manifolds agree (M_combined >= 1), -0.3 when the combined space
    # disagrees (M_combined < 0.5) and there is no large photometric change to fall back on (P < 1). Only enriched candidates
    # have ZTF light curves and therefore a combined score; others are unaffected.
    cs = load('combined_scores.csv')
    m['M_combined'] = np.nan
    if cs is not None and len(cs):
        cs = cs.drop_duplicates('name').set_index('name')
        mc = np.clip(np.fmax(cs.zeltyn_density / 3.0, cs.knn_lit_frac / 0.15), 0, 2)
        m['M_combined'] = m.name.map(mc)
        # discovery tiers only: T3 are confirmed CLAGNs (manifold vote irrelevant), and for controls a high combined
        # score is a contamination flag, not a merit
        agree = (m.M_combined >= 1) & m.tier.isin(['T1', 'T2'])
        disagree = (m.M_combined < 0.5) & (m.P < 1) & (m.tier == 'T1')
        bad_ctrl = (m.M_combined >= 1) & (m.tier == 'T4')
        m.loc[agree, 'priority'] = (m.loc[agree, 'priority'] + 0.5).round(3)
        m.loc[disagree, 'priority'] = (m.loc[disagree, 'priority'] - 0.3).clip(lower=0).round(3)
        m.loc[bad_ctrl, 'priority'] = (m.loc[bad_ctrl, 'priority'] * 0.3).round(3)
        m.loc[agree, 'notes'] = m.loc[agree, 'notes'].astype(str) + '; also in the CLAGN region of the ZTF+W1 manifold'
        m.loc[bad_ctrl, 'notes'] = m.loc[bad_ctrl, 'notes'].astype(str) + '; CLAGN-like in the ZTF+W1 manifold: weak control'
        print(f'   combined-manifold score for {int(m.M_combined.notna().sum())} objects: {int(agree.sum())} T1/T2 agree (+0.5), '
              f'{int(disagree.sum())} T1 disagree without photometric change (-0.3), {int(bad_ctrl.sum())} controls demoted')
    # W1 blending (03e_blending.py): if the WISE beam holds a substantial neighbour, the W1-based change and manifold
    # position are less trustworthy; keep the target but scale its priority by the clean fraction (floor 0.5) unless the
    # optical (ZTF) change confirms it, and say so in the notes.
    bl = [x for x in (load('blending_pool.csv'), load('blending_zeltyn.csv'), load('blending_targets.csv')) if x is not None]
    m['fracflux_w1'] = np.nan; m['n_ps1_8as'] = np.nan; m['nbr_min_dz'] = np.nan; m['w1_contam_est'] = np.nan; m['blend_flag'] = False; m['blend_kind'] = ''
    if bl:
        b = pd.concat(bl).drop_duplicates('name').set_index('name')
        for c in ['fracflux_w1', 'n_ps1_8as', 'nbr_min_dz', 'w1_contam_est']:
            m[c] = m.name.map(b[c])
        m['blend_kind'] = m.name.map(b['blend_kind']).fillna('')
        m['blend_flag'] = m.blend_kind == 'neighbour'
        # neighbour blends: the W1 change may belong to the neighbour -> scale by the clean fraction (floor 0.5) unless the
        # optical change (ZTF, sub-arcsecond) confirms it. Extended hosts only dilute the amplitude: annotate, no penalty.
        clean = (1 - m.w1_contam_est.fillna(0)).where(m.fracflux_w1.isna() | (m.fracflux_w1 > 1 - m.w1_contam_est.fillna(0)), m.fracflux_w1).clip(0.5, 1.0)
        # the blend penalty is waived when the change is confirmed independently of WISE: by the ZTF r change (sub-arcsecond
        # resolution) or by the archival spectra themselves (3" fibre: Hβ EW direction or a pipeline class change)
        ztf_ok = (np.abs(m.get('dr_since_ref', nanS)) > 0.3).fillna(False)
        spec_ok = (m.spec_dir != '') | m.class_change_flag.fillna(False).astype(bool)
        scale = np.where(m.blend_flag & ~ztf_ok & ~spec_ok, clean, 1.0)
        m['priority'] = (m.priority * scale).round(3)
        nb = m.blend_flag
        m.loc[nb, 'notes'] = m.loc[nb, 'notes'].astype(str) + '; W1 NEIGHBOUR BLEND: ~' + (100 * (1 - clean[nb])).round(0).astype(int).astype(str) + '% of the WISE flux from a neighbour within 8"'
        eh = m.blend_kind == 'extended host'
        m.loc[eh, 'notes'] = m.loc[eh, 'notes'].astype(str) + '; extended host in the WISE beam (unWISE own-flux fraction ' + (100 * m.fracflux_w1[eh]).round(0).astype(int).astype(str) + '%)'
        mn = m.blend_kind == 'minor neighbour'
        m.loc[mn, 'notes'] = m.loc[mn, 'notes'].astype(str) + '; minor neighbour within 8" (< 15% of the W1 flux)'
        print(f'   blending: {int(m.fracflux_w1.notna().sum())} objects checked, {int(nb.sum())} neighbour blends (priority scaled unless ZTF confirms), {int(eh.sum())} extended hosts (annotated)')
    rev = ((m.spec_dir == 'turned off') & brightening & ~fading) | ((m.spec_dir == 'turned on') & fading & ~brightening)
    m['reversal_candidate'] = rev.fillna(False)
    m.loc[m.reversal_candidate, 'priority'] = (m.loc[m.reversal_candidate, 'priority'] + 0.5).round(3)
    m.loc[m.spec_dir != '', 'notes'] = m.loc[m.spec_dir != '', 'notes'].astype(str) + '; archival Hβ: ' + m.loc[m.spec_dir != '', 'spec_dir']
    m.loc[m.reversal_candidate, 'notes'] = m.loc[m.reversal_candidate, 'notes'].astype(str) + ' -> photometry now points the other way: STATE-REVERSAL CANDIDATE'
    print(f'   archival Hβ direction known for {(m.spec_dir != "").sum()} objects; state-reversal candidates: {int(m.reversal_candidate.sum())}')
    return m


def allocate(m):
    picked = set(); lists = {}
    for night, nslots in NIGHTS.items():
        hrs, sep = m.get(f'hrs_{night}'), m.get(f'moonsep_{night}')
        if hrs is None:
            print(f'   no observability columns for {night}'); continue
        elig = (hrs >= MIN_HRS) & (sep >= MIN_MOONSEP) & ((m.z <= ZMAX) | (m.tier == 'T3')) \
               & m.tier.isin(['T1', 'T2', 'T3', 'T4']) & ~m.name.isin(picked)
        cand = m[elig].copy()
        cand['priority_night'] = (cand.priority * moon_weight(sep[elig].values)).round(3)   # moon distance weights, not cuts
        cand['t_exp_min'] = exposure_minutes(cand.r_mag.values, cand.z.values, sep[elig].values).round(0)
        cand['t_total_min'] = cand.t_exp_min + OVERHEAD_MIN
        # exposure plan: >= 2 sub-exposures for cosmic-ray rejection, each <= ~10 min
        nexp = np.maximum(2, np.ceil(cand.t_exp_min / 10.0)).astype(int)
        cand['exp_plan'] = [f'{n} x {t/n:.0f} min' for n, t in zip(nexp, cand.t_exp_min)]
        # rank by science return per minute of telescope time, so a bright target is not out-competed by a faint one
        cand['prio_per_hour'] = (60.0 * cand.priority_night / cand.t_total_min).round(3)
        cand = cand.sort_values(['prio_per_hour', 'r_mag'], ascending=[False, True])
        budget = HOURS[night] * 60.0
        chosen, used = [], 0.0
        counts = {'T1': 0, 'T2': 0, 'T3': 0, 'T4': 0}
        # pass 1: tier floors (only targets that cost <= 30 min on source; expensive ones stay as backups)
        for tier, nmin in FLOOR.get(night, {}).items():
            for ix in cand.index[(cand.tier == tier) & (cand.t_exp_min <= 30)][:nmin]:
                chosen.append(ix); counts[tier] += 1; used += cand.at[ix, 't_total_min']
        # pass 2: best return per hour until the night's minutes are spent (caps respected)
        for ix, row in cand.iterrows():
            if len(chosen) >= nslots or used + row.t_total_min > budget:
                if used + 6 + OVERHEAD_MIN > budget or len(chosen) >= nslots:
                    break
                continue
            if ix in chosen:
                continue
            if row.tier in QUOTA and counts[row.tier] >= QUOTA[row.tier]:
                continue
            chosen.append(ix); counts[row.tier] += 1; used += row.t_total_min
        sel = cand.loc[chosen].sort_values('prio_per_hour', ascending=False).copy(); sel['night'] = night
        sel['rank'] = np.arange(1, len(sel) + 1)                       # rank = return-per-hour order within the selected set
        backups = cand[~cand.index.isin(chosen)].head(nslots).copy(); backups['night'] = night; backups['rank'] = 0
        lists[night] = pd.concat([sel, backups])
        picked |= set(sel.name)
        print(f'   {night}: {len(sel)} targets ({counts}) using {used/60:.1f} of {HOURS[night]} h '
              f'(median exposure {sel.t_exp_min.median():.0f} min), {len(backups)} backups; eligible pool {elig.sum()}')
    return lists


if __name__ == '__main__':
    m = build_master()
    m.to_csv(os.path.join(DATA, 'master_list_scored.csv'), index=False)
    print(f'master list: {len(m)} rows; tiers: {m.tier.value_counts().to_dict()}')
    lists = allocate(m)
    cols = ['rank', 'night', 'tier', 'name', 'ra', 'dec', 'z', 'r_mag', 't_exp_min', 'exp_plan', 'prio_per_hour', 'priority_night', 'priority', 'M', 'P', 'trend',
            'years_since_last_spec', 'n_spec', 'last_class', 'clagn_score', 'zeltyn_density_ratio', 'in_region_zeltyn',
            'M_combined', 'fracflux_w1', 'n_ps1_8as', 'blend_flag', 'blend_kind', 'lines_in_ngps', 'spec_dir', 'reversal_candidate', 'notes']
    for night, df in lists.items():
        cols_n = cols + [f'hrs_{night}', f'minX_{night}', f'moonsep_{night}']
        df[[c for c in cols_n if c in df.columns]].to_csv(os.path.join(DATA, f'targets_{night}.csv'), index=False)
        print(f'\n=== {night}: top targets ===')
        print(df[df['rank'] > 0][[c for c in cols_n if c in df.columns]].head(40).to_string(index=False))
