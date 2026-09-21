"""
14_select_targets.py  --  target selection v2 (2026-09-11): calibrated probabilities instead of the M + P score.

Strata (time fractions per night in FRAC):
  D1  turn-off discovery: DR16 pool quasars, no known CL history, no SDSS-V spectrum since 2024, ranked by P_hb_dim
      (probability that broad Hbeta has fallen > 2x since the archival spectrum; 14a model: optical change since the spectrum from
      ZTF vs synthetic photometry, W1 fade/amplitude, W1-W2 colour, Eddington ratio, manifold prior, z)
  D2  turn-on / brightening discovery: same pool ranked by P_hb_bright (the galaxy parent joins when its data are in: 14b)
  K   known changers due for a new epoch: Zeltyn+2024 CL-AGN/EVQ and Camus & Panda 2026 CLAGNs whose photometry moved since their
      last spectrum (recurrence science)
  B   blind manifold stratum: literature-CLAGN region of the W1 manifold AND kNN literature-CLAGN fraction >= 0.15, ranked by that
      fraction, archival spectrum as baseline, WITHOUT any photometric trigger (the paper's selection, tested on its own)
  (C  controls were dropped on 2026-09-13: the SDSS two-epoch calibration set plays that role)
Ranking inside a stratum = probability (or manifold prior for B, quietness for C) per hour of telescope time x moon weight, using the
exposure model and moon weights of 04_score_tiers.py.  Cuts are geometric only (hrs >= 1.5, moon >= 30 deg, z <= 0.8); brightness weights.
Outputs data/targets_<night>_v2.csv (+ backups) and data/selection_v2_summary.txt.  The v1 lists are untouched.
"""
import os, importlib.util, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(HERE, 'data')
spec = importlib.util.spec_from_file_location('f04', os.path.join(HERE, 'scripts', '04_score_tiers.py')); f04 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f04)
NIGHTS, HOURS = ['sep23', 'oct26', 'oct27'], f04.HOURS
FRAC = {'D1': 0.60, 'D2': 0.22, 'K': 0.10, 'B': 0.08}          # controls dropped 2026-09-13: the 5,578-object SDSS two-epoch set is the control
MJD_RECENT_SDSSV = 60300.0          # 2024-01: an SDSS-V spectrum after this already tells us the current state
LIST_LEN = 60


def brightness_now(df):
    r = df.ztf_r_medianmag.where(df.ztf_r_medianmag.notna(), df.get('r_mag', df.get('psfmag_r')))
    return r.fillna(df.get('psfmag_r', pd.Series(19.5, index=df.index)))


def main():
    P = pd.read_csv(os.path.join(DATA, 'pool_probabilities.csv')).set_index('name')
    pool = pd.read_csv(os.path.join(DATA, 'parent_pool_scored.csv')); pool = pool[pool.projected.fillna(False)].copy(); pool['name'] = 'P' + pool.poolid.astype(str)
    obs = pd.read_csv(os.path.join(DATA, 'parent_pool_obs.csv')); obs['name'] = 'P' + obs.poolid.astype(str)
    ocols = [c for c in obs.columns if c.startswith(('hrs_', 'minX_', 'moonsep_'))]
    P = P.join(obs.set_index('name')[ocols], how='left').join(pool.set_index('name')[['sdss_plate_mjd_fiber']] if 'sdss_plate_mjd_fiber' in pool else pool.set_index('name')[[]])
    P['plate_mjd_fiber'] = pool.set_index('name').plate.astype(str) + '-' + pool.set_index('name').mjd.astype(str) + '-' + pool.set_index('name').fiberid.astype(str)
    sv = pd.read_csv(os.path.join(DATA, 'sdssv_internal_epochs.csv')); last_sv = sv.groupby('name').mjd.max()
    P['mjd_last_sdssv'] = last_sv.reindex(P.index); P['recent_sdssv'] = P.mjd_last_sdssv >= MJD_RECENT_SDSSV
    P['r_now'] = brightness_now(P)
    P['clagn_score'] = pd.read_csv(os.path.join(DATA, 'manifold_prior_pool.csv')).set_index('name').clagn_score.reindex(P.index)   # kNN fraction of literature CLAGNs (Hemmati+2026)
    trig_off = (P.d_opt > 0.3) | (P.dw1_neowise > 0.3) | ((P.dW1_full > 0.3) & (P.dcolor < -0.1))
    trig_on = (P.d_opt < -0.3) | (P.dw1_neowise < -0.3)
    quiet = (P.d_opt.abs() < 0.2) & (P.dw1_neowise.abs() < 0.15) & (P.ampW1 < 0.4) & (P.dW1_full.abs() < 0.2)
    known = P.known_clagn.fillna(False).astype(bool)
    # ---------------- strata
    P['stratum'] = ''; P['score'] = np.nan; P['why'] = ''
    d1 = ~known & ~P.recent_sdssv & P.d_opt.notna() & (P.P_hb_dim >= 0.1)
    P.loc[d1, 'stratum'] = 'D1'; P.loc[d1, 'score'] = P.P_hb_dim[d1]
    d2 = ~known & ~P.recent_sdssv & P.d_opt.notna() & (P.P_hb_bright >= 0.1) & (P.stratum == '')
    P.loc[d2, 'stratum'] = 'D2'; P.loc[d2, 'score'] = P.P_hb_bright[d2]
    # blind manifold stratum (user 2026-09-13): the paper's own selection, W1-only manifold, literature-CLAGN region AND kNN literature
    # fraction >= 0.15 (the old threshold), an archival SDSS spectrum as baseline (every pool quasar has one), no photometric trigger;
    # ranked by the kNN literature-CLAGN fraction itself
    b = ~known & ~P.recent_sdssv & P.lit_corner.fillna(False).astype(bool) & (P.clagn_score >= 0.15) & ~trig_off.fillna(False) & ~trig_on.fillna(False) & P.d_opt.notna() & (P.stratum == '')
    P.loc[b, 'stratum'] = 'B'; P.loc[b, 'score'] = P.clagn_score[b]
    lowq = P.P_dim <= P.P_dim.quantile(0.25)
    c = ~known & ~P.recent_sdssv & quiet.fillna(False) & lowq & ~P.lit_corner.fillna(False).astype(bool) & (P.stratum == '') & (P.d_opt.notna())
    # controls (C) are no longer observed (user decision 2026-09-13); the definition is kept for reference:  P.loc[c, 'stratum'] = 'C'
    kq = known & ((P.dw1_neowise.abs() > 0.3) | (P.w1_slope_2yr.abs() > 0.15) | (P.d_opt.abs() > 0.5))
    P.loc[kq, 'stratum'] = 'K'; P.loc[kq, 'score'] = 0.5 + 0.5 * np.clip(P.dw1_neowise.abs()[kq].fillna(0) / 0.5, 0, 1)
    # plain-language reasons
    def why(r):
        parts = []
        if pd.notna(r.d_opt): parts.append(f'optical {"fainter" if r.d_opt > 0 else "brighter"} by {abs(r.d_opt):.2f} mag since the {r.plate_mjd_fiber} spectrum (ZTF vs synthetic)')
        if pd.notna(r.dw1_neowise) and abs(r.dw1_neowise) > 0.1: parts.append(f'W1 {"faded" if r.dw1_neowise > 0 else "brightened"} {abs(r.dw1_neowise):.2f} mag 2014-24')
        if pd.notna(r.dcolor) and abs(r.dcolor) > 0.08: parts.append(f'W1-W2 {"bluer" if r.dcolor < 0 else "redder"} by {abs(r.dcolor):.2f}')
        if pd.notna(r.ampW1) and r.ampW1 > 0.4: parts.append(f'W1 amplitude {r.ampW1:.2f} mag')
        if pd.notna(r.LOGLEDD_RATIO): parts.append(f'log Edd {r.LOGLEDD_RATIO:.2f}')
        parts.append(f'manifold prior P_hb {r.P_hb:.2f}' + (' (literature corner)' if r.lit_corner else '') + f', direction {"turn-on-like" if r.dir > 0.03 else ("turn-off-like" if r.dir < -0.03 else "neutral")}')
        if r.stratum == 'D1': head = f'P(broad Hbeta fell > 2x) = {r.P_hb_dim:.2f}, P(continuum dimmed) = {r.P_cont_dim:.2f}'
        elif r.stratum == 'D2': head = f'P(broad Hbeta rose > 2x) = {r.P_hb_bright:.2f}, P(continuum brightened) = {r.P_cont_bright:.2f}'
        elif r.stratum == 'B': head = f'BLIND manifold stratum (W1 manifold, Hemmati+2026): literature-CLAGN region, {100*r.clagn_score:.0f}% of the 50 nearest training objects are literature CLAGNs, archival SDSS spectrum as baseline, no photometric trigger'
        elif r.stratum == 'C': head = 'CONTROL: quiet manifold region, no photometric trigger'
        else: head = f'KNOWN changer ({r.known_type}) with recent photometric motion'
        return head + '; ' + '; '.join(parts)
    sel = P[P.stratum != ''].copy(); sel['why'] = [why(r) for r in sel.itertuples()]
    # ---------------- known changers from the Zeltyn sample (master list rows), with a trigger since their last spectrum
    m = pd.read_csv(os.path.join(DATA, 'master_list_private.csv'), low_memory=False); zz = m[(m.source_catalog == 'Zeltyn24') & m.tier.isin(['T3', 'T2'])].copy()
    ztrig = (zz.dr_since_ref.abs() > 0.3) | (zz.w1_slope_2yr.abs() > 0.15) | zz.reversal_candidate.fillna(False).astype(bool)
    zz = zz[ztrig].copy(); zz['stratum'] = 'K'; zz['score'] = 0.6 + 0.4 * np.clip(zz.dr_since_ref.abs().fillna(0) / 0.6, 0, 1); zz['r_now'] = zz.r_mag
    zz['why'] = ['KNOWN changer (Zeltyn ' + str(r.zeltyn_class) + '): ' + ('ZTF r moved ' + f'{r.dr_since_ref:+.2f}' + ' mag since the last spectrum; ' if pd.notna(r.dr_since_ref) else '') + (f'W1 slope {r.w1_slope_2yr:+.2f} mag/yr; ' if pd.notna(r.w1_slope_2yr) else '') + ('state-reversal candidate; ' if r.reversal_candidate else '') + f'last spectrum {r.years_since_last_spec:.1f} yr ago' for r in zz.itertuples()]
    zz['known_type'] = zz.zeltyn_class; zz.index = zz.name
    cols_common = ['ra', 'dec', 'z', 'r_now', 'stratum', 'score', 'why', 'known_type'] + ocols
    sel['arm'] = 'quasar'; zz['arm'] = 'zeltyn'
    frames = [sel[cols_common + ['arm', 'P_hb_dim', 'P_cont_dim', 'P_hb_bright', 'P_cont_bright', 'P_hb', 'P_dim', 'lit_corner', 'dir', 'd_opt', 'dw1_neowise', 'dcolor', 'ampW1', 'LOGLEDD_RATIO', 'mjd_last_sdssv']], zz[cols_common + ['arm']]]
    # turn-on arm from the galaxy parent (14b_turnon_arm.py), if built: joins stratum D2 with its literature-calibrated score
    tg = os.path.join(DATA, 'turnon_candidates.csv')
    if os.path.exists(tg):
        g = pd.read_csv(tg); g = g[~g.known_clagn.fillna(False).astype(bool)].copy(); g['stratum'] = 'D2'; g['arm'] = 'galaxy'; g['known_type'] = ''
        g = g.rename(columns={'dw1': 'dw1_neowise'}); g.index = g.name; frames.append(g[[c for c in cols_common + ['arm', 'd_opt', 'dw1_neowise', 'dcolor'] if c in g.columns]])
        print(f'turn-on arm: {len(g)} galaxy candidates added to D2')
    allc = pd.concat(frames, sort=False)
    allc.index.name = 'name'; allc = allc.reset_index()
    # v1-compatible aliases so 06 / 06b / 11 / 07 can consume the v2 tables: tier = stratum, r_mag = r_now, priority = score
    allc['tier'] = allc.stratum; allc['r_mag'] = allc.r_now; allc['priority'] = allc.score.round(3)
    allc['trend'] = np.select([allc.d_opt > 0.2, allc.d_opt < -0.2, allc.dw1_neowise > 0.2, allc.dw1_neowise < -0.2], ['fading', 'brightening', 'fading (IR)', 'brightening (IR)'], 'flat/unknown')
    allc['source_catalog'] = np.where(allc.arm == 'galaxy', 'SDSS_GALAXY', np.where(allc.arm == 'zeltyn', 'Zeltyn24', 'DR16_QSO'))
    allc['notes'] = np.where(allc.arm == 'galaxy', 'turn-on arm: SDSS narrow-line AGN galaxy', np.where(allc.arm == 'zeltyn', 'Zeltyn+2024 object', 'DR16 quasar'))
    mjd_dr16 = pool.set_index('name').mjd.reindex(allc.name).values; last_any = np.fmax(mjd_dr16.astype(float), allc.get('mjd_last_sdssv', pd.Series(np.nan, index=allc.index)).values.astype(float))
    allc['years_since_last_spec'] = np.where(np.isfinite(last_any), (61294.0 - last_any) / 365.25, np.nan)      # 61294 = 2026-09-11; Zeltyn/galaxy rows keep NaN here
    pm = pool.set_index('name'); mp = pd.read_csv(os.path.join(DATA, 'manifold_prior_pool.csv')).set_index('name')
    for col, src in [('umap_x', pm.umap_x), ('umap_y', pm.umap_y), ('clagn_score', pm.clagn_score), ('zeltyn_density_ratio', pm.zeltyn_density_ratio), ('in_region_zeltyn', pm.in_region_zeltyn), ('in_region_clagn', pm.in_region_clagn)]:
        allc[col] = src.reindex(allc.name).values
    allc['last_class'] = np.where(allc.arm == 'quasar', 'QSO', np.where(allc.arm == 'galaxy', 'GALAXY', '')); allc['n_spec'] = np.where(allc.arm == 'quasar', 1 + allc.get('mjd_last_sdssv', pd.Series(np.nan, index=allc.index)).notna().astype(int), np.nan)
    zm = m.set_index('name')
    for col in ['umap_x', 'umap_y', 'clagn_score', 'in_region_zeltyn', 'in_region_clagn', 'n_spec', 'mjd_last_spec', 'last_class', 'years_since_last_spec']:
        if col in zm: allc.loc[allc.arm == 'zeltyn', col] = zm[col].reindex(allc.name[allc.arm == 'zeltyn']).values
    allc.to_csv(os.path.join(DATA, 'candidates_v2.csv'), index=False)
    lines = [f'candidates by stratum: {allc.stratum.value_counts().to_dict()}']
    # ---------------- per-night allocation by expected yield per hour within each stratum's time fraction
    picked = set(); out = {}
    for night in NIGHTS:
        hrs, sep = allc[f'hrs_{night}'], allc[f'moonsep_{night}']
        elig = (hrs >= f04.MIN_HRS) & (sep >= f04.MIN_MOONSEP) & (allc.z <= f04.ZMAX) & ~allc.name.isin(picked)
        cand = allc[elig].copy(); cand['t_exp_min'] = f04.exposure_minutes(cand.r_now.values, cand.z.values, sep[elig].values).round(0); cand['t_total_min'] = cand.t_exp_min + f04.OVERHEAD_MIN
        cand['moon_w'] = f04.moon_weight(sep[elig].values); cand['yield_per_hour'] = (60.0 * cand.score * cand.moon_w / cand.t_total_min).round(3)
        cand['priority_night'] = (cand.score * cand.moon_w).round(3); cand['prio_per_hour'] = cand.yield_per_hour
        nexp = np.maximum(2, np.ceil(cand.t_exp_min / 10.0)).astype(int); cand['exp_plan'] = [f'{n} x {t/n:.0f} min' for n, t in zip(nexp, cand.t_exp_min)]
        cand = cand.sort_values('yield_per_hour', ascending=False); budget = HOURS[night] * 60.0; chosen = []; used = {s: 0.0 for s in FRAC}
        for s, frac in FRAC.items():
            for ix, row in cand[cand.stratum == s].iterrows():
                if used[s] + row.t_total_min > frac * budget: continue
                chosen.append(ix); used[s] += row.t_total_min
        # leftover time (strata that could not fill their share) goes to the best remaining D1/D2 candidates
        left = budget - sum(used.values())
        for ix, row in cand[cand.stratum.isin(['D1', 'D2']) & ~cand.index.isin(chosen)].iterrows():
            if row.t_total_min <= left: chosen.append(ix); left -= row.t_total_min
        prim = cand.loc[chosen].sort_values('yield_per_hour', ascending=False).copy(); prim['rank'] = np.arange(1, len(prim) + 1); prim['night'] = night
        back = cand[~cand.index.isin(chosen)].head(max(LIST_LEN - len(prim), 0)).copy(); back['rank'] = 0; back['night'] = night
        for d in (prim, back):
            d['why_night'] = d.why + f'; {night}: x{{}} moon, {{}} min + 5 overhead -> yield/h {{}}'
            d['why_night'] = [w.format(mw, int(t), y) for w, mw, t, y in zip(d.why_night, d.moon_w, d.t_exp_min, d.yield_per_hour)]
        out[night] = pd.concat([prim, back]); picked |= set(prim.name)
        exp = {s: round(float(prim.score[prim.stratum == s].sum()), 1) for s in ['D1', 'D2']}
        lines.append(f'{night}: {len(prim)} targets, {sum(used.values())/60:.1f} of {HOURS[night]} h; by stratum {prim.stratum.value_counts().to_dict()}; '
                     f'expected broad-line changes D1 {exp["D1"]} / D2 {exp["D2"]}; median r_now {prim.r_now.median():.1f}, median P_hb_dim (D1) {prim.P_hb_dim[prim.stratum == "D1"].median():.2f}')
    cols = ['rank', 'night', 'stratum', 'tier', 'arm', 'name', 'ra', 'dec', 'z', 'r_now', 'r_mag', 't_exp_min', 'exp_plan', 'yield_per_hour', 'prio_per_hour', 'priority_night', 'priority', 'score', 'trend', 'source_catalog', 'notes', 'years_since_last_spec', 'P_hb_dim', 'P_cont_dim', 'P_hb_bright', 'P_cont_bright', 'P_hb', 'lit_corner', 'dir',
            'd_opt', 'dw1_neowise', 'dcolor', 'ampW1', 'LOGLEDD_RATIO', 'known_type', 'mjd_last_sdssv', 'why', 'why_night']
    for night, df in out.items():
        df[[c for c in cols if c in df.columns] + [f'hrs_{night}', f'minX_{night}', f'moonsep_{night}']].to_csv(os.path.join(DATA, f'targets_{night}_v2.csv'), index=False)
    open(os.path.join(DATA, 'selection_v2_summary.txt'), 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
    for night, df in out.items():
        p = df[df['rank'] > 0]; print(f'\n=== {night} ===')
        print(p[['rank', 'stratum', 'name', 'z', 'r_now', 't_exp_min', 'yield_per_hour', 'score', 'd_opt', 'dw1_neowise', 'dcolor', 'LOGLEDD_RATIO', 'P_hb', 'lit_corner']].round(2).head(45).to_string(index=False))


if __name__ == '__main__':
    main()
