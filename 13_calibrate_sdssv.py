"""13_calibrate_sdssv.py  --  calibration test: do manifold position (M) and W1 fade predict a spectral change between the DR16 spectrum and the
2023-2026 SDSS-V epoch?  Uses only pool quasars with a good post-DR19 SDSS-V epoch (proprietary: outputs data/sdssv_calibration_*.csv are git-ignored)."""
import os, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
os.chdir('/Users/shemmati/Desktop/CLAGN')
OUT = 'data'

pool = pd.read_csv('data/parent_pool_scored.csv'); pool = pool[pool.projected.fillna(False)].copy()
pool['name'] = 'P' + pool.poolid.astype(str)
pool['M'] = np.clip(np.fmax(pool.zeltyn_density_ratio.fillna(0) / 3, pool.clagn_score.fillna(0) / 0.15), 0, 2)
print(f'pool projected: {len(pool)}  M>=1: {(pool.M>=1).sum()}  in_region_zeltyn: {pool.in_region_zeltyn.sum()}  in_region_clagn: {pool.in_region_clagn.sum()}')

# ---- unWISE W1 light-curve features for every pool object (2010-2020)
cd = 'data/wise_cache/pool'
lc = pd.concat([pd.read_parquet(os.path.join(cd, f)) for f in sorted(os.listdir(cd)) if f.endswith('.parquet')])
lc = lc[lc.index.get_level_values('band') == 'WISE_W1'].reset_index()[['objectid', 'time', 'flux', 'err']].sort_values(['objectid', 'time'])
g = lc.groupby('objectid')
feat = pd.DataFrame({'n_w1': g.size(), 'w1_med': g.flux.median(), 'w1_first': g.flux.apply(lambda s: s.head(2).median()),
                     'w1_last': g.flux.apply(lambda s: s.tail(3).median()), 'w1_max': g.flux.max(), 'w1_min': g.flux.min(),
                     'w1_std': g.flux.std(), 't_last': g.time.max()})
# 2014-15 level (start of NEOWISE-R) so we can separate the 2010-14 gap from the 2014-20 trend
lc14 = lc[(lc.time > 56600) & (lc.time < 57300)].groupby('objectid').flux.median().rename('w1_2014')
lc18 = lc[(lc.time > 58000) & (lc.time < 58700)].groupby('objectid').flux.median().rename('w1_2018')
feat = feat.join(lc14).join(lc18)
feat['dW1_full'] = -2.5 * np.log10(feat.w1_last / feat.w1_first)       # >0 = faded 2010 -> 2020
feat['dW1_2014_20'] = -2.5 * np.log10(feat.w1_last / feat.w1_2014)
feat['dW1_2018_20'] = -2.5 * np.log10(feat.w1_last / feat.w1_2018)
feat['ampW1'] = 2.5 * np.log10(feat.w1_max / feat.w1_min)
feat['fvar'] = feat.w1_std / feat.w1_med
feat.index.name = 'poolid'
pool = pool.merge(feat.reset_index(), on='poolid', how='left')
print(f'W1 light curves: {pool.n_w1.notna().sum()} objects, median epochs {pool.n_w1.median():.0f}, last epoch MJD {pool.t_last.median():.0f}')

# ---- SDSS-V internal epochs (post-DR19): last good epoch per object
sv = pd.read_csv('data/sdssv_internal_epochs.csv')
sv = sv[sv.name.str.startswith('P')]
good = sv[(sv.zwarning == 0) & (sv.sn_median_all >= 2) & (sv.sep_arcsec < 1.0) & sv['class'].isin(['QSO', 'GALAXY']) & (sv.objtype == 'science')].copy()
good['r_spec'] = 22.5 - 2.5 * np.log10(good.spectroflux_r.where(good.spectroflux_r > 0))
good['g_spec'] = 22.5 - 2.5 * np.log10(good.spectroflux_g.where(good.spectroflux_g > 0))
good = good.sort_values('mjd')
last = good.groupby('name').last()[['mjd', 'class', 'z', 'sn_median_all', 'r_spec', 'g_spec', 'firstcarton']]
last['n_good'] = good.groupby('name').size()
last['n_gal'] = good.groupby('name')['class'].apply(lambda s: (s == 'GALAXY').sum())
last['sn_max'] = good.groupby('name').sn_median_all.max()
best = good.loc[good.groupby('name').sn_median_all.idxmax()].set_index('name')
last['class_bestsn'] = best['class']; last['r_spec_bestsn'] = best.r_spec
c = pool.merge(last.add_prefix('sv_').reset_index(), on='name', how='inner')
print(f'\npool objects with a good post-DR19 SDSS-V epoch: {len(c)} (of {sv.name.nunique()} matched); median SDSS-V MJD {c.sv_mjd.median():.0f}; '
      f'baseline DR16 spectrum -> SDSS-V: median {((c.sv_mjd - c.mjd)/365.25).median():.1f} yr (10-90%: {((c.sv_mjd - c.mjd)/365.25).quantile(.1):.1f}-{((c.sv_mjd - c.mjd)/365.25).quantile(.9):.1f})')
print('DR16 spectrum year distribution of the calibration set:'); print((c.mjd // 365.25 + 1858.9).astype(int).value_counts().sort_index().to_string())

# brightness change: SDSS-V spectrophotometric r vs DR16 PSF r, offset-corrected (fibre vs PSF, calibration)
c['dr_raw'] = c.sv_r_spec - c.psfmag_r
off = c.dr_raw.median(); c['dr'] = c.dr_raw - off
c['dg'] = (c.sv_g_spec - c.psfmag_g) - (c.sv_g_spec - c.psfmag_g).median()
print(f'median r offset (spec - PSF) = {off:+.2f} mag, robust scatter {1.4826*np.median(np.abs(c.dr)):.2f} mag')
c['flip_raw'] = c.sv_class == 'GALAXY'
c['flip_corr'] = c.flip_raw & (c.dr > 0.3)
c['dim07'] = c.dr > 0.7; c['dim05'] = c.dr > 0.5; c['bright07'] = c.dr < -0.7
c['changed'] = c.flip_corr | c.dim07 | c.bright07
c['dimmed_any'] = c.flip_corr | c.dim07
for k in ['flip_raw', 'flip_corr', 'dim05', 'dim07', 'bright07', 'changed']:
    print(f'   {k:10s} {c[k].sum():5d}  ({100*c[k].mean():.1f}%)')
print('   by carton (changed rate):'); print(c.groupby(c.sv_firstcarton.str.replace(r'_(wide|med|bonus).*', '', regex=True)).changed.agg(['size', 'mean']).sort_values('size', ascending=False).head(8).round(3).to_string())

def rate_table(df, col, bins, labels, outs=('flip_corr', 'dim07', 'bright07', 'changed')):
    b = pd.cut(df[col], bins, labels=labels, include_lowest=True)
    t = df.groupby(b, observed=False)[list(outs)].agg(['sum', 'mean'])
    t.columns = [f'{a}_{"n" if b_ == "sum" else "%"}' for a, b_ in t.columns]
    for o in outs:
        t[f'{o}_%'] = (100 * t[f'{o}_%']).round(1)
    t.insert(0, 'N', df.groupby(b, observed=False).size())
    return t

print('\n=== outcome rate vs manifold score M (current definition) ===')
print(rate_table(c, 'M', [-0.01, 0.001, 0.5, 1.0, 1.5, 2.01], ['M=0', '0-0.5', '0.5-1', '1-1.5', '1.5-2']).to_string())
print('\n=== by region flag ===')
for f in ['in_region_zeltyn', 'in_region_clagn']:
    t = c.groupby(f)[['flip_corr', 'dim07', 'changed']].agg(['size', 'mean']); print(f); print((t).round(3).to_string())
print('\n=== outcome rate vs kNN clagn_score (literature CLAGN neighbours) ===')
print(rate_table(c, 'clagn_score', [-0.01, 0.001, 0.08, 0.15, 0.25, 1.01], ['0', '0-0.08', '0.08-0.15', '0.15-0.25', '>0.25']).to_string())
print('\n=== outcome rate vs Zeltyn density ratio ===')
print(rate_table(c, 'zeltyn_density_ratio', [-0.01, 0.3, 1, 3, 6, 100], ['<0.3', '0.3-1', '1-3', '3-6', '>6']).to_string())
print('\n=== outcome rate vs W1 fade 2010 -> 2020 (unWISE, mag; >0 = fainter now) ===')
print(rate_table(c, 'dW1_full', [-9, -0.3, -0.1, 0.1, 0.3, 0.6, 9], ['brightened >0.3', '-0.3..-0.1', 'flat', '0.1-0.3', '0.3-0.6', 'faded >0.6']).to_string())
print('\n=== outcome rate vs W1 change 2014 -> 2020 ===')
print(rate_table(c, 'dW1_2014_20', [-9, -0.3, -0.1, 0.1, 0.3, 0.6, 9], ['brightened >0.3', '-0.3..-0.1', 'flat', '0.1-0.3', '0.3-0.6', 'faded >0.6']).to_string())
print('\n=== outcome rate vs W1 peak-to-peak amplitude ===')
print(rate_table(c, 'ampW1', [0, 0.2, 0.4, 0.7, 1.0, 9], ['<0.2', '0.2-0.4', '0.4-0.7', '0.7-1', '>1']).to_string())

print('\n=== AUC for predicting each outcome (higher = better ranker; 0.5 = chance) ===')
preds = {'M (current)': c.M, 'clagn_score': c.clagn_score, 'zeltyn_density_ratio': c.zeltyn_density_ratio.fillna(0),
         'W1 fade 2010-20 (dW1_full)': c.dW1_full, '|dW1_full|': c.dW1_full.abs(), 'W1 fade 2014-20': c.dW1_2014_20,
         'W1 fade 2018-20': c.dW1_2018_20, 'W1 amplitude': c.ampW1, 'W1 fvar': c.fvar, 'fvar_w1 (pool col)': c.fvar_w1,
         'M + P_unwise (current-style)': c.M + np.clip(c.dW1_full.abs() / (2.5*np.log10(1.5)), 0, 2),
         'DR16 psfmag_r (fainter first)': c.psfmag_r, 'redshift z': c.z, 'baseline yr': (c.sv_mjd - c.mjd)}
rows = []
for nm, p in preds.items():
    ok = p.notna()
    rows.append({'predictor': nm, 'N': int(ok.sum()), **{o: round(roc_auc_score(c[o][ok], p[ok]), 3) for o in ['flip_corr', 'dim07', 'dimmed_any', 'bright07', 'changed']}})
print(pd.DataFrame(rows).to_string(index=False))

print('\n=== does M add anything on top of the W1 fade? (2x2, dimmed_any = corroborated flip or >0.7 mag dimming) ===')
hiM = c.M >= 1; fade = c.dW1_full > 0.3
for a, la in [(~hiM, 'M<1'), (hiM, 'M>=1')]:
    for b, lb in [(~fade, 'W1 fade <=0.3'), (fade, 'W1 fade >0.3')]:
        s = c[a & b]; print(f'   {la:6s} & {lb:14s}: N={len(s):5d}  dimmed {100*s.dimmed_any.mean():5.1f}%  flip {100*s.flip_corr.mean():4.1f}%  dim>0.7 {100*s.dim07.mean():4.1f}%')
ok = c[['M', 'dW1_full', 'ampW1', 'psfmag_r', 'z']].notna().all(axis=1) & c.dimmed_any.notna()
X = c.loc[ok, ['M', 'dW1_full', 'ampW1', 'psfmag_r', 'z']].values; y = c.loc[ok, 'dimmed_any'].values
lr = LogisticRegression(max_iter=2000).fit((X - X.mean(0)) / X.std(0), y)
print('   logistic regression (standardised coefficients) for dimmed_any:', dict(zip(['M', 'dW1_full', 'ampW1', 'psfmag_r', 'z'], lr.coef_[0].round(2))))
X2 = c.loc[ok, ['dW1_full', 'ampW1']].values; lr2 = LogisticRegression(max_iter=2000).fit((X2 - X2.mean(0)) / X2.std(0), y)
print(f'   AUC full model {roc_auc_score(y, lr.decision_function((X - X.mean(0)) / X.std(0))):.3f} vs W1-only (fade + amplitude) {roc_auc_score(y, lr2.decision_function((X2 - X2.mean(0)) / X2.std(0))):.3f}')

print('\n=== redshift / magnitude dependence of the dimming rate ===')
print(rate_table(c, 'z', [0, 0.2, 0.4, 0.6, 0.81], ['0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8'], outs=('flip_corr', 'dim07', 'dimmed_any')).to_string())
print(rate_table(c, 'psfmag_r', [14, 18, 18.5, 19, 19.5, 21], ['<18', '18-18.5', '18.5-19', '19-19.5', '>19.5'], outs=('flip_corr', 'dim07', 'dimmed_any')).to_string())
print(rate_table(c.assign(bl=(c.sv_mjd - c.mjd)/365.25), 'bl', [0, 8, 14, 20, 30], ['<8 yr', '8-14', '14-20', '>20'], outs=('flip_corr', 'dim07', 'dimmed_any')).to_string())

# NEOWISE (enriched subset only, biased towards M>=1) for reference
neo = pd.read_csv('data/neowise_now_pool.csv')
cn = c.merge(neo[['name', 'dw1_neowise', 'w1_slope_2yr']], on='name', how='inner')
print(f'\n=== NEOWISE 2014-2024 W1 change, enriched subset with SDSS-V epoch: N={len(cn)} ===')
print(rate_table(cn, 'dw1_neowise', [-9, -0.3, 0.1, 0.3, 0.6, 9], ['brightened', 'flat', '0.1-0.3', '0.3-0.6', '>0.6'], outs=('flip_corr', 'dim07', 'dimmed_any')).to_string())
print('   AUC dw1_neowise for dimmed_any:', round(roc_auc_score(cn.dimmed_any, cn.dw1_neowise.fillna(0)), 3), ' vs dW1_full on the same objects:', round(roc_auc_score(cn.dimmed_any, cn.dW1_full.fillna(0)), 3))

# what remains for NGPS: pool objects WITHOUT any SDSS-V epoch
nosv = pool[~pool.name.isin(sv.name)]
print(f'\n=== discovery pool without any SDSS-V epoch: {len(nosv)} of {len(pool)}; with M>=1: {(nosv.M>=1).sum()}; with W1 fade>0.3: {(nosv.dW1_full>0.3).sum()}; both: {((nosv.M>=1)&(nosv.dW1_full>0.3)).sum()}; fade>0.6: {(nosv.dW1_full>0.6).sum()} ===')
c.to_csv(os.path.join(OUT, 'sdssv_calibration_pool.csv'), index=False)
nosv.to_csv(os.path.join(OUT, 'sdssv_calibration_pool_unobserved.csv'), index=False)
print('saved sdssv_calibration_pool.csv, sdssv_calibration_pool_unobserved.csv to scratchpad')
