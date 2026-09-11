"""
14a_calibrated_models.py  --  calibrated change probabilities for every pool quasar, fitted on the SDSS-V two-epoch set.

Features (per object): optical change since the archival spectrum d_g = ZTF DR24 median g - synthetic g of the DR16 spectrum (13f/13g),
W1 fade 2010-20 and amplitude (unWISE), NEOWISE W1 change 2014-24, W1-W2 colour change, log Eddington ratio (Wu & Shen 2022),
the manifold prior surface (13h: P_dim / P_hb / P_bright, leave-one-out), redshift.  Missing values -> training median + indicator.
Models (5-fold CV AUC printed): continuum dimmed > 0.7 mag (all calibration objects); Hbeta fell > 2x at 3 sigma (2"-fibre pairs only,
the aperture-clean subset); continuum brightened; Hbeta rose > 2x.  Output: data/pool_probabilities.csv, data/unwise_features_pool.csv (cache).
"""
import os, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import roc_auc_score
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, 'data'); pd.set_option('display.width', 250)
pool = pd.read_csv(os.path.join(DATA, 'parent_pool_scored.csv')); pool = pool[pool.projected.fillna(False)].copy(); pool['name'] = 'P' + pool.poolid.astype(str); pool = pool.set_index('name')
# ---- unWISE 2010-20 features for the whole pool (cached)
uw = os.path.join(DATA, 'unwise_features_pool.csv')
if not os.path.exists(uw):
    cd = os.path.join(DATA, 'wise_cache', 'pool')
    lc = pd.concat([pd.read_parquet(os.path.join(cd, f)) for f in sorted(os.listdir(cd)) if f.endswith('.parquet')])
    lc = lc[lc.index.get_level_values('band') == 'WISE_W1'].reset_index()[['objectid', 'time', 'flux']].sort_values(['objectid', 'time']); g = lc.groupby('objectid')
    f = pd.DataFrame({'n_w1': g.size(), 'w1_first': g.flux.apply(lambda s: s.head(2).median()), 'w1_last': g.flux.apply(lambda s: s.tail(3).median()), 'w1_max': g.flux.max(), 'w1_min': g.flux.min(), 'w1_med': g.flux.median(), 'w1_std': g.flux.std()})
    f['dW1_full'] = -2.5 * np.log10(f.w1_last / f.w1_first); f['ampW1'] = 2.5 * np.log10(f.w1_max / f.w1_min); f['fvarW1'] = f.w1_std / f.w1_med
    f.index = 'P' + f.index.astype(str); f.index.name = 'name'; f.to_csv(uw)
f = pd.read_csv(uw).set_index('name')
oc = pd.read_csv(os.path.join(DATA, 'optical_change_pool.csv')).set_index('name')
neo = pd.concat([pd.read_csv(os.path.join(DATA, 'neowise_now_pool.csv')), pd.read_csv(os.path.join(DATA, 'neowise_now_poolall.csv'))]).drop_duplicates('name').set_index('name')
col = pd.read_csv(os.path.join(DATA, 'neowise_color_pool.csv'), index_col=0); ed = pd.read_csv(os.path.join(DATA, 'external', 'dr16q_prop_pool.csv')).set_index('name')
mp = pd.read_csv(os.path.join(DATA, 'manifold_prior_pool.csv')).set_index('name')
X = pool[['ra', 'dec', 'z', 'psfmag_r', 'plate', 'mjd']].join(f[['dW1_full', 'ampW1', 'fvarW1']]).join(oc[['d_g', 'd_r', 'd_opt', 'd_opt_flag', 'ztf_g_amp', 'ztf_r_amp', 'ztf_r_medianmag', 'ztf_g_medianmag', 'syn_g', 'syn_r']]) \
    .join(neo[['dw1_neowise', 'w1_slope_2yr', 'w1_amp_visits', 'w1_flux_last_mjy', 'mjd_last_neo']]).join(col[['dcolor', 'c_slope']]).join(ed[['LOGLEDD_RATIO', 'LOGMBH', 'LOGLBOL']]) \
    .join(mp[['P_dim', 'P_bright', 'P_hb', 'lit_corner', 'dir', 'f_on', 'f_off', 'known_clagn', 'known_type']])
X['dw1_best'] = X.dw1_neowise.where(X.dw1_neowise.notna(), X.dW1_full)
X['ztf_amp'] = X.ztf_g_amp.where(X.ztf_g_amp.notna(), X.ztf_r_amp).clip(upper=4)      # ZTF DR24 max - min over 2018-25 (g, else r); adds ~0.01 AUC for dimming
FEATS = ['d_opt', 'dW1_full', 'ampW1', 'dw1_neowise', 'dcolor', 'LOGLEDD_RATIO', 'z', 'ztf_amp']
def design(df, prior_col, med=None):
    D = df[FEATS + [prior_col]].copy().rename(columns={prior_col: 'P_manifold'})
    ind = D[['d_opt', 'dw1_neowise', 'LOGLEDD_RATIO', 'ztf_amp']].isna().astype(float).add_prefix('miss_')
    med = D.median() if med is None else med
    D = D.fillna(med); D = pd.concat([D, ind], axis=1); return D, med
c = pd.read_csv(os.path.join(DATA, 'sdssv_calibration_pool.csv')).set_index('name'); ls = pd.read_csv(os.path.join(DATA, 'sdssv_calibration_lines_summary.csv'), index_col=0)
cal = X.join(c[['dimmed_any', 'bright07']], how='inner').join(ls[['hb_dim', 'hb_bright', 'hb_dim3']], how='left')
for k in ['hb_dim', 'hb_bright', 'hb_dim3']: cal[k] = cal[k].fillna(False).astype(bool)
clean = cal.plate >= 3500
print(f'pool {len(X)}; calibration {len(cal)} (clean 2\" pairs {clean.sum()}); feature coverage: ' + ', '.join(f'{k} {100*X[k].notna().mean():.0f}%' for k in FEATS))
models = {}
for name, target, subset, prior in [('P_cont_dim', 'dimmed_any', cal, 'P_dim'), ('P_hb_dim', 'hb_dim', cal[clean], 'P_hb'), ('P_cont_bright', 'bright07', cal, 'P_bright'), ('P_hb_bright', 'hb_bright', cal[clean], 'P_bright')]:
    D, med = design(subset, prior); y = subset[target].astype(int).values; mu, sd = D.mean(), D.std().replace(0, 1); Z = ((D - mu) / sd).values
    lr = LogisticRegression(max_iter=5000, C=0.5); cvp = cross_val_predict(lr, Z, y, cv=StratifiedKFold(5, shuffle=True, random_state=0), method='predict_proba')[:, 1]
    lr.fit(Z, y); models[name] = (lr, med, mu, sd, prior)
    coef = pd.Series(lr.coef_[0], index=D.columns).round(2)
    print(f'\n{name}: N={len(y)}, events={y.sum()}, CV AUC {roc_auc_score(y, cvp):.3f}; standardised coefficients: {coef.to_dict()}')
    q = pd.qcut(cvp, 10, labels=False, duplicates='drop'); t = pd.DataFrame({'pred': cvp, 'obs': y}).groupby(q).agg(N=('obs', 'size'), predicted=('pred', 'mean'), observed=('obs', 'mean'))
    print('   calibration by decile (top 3): ' + '; '.join(f'{100*r.predicted:.1f}% pred vs {100*r.observed:.1f}% obs (N={int(r.N)})' for _, r in t.tail(3).iterrows()))
    # apply to the whole pool
    Dp, _ = design(X, prior, med); Zp = ((Dp - mu) / sd).values; X[name] = lr.predict_proba(Zp)[:, 1]
out_cols = ['ra', 'dec', 'z', 'psfmag_r', 'ztf_r_medianmag', 'ztf_g_medianmag', 'd_g', 'd_r', 'd_opt', 'd_opt_flag', 'ztf_amp', 'syn_g', 'syn_r', 'dW1_full', 'ampW1', 'dw1_neowise', 'w1_slope_2yr', 'w1_amp_visits', 'dcolor', 'LOGLEDD_RATIO', 'LOGMBH',
            'P_dim', 'P_bright', 'P_hb', 'lit_corner', 'dir', 'f_on', 'f_off', 'known_clagn', 'known_type', 'P_cont_dim', 'P_hb_dim', 'P_cont_bright', 'P_hb_bright']
X[out_cols].to_csv(os.path.join(DATA, 'pool_probabilities.csv'))
disc = X[~X.index.isin(pd.read_csv(os.path.join(DATA, 'sdssv_internal_epochs.csv')).name) & ~X.known_clagn.fillna(False).astype(bool)]
print(f'\ndiscovery pool (no SDSS-V epoch, not a known CLAGN): {len(disc)}; P_hb_dim > 0.3: {(disc.P_hb_dim > 0.3).sum()}, > 0.2: {(disc.P_hb_dim > 0.2).sum()}, > 0.1: {(disc.P_hb_dim > 0.1).sum()}; P_hb_bright > 0.2: {(disc.P_hb_bright > 0.2).sum()}, > 0.1: {(disc.P_hb_bright > 0.1).sum()}')
print('top 10 turn-off candidates:'); print(disc.sort_values('P_hb_dim', ascending=False).head(10)[['z', 'psfmag_r', 'ztf_r_medianmag', 'd_g', 'dW1_full', 'dw1_neowise', 'dcolor', 'LOGLEDD_RATIO', 'P_hb', 'lit_corner', 'P_cont_dim', 'P_hb_dim']].round(2).to_string())
print('top 10 brightening candidates:'); print(disc.sort_values('P_hb_bright', ascending=False).head(10)[['z', 'psfmag_r', 'ztf_r_medianmag', 'd_g', 'dW1_full', 'dw1_neowise', 'dcolor', 'LOGLEDD_RATIO', 'P_bright', 'dir', 'P_cont_bright', 'P_hb_bright']].round(2).to_string())
