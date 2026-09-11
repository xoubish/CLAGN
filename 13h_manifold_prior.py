"""
13h_manifold_prior.py  --  the Hemmati+2026 manifold as a CALIBRATED prior and DIRECTION predictor for every pool object.
  P_dim, P_bright, P_hb   local (k nearest calibration objects, k=150) rates of continuum dimming > 0.7 mag, brightening > 0.7 mag and
                          significant Hbeta halving between the DR16 and SDSS-V epochs, at the object's manifold position
  lit_corner              inside the literature-CLAGN region (4-5x enrichment of independently confirmed CLAGNs)
  f_on, f_off, dir        fraction of Sample A turn-on / turn-off objects among the 50 nearest training objects; dir = f_on - f_off
  known_clagn             matched (2") to a spectroscopically confirmed CLAGN in the Camus & Panda 2026 database
Leave-one-out for calibration objects.  Output: data/manifold_prior_pool.csv
"""
import os, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from astropy.coordinates import SkyCoord; import astropy.units as u
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, 'data'); K_CAL, K_DIR = 150, 50
pool = pd.read_csv(os.path.join(DATA, 'parent_pool_scored.csv')); pool = pool[pool.projected.fillna(False)].copy(); pool['name'] = 'P' + pool.poolid.astype(str)
c = pd.read_csv(os.path.join(DATA, 'sdssv_calibration_pool.csv')).set_index('name'); ls = pd.read_csv(os.path.join(DATA, 'sdssv_calibration_lines_summary.csv'), index_col=0)
c = c.join(ls[['hb_dim']], how='left'); c['hb_dim'] = c.hb_dim.fillna(False).astype(bool)
XY = c[['umap_x', 'umap_y']].values; tree = cKDTree(XY); _, idx = tree.query(pool[['umap_x', 'umap_y']].values, k=K_CAL + 1)
is_cal = pool.name.isin(c.index).values
# leave-one-out: for calibration objects drop the neighbour that is the object itself (distance 0 -> first column)
sel = np.where(is_cal[:, None], idx[:, 1:], idx[:, :K_CAL])
for col, out in [('dimmed_any', 'P_dim'), ('bright07', 'P_bright'), ('hb_dim', 'P_hb')]:
    y = c[col].astype(float).values; pool[out] = y[sel].mean(axis=1)
A = pd.read_csv(os.path.join(DATA, 'sampleA_embedding_objectid.csv')); on = A.labels.str.contains('Turn-on', na=False).values; off = A.labels.str.contains('Turn-off', na=False).values
_, ia = cKDTree(A[['umap_x', 'umap_y']].values).query(pool[['umap_x', 'umap_y']].values, k=K_DIR)
pool['f_on'] = on[ia].mean(1); pool['f_off'] = off[ia].mean(1); pool['dir'] = pool.f_on - pool.f_off; pool['lit_corner'] = pool.in_region_clagn.astype(bool)
cat = pd.read_csv(os.path.join(DATA, 'external', 'clagn_catalog_camus_panda2026.csv'), low_memory=False); cat = cat[np.isfinite(cat.ra_deg) & np.isfinite(cat.dec_deg) & (cat.confirmation_status == 'spectroscopic_confirmed')]
i, s, _ = SkyCoord(pool.ra.values * u.deg, pool.dec.values * u.deg).match_to_catalog_sky(SkyCoord(cat.ra_deg.values * u.deg, cat.dec_deg.values * u.deg))
pool['known_clagn'] = s.arcsec < 2; pool['known_type'] = np.where(pool.known_clagn, cat.transition_type.values[i], '')
out = pool[['name', 'poolid', 'umap_x', 'umap_y', 'P_dim', 'P_bright', 'P_hb', 'lit_corner', 'in_region_zeltyn', 'clagn_score', 'f_on', 'f_off', 'dir', 'known_clagn', 'known_type']]
out.to_csv(os.path.join(DATA, 'manifold_prior_pool.csv'), index=False)
print(f'{len(out)} pool objects; P_dim 10/50/90%: {out.P_dim.quantile(.1):.3f}/{out.P_dim.median():.3f}/{out.P_dim.quantile(.9):.3f}; P_hb: {out.P_hb.quantile(.1):.3f}/{out.P_hb.median():.3f}/{out.P_hb.quantile(.9):.3f}; lit_corner {out.lit_corner.sum()}; known CLAGN {out.known_clagn.sum()}')
# sanity: leave-one-out P_dim vs actual outcome on the calibration set
from sklearn.metrics import roc_auc_score
cc = out[out.name.isin(c.index)].set_index('name').join(c[['dimmed_any', 'bright07', 'hb_dim']])
print(f'leave-one-out check on the calibration set: AUC(dimmed | P_dim) {roc_auc_score(cc.dimmed_any, cc.P_dim):.3f}, AUC(brightened | P_bright) {roc_auc_score(cc.bright07, cc.P_bright):.3f}, AUC(Hb fell | P_hb) {roc_auc_score(cc.hb_dim, cc.P_hb):.3f}')
