"""
10_combined_rescore.py  --  would a ZTF + W1 manifold change the target lists?

Builds the combined manifold (ZTF g,r binned to 3 days + W1; manhattan metric) on Sample A, scores it on the held-out
Zeltyn objects exactly as 09_multiband_test.py does, then projects every enriched candidate that has both ZTF and W1
light curves (the 850-object pool subset + the 204 Zeltyn objects) and computes, in the combined space:
  knn_lit_frac      fraction of literature Turn-on/off among the 50 nearest Sample A neighbours
  zeltyn_density    local Zeltyn CL-AGN density relative to the whole manifold (1 = average)
  in_litregion      inside the 3x-enriched bins of the literature CLAGNs
  in_zeltynregion   inside the 3x-enriched bins of the Zeltyn CL-AGNs
Writes data/combined_scores.csv and prints how the current primary lists would move.
"""
import os, sys, glob, time, warnings
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE); sys.path.insert(0, os.path.join(HERE, 'code_src'))
warnings.filterwarnings('ignore')
import umap
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import roc_auc_score
from AGNzoo_functions import unify_lc_gp, stat_bands, combine_bands, normalize_clipmax_objects, update_bitsums
DATA = os.path.join(HERE, 'data'); t0 = time.time()
def say(*a): print(f'[{time.time()-t0:5.0f}s]', *a, flush=True)
BANDS = ['zg', 'zr', 'W1']; XRES = 80; ZTF_END = 60067.0


def bin_ztf(frame, days=3.0):
    d = frame.reset_index(); isz = d.band.isin(['zg', 'zr', 'zi'])
    z = d[isz].copy(); z['tb'] = np.floor(z.time / days) * days + days / 2
    zb = z.groupby(['objectid', 'label', 'band', 'tb']).agg(flux=('flux', 'median'), err=('err', lambda x: np.median(x) / np.sqrt(max(len(x), 1)))).reset_index().rename(columns={'tb': 'time'})
    out = pd.concat([d[~isz][['objectid', 'label', 'band', 'time', 'flux', 'err']], zb[['objectid', 'label', 'band', 'time', 'flux', 'err']]], ignore_index=True)
    return out.set_index(['objectid', 'label', 'band', 'time']).sort_index()


def features(frame):
    objs, dobjs, flabels, keeps = unify_lc_gp(frame, BANDS, xres=XRES, numplots=0, low_limit_size=5)
    fvar, maxarr, meanarr = stat_bands(objs, dobjs, BANDS, sigmacl=5)
    dat = normalize_clipmax_objects(combine_bands(objs, BANDS), maxarr, band=0)
    ok = np.isfinite(dat).all(axis=1)
    ids = np.asarray(frame.index.get_level_values('objectid').unique())[np.asarray(keeps)][ok]
    return dat[ok], np.asarray(flabels)[ok], ids


# ---------------- Sample A
df = update_bitsums(pd.read_parquet(os.path.join(DATA, 'df_lc_020724.parquet.gzip')).pipe(lambda d: d[d.index.get_level_values('label') != '64']))
dfA = bin_ztf(df); say('Sample A binned')
datA, labA, idsA = features(dfA); labA = labA.astype(int)
is_cl = ((labA & 16) > 0) | ((labA & 32) > 0); is_qso = labA == 1
mapp = umap.UMAP(n_neighbors=50, min_dist=0.99, metric='manhattan', random_state=3).fit(datA); emb = mapp.embedding_
say(f'combined manifold fitted on {len(emb)} Sample A objects')

# ---------------- candidates: pool subset (ZTF cache + W1 cache) and Zeltyn (ZTF cache + wise cache)
def ztf_rows(name, oid):
    p = [q for q in glob.glob(os.path.join(DATA, 'ztf_cache', '*', f'{name}.csv')) if os.path.getsize(q) > 10]
    if not p:
        return None
    z = pd.read_csv(p[0])
    if 'mjd' not in z or not len(z):
        return None
    z = z[(z.mjd <= ZTF_END) & np.isfinite(z.mag) & z.filtercode.isin(['zg', 'zr'])]
    fl = 3631e3 * 10 ** (-0.4 * z.mag.values)
    return pd.DataFrame({'objectid': oid, 'band': z.filtercode.values, 'time': z.mjd.values, 'flux': fl, 'err': fl * 0.921 * z.magerr.values})

zc = pd.read_csv(os.path.join(DATA, 'zeltyn_coords.csv')); sub = pd.read_csv(os.path.join(DATA, 'pool_subset_for_enrich.csv'))
names = {}; rows = []
wz = pd.concat([pd.read_parquet(f) for f in glob.glob(os.path.join(DATA, 'wise_cache', 'zeltyn', '*.parquet'))]).reset_index()
wz = wz[wz.band == 'WISE_W1']
for i, nm in enumerate(zc.name):
    oid = 100000 + i + 1; names[oid] = nm
    w = wz[wz.objectid == i + 1]
    if len(w):
        rows.append(pd.DataFrame({'objectid': oid, 'band': 'W1', 'time': w.time.values, 'flux': w.flux.values, 'err': w.err.values}))
    r = ztf_rows(nm, oid)
    if r is not None:
        rows.append(r)
wp = pd.concat([pd.read_parquet(f) for f in glob.glob(os.path.join(DATA, 'wise_cache', 'pool', '*.parquet'))]).reset_index()
wp = wp[wp.band == 'WISE_W1']
for r_ in sub.itertuples():
    oid = int(r_.poolid); names[oid] = r_.name
    w = wp[wp.objectid == oid]
    if len(w):
        rows.append(pd.DataFrame({'objectid': oid, 'band': 'W1', 'time': w.time.values, 'flux': w.flux.values, 'err': w.err.values}))
    r = ztf_rows(r_.name, oid)
    if r is not None:
        rows.append(r)
C = pd.concat(rows, ignore_index=True); C['label'] = 'cand'
C = bin_ztf(C.set_index(['objectid', 'label', 'band', 'time']).sort_index())
say(f'candidate frame: {C.index.get_level_values("objectid").nunique()} objects')
datC, _, idsC = features(C)
embC = mapp.transform(datC); say(f'projected {len(embC)} candidates')

# ---------------- scores in the combined space
nn = NearestNeighbors(n_neighbors=51).fit(emb)
knnA = is_cl[nn.kneighbors(emb)[1][:, 1:]].mean(axis=1); knnC = is_cl[nn.kneighbors(embC)[1][:, :50]].mean(axis=1)
zi = pd.read_csv(os.path.join(DATA, 'zeltyn2024_table5.csv'))
cand = pd.DataFrame({'objectid': idsC, 'name': [names[i] for i in idsC], 'ux': embC[:, 0], 'uy': embC[:, 1], 'knn_lit_frac': knnC})
cand['is_zeltyn'] = cand.objectid >= 100000
zcl_mask = cand.is_zeltyn.values & np.array([zi.loc[zi.Name == n, 'Class'].astype(str).str.startswith('CL-AGN').any() and float(zi.loc[zi.Name == n, 'z'].iloc[0]) < 1 for n in cand.name])
zev_mask = cand.is_zeltyn.values & ~zcl_mask & np.array([float(zi.loc[zi.Name == n, 'z'].iloc[0]) < 1 if (zi.Name == n).any() else False for n in cand.name])
Zc = embC[zcl_mask]
U = np.vstack([emb, Zc]); isZ = np.r_[np.zeros(len(emb), bool), np.ones(len(Zc), bool)]
nnU = NearestNeighbors(n_neighbors=50).fit(U)
cand['zeltyn_density'] = isZ[nnU.kneighbors(embC)[1]].mean(axis=1) / (len(Zc) / len(U))
hist_all, xe, ye = np.histogram2d(emb[:, 0], emb[:, 1], bins=10)
def bins(pts): return (np.clip(np.searchsorted(xe[1:], pts[:, 0]), 0, 9), np.clip(np.searchsorted(ye[1:], pts[:, 1]), 0, 9))
def region_from(pts):
    h, _, _ = np.histogram2d(pts[:, 0], pts[:, 1], bins=(xe, ye))
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(hist_all > 0, h / hist_all, 0) > 3 * len(pts) / len(emb)
def inside(reg, pts): ix, iy = bins(pts); return reg[ix, iy]
reg_lit = region_from(emb[is_cl]); reg_zel = region_from(Zc)
cand['in_litregion'] = inside(reg_lit, embC); cand['in_zeltynregion'] = inside(reg_zel, embC)
# quality of this space on the held-out Zeltyn objects (same metrics as 09)
auc_cl = roc_auc_score(np.r_[np.ones(zcl_mask.sum()), np.zeros(is_qso.sum())], np.r_[knnC[zcl_mask], knnA[is_qso]])
auc_ev = roc_auc_score(np.r_[np.ones(zev_mask.sum()), np.zeros(is_qso.sum())], np.r_[knnC[zev_mask], knnA[is_qso]])
rng = np.random.default_rng(0); idx = np.where(zcl_mask)[0]; enr = []
for _ in range(20):
    rng.shuffle(idx); a, b = idx[:len(idx) // 2], idx[len(idx) // 2:]
    reg = region_from(embC[a]); enr.append(inside(reg, embC[b]).mean() / max(inside(reg, emb[is_qso]).mean(), 1e-9))
say(f'ZTF g,r + W1 (manhattan): AUC CL-AGN vs QSO {auc_cl:.3f}, EVQ {auc_ev:.3f}; held-out Zeltyn-region enrichment {np.mean(enr):.1f}x; '
    f'lit-region: Zeltyn inside {inside(reg_lit, Zc).mean():.3f} vs QSO {inside(reg_lit, emb[is_qso]).mean():.3f}')
cand.to_csv(os.path.join(DATA, 'combined_scores.csv'), index=False)

# ---------------- what would move in the current lists
m = pd.read_csv(os.path.join(DATA, 'master_list_scored.csv'), low_memory=False).set_index('name')
cand = cand.set_index('name').join(m[['tier', 'M', 'P', 'priority', 'r_mag', 'z', 'clagn_score', 'zeltyn_density_ratio', 'in_region_zeltyn', 'in_region_clagn']], how='left')
cand['M_combined'] = np.clip(np.fmax(cand.zeltyn_density / 3.0, cand.knn_lit_frac / 0.15), 0, 2)
prim = pd.concat([pd.read_csv(os.path.join(DATA, f'targets_{n}.csv')) for n in ['sep23', 'oct26', 'oct27']]); prim = prim[prim['rank'] > 0]
p = cand[cand.index.isin(prim.name)]; np_ = cand[~cand.index.isin(prim.name) & (cand.tier == 'T1')]
say(f'current primaries scored in the combined space: {len(p)} | with M_combined >= 1: {(p.M_combined >= 1).sum()} | in combined lit-region: {p.in_litregion.sum()} | in combined Zeltyn-region: {p.in_zeltynregion.sum()}')
say(f'Tier-1 primaries: median M (W1 space) {p[p.tier=="T1"].M.median():.2f} vs M_combined {p[p.tier=="T1"].M_combined.median():.2f}; '
    f'T1 primaries with M_combined < 0.5 (would be dropped by a combined-space cut): {int(((p.tier=="T1") & (p.M_combined < 0.5)).sum())}')
say(f'non-selected enriched T1 candidates with M_combined >= 1 and W1-space M >= 1: {int(((np_.M_combined >= 1) & (np_.M >= 1)).sum())} of {len(np_)}; '
    f'with M_combined >= 1.5 & r < 18.5: {int(((np_.M_combined >= 1.5) & (np_.r_mag < 18.5)).sum())}')
say('rank correlation (Spearman) between W1-space M and combined-space M over all scored candidates: %.2f' % cand[['M', 'M_combined']].corr(method='spearman').iloc[0, 1])
print(p[p.tier == 'T1'].sort_values('M_combined')[['tier', 'z', 'r_mag', 'M', 'M_combined', 'knn_lit_frac', 'zeltyn_density', 'in_litregion', 'in_zeltynregion']].round(2).head(12).to_string())
