"""
09_multiband_test.py  --  is a multi-band (ZTF + WISE) manifold a better CLAGN selector than W1 alone?

For each band set, rebuild the Sample A manifold exactly as the pipeline does (GP on a common grid, normalise by the
sigma-clipped max of the first band, UMAP), project the held-out Zeltyn+2024 CL-AGNs and EVQs, and score:
  * AUC of the kNN CLAGN score (fraction of literature Turn-on/off among 50 nearest Sample A neighbours) for
    Zeltyn CL-AGN vs Sample A SDSS_QSO-only objects, and EVQ vs SDSS_QSO
  * held-out region test: define the 3x-enriched bins from a random half of the Zeltyn CL-AGNs, measure the fraction
    of the other half inside vs the fraction of SDSS_QSO inside (enrichment), averaged over 20 splits
  * literature-CLAGN-defined region: fraction of Zeltyn CL-AGN inside vs SDSS_QSO inside
Band sets: W1 (DTW, the pipeline baseline); W1+W2 (DTW on the concatenation); ZTF g,r (manhattan);
ZTF g,r + W1,W2 (manhattan, as in the paper's combined manifold). ZTF is cut at MJD 60067 to match Sample A.
Writes data/multiband_test.csv and prints the table.  ~10-15 min.
"""
import os, sys, glob, time, warnings
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE); sys.path.insert(0, os.path.join(HERE, 'code_src'))
warnings.filterwarnings('ignore')
import umap
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import roc_auc_score
from AGNzoo_functions import unify_lc_gp, stat_bands, combine_bands, normalize_clipmax_objects, dtw_distance, update_bitsums

DATA = os.path.join(HERE, 'data'); t0 = time.time()
def say(*a): print(f'[{time.time()-t0:5.0f}s]', *a, flush=True)

# ---------------- Sample A
df = pd.read_parquet(os.path.join(DATA, 'df_lc_020724.parquet.gzip'))
df = update_bitsums(df[df.index.get_level_values('label') != '64'])
ZTF_END = 60067.0
# ---------------- Zeltyn multi-band frame in the same MultiIndex format (objectid, label, band, time) with flux in mJy
zc = pd.read_csv(os.path.join(DATA, 'zeltyn_coords.csv'))
wise = pd.concat([pd.read_parquet(f) for f in glob.glob(os.path.join(DATA, 'wise_cache', 'zeltyn', '*.parquet'))]).reset_index()
wise['band'] = wise['band'].str.replace('WISE_', '', regex=False)
rows = [wise[['objectid', 'band', 'time', 'flux', 'err']]]
for i, nm in enumerate(zc.name):
    p = os.path.join(DATA, 'ztf_cache', 'zeltyn', f'{nm}.csv')
    if not os.path.exists(p) or os.path.getsize(p) < 10:
        continue
    z = pd.read_csv(p)
    if 'mjd' not in z or not len(z):
        continue
    z = z[(z.mjd <= ZTF_END) & np.isfinite(z.mag)]
    fl = 3631e3 * 10 ** (-0.4 * z.mag.values); er = fl * 0.921 * z.magerr.values
    rows.append(pd.DataFrame({'objectid': i + 1, 'band': z.filtercode.values, 'time': z.mjd.values, 'flux': fl, 'err': er}))
Z = pd.concat(rows, ignore_index=True); Z['label'] = 'Zeltyn'
Z = Z.set_index(['objectid', 'label', 'band', 'time']).sort_index()
zinfo = pd.read_csv(os.path.join(DATA, 'zeltyn2024_table5.csv'))
is_cl_z = zinfo['Class'].str.startswith('CL-AGN').values; zz = zinfo['z'].values
say('Sample A objects', df.index.get_level_values('objectid').nunique(), '| Zeltyn rows', len(Z))

def bin_ztf(frame, days=3.0):
    """Median-bin the ZTF bands (zg, zr, zi) to `days`-wide bins; other bands untouched. Cuts the GP cost ~100x."""
    d = frame.reset_index()
    isz = d.band.isin(['zg', 'zr', 'zi'])
    z = d[isz].copy(); z['tb'] = np.floor(z.time / days) * days + days / 2
    zb = z.groupby(['objectid', 'label', 'band', 'tb']).agg(flux=('flux', 'median'), err=('err', lambda x: np.median(x) / np.sqrt(max(len(x), 1)))).reset_index().rename(columns={'tb': 'time'})
    out = pd.concat([d[~isz][['objectid', 'label', 'band', 'time', 'flux', 'err']], zb[['objectid', 'label', 'band', 'time', 'flux', 'err']]], ignore_index=True)
    return out.set_index(['objectid', 'label', 'band', 'time']).sort_index()

_binned = {}
def build(bands, metric, xres):
    global df, Z
    if any(b.startswith('z') for b in bands):
        if 'A' not in _binned:
            _binned['A'] = bin_ztf(df); _binned['Z'] = bin_ztf(Z); say('ZTF bands binned to 3 days:', len(_binned['A']), 'Sample A rows,', len(_binned['Z']), 'Zeltyn rows')
        dfA, dfZ = _binned['A'], _binned['Z']
    else:
        dfA, dfZ = df, Z
    objs, dobjs, flabels, keeps = unify_lc_gp(dfA, bands, xres=xres, numplots=0, low_limit_size=5)
    fvar, maxarr, meanarr = stat_bands(objs, dobjs, bands, sigmacl=5)
    dat = normalize_clipmax_objects(combine_bands(objs, bands), maxarr, band=0)
    ok = np.isfinite(dat).all(axis=1)
    lab = np.asarray(flabels)[ok].astype(int); dat = dat[ok]
    mapp = umap.UMAP(n_neighbors=100 if metric == 'dtw' else 50, min_dist=0.99, metric=(dtw_distance if metric == 'dtw' else metric), random_state=3).fit(dat)
    emb = mapp.embedding_
    # Zeltyn
    zo, zd, zl, zk = unify_lc_gp(dfZ, bands, xres=xres, numplots=0, low_limit_size=5)
    zf, zmax, zmean = stat_bands(zo, zd, bands, sigmacl=5)
    zdat = normalize_clipmax_objects(combine_bands(zo, bands), zmax, band=0)
    zok = np.isfinite(zdat).all(axis=1)
    zids = np.asarray(dfZ.index.get_level_values('objectid').unique())[np.asarray(zk)][zok] - 1
    zemb = mapp.transform(zdat[zok])
    return emb, lab, zemb, zids

def evaluate(emb, lab, zemb, zids, tag):
    is_cl = ((lab & 16) > 0) | ((lab & 32) > 0); is_qso = lab == 1
    nn = NearestNeighbors(n_neighbors=51).fit(emb)
    sA = is_cl[nn.kneighbors(emb)[1][:, 1:]].mean(axis=1)
    sZ = is_cl[nn.kneighbors(zemb)[1][:, :50]].mean(axis=1)
    zcl = is_cl_z[zids] & (zz[zids] < 1); zev = ~is_cl_z[zids] & (zz[zids] < 1)
    auc_cl = roc_auc_score(np.r_[np.ones(zcl.sum()), np.zeros(is_qso.sum())], np.r_[sZ[zcl], sA[is_qso]])
    auc_ev = roc_auc_score(np.r_[np.ones(zev.sum()), np.zeros(is_qso.sum())], np.r_[sZ[zev], sA[is_qso]])
    # held-out region test (3x enriched bins, 10x10 grid) with random halves of the Zeltyn CL-AGNs
    hist_all, xe, ye = np.histogram2d(emb[:, 0], emb[:, 1], bins=10)
    def bins(pts): return (np.clip(np.searchsorted(xe[1:], pts[:, 0]), 0, 9), np.clip(np.searchsorted(ye[1:], pts[:, 1]), 0, 9))
    def region_from(pts):
        h, _, _ = np.histogram2d(pts[:, 0], pts[:, 1], bins=(xe, ye))
        with np.errstate(invalid='ignore', divide='ignore'):
            return np.where(hist_all > 0, h / hist_all, 0) > 3 * len(pts) / len(emb)
    def inside(reg, pts): ix, iy = bins(pts); return reg[ix, iy]
    rng = np.random.default_rng(0); idx = np.where(zcl)[0]; enr, fin, fq = [], [], []
    for _ in range(20):
        rng.shuffle(idx); a, b = idx[:len(idx) // 2], idx[len(idx) // 2:]
        reg = region_from(zemb[a]); f_b = inside(reg, zemb[b]).mean(); f_q = inside(reg, emb[is_qso]).mean()
        fin.append(f_b); fq.append(f_q); enr.append(f_b / max(f_q, 1e-9))
    reg_lit = region_from(emb[is_cl]); f_zl = inside(reg_lit, zemb[zcl]).mean(); f_ql = inside(reg_lit, emb[is_qso]).mean()
    r = dict(bands=tag, nA=len(emb), nZ_cl=int(zcl.sum()), nZ_evq=int(zev.sum()), auc_clagn_vs_qso=round(auc_cl, 3), auc_evq_vs_qso=round(auc_ev, 3),
             heldout_frac_inside=round(np.mean(fin), 3), qso_frac_inside=round(np.mean(fq), 3), heldout_enrichment=round(np.mean(enr), 2),
             litregion_zeltyn_inside=round(f_zl, 3), litregion_qso_inside=round(f_ql, 3), litregion_enrichment=round(f_zl / max(f_ql, 1e-9), 2),
             knn_median_zeltyn_cl=round(float(np.median(sZ[zcl])), 3), knn_median_qso=round(float(np.median(sA[is_qso])), 3))
    say(tag, {k: v for k, v in r.items() if k not in ('bands',)})
    return r

OUT = os.path.join(DATA, 'multiband_test.csv')
results = pd.read_csv(OUT).to_dict('records') if os.path.exists(OUT) else []
done = {r['bands'] for r in results}
for bands, metric, xres, tag in [(['W1'], 'dtw', 160, 'W1 (DTW)'), (['W1', 'W2'], 'dtw', 80, 'W1+W2 (DTW)'),
                                 (['zg', 'zr'], 'manhattan', 80, 'ZTF g,r (manhattan)'), (['zg', 'zr', 'W1', 'W2'], 'manhattan', 80, 'ZTF g,r + W1,W2 (manhattan)'),
                                 (['zg', 'zr', 'W1', 'W2'], 'dtw', 60, 'ZTF g,r + W1,W2 (DTW)'), (['zg', 'zr', 'W1'], 'dtw', 60, 'ZTF g,r + W1 (DTW)')]:
    if tag in done:
        continue
    try:
        emb, lab, zemb, zids = build(bands, metric, xres)
        results.append(evaluate(emb, lab, zemb, zids, tag))
        pd.DataFrame(results).to_csv(OUT, index=False)      # save after every band set
    except Exception as e:
        say(tag, 'FAILED', type(e).__name__, str(e)[:120])
R = pd.DataFrame(results)
print('\n' + R.to_string(index=False))
