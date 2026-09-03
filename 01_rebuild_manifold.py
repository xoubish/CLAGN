"""
01_rebuild_manifold.py  --  Stage 0/1 of FOLLOWUP_PLAN.md

Re-runs the WISE-W1 manifold of wise_manifold_zeltyn.ipynb as a script, but keeps the
objectid attached to every embedding row and saves everything needed downstream:

  data/sampleA_embedding_objectid.csv   objectid, label_bits, labels, umap_x, umap_y, fvar_w1, mean_w1_mjy, max_w1_mjy,
                                        clagn_score (kNN fraction of Turn-on/Turn-off), in_region_clagn, in_region_zeltyn
  data/zeltyn_embedding_full.csv        Zeltyn+2024 CL-AGNs and EVQs (all 204 rows; those with unWISE data projected):
                                        name, class, z, ra, dec, umap_x, umap_y, clagn_score, in_region_*
  data/umap_w1_model.pkl                fitted UMAP (for mapp.transform of any new pool)
  data/region_stats.txt                 counts for both region definitions + held-out test numbers
  data/wise_cache/zeltyn/chunk_*.parquet   cached unWISE W1/W2 light curves (resumable)

Run from the CLAGN folder with /opt/anaconda3/bin/python.  Same preprocessing as the notebook:
W1 only, GP (RBF length_scale=200) on a 160-point grid, normalise by sigma-clipped max, drop NaN rows,
UMAP(n_neighbors=100, min_dist=0.99, metric=dtw_distance, random_state=3).  No shuffling (order = objectid).
"""
import os, re, sys, pickle, time, warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, os.path.join(HERE, 'code_src'))
warnings.filterwarnings('ignore')

import umap
import pyarrow.fs as pafs
from astropy.coordinates import SkyCoord
import astropy.units as u
from sklearn.neighbors import NearestNeighbors

from data_structures import MultiIndexDFObject
from wise_functions import wise_get_lightcurves
from sample_selection import clean_sample
from AGNzoo_functions import (unify_lc_gp, stat_bands, combine_bands, normalize_clipmax_objects,
                              dtw_distance, translate_bitwise_sum_to_labels, update_bitsums)

# wise_functions builds an anonymous S3FileSystem with default (short) timeouts; lengthen them.
_orig_s3fs = pafs.S3FileSystem
def _patient_s3fs(*a, **k):
    k.setdefault('request_timeout', 180)
    k.setdefault('connect_timeout', 60)
    return _orig_s3fs(*a, **k)
pafs.S3FileSystem = _patient_s3fs

DATA = os.path.join(HERE, 'data')
PARQUET = os.path.join(DATA, 'df_lc_020724.parquet.gzip')
ZELTYN_CSV = os.path.join(DATA, 'zeltyn2024_table5.csv')
CACHE = os.path.join(DATA, 'wise_cache', 'zeltyn')
BANDS = ['W1']
XRES = 160
K_NN = 50          # neighbours for clagn_score
ENRICH = 3         # enrichment factor for the binned region definitions
NBINS = 10
t0 = time.time()
log = open(os.path.join(DATA, 'region_stats.txt'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); log.write(s + '\n'); log.flush()


def preprocess(df, tag):
    objs, dobjs, flabels, keeps = unify_lc_gp(df, BANDS, xres=XRES, numplots=0, low_limit_size=5)
    fvar, maxarr, meanarr = stat_bands(objs, dobjs, BANDS, sigmacl=5)
    dat = normalize_clipmax_objects(combine_bands(objs, BANDS), maxarr, band=0)
    valid = ~np.isnan(dat).any(axis=1) & np.isfinite(dat).all(axis=1)
    ids = np.asarray(df.index.get_level_values('objectid').unique())[np.asarray(keeps)]
    say(f'[{time.time()-t0:6.0f}s] {tag}: {len(objs)} survive GP cut, {(~valid).sum()} dropped for NaN')
    return dat[valid], ids[valid], np.asarray(flabels)[valid], fvar[0][valid], meanarr[0][valid], maxarr[0][valid]


def fetch_wise_chunked(sample_table, cache_dir, chunk=20, tries=5):
    """unWISE fetch in cached, retried chunks so a network hiccup never loses progress."""
    os.makedirs(cache_dir, exist_ok=True)
    parts = []
    for i in range(0, len(sample_table), chunk):
        f = os.path.join(cache_dir, f'chunk_{i:04d}.parquet')
        if os.path.exists(f):
            parts.append(pd.read_parquet(f)); continue
        sub = sample_table[i:i + chunk]
        for attempt in range(tries):
            try:
                d = wise_get_lightcurves(sub, radius=1.0, bandlist=['WISE_W1', 'WISE_W2']).data
                break
            except ValueError:            # no detections at all in this chunk
                d = pd.DataFrame(); break
            except OSError as e:
                say(f'   chunk {i}: {str(e)[:80]} ... retry {attempt+1}/{tries}')
                time.sleep(15 * (attempt + 1))
        else:
            raise RuntimeError(f'chunk {i} failed after {tries} tries')
        d.to_parquet(f); parts.append(d)
        say(f'[{time.time()-t0:6.0f}s]   fetched chunk {i//chunk+1}/{(len(sample_table)-1)//chunk+1}')
    parts = [p for p in parts if len(p)]
    return pd.concat(parts) if parts else pd.DataFrame()


# ============================================================================= Sample A
df_lc = pd.read_parquet(PARQUET)
df_lc = df_lc[df_lc.index.get_level_values('label') != '64']
df_lc = update_bitsums(df_lc)
say(f'[{time.time()-t0:6.0f}s] Sample A objects after removing SPIDER-only: {df_lc.index.get_level_values("objectid").nunique()}')

datA, idsA, labA, fvarA, meanA, maxA = preprocess(df_lc, 'Sample A')
labA = labA.astype(int)
is_on = (labA & 16) > 0
is_off = (labA & 32) > 0
is_cl = is_on | is_off
say(f'Turn-on {is_on.sum()}, Turn-off {is_off.sum()}, either {is_cl.sum()}')

mapp = umap.UMAP(n_neighbors=100, min_dist=0.99, metric=dtw_distance, random_state=3).fit(datA)
emb = mapp.embedding_
say(f'[{time.time()-t0:6.0f}s] UMAP fitted on {len(emb)} objects')
with open(os.path.join(DATA, 'umap_w1_model.pkl'), 'wb') as f:
    pickle.dump(mapp, f)

# --- Sample-A-only scores and the Turn-on/off-defined region; write immediately
nn = NearestNeighbors(n_neighbors=K_NN + 1).fit(emb)
def clagn_score(points, exclude_self=False):
    _, idx = nn.kneighbors(points)
    idx = idx[:, 1:] if exclude_self else idx[:, :K_NN]
    return is_cl[idx].mean(axis=1)
scoreA = clagn_score(emb, exclude_self=True)

hist_all, xe, ye = np.histogram2d(emb[:, 0], emb[:, 1], bins=NBINS)
def bin_idx(pts):
    ix = np.clip(np.searchsorted(xe[1:], pts[:, 0]), 0, NBINS - 1)
    iy = np.clip(np.searchsorted(ye[1:], pts[:, 1]), 0, NBINS - 1)
    return ix, iy
def region_from(pts):
    h, _, _ = np.histogram2d(pts[:, 0], pts[:, 1], bins=(xe, ye))
    with np.errstate(invalid='ignore', divide='ignore'):
        frac = np.where(hist_all > 0, h / hist_all, 0)
    return frac > ENRICH * len(pts) / len(emb)
def in_region(reg, pts):
    ix, iy = bin_idx(pts); return reg[ix, iy]
reg_clagn = region_from(emb[is_cl])        # Sample A Turn-on/Turn-off enriched (plan Section 8 alternative)

A = pd.DataFrame({'objectid': idsA, 'label_bits': labA,
                  'labels': ['+'.join(translate_bitwise_sum_to_labels(int(l))) for l in labA],
                  'umap_x': emb[:, 0], 'umap_y': emb[:, 1], 'fvar_w1': fvarA, 'mean_w1_mjy': meanA, 'max_w1_mjy': maxA,
                  'clagn_score': scoreA, 'is_known_clagn': is_cl,
                  'in_region_clagn': in_region(reg_clagn, emb)})
A_PATH = os.path.join(DATA, 'sampleA_embedding_objectid.csv')
A.to_csv(A_PATH, index=False)
say(f'[{time.time()-t0:6.0f}s] wrote {A_PATH} (Sample A part; Zeltyn-defined region column added later)')

# ============================================================================= Zeltyn (all 204 rows)
def jname_to_coord(name):
    j = re.sub(r'^(SDSS\s*)?J', '', name.strip())
    m = re.match(r'(\d{2})(\d{2})(\d{2}\.\d+)([+-])(\d{2})(\d{2})(\d{2}\.\d+)$', j)
    if m is None:
        return None
    return SkyCoord(f'{m.group(1)}:{m.group(2)}:{m.group(3)}', f'{m.group(4)}{m.group(5)}:{m.group(6)}:{m.group(7)}',
                    unit=(u.hourangle, u.deg))

zt = pd.read_csv(ZELTYN_CSV)
zt['coord'] = zt['Name'].map(jname_to_coord)
zt = zt[zt['coord'].notna()].reset_index(drop=True)
zt['ra'] = [c.ra.deg for c in zt['coord']]
zt['dec'] = [c.dec.deg for c in zt['coord']]
zt['is_clagn'] = zt['Class'].str.startswith('CL-AGN')
sample_table_Z = clean_sample(list(zt['coord']), list(zt['Name']), consolidate_nearby_objects=False, verbose=0)
# objectid in sample_table_Z is 1..N in the same order as zt
lcZ = fetch_wise_chunked(sample_table_Z, CACHE)
say(f'[{time.time()-t0:6.0f}s] Zeltyn rows {len(zt)}, with unWISE data: {lcZ.index.get_level_values("objectid").nunique()}')
dfZ = lcZ.copy()
dfZ.index = dfZ.index.set_levels(dfZ.index.levels[dfZ.index.names.index('band')].str.replace('WISE_', '', regex=False), level='band')
datZ, idsZ, _, fvarZ, meanZ, maxZ = preprocess(dfZ, 'Zeltyn all')
embZ = mapp.transform(datZ)
zt_emb = pd.DataFrame({'umap_x': embZ[:, 0], 'umap_y': embZ[:, 1], 'fvar_w1': fvarZ, 'mean_w1_mjy': meanZ},
                      index=pd.Index(idsZ - 1, name='row'))   # objectid 1..N -> zt row 0..N-1
ztf = zt.drop(columns=['coord']).join(zt_emb)
ztf['projected'] = ztf['umap_x'].notna()
maskZ = ztf['projected'].values
scoreZ = np.full(len(ztf), np.nan); scoreZ[maskZ] = clagn_score(ztf.loc[maskZ, ['umap_x', 'umap_y']].values)
ztf['clagn_score'] = scoreZ

zcl_pts = ztf.loc[maskZ & ztf['is_clagn'].values & (ztf['z'] < 1).values, ['umap_x', 'umap_y']].values
zev_pts = ztf.loc[maskZ & ~ztf['is_clagn'].values & (ztf['z'] < 1).values, ['umap_x', 'umap_y']].values
reg_zeltyn = region_from(zcl_pts)          # notebook definition (Zeltyn CL-AGN enriched)

A['in_region_zeltyn'] = in_region(reg_zeltyn, emb)
A.to_csv(A_PATH, index=False)
for col, reg in [('in_region_zeltyn', reg_zeltyn), ('in_region_clagn', reg_clagn)]:
    v = np.zeros(len(ztf), bool); v[maskZ] = in_region(reg, ztf.loc[maskZ, ['umap_x', 'umap_y']].values); ztf[col] = v
ztf.to_csv(os.path.join(DATA, 'zeltyn_embedding_full.csv'), index=False)

# ============================================================================= report
say('\n=== Region definitions (%dx%d bins, enrichment > %dx) ===' % (NBINS, NBINS, ENRICH))
for name, reg in [('Zeltyn-defined', reg_zeltyn), ('SampleA Turn-on/off-defined', reg_clagn)]:
    inA = in_region(reg, emb)
    say(f'{name}: {reg.sum()} bins; Sample A inside {inA.sum()} ({100*inA.mean():.1f}%); known CLAGN inside {(inA & is_cl).sum()}/{is_cl.sum()};'
        f' non-CLAGN inside {(inA & ~is_cl).sum()}')
    inZ = in_region(reg, zcl_pts); inE = in_region(reg, zev_pts)
    say(f'   held-out Zeltyn CL-AGN (z<1) inside {inZ.sum()}/{len(zcl_pts)} ({100*inZ.mean():.0f}%),  EVQ inside {inE.sum()}/{len(zev_pts)}'
        f'  -> enrichment vs Sample A fraction: CL-AGN {inZ.mean()/max(inA.mean(),1e-9):.1f}x, EVQ {inE.mean()/max(inA.mean(),1e-9):.1f}x')
    cnt = A.loc[inA & ~is_cl, 'labels'].value_counts()
    say('   non-CLAGN Sample A inside, by label: ' + ', '.join(f'{k} {v}' for k, v in cnt.items()))
say('\n=== kNN clagn_score (k=%d): fraction of Turn-on/off among nearest Sample A neighbours ===' % K_NN)
say(f'Sample A overall CLAGN fraction: {is_cl.mean():.3f}')
for name, m in [('Sample A SDSS_QSO-only', (labA == 1)), ('Sample A WISE_Variable-only', (labA == 2)),
                ('Sample A Turn-on', is_on), ('Sample A Turn-off', is_off)]:
    say(f'  {name:32s} median {np.median(scoreA[m]):.3f}  frac>0.15: {(scoreA[m] > 0.15).mean():.2f}  n={m.sum()}')
for name, m in [('Zeltyn CL-AGN z<1', maskZ & ztf.is_clagn.values & (ztf.z < 1).values),
                ('Zeltyn EVQ z<1', maskZ & ~ztf.is_clagn.values & (ztf.z < 1).values)]:
    s = ztf.loc[m, 'clagn_score']
    say(f'  {name:32s} median {s.median():.3f}  frac>0.15: {(s > 0.15).mean():.2f}  n={m.sum()}')
say(f'\n[{time.time()-t0:6.0f}s] done. Files written to data/.')
log.close()
