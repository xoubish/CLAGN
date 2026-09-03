"""
02_parent_pool.py  --  Stage 1 option of FOLLOWUP_PLAN.md: a bright DR16 quasar parent pool inside the RA windows
of the 2026B NGPS runs, then WISE light curves + projection onto the saved W1 manifold.

Stage 'query'  : SDSS DR16 SpecObj (sciencePrimary) class='QSO', zWarning=0, 0.02<z<0.8, Dec>-15, psfMag_r<20,
                 RA in the observable windows -> data/parent_pool_dr16qso.csv (with plate/mjd/fiber of the
                 archival SDSS spectrum and PSF g,r,i).
Stage 'wise'   : unWISE W1/W2 light curves for the pool in cached, retried chunks -> data/wise_cache/pool/
Stage 'project': GP + normalise exactly as 01_rebuild_manifold.py, mapp.transform with data/umap_w1_model.pkl,
                 kNN clagn_score against Sample A, in_region flags -> data/parent_pool_scored.csv

Usage: /opt/anaconda3/bin/python 02_parent_pool.py query|wise|project
"""
import io, os, sys, time, pickle, warnings
import numpy as np
import pandas as pd
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, os.path.join(HERE, 'code_src'))
warnings.filterwarnings('ignore')
DATA = os.path.join(HERE, 'data')
POOL = os.path.join(DATA, 'parent_pool_dr16qso.csv')
CACHE = os.path.join(DATA, 'wise_cache', 'pool')
SKY16 = 'https://skyserver.sdss.org/dr16/SkyServerWS/SearchTools/SqlSearch'

# RA windows (hours) reachable at airmass<2 on 2026-09-23 (first half) and 2026-10-26/27 (full nights),
# see 05_observability.py output. Dec > -15 for Palomar in bright time.
RA_WINDOWS_H = [(15.5, 24.0), (0.0, 4.5), (6.5, 11.0)]
DEC_MIN, ZMIN, ZMAX, RMAX = -15.0, 0.02, 0.8, 20.0


def skyserver(sql, url=SKY16, tries=4):
    for attempt in range(tries):
        try:
            r = requests.get(url, params={'cmd': sql, 'format': 'csv'}, timeout=600)
            r.raise_for_status()
            txt = r.text
            if txt.startswith('#Table1'):
                txt = txt.split('\n', 1)[1]
            if txt.lstrip().startswith('<') or txt.lstrip().startswith('{'):
                raise RuntimeError(txt[:200])
            return pd.read_csv(io.StringIO(txt)) if txt.strip() else pd.DataFrame()
        except Exception as e:
            print(f'   skyserver attempt {attempt+1}: {str(e)[:150]}', flush=True); time.sleep(10 * (attempt + 1))
    raise RuntimeError('SkyServer query failed')


def stage_query():
    cols = ('s.specObjID AS specobjid, s.bestObjID AS objid, s.plate, s.mjd, s.fiberid, s.ra, s.dec, s.z, s.zErr AS zerr, '
            's.subclass, s.survey, s.programname, s.snMedian AS sn_median, '
            'p.psfMag_g AS psfmag_g, p.psfMag_r AS psfmag_r, p.psfMag_i AS psfmag_i, p.extinction_r, p.type AS phototype')
    parts = []
    step = 7.5  # degrees of RA per query
    for h0, h1 in RA_WINDOWS_H:
        for a in np.arange(h0 * 15, h1 * 15, step):
            b = min(a + step, h1 * 15)
            q = (f"SELECT {cols} FROM SpecObj s JOIN PhotoObjAll p ON s.bestObjID = p.objID "
                 f"WHERE s.class='QSO' AND s.zWarning=0 AND s.z BETWEEN {ZMIN} AND {ZMAX} AND s.dec > {DEC_MIN} "
                 f"AND p.psfMag_r < {RMAX} AND p.psfMag_r > 0 AND s.ra >= {a} AND s.ra < {b}")
            df = skyserver(q)
            parts.append(df)
            print(f'   RA {a/15:5.2f}h-{b/15:5.2f}h: {len(df)} quasars', flush=True)
    pool = pd.concat(parts, ignore_index=True).drop_duplicates('specobjid')
    pool.insert(0, 'poolid', np.arange(1, len(pool) + 1))
    pool.to_csv(POOL, index=False)
    print(f'wrote {POOL}: {len(pool)} quasars; r<19: {(pool.psfmag_r < 19).sum()}, r<19.5: {(pool.psfmag_r < 19.5).sum()}')
    print(pool.groupby(pd.cut(pool.ra / 15, [0, 4.5, 6.5, 11, 15.5, 24])).size().to_string())


def stage_wise():
    """unWISE W1 fetch for the pool. Fixes vs. the naive loop: the 140 MB _metadata file is read once per
    process (memoised parquet_dataset), objects already present in any cached chunk are skipped (so chunk size
    can change between runs), chunks are large (150) and sky-sorted so each S3 partition is read as few times
    as possible, and only W1 is requested.  Stride over chunks with `wise <rmax> <worker> <nworkers>`."""
    import functools
    import pyarrow.fs as pafs
    import pyarrow.dataset as pads
    import hpgeom
    _orig = pafs.S3FileSystem
    def _patient(*a, **k):
        k.setdefault('request_timeout', 300); k.setdefault('connect_timeout', 60); return _orig(*a, **k)
    pafs.S3FileSystem = _patient
    _orig_pd = pads.parquet_dataset
    _cache = {}
    def _memo_parquet_dataset(path, filesystem=None, **k):
        key = (path, tuple(sorted(k.items())))
        if key not in _cache:
            _cache[key] = _orig_pd(path, filesystem=filesystem, **k)
        return _cache[key]
    pads.parquet_dataset = _memo_parquet_dataset
    import wise_functions
    wise_functions.pyarrow.dataset.parquet_dataset = _memo_parquet_dataset
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from sample_selection import clean_sample
    from wise_functions import wise_get_lightcurves

    rmax = float(sys.argv[2]) if len(sys.argv) > 2 else 19.0
    worker, nworkers = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (0, 1)
    pool = pd.read_csv(POOL)
    pix = hpgeom.angle_to_pixel(2 ** 5, pool.ra.values, pool.dec.values, nest=True, lonlat=True)
    rclass = np.digitize(pool.psfmag_r.values, [19.0, 19.5])
    pool = pool.iloc[np.lexsort((pix, rclass))].reset_index(drop=True)
    pool = pool[pool.psfmag_r < rmax]
    os.makedirs(CACHE, exist_ok=True)
    done = set()
    for f in os.listdir(CACHE):
        if f.endswith('.parquet'):
            try:
                done |= set(pd.read_parquet(os.path.join(CACHE, f)).index.get_level_values('objectid').unique())
            except Exception:
                pass
    todo = pool[~pool.poolid.isin(done)]
    coords = [SkyCoord(r.ra * u.deg, r.dec * u.deg) for r in todo.itertuples()]
    st = clean_sample(coords, [str(p) for p in todo.poolid], consolidate_nearby_objects=False, verbose=0)
    st['objectid'] = todo.poolid.values
    chunk, tries, t0 = 150, 6, time.time()
    nchunks = (len(st) - 1) // chunk + 1
    print(f'{len(pool)} pool objects with r<{rmax}; {len(done)} already cached; {len(st)} to fetch in {nchunks} chunks; '
          f'worker {worker}/{nworkers}', flush=True)
    for i in range(0, len(st), chunk):
        if (i // chunk) % nworkers != worker:
            continue
        sub = st[i:i + chunk]
        f = os.path.join(CACHE, f'chunk_{sub["objectid"].min():05d}_{sub["objectid"].max():05d}_{len(sub)}.parquet')
        if os.path.exists(f):
            continue
        for attempt in range(tries):
            try:
                d = wise_get_lightcurves(sub, radius=1.0, bandlist=['WISE_W1']).data; break
            except ValueError:
                d = pd.DataFrame(); break
            except OSError as e:
                print(f'   chunk {i}: {str(e)[:80]} ... retry {attempt+1}/{tries}', flush=True); time.sleep(20 * (attempt + 1))
        else:
            raise RuntimeError(f'chunk {i} failed')
        d.to_parquet(f)
        print(f'[{time.time()-t0:6.0f}s] chunk {i//chunk+1}/{nchunks} ({len(d)} rows)', flush=True)
    print('WISE stage complete.')


def stage_project():
    from sklearn.neighbors import NearestNeighbors
    from AGNzoo_functions import unify_lc_gp, stat_bands, combine_bands, normalize_clipmax_objects, dtw_distance  # noqa (dtw needed to unpickle)
    pool = pd.read_csv(POOL)
    files = sorted(os.path.join(CACHE, f) for f in os.listdir(CACHE) if f.endswith('.parquet'))
    lc = pd.concat([pd.read_parquet(f) for f in files])
    lc = lc[lc.index.get_level_values('band') == 'WISE_W1']
    lc.index = lc.index.set_levels(lc.index.levels[lc.index.names.index('band')].str.replace('WISE_', '', regex=False), level='band')
    print(f'pool objects with unWISE W1: {lc.index.get_level_values("objectid").nunique()} / {len(pool)}')
    objs, dobjs, flabels, keeps = unify_lc_gp(lc, ['W1'], xres=160, numplots=0, low_limit_size=5)
    fvar, maxarr, meanarr = stat_bands(objs, dobjs, ['W1'], sigmacl=5)
    dat = normalize_clipmax_objects(combine_bands(objs, ['W1']), maxarr, band=0)
    valid = ~np.isnan(dat).any(axis=1) & np.isfinite(dat).all(axis=1)
    ids = np.asarray(lc.index.get_level_values('objectid').unique())[np.asarray(keeps)][valid]   # objectid == poolid
    print(f'after GP cut / NaN drop: {valid.sum()}')
    with open(os.path.join(DATA, 'umap_w1_model.pkl'), 'rb') as f:
        mapp = pickle.load(f)
    t0 = time.time()
    emb = mapp.transform(dat[valid])
    print(f'transform done in {time.time()-t0:.0f}s')
    A = pd.read_csv(os.path.join(DATA, 'sampleA_embedding_objectid.csv'))
    embA = A[['umap_x', 'umap_y']].values
    is_cl = A.is_known_clagn.values
    nn = NearestNeighbors(n_neighbors=50).fit(embA)
    _, idx = nn.kneighbors(emb)
    score = is_cl[idx].mean(axis=1)
    # region flags: reuse the 10x10 grid of Sample A and the bins flagged in sampleA_embedding_objectid.csv
    NB = 10
    hist_all, xe, ye = np.histogram2d(embA[:, 0], embA[:, 1], bins=NB)
    def bin_idx(pts):
        return (np.clip(np.searchsorted(xe[1:], pts[:, 0]), 0, NB - 1), np.clip(np.searchsorted(ye[1:], pts[:, 1]), 0, NB - 1))
    ixA, iyA = bin_idx(embA)
    flags = {}
    for col in ['in_region_zeltyn', 'in_region_clagn']:
        reg = np.zeros((NB, NB), bool)
        reg[ixA[A[col].values], iyA[A[col].values]] = True
        ix, iy = bin_idx(emb); flags[col] = reg[ix, iy]
    # distance in embedding space to nearest Zeltyn CL-AGN (z<1) as a continuous alternative
    Z = pd.read_csv(os.path.join(DATA, 'zeltyn_embedding_full.csv'))
    Zc = Z[Z.projected & Z.is_clagn & (Z.z < 1)][['umap_x', 'umap_y']].values
    nnZ = NearestNeighbors(n_neighbors=5).fit(Zc)
    dZ, _ = nnZ.kneighbors(emb)
    # local density of Zeltyn CL-AGN vs Sample A (kNN ratio): fraction of Zeltyn among 50 nearest of the union
    U = np.vstack([embA, Zc]); isZ = np.r_[np.zeros(len(embA), bool), np.ones(len(Zc), bool)]
    nnU = NearestNeighbors(n_neighbors=50).fit(U)
    _, iu = nnU.kneighbors(emb)
    zeltyn_score = isZ[iu].mean(axis=1) / (len(Zc) / len(U))     # 1 = average density, >3 = enriched
    res = pd.DataFrame({'poolid': ids, 'umap_x': emb[:, 0], 'umap_y': emb[:, 1], 'fvar_w1': fvar[0][valid],
                        'mean_w1_mjy': meanarr[0][valid], 'clagn_score': score, 'zeltyn_density_ratio': zeltyn_score.round(2),
                        'd5_zeltyn': dZ.mean(axis=1).round(3), **flags})
    out = pool.merge(res, on='poolid', how='left')
    out['projected'] = out.umap_x.notna()
    out.to_csv(os.path.join(DATA, 'parent_pool_scored.csv'), index=False)
    p = out[out.projected]
    print(f'wrote data/parent_pool_scored.csv: {len(out)} rows, {len(p)} projected')
    print(f'  in_region_zeltyn: {p.in_region_zeltyn.sum()}  in_region_clagn: {p.in_region_clagn.sum()}  clagn_score>=0.15: {(p.clagn_score>=0.15).sum()}  '
          f'zeltyn_density_ratio>=3: {(p.zeltyn_density_ratio>=3).sum()}')
    print(f'  r<19 & in_region_zeltyn: {((p.psfmag_r<19)&p.in_region_zeltyn).sum()}  r<19 & density>=3: {((p.psfmag_r<19)&(p.zeltyn_density_ratio>=3)).sum()}')


if __name__ == '__main__':
    stage = sys.argv[1] if len(sys.argv) > 1 else 'query'
    {'query': stage_query, 'wise': stage_wise, 'project': stage_project}[stage]()
