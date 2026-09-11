"""
13e_dr16_spectra_pool.py  --  download the DR16 (DR17 spec-lite) spectrum of every pool quasar that has NO SDSS-V epoch, so that
synthetic ZTF g/r magnitudes at the archival epoch can be compared with current ZTF photometry (the CLAS+ test of Nakazono+2026).
Cache: data/spectra_cache/<name>/ (git-ignored), same layout 03d uses.  Resumable.  Usage: 13e_dr16_spectra_pool.py [--threads 8]
"""
import os, sys, time, argparse, importlib.util, pandas as pd
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, 'data')
def load(mod):
    s = importlib.util.spec_from_file_location(mod, os.path.join(HERE, mod + '.py')); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
f03d = load('03d_fetch_spectra'); f13b = load('13b_calibration_lines')
ap = argparse.ArgumentParser(); ap.add_argument('--threads', type=int, default=8); ap.add_argument('--input', default=None, help='alternative position table with name, plate, mjd, fiberid (e.g. the turn-on galaxy parent)'); a = ap.parse_args()
pool = pd.read_csv(os.path.join(DATA, 'parent_pool_scored.csv')); pool = pool[pool.projected.fillna(False)].copy(); pool['name'] = 'P' + pool.poolid.astype(str)
sv = set(pd.read_csv(os.path.join(DATA, 'sdssv_internal_epochs.csv')).name)
todo = pool[~pool.name.isin(sv)]
if a.input:
    todo = pd.read_csv(a.input); todo['survey'] = todo.get('survey', 'sdss')
print(f'{len(todo)} objects to fetch', flush=True)
def one(r):
    d = os.path.join(f03d.CACHE, r.name); os.makedirs(d, exist_ok=True)
    if any(f.endswith('.fits') for f in os.listdir(d)):
        return 'cached'
    for u in f13b.dr16_urls(r.plate, r.mjd, r.fiberid, r.survey):
        if f03d.fetch(u, os.path.join(d, os.path.basename(u))):
            return 'ok'
    return 'missing'
t0 = time.time(); n = {'ok': 0, 'cached': 0, 'missing': 0}
with ThreadPoolExecutor(max_workers=a.threads) as ex:
    for k, s in enumerate(ex.map(one, todo.itertuples(index=False))):
        n[s] += 1
        if (k + 1) % 500 == 0: print(f'[{time.time()-t0:5.0f}s] {k+1}/{len(todo)} {n}', flush=True)
print(f'done in {time.time()-t0:.0f}s: {n}', flush=True)
