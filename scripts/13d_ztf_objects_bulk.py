"""
13d_ztf_objects_bulk.py  --  ZTF DR object-table statistics for a position list via the IRSA Gator bulk upload
(catalog ztf_objects_dr24: per object and filter, medianmag/minmag/maxmag over 2018-2025, refmag, ngoodobs, magrms, medianabsdev).
100 positions per request take ~6 s, so the whole pool runs in ~20 min.  Cached per batch in data/ztf_objects_cache/<tag>/.
Usage: /opt/anaconda3/bin/python 13d_ztf_objects_bulk.py <in.csv> <tag> [ra_col dec_col name_col]
Output: data/ztf_objects_<tag>.csv  (one row per object: for g and r the entry with most good observations)
"""
import io, os, sys, time, requests, numpy as np, pandas as pd
from astropy.table import Table; from astropy.io import ascii
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(HERE, 'data')
GATOR = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query'; BATCH = 100; CAT = 'ztf_objects_dr24'
COLS = ['oid', 'filtercode', 'ngoodobs', 'nobs', 'refmag', 'refmagerr', 'magrms', 'maxmag', 'medianabsdev', 'medianmag', 'medmagerr', 'minmag', 'meanmag', 'lineartrend', 'chisq']


def gator(pos, tries=4):
    up = Table({'ra': pos.ra.values, 'dec': pos.dec.values}); buf = io.StringIO(); ascii.write(up, buf, format='ipac')
    for attempt in range(tries):
        try:
            r = requests.post(GATOR, data={'catalog': CAT, 'spatial': 'Upload', 'uplist': 'pos.tbl', 'radius': 2, 'radunits': 'arcsec', 'outfmt': 1, 'selcols': ','.join(COLS)},
                              files={'filename': ('pos.tbl', buf.getvalue())}, timeout=900)
            r.raise_for_status()
            if 'stat="ERROR"' in r.text[:400]:
                raise RuntimeError(r.text[:200])
            return ascii.read(r.text, format='ipac').to_pandas()
        except Exception as e:
            print(f'   gator attempt {attempt+1}: {str(e)[:120]}', flush=True); time.sleep(20 * (attempt + 1))
    raise RuntimeError('Gator failed')


def main():
    inp, tag = sys.argv[1], sys.argv[2]; ra_col, dec_col, name_col = (sys.argv[3:6] if len(sys.argv) >= 6 else ('ra', 'dec', 'name'))
    pos = pd.read_csv(inp)[[name_col, ra_col, dec_col]].rename(columns={name_col: 'name', ra_col: 'ra', dec_col: 'dec'}).drop_duplicates('name').reset_index(drop=True)
    cache = os.path.join(DATA, 'ztf_objects_cache', tag); os.makedirs(cache, exist_ok=True); t0 = time.time(); frames = []
    for i in range(0, len(pos), BATCH):
        f = os.path.join(cache, f'batch_{i:06d}.csv')
        if os.path.exists(f):
            frames.append(pd.read_csv(f)); continue
        sub = pos.iloc[i:i + BATCH]; g = gator(sub)
        if len(g):
            g['name'] = sub.name.values[g.cntr_01.values.astype(int) - 1]      # cntr_01 = 1-based row of the uploaded list
        g.to_csv(f, index=False); frames.append(g)
        if (i // BATCH) % 10 == 0:
            print(f'[{time.time()-t0:5.0f}s] batch {i//BATCH+1}/{(len(pos)-1)//BATCH+1}: {len(g)} rows', flush=True)
    z = pd.concat([x for x in frames if len(x)], ignore_index=True)
    # nearest ZTF object to the uploaded position (a brighter neighbour within 2" would otherwise win on ngoodobs), within 1.5", >= 20 good epochs
    z = z[z.filtercode.isin(['zg', 'zr']) & (z.dist_x <= 1.5) & (z.ngoodobs >= 20)].sort_values('dist_x').drop_duplicates(['name', 'filtercode'])
    out = pos.set_index('name')[['ra', 'dec']]
    for fc, b in [('zg', 'g'), ('zr', 'r')]:
        s = z[z.filtercode == fc].set_index('name')
        for c in ['ngoodobs', 'medianmag', 'minmag', 'maxmag', 'refmag', 'magrms', 'medianabsdev', 'lineartrend', 'dist_x']:
            out[f'ztf_{b}_{c}'] = s[c].reindex(out.index)
        out[f'ztf_{b}_amp'] = out[f'ztf_{b}_maxmag'] - out[f'ztf_{b}_minmag']
    out.reset_index().to_csv(os.path.join(DATA, f'ztf_objects_{tag}.csv'), index=False)
    print(f'wrote data/ztf_objects_{tag}.csv: {out.ztf_r_medianmag.notna().sum()} of {len(out)} objects with r, {out.ztf_g_medianmag.notna().sum()} with g  [{time.time()-t0:.0f}s]')


if __name__ == '__main__':
    main()
