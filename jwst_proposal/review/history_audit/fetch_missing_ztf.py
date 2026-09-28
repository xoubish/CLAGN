"""Retrieve public DR24 light curves missing from the local sample caches."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from datetime import datetime, timezone
import json
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent / 'raw_ztf'
IDS = ['F01', 'F09', 'F11', 'F12', 'R03', 'R06', 'R07', 'R08', 'R09', 'R10', 'R11', 'R12']

def fetch(row):
    path = OUT / f'{row.id}.csv'
    status = dict(id=row.id, collection='ztf_dr24', ra=row.ra, dec=row.dec,
                  queried_utc=datetime.now(timezone.utc).isoformat())
    if path.exists():
        return dict(status, status='cached', bytes=path.stat().st_size)
    params = dict(POS=f'CIRCLE {row.ra:.8f} {row.dec:.8f} {3/3600:.10f}',
                  BANDNAME='g,r', FORMAT='csv', BAD_CATFLAGS_MASK=32768,
                  COLLECTION='ztf_dr24')
    try:
        response = requests.get('https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves',
                                params=params, timeout=55)
        response.raise_for_status()
        if 'mjd' not in response.text[:2000]:
            raise ValueError('No CSV light curve returned')
        path.write_text(response.text)
        status.update(status='downloaded', bytes=len(response.content))
    except Exception as exc:
        status.update(status='failed', reason=str(exc)[:400])
    path.with_suffix('.json').write_text(json.dumps(status, indent=2)+'\n')
    print(row.id, status['status'], flush=True)
    return status

if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    sample = pd.read_csv(ROOT/'jwst_proposal/inputs/jwst_sample_cycle6.csv')
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(fetch, sample[sample.id.isin(IDS)].itertuples()))
    (OUT/'manifest.json').write_text(json.dumps(results, indent=2)+'\n')
