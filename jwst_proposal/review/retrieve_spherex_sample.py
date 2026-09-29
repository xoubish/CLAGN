"""Retrieve public QR2 cutouts and apply the existing Figure 1 extractor.

Stage new photometry in local_data/spherex/sample; never overwrite supplied spectra.
Keep per-image outcomes and hashes, with no selection on measured source flux.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse, io, json, hashlib, sys, time, gzip, threading, signal
from itertools import zip_longest
import requests
import pandas as pd
from astropy.io.votable import parse_single_table

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT/'local_data/spherex/sample'
sys.path.insert(0, str(ROOT))
from spherex_extract import extract_one
LOCAL = threading.local()
STOP = threading.Event()

def save_json(path, value):
    temporary = path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)

def session():
    if not hasattr(LOCAL, 'session'):
        LOCAL.session = requests.Session()
    return LOCAL.session

def index(row):
    path = OUT/f'{row.id}_index.csv'
    if path.exists():
        cached = pd.read_csv(path)
        if 'access_url' in cached: return cached
    params = dict(COLLECTION='spherex_qr2', POS=f'CIRCLE {row.ra} {row.dec} 0.00001')
    xml = OUT/f'{row.id}_query.xml'
    if not xml.exists():
        r = requests.get('https://irsa.ipac.caltech.edu/SIA', params=params, timeout=90)
        r.raise_for_status()
        xml.write_bytes(r.content)
    df = parse_single_table(io.BytesIO(xml.read_bytes())).to_table(use_names_over_ids=True).to_pandas()
    for col in df.select_dtypes('object'):
        df[col] = df[col].map(lambda v: v.decode() if isinstance(v, bytes) else v)
    df.to_csv(path, index=False)
    return df

def fetch(item):
    row, url = item
    key = hashlib.sha256(url.encode()).hexdigest()[:20]
    folder = OUT/row.id
    folder.mkdir(exist_ok=True)
    record = folder/f'{key}.json'
    if record.exists():
        try:
            cached = json.loads(record.read_text())
            if cached['status'] in ('accepted', 'quality_rejected'): return cached
        except (ValueError, KeyError):
            pass  # Recover records truncated by an earlier interrupted write.
    cuturl = url.split('?')[0]+f'?center={row.ra},{row.dec}&size=100arcsec&gzip=true'
    result = dict(id=row.id, url=cuturl, status='error')
    path = folder/f'{key}.fits'
    for attempt in range(2):
        try:
            r = session().get(cuturl, timeout=(12,30))
            r.raise_for_status()
            raw = gzip.decompress(r.content) if r.content[:2] == b'\x1f\x8b' else r.content
            if not raw.startswith(b'SIMPLE'): raise ValueError('Not a FITS response')
            path.write_bytes(raw)
            result.update(bytes=len(raw), transferred_bytes=len(r.content), sha256=hashlib.sha256(raw).hexdigest())
            measured = extract_one(path, row.ra, row.dec)
            result.update(status='accepted' if measured else 'quality_rejected', measurement=measured)
            result.pop('error', None)
            break
        except Exception as e:
            result['error'] = f'{type(e).__name__}: {e}'
    save_json(record, result)
    # The URL, hash and photometry remain reproducible; do not retain thousands
    # of identical full PSF cubes from the public cutout service.
    if path.exists(): path.unlink()
    return result

def main():
    signal.signal(signal.SIGINT, lambda *_: STOP.set())
    signal.signal(signal.SIGTERM, lambda *_: STOP.set())
    ap = argparse.ArgumentParser()
    ap.add_argument('--index-only', action='store_true')
    ap.add_argument('--ids', nargs='*')
    ap.add_argument('--limit', type=int)
    ap.add_argument('--max-errors', type=int, default=12)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    sample = pd.read_csv(ROOT/'inputs/jwst_sample_cycle6.csv')
    if args.ids: sample = sample[sample.id.isin(args.ids)]
    tasks = []; summary = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(index, row): row for row in sample.itertuples()}
        for future in as_completed(futures):
            row = futures[future]
            try:
                df = future.result()
                if 'dataproduct_subtype' in df: df = df[df.dataproduct_subtype.eq('science')]
                urls = sorted(set(df.access_url.dropna()))
                print(row.id, 'images', len(urls), flush=True)
                summary.append(dict(id=row.id, images=len(urls), status='indexed'))
                if args.limit: urls = urls[:args.limit]
                tasks.extend((row, u) for u in urls)
            except Exception as e:
                summary.append(dict(id=row.id, status='query_failed', error=str(e)))
                print(row.id, 'QUERY FAILED', str(e)[:160], flush=True)
    pd.DataFrame(summary).to_csv(OUT/'coverage_latest.csv', index=False)
    if args.index_only: return
    results = []
    # Round-robin targets so an interrupted run does not leave all other
    # sources untouched. Save each finished image immediately.
    groups=[[t for t in tasks if t[0].id==row.id] for row in sample.itertuples()]
    tasks=[t for batch in zip_longest(*groups) for t in batch if t is not None]
    errors=0
    with ThreadPoolExecutor(max_workers=4) as pool:
        for start in range(0,len(tasks),8):
            futures=[pool.submit(fetch,t) for t in tasks[start:start+8]]
            for future in as_completed(futures):
                result=future.result();results.append(result)
                if result['status']=='error':
                    errors+=1;print('ERROR',result['id'],result['error'][:180],flush=True)
            print('processed',len(results),'/',len(tasks),'errors',errors,flush=True)
            # A compact manifest is readable while the process runs. All
            # detailed per-image results are already durable JSON files.
            (OUT/'progress.json').write_text(json.dumps(dict(
                requested_images=len(tasks),processed_this_run=len(results),
                errors_this_run=errors,last_update_unix=time.time()),indent=2)+'\n')
            if errors>=args.max_errors or STOP.is_set():
                print('Stopped at checkpoint; rerun resumes saved images.',flush=True)
                break
    status=[]
    for row in sample.itertuples():
        records = [json.loads(f.read_text()) for f in (OUT/row.id).glob('*.json')]
        accepted = [dict(**r['measurement'], source_url=r['url']) for r in records if r['status']=='accepted']
        if accepted: pd.DataFrame(accepted).sort_values('mjds').to_csv(OUT/f'{row.id}.csv',index=False)
        requested=sum(t[0].id==row.id for t in tasks)
        bad=sum(r['status']=='error' for r in records)
        status.append(dict(id=row.id,requested=requested,processed=len(records),accepted=len(accepted),
                           errors=bad,complete=len(records)>=requested and bad==0))
        print(row.id, 'accepted', len(accepted), 'of',len(records), flush=True)
    pd.DataFrame(status).to_csv(OUT/'extraction_status.csv',index=False)
    (OUT/'retrieval_latest.json').write_text(json.dumps(dict(scope=__doc__, records=results),indent=2)+'\n')

if __name__ == '__main__': main()
