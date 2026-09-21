"""Recover the September candidates' baseline spectra for scientific review.

Only the public parent-reference spectrum is added to spectra_dl. Internal
spectra remain in the ignored collaboration directory. No publication occurs.
"""
from pathlib import Path
import importlib,json,argparse
from concurrent.futures import ThreadPoolExecutor
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/reselection_2026-09-20'
PARSE=importlib.import_module('03d_fetch_spectra')
URLS=importlib.import_module('13b_calibration_lines')
NATIVE=importlib.import_module('17_fetch_review_spectra').finite

def one(row):
    name='P'+str(int(row.poolid))
    directory=ROOT/'data/spectra_cache'/name
    directory.mkdir(parents=True,exist_ok=True)
    for url in URLS.dr16_urls(row.plate,row.mjd,row.fiberid,row.survey):
        path=directory/url.rsplit('/',1)[-1]
        if not PARSE.fetch(url,str(path),tries=1):continue
        record=PARSE.parse(str(path))
        record.update(mjd=int(row.mjd),phase=None,program=str(row.programname),
                      coadd=False,url=url,source='SDSS',proprietary=False)
        dest=ROOT/'data/spectra_dl'/f'{name}.json'
        records=json.loads(dest.read_text()) if dest.exists() else []
        if not any(r.get('source','SDSS')=='SDSS' and r.get('mjd')==row.mjd
                   and not r.get('coadd') and not r.get('proprietary') for r in records):
            records.append(NATIVE(record))
            dest.write_text(json.dumps(NATIVE(records),separators=(',',':'),allow_nan=False))
        return dict(name=name,status='available',path=str(path),mjd=int(row.mjd))
    return dict(name=name,status='unavailable',mjd=int(row.mjd))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--all-nights',action='store_true');args=ap.parse_args()
    targets=pd.read_csv(OUT/'compact_review_objects.csv')
    names=set(targets.name if args.all_nights else targets.loc[targets.eligible_nights.str.contains('sep23'),'name'])
    parent=pd.read_csv(ROOT/'data/parent_pool_scored.csv')
    parent=parent[('P'+parent.poolid.astype(str)).isin(names)]
    with ThreadPoolExecutor(max_workers=4) as executor:
        rows=list(executor.map(one,parent.itertuples(index=False)))
    represented=set('P'+parent.poolid.astype(str))
    rows.extend(dict(name=name,status='use public archive inventory; not an old DR16 parent identity') for name in sorted(names-represented))
    pd.DataFrame(rows).to_csv(OUT/('three_night_baseline_manifest.csv' if args.all_nights else 'sep23_baseline_manifest.csv'),index=False)
    print(pd.DataFrame(rows).status.value_counts().to_dict(),flush=True)

if __name__=='__main__':main()
