"""Download SDSS-V daily spectra for the compact local review, resumably.

Uses the existing authorized collaboration credentials without logging them.
All FITS, derived spectra and manifests remain in the git-ignored review folder.
No public export, commit, or push is performed.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import importlib
import json
import threading
import time
import numpy as np
import pandas as pd
import requests

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'data/reselection_2026-09-20'
DEST=OUT/'sdssv_spectra'
ARCHIVE=importlib.import_module('03f_sdssv_internal')
PARSE=importlib.import_module('03d_fetch_spectra')
LOCAL=threading.local()


def finite(value):
    if isinstance(value,dict):return {k:finite(v) for k,v in value.items()}
    if isinstance(value,list):return [finite(v) for v in value]
    if isinstance(value,np.generic):value=value.item()
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def session():
    if not hasattr(LOCAL,'session'):
        LOCAL.session=requests.Session()
        LOCAL.session.auth=ARCHIVE.auth()
    return LOCAL.session


def one(row):
    name=row.canonical_name
    version=str(row.run2d_row).strip() or 'v6_2_1'
    url=ARCHIVE.spec_url(version,'daily',int(row.field),int(row.mjd),int(row.catalogid),str(row.spec_file).strip())
    directory=OUT/'spectra_files'/version/name
    directory.mkdir(parents=True,exist_ok=True)
    path=directory/str(row.spec_file).strip()
    metadata=dict(name=name,mjd=int(row.mjd),field=int(row.field),status='failed',version=version)
    for attempt in range(3):
        try:
            if not path.exists():
                response=session().get(url,timeout=90)
                if response.status_code in (401,403,404):
                    return metadata|{'reason':f'HTTP {response.status_code}'},None
                response.raise_for_status()
                if len(response.content)<10000 or not response.content.startswith(b'SIMPLE'):
                    return metadata|{'reason':'not a spectrum FITS file'},None
                temp=path.with_suffix('.tmp');temp.write_bytes(response.content);temp.replace(path)
            rec=PARSE.parse(str(path))
            if sum(v is not None and np.isfinite(v) for v in rec['flux'])<20:
                return metadata|{'reason':'too few finite spectral bins'},None
            rec.update(mjd=int(row.mjd),phase=5,program=str(row.programname),coadd=False,
                       source='SDSS',proprietary=True,run2d=version,
                       metadata_quality_ok=bool(row.metadata_quality_ok),fieldquality=str(row.fieldquality),
                       sn_median_all=float(row.sn_median_all))
            return metadata|{'status':'available','reason':''},finite(rec)
        except Exception as exc:
            # Exception messages can include request internals; record only type.
            if attempt==2:return metadata|{'reason':type(exc).__name__},None
            time.sleep(attempt+1)


def main():
    DEST.mkdir(parents=True,exist_ok=True)
    targets=pd.read_csv(OUT/'compact_review_objects.csv')
    epochs=pd.read_parquet(OUT/'sdssv_epochs_matched.parquet')
    epochs=epochs[epochs.canonical_name.isin(targets.name)].copy()
    epochs=epochs.sort_values(['metadata_quality_ok','sn_median_all'],ascending=False).drop_duplicates(['canonical_name','mjd'])
    cache={}
    for name in epochs.canonical_name.unique():
        path=DEST/f'{name}.json'
        cache[name]=json.loads(path.read_text()) if path.exists() else []
    done={(name,int(r['mjd'])) for name,rows in cache.items() for r in rows}
    todo=[r for r in epochs.itertuples() if (r.canonical_name,int(r.mjd)) not in done]
    manifest=[];start=time.monotonic()
    print(f'{len(epochs)} distinct daily epochs for {epochs.canonical_name.nunique()} targets; {len(todo)} remaining',flush=True)
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures=[executor.submit(one,r) for r in todo]
        for k,f in enumerate(as_completed(futures),1):
            status,record=f.result();manifest.append(status)
            if record is not None:
                name=status['name'];cache[name].append(record)
                path=DEST/f'{name}.json';temp=path.with_suffix('.tmp')
                temp.write_text(json.dumps(sorted(cache[name],key=lambda r:r['mjd']),separators=(',',':'),allow_nan=False));temp.replace(path)
            if k%20==0 or k==len(todo):
                print(f'{k}/{len(todo)} attempted; {sum(s["status"]!="available" for s in manifest)} unavailable; {time.monotonic()-start:.0f}s',flush=True)
    pd.DataFrame(manifest).to_csv(DEST/'fetch_manifest.csv',index=False)
    all_records=[r for rows in cache.values() for r in rows]
    summary=dict(requested_epochs=len(epochs),available_epochs=len(all_records),
                 targets_with_spectra=sum(bool(rows) for rows in cache.values()),
                 failed_this_run=[s for s in manifest if s['status']!='available'])
    (DEST/'status.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main()
