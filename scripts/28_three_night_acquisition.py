"""Resumable acquisition for the three-night review; no publishing or commits.

Each archive attempt has its own status. A failed query is not an empty sky.
Public optical spectra and internal spectra remain in their separate caches.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from types import SimpleNamespace
import argparse, importlib, io, json, time, fcntl, threading
import numpy as np
import pandas as pd
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/reselection_2026-09-20'
CACHE=OUT/'three_night_archive_queries'
RELEASE='dr20'
LOCAL=threading.local()
PARSE=importlib.import_module('03d_fetch_spectra')
NATIVE=importlib.import_module('17_fetch_review_spectra').finite


def sdss(row):
    path=CACHE/f'{row.name}_{RELEASE}_sdss.csv';meta=path.with_suffix('.json')
    status=dict(name=row.name,archive=f'SDSS {RELEASE.upper()} allspec',status='query failed',ra=row.ra,dec=row.dec)
    if path.exists() and meta.exists():
        prior=json.loads(meta.read_text())
        if prior.get('status')=='available' and prior.get('ra')==row.ra and prior.get('dec')==row.dec:
            return pd.read_csv(path,dtype={'specobjid':str,'catalogid':str}),prior
    radius=2/3600
    half=np.rad2deg(np.arcsin(min(1.,np.sin(np.deg2rad(radius))/max(1e-12,np.cos(np.deg2rad(row.dec))))))*1.000001
    lo,hi=row.ra-half,row.ra+half
    if abs(row.dec)+radius>=90:ra_box='1=1'
    elif lo<0:ra_box=f'(ra>={lo+360:.10f} OR ra<={hi:.10f})'
    elif hi>=360:ra_box=f'(ra>={lo:.10f} OR ra<={hi-360:.10f})'
    else:ra_box=f'ra BETWEEN {lo:.10f} AND {hi:.10f}'
    # Cheap coordinate bounds avoid evaluating a spherical-distance UDF over
    # the entire allspec catalogue; the exact angular cut remains unchanged.
    sql=(f'SELECT allspec_id,sdss_phase,instrument,sdss_id,catalogid,fiberid,plate_or_fps_field,mjd,run2d,coadd,programname,survey,ra,dec,specobjid,sas_url '
         f'FROM allspec WHERE {ra_box} AND dec BETWEEN {row.dec-radius*1.000001:.10f} AND {row.dec+radius*1.000001:.10f} '
         f'AND dbo.fDistanceArcMinEq({row.ra:.8f},{row.dec:.8f},ra,dec)<0.0333333333')
    if not hasattr(LOCAL,'session'):LOCAL.session=requests.Session()
    for attempt in range(3):
        try:
            r=LOCAL.session.get(f'https://skyserver.sdss.org/{RELEASE}/SkyServerWS/SearchTools/SqlSearch',params={'cmd':sql,'format':'csv'},timeout=(20,90))
            r.raise_for_status()
            d=pd.read_csv(io.StringIO(r.text),comment='#',dtype={'specobjid':str,'catalogid':str})
            if not {'instrument','mjd','sas_url'}.issubset(d):raise ValueError('Unexpected archive response')
            d=d[d.instrument.str.lower().isin(['boss','sdss'])].copy()
            d['name']=row.name;d['source']='SDSS';d['is_coadd']=d.coadd.fillna('').astype(str).str.lower().isin(['epoch','allepoch'])
            d['proprietary']=False
            d.to_csv(path,index=False);status.update(status='available',reason='',rows=len(d),queried_utc=pd.Timestamp.now(tz='UTC').isoformat())
            meta.write_text(json.dumps(status,indent=2));return d,status
        except Exception as exc:
            status['reason']=type(exc).__name__;time.sleep(attempt+1)
    meta.write_text(json.dumps(status,indent=2));return pd.DataFrame(),status


def public_spectra(row):
    epochs,status=sdss(row)
    dest=ROOT/'data/spectra_dl'/f'{row.name}.json'
    old=json.loads(dest.read_text()) if dest.exists() else []
    fetched=[];fail=[]
    if len(epochs):
        directory=ROOT/'data/spectra_cache'/row.name;directory.mkdir(parents=True,exist_ok=True)
        for e in epochs.drop_duplicates('sas_url').itertuples():
            url=str(e.sas_url)
            if not url.startswith('https://'):continue
            if any(v.get('url')==url and v.get('cache_identity_version')==2 for v in old):continue
            # Daily and epoch coadds can share a basename. Key native files by
            # the complete archive URL, including release/reduction/coadd path.
            import hashlib
            specific=directory/hashlib.sha256(url.encode()).hexdigest()[:20]
            specific.mkdir(exist_ok=True)
            path=specific/url.rsplit('/',1)[-1]
            try:
                if not PARSE.fetch(url,str(path),tries=2):raise ValueError('No FITS returned')
                rec=PARSE.parse(str(path))
                rec.update(mjd=float(e.mjd),phase=int(e.sdss_phase),run2d=str(e.run2d),program=str(e.programname),coadd=bool(e.is_coadd),url=url,source='SDSS',proprietary=False,cache_identity_version=2)
                fetched.append(NATIVE(rec))
            except Exception as exc:fail.append(dict(mjd=float(e.mjd),reason=type(exc).__name__,url=url))
    # Existing records survive any archive failure. Coadds are labeled, not
    # interpreted as another independently observed state.
    records=old+fetched
    if fetched:
        with dest.with_suffix('.lock').open('w') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            latest=json.loads(dest.read_text()) if dest.exists() else []
            urls={r.get('url') for r in fetched}
            records=[r for r in latest if r.get('url') not in urls]+fetched
            temp=dest.with_suffix('.sdss.tmp');temp.write_text(json.dumps(NATIVE(records),separators=(',',':'),allow_nan=False));temp.replace(dest)
    status.update(files_added=len(fetched),files_failed=fail,public_records_cached=sum(not r.get('proprietary') for r in records))
    (CACHE/f'{row.name}_spectra.json').write_text(json.dumps(status,indent=2))
    return status


def main():
    global RELEASE
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['spectra','ztf']);ap.add_argument('--workers',type=int,default=6);ap.add_argument('--exclude-original-ztf',action='store_true');ap.add_argument('--targets');ap.add_argument('--release',choices=['dr19','dr20'],default='dr20');args=ap.parse_args();RELEASE=args.release
    CACHE.mkdir(exist_ok=True)
    targets=pd.read_csv(args.targets or OUT/'compact_review_objects.csv')
    if args.stage=='spectra':function=public_spectra
    else:
        module=importlib.import_module('19_refresh_review_ztf');function=module.one
        (module.CACHE/'raw').mkdir(parents=True,exist_ok=True)
        if args.exclude_original_ztf:
            original=pd.read_csv(OUT/'ztf_cache_audit_before.csv')
            targets=targets[~targets.name.isin(original.name)]
    results=[];start=time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures={ex.submit(function,r):r.name for r in targets.itertuples(index=False)}
        for i,f in enumerate(as_completed(futures),1):
            try:results.append(f.result())
            except Exception as exc:results.append(dict(name=futures[f],status='worker failed',reason=type(exc).__name__))
            if i%20==0 or i==len(futures):
                suffix='_'+RELEASE if args.stage=='spectra' else ''
                (OUT/f'three_night_{args.stage}_acquisition{suffix}.json').write_text(json.dumps(NATIVE(dict(total=len(targets),attempted=i,results=results)),indent=2))
                print(args.stage,i,'/',len(futures),pd.Series([r['status'] for r in results]).value_counts().to_dict(),round(time.monotonic()-start),'s',flush=True)
    if args.stage=='spectra':
        frames=[pd.read_csv(CACHE/f'{n}_{RELEASE}_sdss.csv',dtype={'specobjid':str,'catalogid':str}) for n in targets.name if (CACHE/f'{n}_{RELEASE}_sdss.csv').exists()]
        if frames:
            destination=ROOT/'data'/f'spectra_epochs_three_night_public_{RELEASE}.csv'
            with (OUT/f'sdss_{RELEASE}_inventory.lock').open('w') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX)
                fresh=pd.concat(frames,ignore_index=True)
                if destination.exists():
                    previous=pd.read_csv(destination,dtype={'specobjid':str,'catalogid':str})
                    fresh=pd.concat([previous[~previous.name.isin(targets.name)],fresh],ignore_index=True)
                fresh.drop_duplicates(['name','sas_url']).to_csv(destination,index=False)


if __name__=='__main__':main()
