"""Query DESI DR1 coadds and retrieve every matched SPARCL spectrum.

DESI coadds retain their individual mean/min/max dates and target IDs. They are
not treated as individual exposures or assigned another coadd's median date.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import hashlib,io,json,time,importlib,fcntl,re,threading,argparse
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u
from sparcl.client import SparclClient

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'data/reselection_2026-09-20'
CACHE=OUT/'three_night_desi';CACHE.mkdir(exist_ok=True)
PARSE=importlib.import_module('03d_fetch_spectra');NATIVE=importlib.import_module('17_fetch_review_spectra').finite
LOCAL=threading.local();RATE_LOCK=threading.Lock();NEXT_CALL=0.


def paced(method,**kwargs):
    global NEXT_CALL
    for attempt in range(5):
        with RATE_LOCK:
            time.sleep(max(0,NEXT_CALL-time.monotonic()));NEXT_CALL=time.monotonic()+1.5
        try:return method(**kwargs)
        except Exception as exc:
            if attempt==4:raise
            seconds=re.search(r'in (\d+) seconds',str(exc))
            time.sleep(int(seconds.group(1))+2 if seconds else 3*(attempt+1))


def client():
    if not hasattr(LOCAL,'client'):
        LOCAL.client=paced(SparclClient,connect_timeout=10,read_timeout=90,announcement=False)
    return LOCAL.client


def query(sub):
    key=hashlib.sha256(sub[['name','ra','dec']].to_csv(index=False).encode()).hexdigest()[:20]
    path=CACHE/f'{key}.csv';meta=path.with_suffix('.json')
    if path.exists() and meta.exists() and json.loads(meta.read_text())['status']=='available':return pd.read_csv(path,dtype={'targetid':str})
    conditions=[]
    for r in sub.itertuples():
        dr=2/3600/np.cos(np.deg2rad(r.dec));dd=2/3600
        conditions.append(f'(mean_fiber_ra BETWEEN {r.ra-dr:.9f} AND {r.ra+dr:.9f} AND mean_fiber_dec BETWEEN {r.dec-dd:.9f} AND {r.dec+dd:.9f})')
    sql='SELECT targetid,mean_fiber_ra AS ra,mean_fiber_dec AS dec,z,zwarn AS zwarning,spectype AS class,survey,program,mean_mjd AS mjd,min_mjd,max_mjd,coadd_numexp,coadd_numnight FROM desi_dr1.zpix WHERE '+' OR '.join(conditions)
    status=dict(status='query failed',targets=sub.name.tolist())
    for attempt in range(3):
        try:
            response=requests.post('https://datalab.noirlab.edu/tap/sync',data={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':sql},timeout=90);response.raise_for_status()
            d=pd.read_csv(io.StringIO(response.text),dtype={'targetid':str})
            if not {'targetid','ra','dec','mjd'}.issubset(d):raise ValueError('Unexpected DESI response')
            if len(d):
                idx,sep,_=SkyCoord(d.ra.to_numpy()*u.deg,d.dec.to_numpy()*u.deg).match_to_catalog_sky(SkyCoord(sub.ra.to_numpy()*u.deg,sub.dec.to_numpy()*u.deg))
                d['name']=sub.name.to_numpy()[idx];d['match_arcsec']=sep.arcsec;d=d[d.match_arcsec<=2].copy()
            else:d['name']=pd.Series(dtype=str)
            d['source']='DESI';d['is_coadd']=True;d['proprietary']=False
            d.to_csv(path,index=False);meta.write_text(json.dumps(status|dict(status='available',rows=len(d)),indent=2));return d
        except Exception as exc:status['reason']=type(exc).__name__;time.sleep(attempt+1)
    meta.write_text(json.dumps(status,indent=2));return pd.DataFrame()


def fetch(item):
    name,epochs=item;manifest=CACHE/f'{name}_files.json';records=[];failed=[]
    for e in epochs.drop_duplicates(['targetid','survey','program']).itertuples():
        cache=CACHE/f'{name}_{e.targetid}_{e.survey}_{e.program}.v2.json'
        if cache.exists():records.extend(json.loads(cache.read_text()));continue
        try:
            found=paced(client().find,outfields=['sparcl_id','specid','ra','dec'],constraints={'data_release':['DESI-DR1'],'specid':[str(e.targetid)]},limit=1000)
            if not found.ids:raise ValueError('Catalog match has no SPARCL file')
            got=paced(client().retrieve,uuid_list=found.ids,dataset_list=['DESI-DR1'],include=['sparcl_id','specid','data_release','redshift','spectype','flux','wavelength','ivar','ra','dec','survey','program','mean_mjd','coadd_numnight','file'],limit=1000)
            current=[]
            for r in got.records:
                if str(r.survey)!=str(e.survey) or str(r.program)!=str(e.program):continue
                if abs(float(r.mean_mjd)-float(e.mjd))>.05:raise ValueError('DESI mean date does not agree with catalog row')
                sep=SkyCoord(float(r.ra)*u.deg,float(r.dec)*u.deg).separation(SkyCoord(e.ra*u.deg,e.dec*u.deg)).arcsec
                if sep>2:raise ValueError('Spectrum position mismatch')
                lam=np.asarray(r.wavelength);flux=np.asarray(r.flux,float);iv=np.asarray(r.ivar,float);good=(iv>0)&np.isfinite(flux)
                idx=np.searchsorted(PARSE.GRID,lam)-1;fb=np.full(len(PARSE.GRID),np.nan)
                for i in np.unique(idx[(idx>=0)&(idx<len(fb))]):
                    valid=(idx==i)&good
                    if valid.any():fb[i]=np.median(flux[valid])
                current.append(NATIVE(dict(wave=PARSE.GRID.tolist(),flux=fb.tolist(),meta={'class':str(r.spectype),'z':float(r.redshift),'specid':str(r.specid),'zwarning':int(e.zwarning)},
                    source='DESI',coadd=True,proprietary=False,mjd=float(r.mean_mjd),min_mjd=float(e.min_mjd),max_mjd=float(e.max_mjd),coadd_numnight=int(e.coadd_numnight),survey=str(r.survey),program=str(r.program),archive_file=str(r.file),date_verified=True,
                    metadata_quality_ok=bool(e.zwarning==0),url=f'sparcl:{r.sparcl_id}',lines={},ew={})))
            if not current:raise ValueError('No SPARCL record matching this survey/program coadd')
            cache.write_text(json.dumps(current,separators=(',',':'),allow_nan=False));records.extend(current)
        except Exception as exc:failed.append(dict(targetid=str(e.targetid),reason=type(exc).__name__,detail=str(exc)[:180]))
    if records:
        path=ROOT/'data/spectra_dl'/f'{name}.json'
        with path.with_suffix('.lock').open('w') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            old=json.loads(path.read_text()) if path.exists() else []
            urls={r['url'] for r in records}
            old=[r for r in old if r.get('url') not in urls]
            tmp=path.with_suffix('.desi.tmp');tmp.write_text(json.dumps(NATIVE(old+records),separators=(',',':'),allow_nan=False));tmp.replace(path)
    result=dict(name=name,records=len(records),failed=failed);manifest.write_text(json.dumps(result,indent=2));return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--targets');args=ap.parse_args()
    targets=pd.read_csv(args.targets or OUT/'compact_review_objects.csv');frames=[]
    with ThreadPoolExecutor(max_workers=3) as ex:
        futures=[ex.submit(query,targets.iloc[i:i+20]) for i in range(0,len(targets),20)]
        for i,f in enumerate(as_completed(futures),1):
            frames.append(f.result());print('DESI inventory batch',i,'/',len(futures),flush=True)
    good=[d for d in frames if len(d)]
    if not good:print('No matched DESI records; consult query status for failures.');return
    epochs=pd.concat(good,ignore_index=True).drop_duplicates(['name','targetid','survey','program'])
    epochs.to_csv(ROOT/'data/spectra_epochs_three_night_desi.csv',index=False)
    results=[]
    with ThreadPoolExecutor(max_workers=2) as ex:
        futures={ex.submit(fetch,item):item[0] for item in epochs.groupby('name')}
        for i,f in enumerate(as_completed(futures),1):
            try:results.append(f.result())
            except Exception as exc:results.append(dict(name=futures[f],status='worker failed',reason=type(exc).__name__))
            if i%20==0 or i==len(futures):
                (OUT/'three_night_desi_status.json').write_text(json.dumps(results,indent=2));print('DESI files',i,'/',len(futures),flush=True)


if __name__=='__main__':main()
