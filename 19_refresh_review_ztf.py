"""Refresh public ZTF DR24 g/r curves for the compact review, with provenance.

Raw cone-search responses are retained. Usable cached photometry is never erased
by an HTTP/parsing failure. Empty successful queries and failures are distinct.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import io,json,threading,time,subprocess,sys
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'data/reselection_2026-09-20'
CACHE=ROOT/'data/ztf_cache/review_dr24_20260920'
COLLECTION='ztf_dr24'
URL='https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves'
LOCAL=threading.local()


def select_source(raw,ra,dec):
    if raw.empty:return raw.copy(),dict(n_raw=0,n_kept=0,n_oids=0,association_warning=False)
    needed={'oid','mjd','mag','magerr','ra','dec','filtercode','catflags','expid'}
    if not needed.issubset(raw.columns):raise ValueError('Response does not have the expected photometry columns')
    df=raw.copy()
    for c in ['mjd','mag','magerr','ra','dec','catflags']:df[c]=pd.to_numeric(df[c],errors='coerce')
    valid=(np.isfinite(df.mjd)&np.isfinite(df.mag)&np.isfinite(df.magerr)&(df.magerr>0)&
           np.isfinite(df.ra)&np.isfinite(df.dec)&df.filtercode.isin(['zg','zr'])&df.catflags.notna())
    df=df[valid].copy()
    df=df[(df.catflags.astype(np.int64)&32768)==0].copy()
    groups=['oid','filtercode']
    obj=df.groupby(groups,as_index=False).agg(ra=('ra','median'),dec=('dec','median'),n=('mjd','size'))
    if obj.empty:return df,dict(n_raw=len(raw),n_kept=0,n_oids=0,association_warning=False)
    obj['sep_arcsec']=SkyCoord(obj.ra.to_numpy()*u.deg,obj.dec.to_numpy()*u.deg).separation(SkyCoord(ra*u.deg,dec*u.deg)).arcsec
    obj=obj[obj.sep_arcsec<=1.5].copy()
    # The same source can have independent light curves in overlapping fields.
    # Within each detector/reference field keep the nearest object, not the brightest.
    keys=[k for k in ['field','ccdid','qid'] if k in df]
    obj=obj.merge(df[groups+keys].drop_duplicates(groups),on=groups,how='left')
    warning=bool(obj.duplicated(['filtercode']+keys,keep=False).any()) if keys else False
    if keys:obj=obj.sort_values(['sep_arcsec','n'],ascending=[True,False]).drop_duplicates(['filtercode']+keys)
    df=df.merge(obj[groups+['sep_arcsec']],on=groups,how='inner')
    df=df.sort_values('magerr').drop_duplicates(['filtercode','expid']).sort_values(['mjd','filtercode'])
    return df,dict(n_raw=len(raw),n_kept=len(df),n_oids=int(obj.oid.nunique()),association_warning=warning,
                   max_match_arcsec=float(obj.sep_arcsec.max()) if len(obj) else None)


def one(row):
    path=CACHE/f'{row.name}.csv';meta=path.with_suffix('.json');rawpath=CACHE/'raw'/f'{row.name}.csv'
    if path.exists() and meta.exists():
        status=json.loads(meta.read_text())
        if status.get('collection')==COLLECTION and status.get('status') in ['available','no usable photometry']:
            return status
    status=dict(name=row.name,ra=row.ra,dec=row.dec,collection=COLLECTION,status='query failed',
                queried_utc=datetime.now(timezone.utc).isoformat(),n_kept=0)
    if not hasattr(LOCAL,'session'):LOCAL.session=requests.Session()
    params={'POS':f'CIRCLE {row.ra:.7f} {row.dec:.7f} {3/3600:.9f}','BANDNAME':'g,r',
            'FORMAT':'csv','BAD_CATFLAGS_MASK':32768,'COLLECTION':COLLECTION}
    for attempt in range(3):
        try:
            r=LOCAL.session.get(URL,params=params,timeout=150);r.raise_for_status()
            if not r.text.strip() or r.text.lstrip().startswith('<'):
                raise ValueError('No valid CSV response')
            raw=pd.read_csv(io.StringIO(r.text),dtype={'oid':str,'expid':str})
            if not {'mjd','mag','filtercode'}.issubset(raw.columns):raise ValueError('Unexpected response columns')
            clean,info=select_source(raw,row.ra,row.dec)
            rawtmp=rawpath.with_suffix('.tmp');raw.to_csv(rawtmp,index=False);rawtmp.replace(rawpath)
            tmp=path.with_suffix('.tmp');clean.to_csv(tmp,index=False);tmp.replace(path)
            status.update(info,status='available' if len(clean) else 'no usable photometry',reason='',
                          n_g=int(clean.filtercode.eq('zg').sum()),n_r=int(clean.filtercode.eq('zr').sum()),
                          first_mjd=float(clean.mjd.min()) if len(clean) else None,
                          last_mjd=float(clean.mjd.max()) if len(clean) else None)
            meta.write_text(json.dumps(status,indent=2,allow_nan=False));return status
        except Exception as exc:
            status['reason']=type(exc).__name__
            if attempt<2:time.sleep(2*(attempt+1))
    meta.write_text(json.dumps(status,indent=2));return status


def main():
    (CACHE/'raw').mkdir(parents=True,exist_ok=True)
    targets=pd.read_csv(OUT/'compact_review_objects.csv')
    prior=pd.read_csv(OUT/'ztf_cache_audit_before.csv').set_index('name')
    targets['already_shown']=targets.name.map(prior.has_displayed_ztf)
    targets=targets.sort_values(['already_shown','ra'])
    results=[];start=time.monotonic()
    with ThreadPoolExecutor(max_workers=4) as executor:
        tasks=[executor.submit(one,row) for row in targets.itertuples()]
        for i,f in enumerate(as_completed(tasks),1):
            results.append(f.result())
            print(f'{i}/{len(tasks)}: {results[-1]["name"]} — {results[-1]["status"]}',flush=True)
            if i%10==0 or i==len(tasks):
                pd.DataFrame(results).to_csv(OUT/'ztf_refresh_manifest.csv',index=False)
                print(f'{i}/{len(tasks)} queried; {sum(r["status"]=="available" for r in results)} usable, '
                      f'{sum(r["status"]=="query failed" for r in results)} failed; {time.monotonic()-start:.0f}s',flush=True)
    summary=dict(collection=COLLECTION,targets=len(targets),available=sum(r['status']=='available' for r in results),
                 no_usable=sum(r['status']=='no usable photometry' for r in results),failed=[r for r in results if r['status']=='query failed'],
                 association_flags=sum(r.get('association_warning',False) for r in results))
    (OUT/'ztf_refresh_status.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)
    subprocess.run([sys.executable,str(ROOT/'16_candidate_webpage.py')],cwd=ROOT,check=True)


if __name__=='__main__':main()
