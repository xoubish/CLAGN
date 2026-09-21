"""Cache public ALeRCE/ZTF alerts for a local science-review target list.

Alert detection photometry is selected by variability. It is not a complete
forced-photometry light curve, and missing alerts never mean unchanged flux.
Corrected magnitudes are kept separate from DR24 catalogue photometry.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import argparse, json, threading
import requests
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'observing/science_review'
CACHE=DEST/'alerts'; LOCAL=threading.local()
BASE='https://api.alerce.online/ztf/v1/objects/'


def get(url,params=None):
    if not hasattr(LOCAL,'session'): LOCAL.session=requests.Session()
    r=LOCAL.session.get(url,params=params,timeout=40);r.raise_for_status();return r.json()


def one(t):
    path=CACHE/f'{t.name}.json'
    if path.exists():return json.loads(path.read_text())['summary']
    result=dict(name=t.name,ra=t.ra,dec=t.dec,queried_utc=datetime.now(timezone.utc).isoformat(),
                caveat='Alert-selected detections; corrected flux can be affected by reference/host systematics. No alerts is not a nonvariability constraint.')
    try:
        objects=get(BASE,dict(ra=t.ra,dec=t.dec,radius=2,page_size=100))
        result['objects']=objects;result['detections']={};allrows=[]
        for obj in objects.get('items',[]):
            sep=SkyCoord(t.ra*u.deg,t.dec*u.deg).separation(SkyCoord(obj['meanra']*u.deg,obj['meandec']*u.deg)).arcsec
            if sep>1.5:continue
            oid=obj['oid'];rows=get(BASE+oid+'/detections');result['detections'][oid]=rows
            allrows.extend([dict(oid=oid,match_arcsec=sep,**r) for r in rows])
        d=pd.DataFrame(allrows)
        summary=dict(name=t.name,status='available' if len(d) else 'no matching alerts',n_raw=len(d))
        if len(d):
            d=d.drop_duplicates('candid')
            for col in ['mjd','magpsf_corr','sigmapsf_corr_ext','drb','distnr']:
                if col not in d:d[col]=np.nan
                d[col]=pd.to_numeric(d[col],errors='coerce')
            # Reference-corrected photometry only; never compare difference
            # magnitudes directly to total-light catalogue or spectral magnitudes.
            good=(d.magpsf_corr.notna()&d.sigmapsf_corr_ext.between(0,.2,inclusive='neither')
                  &d.get('corrected',pd.Series(False,index=d.index)).eq(True)
                  &~d.get('dubious',pd.Series(True,index=d.index)).fillna(True)
                  &d.drb.ge(.8)&d.distnr.le(1.5))
            clean=d[good].copy()
            clean.to_csv(CACHE/f'{t.name}_corrected.csv',index=False)
            summary['n_corrected_usable']=len(clean)
            summary['last_alert_date']=Time(d.mjd.max(),format='mjd').strftime('%Y-%m-%d')
            for fid,band in [(1,'g'),(2,'r')]:
                b=clean[clean.fid.eq(fid)].copy()
                if b.empty:continue
                b['day']=np.floor(b.mjd);night=b.groupby('day').magpsf_corr.median()
                recent=night[night.index>=max(Time('2026-01-01').mjd,night.index.max()-180)]
                old=night[(night.index>=Time('2025-01-01').mjd)&(night.index<Time('2026-01-01').mjd)]
                summary[f'{band}_last_corrected_date']=Time(night.index.max(),format='mjd').strftime('%Y-%m-%d')
                summary[f'{band}_recent_nights']=len(recent)
                summary[f'{band}_recent_median']=float(recent.median()) if len(recent)>=3 else None
                summary[f'{band}_2026_minus_2025']=float(recent.median()-old.median()) if len(recent)>=3 and len(old)>=3 else None
        result['summary']=summary
        path.write_text(json.dumps(result,indent=2,allow_nan=False))
        return summary
    except Exception as exc:
        summary=dict(name=t.name,status='query failed',reason=type(exc).__name__)
        # Keep failures separate so they can be retried.
        (CACHE/f'{t.name}_failure.json').write_text(json.dumps(summary));return summary


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('targets');a=ap.parse_args()
    CACHE.mkdir(parents=True,exist_ok=True)
    targets=pd.read_csv(a.targets).drop_duplicates('name');results=[]
    with ThreadPoolExecutor(max_workers=4) as ex:
        futures=[ex.submit(one,t) for t in targets.itertuples(index=False)]
        for i,f in enumerate(as_completed(futures),1):
            results.append(f.result());print(i,'/',len(futures),results[-1]['name'],results[-1]['status'],flush=True)
    # Include earlier successful searches in the combined review manifest.
    allresults={r['name']:r for r in results}
    for path in CACHE.glob('*.json'):
        d=json.loads(path.read_text())
        if 'summary' in d:allresults.setdefault(d['name'],d['summary'])
    pd.DataFrame(allresults.values()).to_csv(DEST/'public_alert_summary.csv',index=False)
