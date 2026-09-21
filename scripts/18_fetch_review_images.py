"""Fetch verified SDSS cutouts for the current compact candidate pool.

An absent file or failed request is never interpreted as missing sky coverage.
Existing thumbnails are preserved if a refresh fails. No commit or publication.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import io
import json
import threading
import time
import argparse
import pandas as pd
import requests
from PIL import Image, ImageStat

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/reselection_2026-09-20'
CACHE=ROOT/'data/cutouts'
URL='https://skyserver.sdss.org/dr18/SkyServerWS/ImgCutout/getjpeg'
LOCAL=threading.local()


def valid_image(data):
    try:
        with Image.open(io.BytesIO(data)) as img:
            img.load()
            return min(img.size)>=100 and max(ImageStat.Stat(img.convert('RGB')).stddev)>.5
    except Exception:
        return False


def one(row,field_arcsec=40,pixels=256):
    kind='sdss_wide' if field_arcsec>40 else 'sdss'
    path=CACHE/f'{row.name}_{kind}.jpg'
    meta_path=path.with_suffix('.json')
    if path.exists() and meta_path.exists() and valid_image(path.read_bytes()):
        m=json.loads(meta_path.read_text())
        if m.get('ra')==row.ra and m.get('dec')==row.dec and m.get('source')=='SDSS DR18' and m.get('field_arcsec')==field_arcsec:return m
    if not hasattr(LOCAL,'session'):LOCAL.session=requests.Session()
    result=dict(name=row.name,ra=row.ra,dec=row.dec,status='request failed',source='SDSS DR18',
                field_arcsec=field_arcsec,pixels=pixels,updated_utc=datetime.now(timezone.utc).isoformat())
    params=dict(ra=row.ra,dec=row.dec,scale=field_arcsec/pixels,width=pixels,height=pixels)
    for attempt in range(3):
        try:
            r=LOCAL.session.get(URL,params=params,timeout=40)
            if r.ok and r.headers.get('content-type','').startswith('image/') and valid_image(r.content):
                tmp=path.with_suffix('.tmp');tmp.write_bytes(r.content);tmp.replace(path)
                result.update(status='available',reason='',bytes=len(r.content))
                meta_path.write_text(json.dumps(result,indent=2));return result
            result['reason']=f'HTTP {r.status_code}; no validated image returned'
        except Exception as exc:
            result['reason']=type(exc).__name__
        if attempt<2:time.sleep(attempt+1)
    if path.exists() and valid_image(path.read_bytes()):
        result.update(status='cached image; refresh failed',source='SDSS archival cache',field_arcsec=None,pixels=None)
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--wide',action='store_true');args=ap.parse_args()
    CACHE.mkdir(parents=True,exist_ok=True)
    targets=pd.read_csv(OUT/'compact_review_objects.csv')
    results=[];start=time.monotonic()
    with ThreadPoolExecutor(max_workers=4) as executor:
        tasks=[executor.submit(one,row,240 if args.wide else 40,384 if args.wide else 256) for row in targets.itertuples()]
        for i,f in enumerate(as_completed(tasks),1):
            results.append(f.result())
            if i%25==0 or i==len(tasks):
                print(f'{i}/{len(tasks)} checked; {sum(r["status"]=="available" for r in results)} verified SDSS images; {time.monotonic()-start:.0f}s',flush=True)
    pd.DataFrame(results).sort_values('name').to_csv(OUT/'image_manifest.csv',index=False)
    summary=dict(targets=len(targets),verified_sdss=sum(r['status']=='available' for r in results),
                 failed=[r for r in results if r['status']!='available'])
    (OUT/'image_fetch_status.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main()
