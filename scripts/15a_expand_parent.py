"""Fetch an all-RA DR16 spectral parent without photometric/Moon/variability cuts.

The Balmer follow-up branch covers 0.001<z<1.14 (Hbeta near the NGPS red edge).
Retains QSO and GALAXY spectra with AGN subclass, including missing photometry.
Original 750,414-object DR16Q catalog is separately retained by 15b as provenance.
No public/private observing products are replaced by this acquisition script.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import io, json, time
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/reselection_2026-09-20'
URL = 'https://skyserver.sdss.org/dr16/SkyServerWS/SearchTools/SqlSearch'


def fetch(a):
    path = OUT / 'catalog_chunks' / f'ra_{a:03d}.csv'
    if path.exists():
        return pd.read_csv(path, dtype={'specobjid':str, 'objid':str})
    sql = f"""SELECT CAST(s.specObjID AS varchar(25)) AS specobjid,
      CAST(s.bestObjID AS varchar(25)) AS objid, s.plate, s.mjd, s.fiberid,
      s.ra, s.dec, s.z, s.zWarning AS zwarning, s.class, s.subclass,
      s.snMedian AS sn_median, p.psfMag_g AS psfmag_g,
      p.psfMag_r AS psfmag_r, p.psfMag_i AS psfmag_i,
      p.modelMag_r AS modelmag_r, p.petroRad_r AS petrorad_r
      FROM SpecObj s LEFT JOIN PhotoObjAll p ON s.bestObjID=p.objID
      WHERE s.z BETWEEN 0.001 AND 1.14
      AND (s.class='QSO' OR (s.class='GALAXY' AND s.subclass LIKE '%AGN%'))
      AND s.ra >= {a} AND s.ra < {a+15}"""
    for attempt in range(3):
        try:
            r = requests.get(URL, params={'cmd':sql, 'format':'csv'}, timeout=180)
            r.raise_for_status()
            txt = r.text
            if txt.startswith('#Table1'): txt = txt.split('\n',1)[1]
            if not txt.startswith('specobjid,'): raise ValueError(txt[:200])
            d = pd.read_csv(io.StringIO(txt), dtype={'specobjid':str,'objid':str})
            d.to_csv(path,index=False)
            print(f'RA {a:03d}-{a+15:03d}: {len(d):,} spectra',flush=True)
            return d
        except Exception:
            if attempt==2: raise
            time.sleep(2*(attempt+1))


if __name__=='__main__':
    (OUT/'catalog_chunks').mkdir(parents=True,exist_ok=True)
    parts=[]; failures=[]
    with ThreadPoolExecutor(max_workers=4) as ex:
        fs={ex.submit(fetch,a):a for a in range(0,360,15)}
        for f in as_completed(fs):
            try: parts.append(f.result())
            except Exception as e:
                failures.append({'ra_start':fs[f],'error':str(e)})
                print('FAILED',fs[f],str(e),flush=True)
    (OUT/'catalog_query_status.json').write_text(json.dumps({'failed_chunks':failures,'completed_chunks':len(parts),'url':URL,'z_range':[.001,1.14],'photometry_cut':None,'ra_cut':None},indent=2))
    if failures: raise SystemExit('Incomplete catalog: failed chunks recorded; rerun to resume.')
    d=pd.concat(parts,ignore_index=True).drop_duplicates('specobjid').sort_values(['ra','dec','mjd'])
    d.to_csv(OUT/'dr16_balmer_parent.csv',index=False)
    print('Complete parent:',len(d),d['class'].value_counts().to_dict(),flush=True)
