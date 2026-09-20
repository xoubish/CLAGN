"""Public-catalogue neighbour screening, with failed queries kept unknown.

SDSS primary detections within 60 arcsec; Gaia DR3 within 120 arcsec,
propagated to the observing epoch when proper motions are available; unWISE
deblending diagnostics. Thresholds are conservative inspection flags, not
predicted slit contamination or certified clean-field classifications.
"""
from pathlib import Path
from datetime import datetime,timezone
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,io,time,threading,hashlib
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'data/reselection_2026-09-20'
CACHE=OUT/'neighbour_cache'
LOCAL=threading.local()

def session():
    if not hasattr(LOCAL,'session'):LOCAL.session=requests.Session()
    return LOCAL.session

def query(name,kind,url,params):
    path=CACHE/f'{name}_{kind}.csv';meta=path.with_suffix('.json')
    if path.exists() and meta.exists():
        m=json.loads(meta.read_text())
        if m['status']=='available':return pd.read_csv(path,dtype={'objid':str,'source_id':str}),m
    status=dict(status='query failed',catalog=kind,queried_utc=datetime.now(timezone.utc).isoformat())
    for attempt in range(2):
        try:
            r=session().get(url,params=params,timeout=45);r.raise_for_status()
            df=pd.read_csv(io.StringIO(r.text),comment='#',dtype={'objid':str,'source_id':str})
            if not {'ra','dec'}.issubset(df.columns):raise ValueError('Unexpected catalogue response')
            df.to_csv(path,index=False);status.update(status='available',rows=len(df))
            meta.write_text(json.dumps(status));return df,status
        except Exception as exc:status['reason']=type(exc).__name__
    meta.write_text(json.dumps(status));return pd.DataFrame(),status

def one(row):
    name=row['name'];ra=float(row['ra']);dec=float(row['dec']);target=SkyCoord(ra*u.deg,dec*u.deg)
    sql=f'SELECT p.objid,p.ra,p.dec,p.type,p.psfMag_r,p.modelMag_r,p.clean,n.distance FROM dbo.fGetNearbyObjEq({ra:.8f},{dec:.8f},1) n JOIN PhotoObj p ON p.objid=n.objid'
    sdss,sm=query(name,'sdss','https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch',{'cmd':sql,'format':'csv'})
    adql=f"SELECT source_id,ra,dec,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,pmra,pmdec,ref_epoch FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.0333333333))"
    gaia,gm=query(name,'gaia','https://gea.esac.esa.int/tap-server/tap/sync',{'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':adql})
    neighbours=[];target_r=np.nan;target_g=np.nan;flags=[]
    if len(sdss):
        c=SkyCoord(sdss.ra.to_numpy()*u.deg,sdss.dec.to_numpy()*u.deg)
        sdss['sep']=target.separation(c).arcsec;sdss['pa']=target.position_angle(c).deg
        sdss['mag']=np.where(sdss.type.eq(6),sdss.psfMag_r,sdss.modelMag_r)
        sdss.loc[~sdss.mag.between(0,30),'mag']=np.nan
        near=sdss[sdss.sep<1].sort_values('sep')
        if len(near):target_r=float(near.iloc[0]['mag'])
        if len(near)>1:flags.append('multiple catalogue detections within 1 arcsec; check deblending')
        for r in sdss[sdss.sep>=1].itertuples():
            dr=float(r.mag-target_r) if np.isfinite(target_r) and np.isfinite(r.mag) else None
            neighbours.append(dict(catalog='SDSS',id=r.objid,ra=r.ra,dec=r.dec,sep_arcsec=float(r.sep),pa_deg=float(r.pa),mag=float(r.mag) if np.isfinite(r.mag) else None,band='r',delta_mag=dr,photometry_clean=bool(r.clean),object_type='star' if r.type==6 else 'extended'))
    if len(gaia):
        # pmra is mu_alpha*cos(dec); ten-year linear propagation suffices for
        # screening. Missing PM stays at the catalogue position and is flagged.
        years=2026.73-gaia.ref_epoch.fillna(2016)
        gaia['screen_ra']=gaia.ra+gaia.pmra.fillna(0)*years/(3.6e6*np.cos(np.deg2rad(gaia.dec)))
        gaia['screen_dec']=gaia.dec+gaia.pmdec.fillna(0)*years/3.6e6
        c=SkyCoord(gaia.screen_ra.to_numpy()*u.deg,gaia.screen_dec.to_numpy()*u.deg)
        gaia['sep']=target.separation(c).arcsec;gaia['pa']=target.position_angle(c).deg
        near=gaia[gaia.sep<1].sort_values('sep')
        if len(near):target_g=float(near.iloc[0].phot_g_mean_mag)
        for r in gaia[gaia.sep>=1].itertuples():
            mag=float(r.phot_g_mean_mag) if np.isfinite(r.phot_g_mean_mag) else None
            neighbours.append(dict(catalog='Gaia',id=r.source_id,ra=r.screen_ra,dec=r.screen_dec,sep_arcsec=float(r.sep),pa_deg=float(r.pa),mag=mag,band='G',delta_mag=float(mag-target_g) if mag is not None and np.isfinite(target_g) else None,proper_motion_available=bool(np.isfinite(r.pmra) and np.isfinite(r.pmdec))))
    # Same-band comparisons only: never subtract target r from neighbour G.
    for n in neighbours:
        sep=n['sep_arcsec'];dm=n['delta_mag'];mag=n['mag']
        reason=[]
        if dm is not None:
            if sep<3 and dm<=3:reason.append('close companion / target-deblend inspection')
            elif sep<8 and dm<=2:reason.append('nearby source: WISE blending and slice contamination review')
            elif sep<30 and dm<=1:reason.append('comparable or brighter source near the slit field; inspect position angle')
            elif sep<60 and dm<=-3:reason.append('much brighter nearby source; inspect halo and acquisition')
            elif sep<120 and dm<=-4:reason.append('much brighter source in wider field; inspect halo and scattered light')
        if mag is not None and mag<=12 and sep<120:reason.append('bright-star halo inspection')
        n['review_reasons']=reason
        flags.extend(reason)
    if not np.isfinite(target_r):flags.append('SDSS target photometry not matched; contrast uncertain')
    if sm['status']!='available' or gm['status']!='available':flags.append('catalogue query incomplete')
    result=dict(name=name,ra=ra,dec=dec,sdss_query=sm['status'],gaia_query=gm['status'],
                status='review' if flags else 'no catalogue flag',flags=sorted(set(flags)),
                target_r=target_r if np.isfinite(target_r) else None,target_g=target_g if np.isfinite(target_g) else None,
                neighbours=sorted(neighbours,key=lambda n:n['sep_arcsec']))
    (CACHE/f'{name}_screen.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    return result

def unwise(targets):
    result=[]
    for i in range(0,len(targets),20):
        sub=targets.iloc[i:i+20]
        # Batch membership is retained, so a rebuilt target list cannot reuse
        # another batch's catalogue result by accident.
        key=';'.join(sub.name);cache=CACHE/('unwise_'+hashlib.sha256(key.encode()).hexdigest()[:16]+'.csv');keyfile=cache.with_suffix('.key')
        if cache.exists() and keyfile.exists() and keyfile.read_text()==key:
            result.extend(pd.read_csv(cache).to_dict('records'));continue
        conditions=' OR '.join(f"(ra BETWEEN {r.ra-10/3600/np.cos(np.deg2rad(r.dec)):.8f} AND {r.ra+10/3600/np.cos(np.deg2rad(r.dec)):.8f} AND dec BETWEEN {r.dec-10/3600:.8f} AND {r.dec+10/3600:.8f})" for r in sub.itertuples())
        q='SELECT ra,dec,flux_w1,fracflux_w1,fracflux_w2,flags_unwise_w1,flags_info_w1 FROM unwise_dr1.object WHERE '+conditions
        rows=[]
        try:
            r=session().post('https://datalab.noirlab.edu/tap/sync',data={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':q},timeout=60);r.raise_for_status()
            d=pd.read_csv(io.StringIO(r.text));assert {'ra','dec','fracflux_w1'}.issubset(d.columns)
            coords=SkyCoord(d.ra.to_numpy()*u.deg,d.dec.to_numpy()*u.deg)
            for t in sub.itertuples():
                sep=coords.separation(SkyCoord(t.ra*u.deg,t.dec*u.deg)).arcsec
                j=int(np.argmin(sep)) if len(d) else -1
                if j>=0 and sep[j]<2:
                    v=d.iloc[j];rows.append(dict(name=t.name,status='matched',sep_arcsec=float(sep[j]),fracflux_w1=v.fracflux_w1,flags_unwise_w1=v.flags_unwise_w1))
                else:rows.append(dict(name=t.name,status='no match within 2 arcsec'))
            pd.DataFrame(rows).to_csv(cache,index=False);keyfile.write_text(key)
        except Exception as e:rows=[dict(name=t.name,status='query failed',reason=type(e).__name__) for t in sub.itertuples()]
        result.extend(rows);print('unWISE',min(i+20,len(targets)),'/',len(targets),flush=True)
    dest=OUT/'neighbour_unwise_manifest.csv';fresh=pd.DataFrame(result)
    if dest.exists():
        old=pd.read_csv(dest);fresh=pd.concat([old[~old.name.isin(fresh.name)],fresh],ignore_index=True)
    fresh.to_csv(dest,index=False)

def main():
    CACHE.mkdir(exist_ok=True)
    targets=pd.read_csv(OUT/'compact_review_objects.csv')
    standby=pd.read_csv(OUT/'bright_quasar_standby_screen.csv')
    targets=pd.concat([targets,standby],ignore_index=True).drop_duplicates('name')
    sep=pd.read_csv(OUT/'compact_review_object_nights.csv');sep=sep[sep.night.eq('sep23')].set_index('name')
    targets['sep_window']=targets.name.map(sep.preferred_longest_minutes).fillna(0)
    targets=targets.sort_values(['sep_window','r_planning'],ascending=[False,True])
    results=[];start=time.monotonic()
    with ThreadPoolExecutor(max_workers=6) as executor:
        futures=[executor.submit(one,row) for row in targets.to_dict('records')]
        for i,f in enumerate(as_completed(futures),1):
            results.append(f.result())
            if i%20==0 or i==len(futures):
                pd.DataFrame([{k:v for k,v in r.items() if k!='neighbours'} for r in results]).to_csv(OUT/'neighbour_screen_manifest.csv',index=False)
                print(f'{i}/{len(futures)} optical fields checked; {time.monotonic()-start:.0f}s',flush=True)
    unwise(targets.reset_index(drop=True))
    print('Screening complete',flush=True)

if __name__=='__main__':main()
