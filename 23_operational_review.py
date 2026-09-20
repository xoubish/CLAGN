"""Rebuild the operational pool: Moon >=40 deg, X<=1.5, explicit reserves.

Retains 60-minute eligibility, prefers 120-minute windows, and keeps known
bright quasars outside the two manifold regions in a separate reserve table.
No probability of a state change is assigned. All outputs are local.
"""
from pathlib import Path
import importlib,json
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord,AltAz,get_body
from astropy.utils import iers

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'data/reselection_2026-09-20'
OBS=importlib.import_module('05_observability')
iers.conf.auto_download=False; iers.conf.auto_max_age=None

def geometry(parent):
    pieces=[]
    for night,(date,part) in OBS.NIGHTS.items():
        t0,t1,_,_=OBS.night_window(date,part)
        duration=(t1-t0).to_value(u.min);edges=np.r_[np.arange(0,duration,5),duration];dt=np.diff(edges)
        times=t0+edges*u.min
        frame=AltAz(obstime=times,location=OBS.PALOMAR.location,pressure=0*u.hPa)
        moon=get_body('moon',times,OBS.PALOMAR.location).transform_to(frame)
        for start in range(0,len(parent),1000):
            p=parent.iloc[start:start+1000]
            aa=SkyCoord(p.ra.to_numpy()*u.deg,p.dec.to_numpy()*u.deg)[:,None].transform_to(frame)
            x=aa.secz.value;sep=aa.separation(moon).deg
            for i,r in enumerate(p.itertuples()):
                good=(aa.alt.deg[i]>0)&(x[i]>0)&(x[i]<=1.5)&(sep[i]>=40)
                accepted=good[:-1]&good[1:]
                change=np.diff(np.r_[False,accepted,False].astype(int))
                runs=list(zip(np.where(change==1)[0],np.where(change==-1)[0]))
                longest=max((edges[b]-edges[a] for a,b in runs),default=0)
                if longest<60:continue
                qualifying=[(a,b) for a,b in runs if edges[b]-edges[a]>=60]
                idx=np.concatenate([np.arange(a,b+1) for a,b in qualifying])
                near=(x[i]<=1.3)&good;near=near[:-1]&near[1:]
                best=idx[np.argmin(x[i,idx])]
                pieces.append(dict(name=r.name,night=night,preferred_longest_minutes=longest,
                    preferred_minutes=float((accepted*dt).sum()),minutes_airmass_le1p3=float((near*dt).sum()),
                    min_airmass=float(x[i,idx].min()),best_airmass_utc=times[best].isot,
                    moon_min_visible=float(sep[i,idx].min()),moon_max_visible=float(sep[i,idx].max()),
                    moon_at_best_airmass=float(sep[i,best]),selection_moon_min=40,
                    ranges_json=json.dumps([dict(start_utc=times[a].isot,end_utc=times[b].isot,minutes=float(edges[b]-edges[a])) for a,b in qualifying])))
    return pd.DataFrame(pieces)

def build_review(p=None,g=None):
    if p is None:p=pd.read_parquet(OUT/'parent_with_manifold_descriptors.parquet')
    in_region=p.in_region_clagn.eq(True)|p.in_region_zeltyn.eq(True)
    # The public DR16 parent was queried with class=QSO and zWarning=0.
    reserves=(~in_region)&p.r_planning.le(18)&p.origin.eq('old_DR16_parent')
    eligible=p[p.has_projection&p.r_planning.le(18.5)&(in_region|reserves)].copy()
    eligible['review_region']=np.select([eligible.in_region_clagn.eq(True)&eligible.in_region_zeltyn.eq(True),eligible.in_region_clagn.eq(True),eligible.in_region_zeltyn.eq(True)],['both','literature turn-on/off','Zeltyn'],default='Other manifold region')
    geom=geometry(eligible);windows=geom.merge(eligible,on='name')
    main=windows[windows.review_region.ne('Other manifold region')].copy()
    reserve=windows[windows.review_region.eq('Other manifold region')&windows.preferred_longest_minutes.ge(120)].copy()
    reserve['catalog_basis']='SDSS DR16 sciencePrimary QSO, zWarning=0'
    reserve.to_csv(OUT/'bright_quasar_reserve_object_nights.csv',index=False)
    standby_pool=reserve[reserve.night.eq('sep23')].copy()
    uw_path=OUT/'neighbour_unwise_manifest.csv'
    uw=pd.read_csv(uw_path).set_index('name') if uw_path.exists() else pd.DataFrame()
    def screen_pass(name):
        path=OUT/'neighbour_cache'/f'{name}_screen.json'
        if not path.exists() or name not in uw.index:return False
        d=json.loads(path.read_text());v=uw.loc[name]
        return (d.get('sdss_query')=='available' and d.get('gaia_query')=='available'
                and not d['flags'] and v['status']=='matched' and v.fracflux_w1>=.8
                and v.flags_unwise_w1==0)
    screened=standby_pool.name.map(screen_pass)
    standby_pool=standby_pool[screened]
    standby=standby_pool.sort_values(
        ['minutes_airmass_le1p3','preferred_longest_minutes','r_planning'],ascending=[False,False,True]).head(8)
    standby.to_csv(OUT/'bright_quasar_standby_screen.csv',index=False)
    main['pool_role']='manifold'
    prepared=reserve[reserve.name.isin(standby.name)].copy();prepared['pool_role']='reserve'
    main=pd.concat([main,prepared],ignore_index=True)
    main.to_csv(OUT/'compact_review_object_nights.csv',index=False)
    objects=main.sort_values(['name','preferred_longest_minutes'],ascending=[True,False]).drop_duplicates('name').rename(columns={'night':'best_night_by_window_length'})
    objects['eligible_nights']=objects.name.map(main.groupby('name').night.agg(lambda x:','.join(sorted(set(x)))))
    objects.sort_values(['ra','name']).to_csv(OUT/'compact_review_objects.csv',index=False)
    # Keep the complete outside-region catalogue; expose only a small screened
    # standby subset after neighbour checking, not an arbitrary science top-N.
    reserve.drop_duplicates('name').to_csv(OUT/'bright_quasar_reserve_objects.csv',index=False)
    config=dict(moon_min_deg=40,airmass_max=1.5,minimum_window_minutes=60,preferred_window_minutes=120,
                intermediate_window_minutes=90,historical_r_limit=18.5,reserve_r_limit=18.,
                reserve_minimum_window_minutes=120,version='moon40-long-window-review-2026-09-20',
                default_night='sep23',default_window_minutes=120,
                primary_gate='Neighbour screening and spectral/science review are required; no automatic observing assignment.')
    (OUT/'review_selection.json').write_text(json.dumps(config,indent=2))
    count=main.groupby(['pool_role','night']).agg(eligible=('name','nunique'),preferred=('preferred_longest_minutes',lambda x:int((x>=120).sum())))
    (OUT/'COMPACT_REVIEW.md').write_text('# Operational candidate review\n\n'+count.to_markdown()+'\n\nMoon >=40 degrees and X<=1.5 simultaneously for at least 60 minutes. Prefer >=120 minutes; 90-minute alternatives and shorter exceptions remain available. Historical r<=18.5. Neighbour screens are review flags, not certified contamination measurements. Science priorities remain pending.\n\nBright SDSS DR16 spectroscopic quasars (r<=18) in other manifold regions with >=120-minute windows are a separate reserve; they do not have the same manifold selection rationale.\n')
    print(count.to_string(),flush=True)
    print('Main objects:',int(objects.pool_role.eq('manifold').sum()),'prepared reserves:',len(standby),'reserve object-night rows:',len(reserve),flush=True)
    return main,objects,reserve

if __name__=='__main__':build_review()
