"""Rebuild the operational pool: Moon >=40 deg, X<=1.5, explicit reserves.

Uses 30-minute eligibility for a 2x600-second visit, and keeps known
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

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/reselection_2026-09-20'
OBS=importlib.import_module('05_observability')
iers.conf.auto_download=False; iers.conf.auto_max_age=None

def geometry(parent, minimum_minutes=30):
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
                if longest<minimum_minutes:continue
                qualifying=[(a,b) for a,b in runs if edges[b]-edges[a]>=minimum_minutes]
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
    eligible=p[p.has_projection & p.r_planning.lt(19) & p.Hb_in_NGPS.eq(True)].copy()
    if 'agn_identity_eligible' in eligible:eligible=eligible[eligible.agn_identity_eligible].copy()
    eligible['review_region']=np.select([eligible.in_region_clagn.eq(True)&eligible.in_region_zeltyn.eq(True),eligible.in_region_clagn.eq(True),eligible.in_region_zeltyn.eq(True)],['both','literature turn-on/off','Zeltyn'],default='Other manifold region')
    if g is None or not {'ranges_json','selection_moon_min'}.issubset(g.columns):
        g=geometry(eligible,minimum_minutes=30)
    else:
        assert g.selection_moon_min.eq(40).all(), 'Geometry must use the active Moon cut.'
    windows=g.merge(eligible,on='name',validate='many_to_one')
    main=windows[windows.review_region.ne('Other manifold region')].copy()
    # Reserve preparation is driven by coverage in time, never by longest window.
    reserve=windows[windows.review_region.eq('Other manifold region') & windows.r_planning.lt(19) & windows.origin.isin(['old_DR16_parent','expanded_DR16_QSO','full_DR16Q_catalog'])].copy()
    reserve['catalog_basis']='Public SDSS spectroscopic quasar; reference metadata quality checked during preparation'
    reserve.to_csv(OUT/'bright_quasar_reserve_object_nights.csv',index=False)
    reserve.drop_duplicates('name').to_csv(OUT/'bright_quasar_reserve_objects.csv',index=False)
    selected=set()
    for night,(date,part) in OBS.NIGHTS.items():
        t0,t1,_,_=OBS.night_window(date,part)
        rows=reserve[reserve.night.eq(night)].copy()
        for minute in np.arange(0,(t1-t0).to_value(u.min)-30+1e-6,30):
            begin=(t0+minute*u.min).isot;end=(t0+(minute+30)*u.min).isot
            available=rows[rows.ranges_json.map(lambda value:any(r['start_utc']<=begin and r['end_utc']>=end for r in json.loads(value)))]
            # Screen a wider reserve queue to obtain three clean alternatives
            # after neighbour checks. Inspection remains a
            # separate status: this is not a claim that these are clean fields.
            selected.update(available.sort_values(['r_planning','min_airmass','name']).head(30).name)
    prepared=reserve[reserve.name.isin(selected)].copy()
    prepared.to_csv(OUT/'bright_quasar_standby_screen.csv',index=False)
    main['pool_role']='manifold';prepared['pool_role']='reserve'
    main=pd.concat([main,prepared],ignore_index=True)
    main.to_csv(OUT/'compact_review_object_nights.csv',index=False)
    objects=main.sort_values(['name','preferred_longest_minutes'],ascending=[True,False]).drop_duplicates('name').rename(columns={'night':'best_night_by_window_length'})
    objects['eligible_nights']=objects.name.map(main.groupby('name').night.agg(lambda x:','.join(sorted(set(x)))))
    objects.sort_values(['ra','name']).to_csv(OUT/'compact_review_objects.csv',index=False)
    config=dict(moon_min_deg=40,airmass_max=1.5,minimum_window_minutes=30,preferred_window_minutes=30,
                historical_r_limit=19,brightness_comparison='strictly less than',reserve_r_limit=19.,
                reserve_minimum_window_minutes=30,version='three-nights-2x600-review-2026-09-20',
                default_night='sep23',default_window_minutes=30,exposures=2,exposure_seconds=600,
                bin_spatial=2,bin_spectral=2,slice_arcsec=1.,planning_block_minutes=30,
                overhead_note='Ten minutes provisionally allowed for acquisition, reads and margin; verify on site.',
                primary_gate='Candidate review pool. Require a usable spectral baseline, field inspection and H-beta sensitivity review before assigning an observation.')
    (OUT/'review_selection.json').write_text(json.dumps(config,indent=2))
    count=main.groupby(['pool_role','night']).agg(eligible=('name','nunique'))
    (OUT/'COMPACT_REVIEW.md').write_text('# Three-night candidate review\n\n'+count.to_markdown()+'\n\nHistorical r<19; X<=1.5 and topocentric Moon separation>=40 degrees simultaneously for at least 30 minutes. Planned 2x600 seconds, 2x2 binning and 1 arcsec slice. Thirty minutes includes a provisional ten-minute overhead/margin budget. All projected candidates in either manifold region are retained; no ranking by spectral count or two-hour cut. Bright outside-region quasars are separately labeled reserves. Field flags, archival baseline quality and H-beta sensitivity still require review. The larger parent and pending W1 projections remain in separate audit tables.\n')
    print(count.to_string(),flush=True)
    print('Main objects:',int(objects.pool_role.eq('manifold').sum()),'prepared reserves:',len(selected),flush=True)
    return main,objects,reserve

if __name__=='__main__':build_review()
