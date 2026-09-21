"""Independent science experiment review; never edits the telescope packet.

The observable parent includes all already projected AGNs with historical r<19.5,
not just the prepared lists. Missing evidence is explicitly distinct from a null
result. Cached spectra supply screening indices, NOT changing-look classifications.
All outputs can contain collaboration data and remain in the ignored observing tree.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, importlib, json, warnings
from functools import lru_cache
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz, get_body, search_around_sky
from astropy.time import Time
from astropy.utils import iers

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data/reselection_2026-09-20'
DEST = ROOT/'observing/science_review'
OBS = importlib.import_module('05_observability')
WEB = importlib.import_module('16_candidate_webpage')
WEB.review_ztf = lru_cache(maxsize=1)(WEB.review_ztf)
TRIAGE = importlib.import_module('31_three_night_review')
iers.conf.auto_download = False
iers.conf.auto_max_age = None


def finite(value):
    if isinstance(value, dict): return {k: finite(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [finite(v) for v in value]
    if isinstance(value, np.generic): value=value.item()
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def inventory():
    DEST.mkdir(parents=True, exist_ok=True)
    p=pd.read_parquet(OUT/'parent_with_manifold_descriptors.parquet')
    gate=p.agn_identity_eligible.fillna(False)&p.r_planning.lt(19.5)&p.z.between(.001,1.02)
    parent=p[gate&p.has_projection.fillna(False)].copy()
    parent.to_parquet(DEST/'projected_parent_snapshot.parquet',index=False)
    # Deliberately do not require membership of either hand-drawn UMAP region.
    rows=[]
    for night,(date,part) in OBS.NIGHTS.items():
        t0,t1,_,_=OBS.night_window(date,part)
        duration=(t1-t0).to_value(u.min)
        edges=np.r_[np.arange(0,duration,5),duration]; times=t0+edges*u.min
        frame=AltAz(obstime=times,location=OBS.PALOMAR.location,pressure=0*u.hPa)
        moon=get_body('moon',times,OBS.PALOMAR.location).transform_to(frame)
        for start in range(0,len(parent),1000):
            sub=parent.iloc[start:start+1000]
            aa=SkyCoord(sub.ra.to_numpy()*u.deg,sub.dec.to_numpy()*u.deg)[:,None].transform_to(frame)
            x=aa.secz.value; sep=aa.separation(moon).deg
            for i,name in enumerate(sub.name):
                for tier,limit in [('preferred',1.5),('extended',1.8)]:
                    good=(aa.alt.deg[i]>0)&(x[i]>=1)&(x[i]<=limit)&(sep[i]>=40)
                    accepted=good[:-1]&good[1:]
                    changes=np.diff(np.r_[False,accepted,False].astype(int))
                    runs=[(a,b) for a,b in zip(np.where(changes==1)[0],np.where(changes==-1)[0]) if edges[b]-edges[a]>=20]
                    for a,b in runs:
                        rows.append(dict(name=name,night=night,tier=tier,start_utc=times[a].isot,
                                         end_utc=times[b].isot,minutes=float(edges[b]-edges[a]),
                                         min_airmass=float(x[i,a:b+1].min()),moon_min=float(sep[i,a:b+1].min())))
        print('Geometry',night,len(rows),flush=True)
    windows=pd.DataFrame(rows);windows.to_csv(DEST/'windows.csv',index=False)
    parent=parent[parent.name.isin(windows.name)].copy()
    parent.to_parquet(DEST/'observable_parent.parquet',index=False)
    meta=dict(created_utc=datetime.now(timezone.utc).isoformat(),parent_rows=len(p),
              projected_brightness_redshift_eligible=int(gate[p.has_projection.fillna(False)].sum()),
              missing_projection_before_geometry=int((gate&~p.has_projection.fillna(False)).sum()),
              observable_projected=len(parent),night_counts=windows.groupby('night').name.nunique().to_dict(),
              r_limit=19.5,redshift_range=[.001,1.02],minimum_window_minutes=20,moon_min_deg=40,
              airmass_tiers=[1.5,1.8],geometry_grid_minutes=5,
              caveat='Five-minute conservative windows are screening only; recheck any final visit at finer cadence. Unprojected identities remain outside this first-pass manifold experiment, not rejected scientifically.',
              broad_parent_pipeline=json.loads((OUT/'completion_pipeline_status.json').read_text()))
    (DEST/'scope.json').write_text(json.dumps(meta,indent=2))


def records(name):
    result=[]
    for path in [ROOT/'data/spectra_dl'/f'{name}.json',OUT/'sdssv_spectra'/f'{name}.json']:
        if path.exists(): result.extend(json.loads(path.read_text()))
    return result


def measure(record,z):
    """Conservative metadata/coverage screening and transparent line-wing index.

    Rebinned plotting caches have no reliable propagated pixel variance. Therefore
    indices carry no fitted significance and cannot certify a spectral state.
    """
    meta=record.get('meta',{}); reasons=[]
    if record.get('metadata_quality_ok') is False: reasons.append('metadata quality flag')
    if str(record.get('fieldquality','')).lower()=='bad': reasons.append('bad field')
    if meta.get('zwarning') not in [None,0,0.]: reasons.append('redshift warning')
    sn=record.get('sn_median_all',meta.get('sn_median_all'))
    if sn is not None and np.isfinite(sn) and sn<5: reasons.append('catalogue S/N below 5')
    lo=record.get('min_mjd');hi=record.get('max_mjd');coadd=bool(record.get('coadd'))
    lo=lo if lo is not None and np.isfinite(lo) else None
    hi=hi if hi is not None and np.isfinite(hi) else None
    # A coadd can be useful, but its filename date is not an independent epoch.
    if coadd and (lo is None or hi is None): reasons.append('coadd time interval unknown')
    if lo is not None and hi is not None and hi-lo>30: reasons.append('coadd spans over 30 days')
    mjd=hi if coadd and hi is not None else record.get('mjd')
    if mjd is None or not np.isfinite(mjd) or not 40000<mjd<70000:
        mjd=None;reasons.append('unknown observation date')
    w=np.asarray(record['wave'],float);f=np.asarray(record['flux'],float);rest=w/(1+z)
    good=np.isfinite(f)&np.isfinite(w)
    a=f[good&(rest>=4750)&(rest<=4790)];b=f[good&(rest>=5100)&(rest<=5140)]
    info=dict(mjd=mjd,date=Time(mjd,format='mjd').strftime('%Y-%m-%d') if mjd else '',
              proprietary=bool(record.get('proprietary')),source=record.get('source','unknown'),
              coadd=coadd,min_mjd=lo,max_mjd=hi,catalog_sn=sn,reference_url=record.get('url',''))
    if len(a)<4 or len(b)<4:
        reasons.append('Hbeta continuum coverage insufficient')
    else:
        c=np.interp(rest,[4770,5120],[np.median(a),np.median(b)])
        local=float(np.interp(4862.7,[4770,5120],[np.median(a),np.median(b)]))
        if local<=0: reasons.append('nonpositive continuum')
        else:
            info['continuum_flux']=local
            info['continuum_AB']=-2.5*np.log10(local*1e-17*(4862.7*(1+z))**2/2.99792458e18)-48.6
            # Omit the narrow Hbeta core (+/-~600 km/s); FeII and host absorption
            # remain possible contaminants, so this is only a visual-review aid.
            mask=good&(((rest>=4800)&(rest<=4853))|((rest>=4873)&(rest<=4930)))
            info['hbeta_wing_EW_index']=float(np.sum((f[mask]-c[mask])/local*np.gradient(rest)[mask])) if mask.sum()>=8 else np.nan
            core=good&(rest>=4994)&(rest<=5020)
            info['oiii_flux_index']=float(np.sum((f[core]-c[core])*np.gradient(w)[core])) if core.sum()>=3 else np.nan
    info['accepted']=not reasons
    info['rejection_reasons']='; '.join(reasons)
    return info


def evidence():
    parent=pd.read_parquet(DEST/'observable_parent.parquet')
    windows=pd.read_csv(DEST/'windows.csv')
    nights=windows.groupby('name').night.agg(lambda v:','.join(sorted(set(v)))).to_dict()
    sv=pd.read_parquet(OUT/'sdssv_epochs_matched.parquet')
    svgood=sv[sv.metadata_quality_ok.fillna(False)]
    latest_catalog=svgood.groupby('canonical_name').mjd.max().to_dict()
    svquality=sv.groupby(['canonical_name','mjd']).metadata_quality_ok.any().to_dict()
    field=pd.read_csv(OUT/'all_three_night_candidates_screened.csv').set_index('name')
    known={};matches=[]
    coords=SkyCoord(parent.ra.to_numpy()*u.deg,parent.dec.to_numpy()*u.deg)
    for label,path,ra,dec in [('compilation',ROOT/'data/external/clagn_catalog_camus_panda2026.csv','ra_deg','dec_deg'),
                              ('literature',ROOT/'data/literature_clagn.csv','ra','dec'),
                              ('Zeltyn',ROOT/'data/zeltyn_coords.csv','ra','dec')]:
        cat=pd.read_csv(path);cat=cat[cat[ra].notna()&cat[dec].notna()].reset_index(drop=True)
        ii,jj,sep,_=search_around_sky(coords,SkyCoord(cat[ra].to_numpy()*u.deg,cat[dec].to_numpy()*u.deg),2*u.arcsec)
        for i,j,dist in zip(ii,jj,sep.arcsec):
            r=cat.iloc[j];n=parent.iloc[i]['name'];status=r.get('confirmation_status',r.get('class_zeltyn','unverified literature match'))
            matches.append(dict(name=n,catalog=label,status=status,sep_arcsec=dist,transition=r.get('transition_type',''),reference=r.get('reference',r.get('ref','Zeltyn 2024'))))
            value='catalog-confirmed CLAGN' if status in ['spectroscopic_confirmed','CL-AGN'] else 'reported candidate' if 'candidate' in str(status).lower() else 'literature match; confirmation needs audit'
            if known.get(n)!='catalog-confirmed CLAGN':known[n]=value
    pd.DataFrame(matches).to_csv(DEST/'literature_matches.csv',index=False)
    frames=[]
    for path in list((ROOT/'data').glob('neowise_visits_*.csv'))+list(OUT.glob('neowise_visits_prepared*.csv')):
        try:
            x=pd.read_csv(path)
            if {'name','mjd','w1'}.issubset(x): frames.append(x)
        except pd.errors.EmptyDataError: pass
    neo=pd.concat(frames,ignore_index=True).drop_duplicates(['name','mjd'],keep='last')
    neowise=dict(tuple(neo.groupby('name')))
    epochs=[];rows=[]
    old={v['name'] for v in json.loads((ROOT/'observing/sep23/packet.json').read_text())['primaries']}
    for i,t in enumerate(parent.to_dict('records')):
        name=t['name']; rr=records(name)
        ee=[]
        for r in rr:
            r=dict(r)
            if r.get('phase')==5 and r.get('mjd') is not None:
                q=svquality.get((name,r['mjd']))
                if q is not None and not q:r['metadata_quality_ok']=False
            ee.append(measure(r,t['z']))
        epochs.extend([dict(name=name,**e) for e in ee])
        accepted=[e for e in ee if e['accepted'] and e['mjd']]
        latest=max(accepted,key=lambda e:e['mjd']) if accepted else {}
        # Independent dates, not reductions/coadds, are the displayed count.
        ndays=len({int(e['mjd']) for e in accepted})
        region=bool(t['in_region_clagn'] or t['in_region_zeltyn'])
        t['pool_role']='manifold' if region else 'comparison'
        ev=TRIAGE.science(t,{name:{'reference_mjd':latest.get('mjd',np.nan)}},known,neowise)
        ev.update({k:t[k] for k in ['ra','dec','z','r_planning','r_source','manifold_cl_neighbor_fraction','manifold_on_neighbor_fraction','manifold_off_neighbor_fraction','umap_x','umap_y']})
        ev.update(manifold_region=region,old_primary=name in old,raw_spectra=len(rr),accepted_dates=ndays,
                  baseline_status='available for line review' if latest else 'missing or rejected',
                  last_reference_private=latest.get('proprietary',False),continuum_AB=latest.get('continuum_AB',np.nan),
                  latest_hbeta_wing_index=latest.get('hbeta_wing_EW_index',np.nan),
                  latest_good_catalog_mjd=latest_catalog.get(name,np.nan),
                  newer_good_spectrum_not_reviewed=bool(latest_catalog.get(name,0)>latest.get('mjd',0)+1),
                  field_status=field.loc[name,'field_status'] if name in field.index else 'pending',
                  field_notes=field.loc[name,'field_notes'] if name in field.index else 'Needs optical and WISE field screening',
                  eligible_nights=nights[name],
                  broad_line_state='not classified; inspect profiles and native errors',
                  repeat_transition_status='not established')
        series,info=WEB.review_ztf(name)
        ev['optical_association_warning']=bool(info.get('association_warning',False))
        measurable=any(np.isfinite(ev.get(f'ztf_{b}_change_after_reference',np.nan)) for b in ['g','r'])
        ev['optical_timing_test']='flagged change' if ev['post_spectrum_optical_trigger'] else 'measurable; below screening threshold' if measurable else 'not testable with cached baseline photometry'
        ev['photometry_age_days']=(Time('2026-09-23').mjd-Time(ev['ztf_last_date']).mjd) if ev.get('ztf_last_date') else np.nan
        changes=[ev.get(f'ztf_{b}_change_after_reference',np.nan) for b in ['g','r']]
        ev['optical_change_abs']=float(np.nanmax(np.abs(changes))) if np.isfinite(changes).any() else np.nan
        ev['balmer_pair']=6564.61*(1+t['z'])<10350
        ev['hbeta_red_difficult']=4862.7*(1+t['z'])>9000
        ev['historical_faint_extension']=bool(t['r_planning']>=19)
        ev['current_faint_exception']=bool(ev.get('ztf_r_latest180_mag',0)>=19)
        dated=bool(ev['post_spectrum_optical_trigger'] or ev['post_spectrum_ir_flag'])
        ev['experiment_group']=('A: manifold + dated change' if dated else 'B: manifold without dated trigger') if region else ('C: variability comparison' if dated else 'D: comparison pool without dated trigger')
        # This score orders evidence review, not CL probability or telescope time.
        # Missing spectra have a separate acquisition queue, not a zero science value.
        ev['review_priority']= (100*bool(ev['post_spectrum_optical_trigger'])+50*bool(ev['post_spectrum_ir_flag'])
                               +10*min(1,float(t.get('manifold_cl_neighbor_fraction') or 0))
                               +5*ev['balmer_pair'])
        if ev['post_spectrum_optical_trigger']:
            ev['palomar_test']='Measure broad Balmer response to '+ev['trigger_direction']+' since the accepted reference; establish prior broad-line state before assigning on/off.'
        elif ev['post_spectrum_ir_flag']:
            ev['palomar_test']='Test broad Balmer change across the dated IR evolution; dust lag and post-2024 state remain unresolved.'
        elif region:
            ev['palomar_test']='Blind manifold validation against a measured spectral baseline; no direction forecast from coordinates alone.'
        else:
            ev['palomar_test']='Potential matched comparison; no control status until matched on spectral state and observable properties.'
        rows.append(ev)
        if (i+1)%1000==0:
            print('Evidence',i+1,'/',len(parent),flush=True)
            pd.DataFrame(rows).to_csv(DEST/'partial_evidence.csv',index=False)
    d=pd.DataFrame(rows).sort_values(['review_priority','name'],ascending=[False,True])
    d.to_csv(DEST/'all_candidate_evidence.csv',index=False)
    pd.DataFrame(epochs).to_csv(DEST/'spectral_epoch_audit.csv',index=False)
    usable=d[d.accepted_dates.gt(0)&~d.optical_association_warning]
    for night in OBS.NIGHTS:
        q=usable[usable.eligible_nights.str.contains(night)].copy()
        q.to_csv(DEST/f'{night}_science_ranking.csv',index=False)
    # Balance review effort by experiment group and night, without locking visits.
    proposals=[]
    for night in OBS.NIGHTS:
        q=usable[usable.eligible_nights.str.contains(night)]
        for group,n in [('A: manifold + dated change',10),('B: manifold without dated trigger',5),('C: variability comparison',5)]:
            g=q[q.experiment_group.eq(group)].sort_values(['field_status','review_priority'],ascending=[True,False])
            for r in g.head(n).to_dict('records'): proposals.append(dict(review_night=night,**r))
    pd.DataFrame(proposals).to_csv(DEST/'discussion_shortlist.csv',index=False)
    print('Observable',len(d),'cached baselines',len(usable),'groups',usable.experiment_group.value_counts().to_dict(),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['inventory','evidence']);args=ap.parse_args()
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Tried to get polar motions')
        globals()[args.stage]()
