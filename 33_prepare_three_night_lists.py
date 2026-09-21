"""Prepare manageable, time-balanced lists from the full screened candidate pool.

The full qualifying catalog remains saved. Preparation priorities use public
WISE changes and public baseline dates where available; spectral interpretation
and the SDSS-V timing review remain a separate local scientific step.
"""
from pathlib import Path
from datetime import timezone
import importlib,json
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
import astropy.units as u

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'data/reselection_2026-09-20'
OPS=importlib.import_module('23_operational_review');REVIEW=importlib.import_module('31_three_night_review')
OLD=importlib.import_module('07_make_webpage')
EXPOSURE=importlib.import_module('22_september_etc')


def main():
    config=json.loads((OUT/'review_selection.json').read_text())
    # Re-running this stage uses the saved full catalog, not its own short list.
    full=config.get('selection_phase')=='prepared'
    objects=pd.read_csv(OUT/('all_three_night_candidates.csv' if full else 'compact_review_objects.csv'))
    windows=pd.read_csv(OUT/('all_three_night_candidate_windows.csv' if full else 'compact_review_object_nights.csv'))
    if not full:
        objects.to_csv(OUT/'all_three_night_candidates.csv',index=False);windows.to_csv(OUT/'all_three_night_candidate_windows.csv',index=False)
    uw=pd.read_csv(OUT/'neighbour_unwise_manifest.csv').set_index('name').to_dict('index')
    fields=[REVIEW.field(n,uw) for n in objects.name]
    objects['field_status']=[f[0] for f in fields];objects['field_notes']=[f[1] for f in fields]
    aliases=pd.read_parquet(OUT/'parent_aliases.parquet');amap=aliases.drop_duplicates('name').set_index('name').canonical_name
    sdss=pd.read_csv(OUT/'dr16_balmer_parent.csv',dtype={'specobjid':str});sdss['name']=('S'+sdss.specobjid).map(amap)
    old=pd.read_csv(ROOT/'data/parent_pool_scored.csv');old['name']=('P'+old.poolid.astype(str)).map(amap);old['zwarning']=0
    public=pd.concat([sdss,old],ignore_index=True);public=public[public.sn_median.ge(5)&public.zwarning.eq(0)]
    quality=set(public.name);last=public.groupby('name').mjd.max()
    objects['public_baseline_quality']=objects.name.isin(quality)
    objects['baseline_quality']=objects.public_baseline_quality|objects.sdssv_n_quality_epochs.ge(1)
    objects['public_reference_mjd']=objects.name.map(last)
    # A generic optical spectrum can cover MgII but miss Hbeta at high z.
    # Keep such objects in the full catalogue, outside the Hbeta primary list.
    # Before a file is loaded, use conservative wavelength limits for the
    # identified reference instrument; validate against actual files afterward.
    public['comparison_red_limit']=np.where(public.plate.lt(3500),9200.,10300.)
    comparable=public[5140*(1+public.z)<=public.comparison_red_limit]
    expected=set(comparable.name)
    objects['hbeta_comparison_possible']=objects.name.isin(expected)|objects.sdssv_n_quality_epochs.ge(1)
    objects['hbeta_comparison_possible']&=(5140*(1+objects.z)<=10300)
    checked=[];accepted=[]
    for target in objects.to_dict('records'):
        paths=[ROOT/'data/spectra_dl'/f"{target['name']}.json",OUT/'sdssv_spectra'/f"{target['name']}.json"]
        exists=any(path.exists() for path in paths)
        # A DESI-only cache does not mean the identified SDSS baseline has
        # already been retrieved, so don't reject it on that basis alone.
        daily=any(any(not r.get('coadd') for r in json.loads(path.read_text())) for path in paths if path.exists())
        checked.append(exists and daily)
        accepted.append(EXPOSURE.reference(target) is not None if exists and daily else False)
    objects['hbeta_baseline_file_checked']=checked
    objects['hbeta_baseline_file_accepted']=accepted
    mask=objects.hbeta_baseline_file_checked
    objects.loc[mask,'hbeta_comparison_possible']&=objects.loc[mask,'hbeta_baseline_file_accepted']
    wise=OLD.load_wise(str(ROOT/'data/wise_cache/pool'))
    expansion=OUT/'three_night_w1.parquet'
    if expansion.exists():
        for name,g in pd.read_parquet(expansion).groupby('name'):
            wise[name]={'W1':g.sort_values('time')[['time','flux','err']].values.tolist()}
    changes=[]
    for r in objects.itertuples():
        a=np.asarray(wise.get(r.name,{}).get('W1',[]),float);value=np.nan;significant=False;last_mjd=np.nan
        if len(a)>=6:
            a=a[np.isfinite(a).all(axis=1)&(a[:,1]>0)&(a[:,2]>0)];a=a[np.argsort(a[:,0])]
            if len(a)>=6:
                first,last_flux=np.median(a[:3,1]),np.median(a[-3:,1]);last_mjd=float(a[-1,0])
                value=float(-2.5*np.log10(last_flux/first))
                error=1.0857*np.hypot(np.sqrt(np.sum(a[:3,2]**2))/3/first,np.sqrt(np.sum(a[-3:,2]**2))/3/last_flux)
                if np.isfinite(r.public_reference_mjd):
                    baseline=a[:3] if r.public_reference_mjd<a[0,0] else a[np.abs(a[:,0]-r.public_reference_mjd)<=365]
                    after=a[-3:][a[-3:,0]>=r.public_reference_mjd+180]
                    if len(baseline)>=2 and len(after)>=2:
                        f0,f1=np.median(baseline[:,1]),np.median(after[:,1])
                        delta=-2.5*np.log10(f1/f0)
                        uncertainty=1.0857*np.hypot(np.sqrt(np.sum(baseline[:,2]**2))/len(baseline)/f0,np.sqrt(np.sum(after[:,2]**2))/len(after)/f1)
                        significant=abs(delta)>=.3 and abs(delta)>=3*uncertainty
        changes.append(dict(name=r.name,unwise_end_minus_start_mag=value,unwise_last_mjd=last_mjd,dated_w1_preparation_flag=significant))
    objects=objects.merge(pd.DataFrame(changes),on='name',validate='one_to_one')
    objects['balmer_pair_in_range']=objects.Ha_A.le(10350)&objects.Hb_A.ge(3150)
    objects['preparation_score']=20*objects.dated_w1_preparation_flag.astype(int)+2*objects.balmer_pair_in_range.astype(int)+(19-objects.r_planning).clip(0,6)
    objects.to_csv(OUT/'all_three_night_candidates_screened.csv',index=False)
    usable=objects[objects.field_status.eq('clear')&objects.baseline_quality&objects.hbeta_comparison_possible].copy()
    assigned={};bundle=[];summary=[];used=set()
    for night,(date,part) in OPS.OBS.NIGHTS.items():
        t0,t1,_,_=OPS.OBS.night_window(date,part);duration=(t1-t0).to_value(u.min)
        slots=[((t0+m*u.min).to_datetime(timezone=timezone.utc),(t0+(m+30)*u.min).to_datetime(timezone=timezone.utc)) for m in np.arange(0,duration-30+1e-6,30)]
        w=windows[windows.night.eq(night)][['name','ranges_json']].merge(usable,on='name')
        candidates=w[~w.name.isin(used)].reset_index(drop=True)
        nslots=len(slots)*3;cost=np.full((nslots,len(candidates)+nslots),1e6)
        for i,(start,end) in enumerate(slots):
            ok=candidates.ranges_json_x.map(lambda value:REVIEW.covers(value,start,end)) if 'ranges_json_x' in candidates else candidates.ranges_json.map(lambda value:REVIEW.covers(value,start,end))
            for j in np.where(ok.to_numpy())[0]:
                row=candidates.iloc[j];cost[3*i:3*i+3,j]=(1000 if row.pool_role=='reserve' else 0)-row.preparation_score
        rr,cc=linear_sum_assignment(cost);selected=set();gaps=0
        for r,c in zip(rr,cc):
            good=c<len(candidates) and cost[r,c]<1e6;name=candidates.iloc[c]['name'] if good else ''
            if name:selected.add(name)
            else:gaps+=1
            start,end=slots[r//3]
            bundle.append(dict(night=night,start_pdt=start.astimezone(REVIEW.TZ).strftime('%Y-%m-%d %H:%M'),end_pdt=end.astimezone(REVIEW.TZ).strftime('%H:%M'),choice=r%3+1,name=name,
                               role=candidates.iloc[c].pool_role if good else 'gap'))
        desired=max(45,len(slots)*3)
        main=candidates[candidates.pool_role.eq('manifold')].sort_values(['preparation_score','r_planning'],ascending=[False,True])
        main_names=set(main.name)
        for name in main.name:
            if len(selected&main_names)>=desired:break
            selected.add(name)
        # A few screened outside-region backups are explicitly separate.
        reserve=candidates[candidates.pool_role.eq('reserve')].sort_values('r_planning')
        for name in reserve.name:
            if len(selected&set(reserve.name))>=6:break
            selected.add(name)
        for name in selected:assigned.setdefault(name,[]).append(night)
        used.update(selected)
        summary.append(dict(night=night,main_prepared=len(selected&main_names),reserves_prepared=len(selected&set(reserve.name)),unfilled_distinct_alternative_positions=gaps,
                            minimum_goal=desired,available_screened_main=int(w.pool_role.eq('manifold').sum())))
    prepared=objects[objects.name.isin(assigned)].copy();prepared['prepared_for_nights']=prepared.name.map(lambda n:','.join(assigned[n]))
    prepared.to_csv(OUT/'compact_review_objects.csv',index=False)
    windows[windows.name.isin(assigned)].to_csv(OUT/'compact_review_object_nights.csv',index=False)
    pd.DataFrame(bundle).to_csv(OUT/'prepared_distinct_alternatives.csv',index=False)
    pd.DataFrame(summary).to_csv(OUT/'prepared_night_summary.csv',index=False)
    config.update(selection_phase='prepared',prepared_objects=len(prepared),full_candidate_objects=len(objects),
                  hbeta_baseline_requirement='An accepted prior spectrum covering Hbeta and both comparison continua (rest 4750-4790 and 5100-5140 Angstrom), with the red continuum inside the planning bandpass. Instrument metadata screen first; loaded-file verification before final review.',
                  preparation_note='Three distinct alternatives per 30-minute block, where feasible; at least 45 manifold candidates per night, or three per block for longer nights. Lists use different objects across nights. Full qualifying catalog retained separately.',
                  priority_note='Preparation heuristic: dated significant W1 change after a public reference, Balmer-pair coverage, then brightness. This does not classify a spectral transition or forecast the present state.')
    (OUT/'review_selection.json').write_text(json.dumps(config,indent=2));print(pd.DataFrame(summary).to_string(index=False),flush=True);print('Prepared distinct objects',len(prepared),flush=True)


if __name__=='__main__':main()
