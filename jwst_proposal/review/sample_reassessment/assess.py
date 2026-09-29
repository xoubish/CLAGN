"""Attach visual assessments and evaluate conditional timing sensitivity.

Assessments reflect inspection of the saved four-band plots and five available
2026 alert extensions. These are recommendations for further modelling, not a
validated replacement observing list. No proposal or input sample is changed.
"""
from pathlib import Path
import os, sys, json, hashlib
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/clagn-matplotlib')
import numpy as np
import pandas as pd
from astropy.time import Time
OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
sys.path.insert(0,str(OUT.parent/'history_audit'))
from response_stress_test import exponential, uniform, weights

def recent_updates():
    """Reuse Figure 1's documented IRSA/alert matching for five local targets."""
    sys.path.insert(0,str(ROOT/'jwst_proposal'))
    import make_fig1_five as fig1
    import matplotlib.pyplot as plt
    ids=['F02','F04','F06','F07','R01']
    sample=pd.read_csv(ROOT/'jwst_proposal/inputs/jwst_sample_cycle6.csv').set_index('id')
    fig,axes=plt.subplots(3,2,figsize=(12,10),sharex=True)
    inventory=[];binned=[];hashes={}
    for ax,ident in zip(axes.flat,ids):
        for suffix in ['ztf_irsa','alerce']:
            p=ROOT/f'jwst_proposal/inputs/lightcurves/{ident}_{suffix}.csv'
            hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
        for band,(d,a) in fig1.ztf(ident).items():
            color={'g':'#008562','r':'#c34054'}[band]
            ax.scatter(Time(d[:,0],format='mjd').decimalyear,d[:,1],s=2,c=color,alpha=.25)
            if len(a):ax.scatter(Time(a[:,0],format='mjd').decimalyear,a[:,1],s=8,facecolors='none',edgecolors=color,label=band+' alerts')
            for kind,arr in [('IRSA',d),('alerts',a)]:
                if len(arr):inventory.append(dict(id=ident,band=band,source=kind,n_nights=len(arr),last_year=float(Time(arr[-1,0],format='mjd').decimalyear),last_180d_median_mag=float(np.median(arr[arr[:,0]>=arr[-1,0]-180,1]))))
            combined=np.vstack([d,a]) if len(a) else d
            c=pd.DataFrame(combined,columns=['mjd','mag']);c['bin']=np.floor(c.mjd/90)
            c=c.groupby('bin').agg(mjd=('mjd','median'),mag=('mag','median'),n=('mag','size'))
            c=c[c.n>=3].copy();c['name']=sample.loc[ident,'internal_id'];c['band']=band
            binned.append(c.reset_index(drop=True))
        ax.invert_yaxis();ax.set_title(ident);ax.legend(fontsize=8);ax.grid(alpha=.15);ax.set_xlim(2018,2026.9)
    axes.flat[-1].set_visible(False);fig.suptitle('Available alert extensions: photometry to September 2026')
    fig.tight_layout();fig.savefig(OUT/'current_recent_updates.png',dpi=130);plt.close(fig)
    pd.DataFrame(inventory).to_csv(OUT/'recent_alert_inventory.csv',index=False)
    updates=pd.concat(binned,ignore_index=True);updates.to_csv(OUT/'recent_alert_binned.csv',index=False)
    (OUT/'recent_alert_provenance.json').write_text(json.dumps(dict(inputs=hashes,
        method='Figure 1 IRSA nightly medians and post-IRSA corrected alerts, offset matching where >=5 pairs; 90-day medians with >=3 nights. Different bin weights/object handling from original single-OID audit; for sensitivity only.'),indent=2)+'\n')
    return updates

# Priority means inspect/model first, not demonstrated detectability.
CURRENT={
'F01':('priority','Sustained optical decline since 2018-2019; W1/W2 decline after a broad high state.','Sparse optical bins; refine event timing and current state.'),
'F02':('hold','Long infrared decline but renewed optical brightening; 2026 alerts support the recovery.','Do not count as a clean continuing fade without updated joint modelling.'),
'F03':('exclude_from_clean_core','Recurrent infrared peaks and 2021-2022 optical brightening followed by fading.','Endpoint decline hides a recurrent history.'),
'F04':('exclude_from_clean_core','Infrared decline followed by a large optical recovery persisting into 2026.','Opposite recent optical and historical infrared directions.'),
'F05':('hold','Strong historical decline, with renewed optical brightening at the end of the saved record.','Establish whether the 2025 recovery persists in 2026.'),
'F06':('exclude_from_clean_core','2018-2019 flare and 2022-2023 recovery; September 2026 data remain at a lower state.','Useful complex-history object, not a single sustained decline.'),
'F07':('exclude_from_clean_core','Largest infrared decline precedes the optical record; later optical variability is recurrent and smaller.','The 2026 alerts are sparse and do not reconstruct the early transition.'),
'F08':('exclude_from_clean_core','Infrared decline followed by recovery and renewed fading.','Multiple phases weaken a simple sustained-fade interpretation.'),
'F09':('priority','NGC 2617: optical and infrared decline into a low state around 2022-2025.','Earlier fluctuations remain; old shell forecast is strongly host-limited.'),
'F10':('exclude_from_clean_core','Early infrared flare and later recurrence; optical record misses the main event.','No well-observed single driver transition.'),
'F11':('exclude_from_clean_core','Recurrent optical and infrared episodes with a later optical decline.','Complex response history, not a clean persistent fading trajectory.'),
'F12':('exclude_from_clean_core','Early decline, recurrent infrared excursions and optical recovery.','Current direction and original class differ.'),
'R01':('priority_with_downturn','P9694: retained optical/infrared rise, but optical emission declines from its 2025 high into 2026.','Keep for modelling; do not describe as monotonically or currently rising.'),
'R02':('hold','Rise to a 2020-2021 infrared high state followed by a decline; optical plateau with fluctuations.','May retain a brighter state, but not a clean continuing rise.'),
'R03':('hold','Strong infrared rise from the 2020 trough and a supporting but sparse optical rise.','Earlier decline; retained optical bins stop in 2023. Update current state.'),
'R04':('exclude_from_clean_core','Optical high state around 2019-2022 followed by substantial fading.','Do not retain as a rising target or relabel without event-response modelling.'),
'R05':('hold','Overall infrared rise with large optical dips and recovery.','Sustained mean brightening possible; inspect recurring phases before inclusion.'),
'R06':('exclude_from_clean_core','Infrared rise and renewed decline with very sparse optical coverage.','Both classification and recent driving history are uncertain.'),
'R07':('priority','J101152.98+544206.4: optical rise to a sustained high plateau; delayed W1/W2 brightening.','A clear transition is not a guarantee of surviving warm memory by JWST.'),
'R08':('hold','Modest infrared rise and sparse, mixed optical changes.','Needs denser current-state coverage and a host-light assessment.'),
'R09':('exclude_from_clean_core','Rise toward 2020 followed by fluctuations and a late infrared decline.','Original rising label does not describe the latest infrared direction.'),
'R10':('exclude_from_clean_core','Recovery from an infrared trough with little net long-baseline change.','A recovery experiment rather than a sustained-rise anchor.'),
'R11':('exclude_from_clean_core','Major infrared peak precedes ZTF; subsequent decline and optical recovery.','Main driver episode is unobserved optically.'),
'R12':('hold','Large retained infrared rise; optical variability is modest, with dips and a later downturn.','Promising infrared history, but establish the nuclear driving history.'),
}

REPLACEMENTS={
'P20778':('priority','Large prolonged decline, broadly supported by both optical bands; infrared fading continues through 2024.'),
'P16919':('reserve','Early infrared decline and later broad hump; optical settles into a low plateau. Only one catalogue spectrum.'),
'P11719':('reserve','Most decline is old; recent optical upturn warrants a current-state check.'),
'P14318':('priority','Both optical bands and infrared decline to a retained low state; superposed fluctuations.'),
'P18762':('priority','Broad multi-year decline in optical and infrared; intermediate rebounds need modelling.'),
'P22405':('reserve','Long infrared decline, but a strong 2019 optical flare dominates the observed optical history.'),
'P24210':('reserve','Infrared decline with sparse and uneven optical coverage; g and r not equally constraining.'),
'P13653':('exclude_from_clean_core','Large optical flare and infrared recovery/decline; not a clean single-direction history.'),
'P16096':('priority','Long infrared decline with a further 2022-2023 drop and broadly supportive optical fading.'),
'P14398':('reserve','Clear drop to a low state around 2019, then near-constant infrared; timing leverage may be small.'),
'P13844':('reserve','Gradual modest decline with appreciable fluctuations; establish nuclear contrast.'),
'P8748':('reserve','Infrared decline with fluctuations; optical late record nearly flat.'),
'P25402':('exclude_from_clean_core','Strong 2018-2019 flare, later recurrence and late infrared upturn.'),
'P16524':('priority','Gradual infrared fading with a recent additional decline; supporting optical trend. One catalogue spectrum.'),
'P20262':('reserve','Decline into a plateau after 2021; old high state and late fluctuations require timing check.'),
'P14736':('reserve','Main decline around 2019; later optical variability and small infrared changes.'),
'P22654':('reserve','Net decline with a substantial 2022 recovery before renewed fading.'),
'P5293':('reserve','Long decline but recurrent optical episodes and late upturn.'),
'P17336':('reserve','Overall decline with large fluctuations and earlier infrared recovery.'),
'P21253':('reserve','Abrupt 2019 drop followed by fluctuations and optical recovery; assess event age.'),
'P14641':('reserve','Strong net infrared rise, but optical dips return near the earlier level before renewed brightening.'),
'P12673':('priority','Optical and infrared rise into a retained high state around 2021-2022; minor late infrared retreat.'),
'P21702':('reserve','Earlier rise with recurrent dips; recent record largely a plateau.'),
'P12180':('reserve','Rise retained overall, but substantial optical fluctuations and recent downturn.'),
'P20003':('reserve','Retained brighter state with repeated optical dips; only one catalogue spectrum.'),
'P7750':('priority','Large optical rise around 2020 followed by infrared rise; high state retained despite a later dip. One catalogue spectrum.'),
'P5609':('reserve','Optical and infrared rise then recent strong optical downturn; update before using as rising.'),
'P635':('reserve','Rise followed by a retreat from the optical peak; current state must be checked.'),
'P21631':('reserve','Overall rise with large fluctuations and a recent optical decline.'),
'P20180':('reserve','Modest infrared rise and recurrent optical variability; weak single-event interpretation.'),
'P19631':('priority','Strong optical rise since 2019-2020 and continuing infrared rise; small late optical retreat.'),
'P24649':('reserve','Long infrared rise, but large optical fluctuations and gaps.'),
'P22507':('exclude_from_clean_core','2023 optical flare followed by decline; would reintroduce the flare concern.'),
'P14144':('reserve','Retained optical step up around 2021; modest infrared amplitude and sparse latest optical coverage.'),
}

def main():
    sample=pd.read_csv(ROOT/'jwst_proposal/inputs/jwst_sample_cycle6.csv')
    metrics=pd.read_csv(OUT/'current_metrics.csv')
    forecast=pd.read_csv(ROOT/'jwst_proposal/inputs/warm_response_forecast.csv')
    rows=[]
    for r in metrics.to_dict('records'):
        status,reason,need=CURRENT[r['id']]
        r.update(review_status=status,assessment=reason,required_check=need)
        f=forecast[forecast.id==r['id']].iloc[0]
        vals=[f['delayed_minus_instant'+suffix] for suffix in ['', '_short','_long']]
        r.update(old_shell_memory_min_dex=min(vals),old_shell_memory_max_dex=max(vals),old_shell_sigma_dex=f.sigma_dex)
        rows.append(r)
    current=pd.DataFrame(rows);current.to_csv(OUT/'current_target_assessment.csv',index=False)
    replacements=pd.read_csv(OUT/'replacement_candidates.csv')
    replacements=replacements[replacements.optical_status=='supports_direction'].copy()
    assert set(replacements.name)==set(REPLACEMENTS)
    replacements['review_status']=replacements.name.map(lambda n:REPLACEMENTS[n][0])
    replacements['assessment']=replacements.name.map(lambda n:REPLACEMENTS[n][1])
    replacements['spherex_status']='not_verified_for_replacement'
    replacements.to_csv(OUT/'replacement_assessment.csv',index=False)
    priority=replacements[replacements.review_status=='priority'].sort_values(['direction','history_rank'],ascending=[True,False])
    priority.to_csv(OUT/'priority_replacements.csv',index=False)
    import screen
    allcurves=pd.read_csv(OUT/'replacement_binned_measurements.csv')
    grouped={k:g for k,g in allcurves.groupby('name')}
    shown=pd.read_csv(OUT/'replacement_visual_review.csv')
    remaining=replacements[~replacements.name.isin(shown.name)].sort_values(['direction','history_rank'],ascending=[True,False])
    screen.atlas(remaining,grouped,'additional_candidates')
    screen.atlas(priority,grouped,'priority_histories')

    # Same transparent linear-response stress test for current and candidate histories.
    # No host subtraction, no UV conversion, and no source-specific lag assigned.
    old=pd.read_csv(OUT.parent/'history_audit/binned_measurements.csv')
    old['name']=old.internal_id
    updates=recent_updates()
    old=old[~(old.name.isin(updates.name)&old.band.isin(['g','r']))]
    new=pd.read_csv(OUT/'replacement_binned_measurements.csv')
    new=new[new.name.isin(replacements.name)]
    curves=pd.concat([old[['name','band','mjd','mag']],new[['name','band','mjd','mag']],updates[['name','band','mjd','mag']]],ignore_index=True)
    meta=pd.concat([sample.rename(columns={'internal_id':'name','family':'direction'})[['name','z','direction']],replacements[['name','z','direction']]]).drop_duplicates('name').set_index('name')
    rows=[]
    for (name,band),d in curves.groupby(['name','band']):
        if band not in ['g','r'] or len(d)<4:continue
        d=d.sort_values('mjd');z=meta.loc[name,'z'];sign=1 if meta.loc[name,'direction']=='fade' else -1
        t=(d.mjd.to_numpy()-d.mjd.min())/(365.25*(1+z))
        f=10**(-.4*(d.mag.to_numpy()-d.mag.median()))
        f[-1]=np.median(f[d.mjd>=d.mjd.max()-400])
        for date in ['2027-07-01','2028-01-01','2028-06-30']:
            now=(Time(date).mjd-d.mjd.min())/(365.25*(1+z))
            for tau in [1.,3.,5.,10.]:
                for kernel,fn in [('exponential',exponential),('uniform',uniform)]:
                    value=float(np.log10(fn(t,f,now,tau)/f[-1]));pre,inside,future=weights(t,now,tau,kernel)
                    rows.append(dict(name=name,band=band,date=date,tau_rest_yr=tau,kernel=kernel,
                        history_contrast_dex=value,sign_for_selected_direction=sign,signed_contrast_dex=sign*value,
                        prehistory_weight=pre,record_window_weight=inside,unobserved_late_weight=future,
                        n_bins=len(d),last_observed_year=float(Time(d.mjd.max(),format='mjd').decimalyear)))
    timing=pd.DataFrame(rows);timing.to_csv(OUT/'timing_sensitivity.csv',index=False)
    selected=set(priority.name)|set(current[current.review_status.str.startswith('priority')].name)
    summary=[]
    for name,d in timing[(timing.name.isin(selected))&(timing.date=='2028-01-01')].groupby('name'):
        # Prefer r at comparable coverage; use g if it extends >=0.5 yr later.
        # This retains the 2026 g-band update for R01 rather than stale r alone.
        r=d[d.band=='r'];g=d[d.band=='g']
        band='r' if len(r) and r.n_bins.iloc[0]>=10 else 'g' if len(g) and g.n_bins.iloc[0]>=10 else d.sort_values(['n_bins','last_observed_year'],ascending=False).band.iloc[0]
        if len(g) and g.n_bins.iloc[0]>=10 and (not len(r) or g.last_observed_year.iloc[0]>r.last_observed_year.iloc[0]+.5):band='g'
        r=d[d.band==band]
        if r.empty:continue
        row=dict(name=name,band=band,n_bins=int(r.n_bins.iloc[0]),last_observed_year=r.last_observed_year.iloc[0])
        for tau in [1.,3.,5.,10.]:
            g=r[r.tau_rest_yr==tau]
            row[f'contrast_tau{int(tau)}_min_dex']=g.signed_contrast_dex.min()
            row[f'contrast_tau{int(tau)}_max_dex']=g.signed_contrast_dex.max()
            row[f'record_weight_tau{int(tau)}_min']=g.record_window_weight.min()
        summary.append(row)
    pd.DataFrame(summary).to_csv(OUT/'priority_timing_summary.csv',index=False)

    # Sanity checks on the analysis products, not a scientific validation claim.
    assert len(current)==24 and current.id.is_unique
    assert len(replacements)==34 and replacements.name.is_unique
    assert np.isfinite(timing.history_contrast_dex).all()
    assert np.allclose(timing.prehistory_weight+timing.record_window_weight+timing.unobserved_late_weight,1)
    assert set(priority.name).isdisjoint(set(sample.internal_id))
    checks=dict(current_targets=len(current),replacement_candidates=len(replacements),
                current_review_counts=current.review_status.value_counts().to_dict(),
                replacement_review_counts=replacements.groupby(['direction','review_status']).size().to_dict().__str__(),
                priority_replacements=len(priority),timing_rows=len(timing),passed=True,
                meaning='Bookkeeping and numerical checks only; not validated MIRI detectability.',
                source_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'assessment_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(current.review_status.value_counts().to_string())
    print(priority[['name','direction','ra','dec','z','w3_mjy','n_spec']].round(4).to_string(index=False))
    print(pd.DataFrame(summary).round(3).to_string(index=False))

if __name__=='__main__':main()
