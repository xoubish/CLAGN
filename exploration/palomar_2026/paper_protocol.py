"""Palomar-only analysis protocol, timing audit and outcome-blind pilot comparisons.

This is a post-first-look protocol, not a retroactive preregistration. The 11
science cases are outcome-informed; matching uses only pre-run covariates.
"""
from pathlib import Path
import hashlib,json,sys
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from astropy.time import Time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=HERE/'outputs/palomar_paper'
SCIENCE=['P21422','P24930','P7947','P17881','P7837','P8544','P9506','P9694','P10381','P11530','P1823']
INPUTS={}


def track(p):
    p=Path(p);INPUTS[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest();return p


def frame(p,**kw):return pd.read_csv(track(p),low_memory=False,**kw)


def build():
    OUT.mkdir(parents=True,exist_ok=True)
    objects=frame(HERE/'outputs/objects.csv')
    epochs=frame(HERE/'outputs/epochs.csv')
    compact=frame(ROOT/'data/reselection_2026-09-20/compact_review_objects.csv')
    cols=['name','poolid','origin','r_planning','manifold_direction','manifold_on_neighbor_fraction',
          'manifold_off_neighbor_fraction','manifold_cl_neighbor_fraction','in_region_clagn','in_region_zeltyn']
    sample=objects.merge(compact[cols],on='name',how='left')
    proposal=frame(ROOT/'jwst_proposal/inputs/jwst_sample_cycle6.csv').set_index('internal_id')
    # Only historical optical brightness is used to match the pilot.
    sample['r_pre_run']=sample.r_planning.fillna(sample.name.map(proposal.r_mag))
    sample['r_provenance']=np.where(sample.r_planning.notna(),'compact pre-run planning','proposal historical r')
    sample['latest_archive_mjd']=sample.name.map(epochs[epochs.source!='NGPS'].groupby('name').max_mjd.max())
    sample['first_ngps_mjd']=sample.name.map(epochs[epochs.source=='NGPS'].groupby('name').min_mjd.min())
    sample['comparison_baseline_years']=(sample.first_ngps_mjd-sample.latest_archive_mjd)/365.25/(1+sample.z)
    # Recover the time support of the original pool W1 series, not the later NEOWISE curve.
    wanted=set(sample.poolid.dropna().astype(int));cache=[]
    for p in sorted((ROOT/'data/wise_cache/pool').glob('*.parquet')):
        d=pd.read_parquet(p).reset_index()
        m=d[d.objectid.isin(wanted)&d.band.eq('WISE_W1')]
        if len(m):track(p);cache.append(m.assign(input_file=str(p.relative_to(ROOT))))
    raw=pd.concat(cache,ignore_index=True).drop_duplicates(['objectid','time'])
    histories=[]
    for oid,g in raw.groupby('objectid'):
        good=np.isfinite(g.flux)&np.isfinite(g.err)&g.err.gt(0)&g.flux.gt(0)
        g=g[good].sort_values('time')
        if len(g)<4:continue
        mag=-2.5*np.log10(g.flux.to_numpy())
        histories.append(dict(poolid=oid,w1_input_first_mjd=g.time.min(),w1_input_last_mjd=g.time.max(),
            w1_n=len(g),w1_p95_p05_mag=np.ptp(np.quantile(mag,[.05,.95])),
            w1_endpoint_delta_mag=float(np.median(mag[-3:])-np.median(mag[:3])),
            w1_input_files=';'.join(sorted(set(g.input_file)))))
    sample=sample.merge(pd.DataFrame(histories),on='poolid',how='left')
    sample['feature_support_status']=np.where(sample.origin.eq('old_DR16_parent')&sample.w1_input_last_mjd.notna(),
        'original pool W1 cache linked by poolid; historical snapshot still to verify','not established for this projection')
    established=sample.feature_support_status.str.startswith('original')
    sample['temporal_forecast_candidate']=established & sample.latest_archive_mjd.ge(sample.w1_input_last_mjd) & sample.latest_archive_mjd.lt(sample.first_ngps_mjd)
    sample['timing_interpretation']=np.select([
        ~established, sample.latest_archive_mjd.isna(), sample.temporal_forecast_candidate],
        ['feature cutoff unresolved','no accepted spectral baseline','baseline follows all linked W1 measurements'],
        default='baseline precedes some W1 inputs; change may predate selection')
    sample['predicted_direction']=np.where(sample.manifold_direction>0,'strengthening-associated neighborhood',
        np.where(sample.manifold_direction<0,'weakening-associated neighborhood',
        np.where(sample.manifold_direction.eq(0),'tied neighborhood','unavailable')))
    # These fractions describe nearby labeled training sources; they are not event probabilities.
    sample['prediction_is_probability']=False
    candidates=sample[sample.experiment_group.str.startswith('D:') & ~sample.name.isin(SCIENCE)
        & sample.known_state_status.ne('catalog-confirmed CLAGN') & sample.n_archive.gt(0)].dropna(subset=['z','r_pre_run','comparison_baseline_years'])
    anchors=sample[sample.name.isin(SCIENCE)].dropna(subset=['z','r_pre_run','comparison_baseline_years'])
    cost=np.full((len(candidates),len(anchors)),1e6)
    for i,c in enumerate(candidates.itertuples(index=False)):
        for j,a in enumerate(anchors.itertuples(index=False)):
            dz=abs(c.z-a.z);dm=abs(c.r_pre_run-a.r_pre_run)
            dt=abs(np.log10((c.comparison_baseline_years+.1)/(a.comparison_baseline_years+.1)))
            if dz<=.15 and dm<=1.5 and dt<=1:
                cost[i,j]=(dz/.1)**2+(dm/1.)**2+(dt/.6)**2
    ii,jj=linear_sum_assignment(cost)
    matched=[]
    for i,j in sorted(zip(ii,jj),key=lambda p:cost[p]):
        if cost[i,j]>=1e6 or len(matched)==4:continue
        c,a=candidates.iloc[i],anchors.iloc[j]
        matched.append(dict(comparison=c['name'],anchor=a['name'],distance=cost[i,j],
            delta_z=c.z-a.z,delta_r=c.r_pre_run-a.r_pre_run,
            baseline_ratio=(c.comparison_baseline_years+.1)/(a.comparison_baseline_years+.1)))
    matches=pd.DataFrame(matched)
    pilot=sample[sample.name.isin(SCIENCE+matches.comparison.tolist())].copy()
    pilot['pilot_role']=np.where(pilot.name.isin(SCIENCE),'outcome-informed science pilot','pre-covariate matched comparison')
    # Every NGPS visit gets a reference chosen by date/quality, not the largest change.
    ep=[]
    for name in pilot.name:
        g=epochs[epochs.name==name];old=g[g.source!='NGPS'].sort_values('max_mjd')
        chosen=set(g.loc[g.source=='NGPS','epoch_id'])
        if len(old):chosen.add(old.iloc[0].epoch_id)
        for n in g[g.source=='NGPS'].itertuples(index=False):
            before=old[old.max_mjd<n.min_mjd]
            if len(before):chosen.add(before.iloc[-1].epoch_id)
        ep.append(g[g.epoch_id.isin(chosen)])
    selected=pd.concat(ep,ignore_index=True)
    protocol=dict(date='2026-11-03',scope='Palomar spectral evolution; dust-physics work deferred to JWST',
        status='post-first-look analysis protocol; not preregistered',
        science_cases=SCIENCE,comparison_selection='D group; excluding catalog-confirmed and science cases; optimal distinct matches using z, historical r, and baseline only',
        calipers=dict(abs_z=.15,abs_r_mag=1.5,abs_log10_baseline=1.),
        primary_question='Does frozen manifold information improve prediction of broad-Hbeta evolution beyond simple pre-run W1 variability?',
        primary_outcome='Continuous broad Hbeta evolution, with host/FeII/calibration-aware uncertainty; no CL label inferred from a wing index',
        reference_rule='latest accepted archive ending before first relevant NGPS visit; first archive secondary',
        covariates=['redshift','historical optical brightness','spectral baseline','W1 amplitude','W1 direction','line detectability'],
        restrictions=['pilot is outcome-informed and cannot estimate discovery yield',
          'post-run locked snapshot cannot establish genuine pre-run prediction provenance',
          'known CLAGN separate; no extrapolation of targeted-sample fractions to all AGN',
          'forecast eligibility is only a timing screen until original feature/training overlap audit passes'],
        inputs=INPUTS)
    lock=OUT/'protocol_snapshot.json'
    if lock.exists():
        prior=json.loads(lock.read_text())
        if prior['inputs']!=INPUTS:raise RuntimeError('Inputs changed: preserve prior snapshot and create a versioned protocol.')
    sample.to_csv(OUT/'analysis_sample.csv',index=False)
    pilot.to_csv(OUT/'pilot_sample.csv',index=False)
    matches.to_csv(OUT/'pilot_matches.csv',index=False)
    selected.to_csv(OUT/'pilot_epochs.csv',index=False)
    lock.write_text(json.dumps(protocol,indent=2)+'\n')
    print('Pilot sources',len(pilot),'epochs',len(selected))
    print(matches.to_string(index=False))
    print(sample.timing_interpretation.value_counts().to_string())
    return sample,pilot,selected


if __name__=='__main__':build()
