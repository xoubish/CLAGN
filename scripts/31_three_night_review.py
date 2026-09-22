"""Dated science triage and time-balanced alternatives for all three nights.

No automated changing-look classifications or calibrated transition probabilities.
Exposure sensitivity is reported, not silently used as a broad-line S/N cut.
All products stay in the ignored local research directory.
"""
from pathlib import Path
from datetime import datetime,timezone,timedelta
from zoneinfo import ZoneInfo
import importlib,json
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from astropy.time import Time
import astropy.units as u
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/reselection_2026-09-20'
DEST=OUT/'three_night_review';DEST.mkdir(exist_ok=True)
OBS=importlib.import_module('05_observability');WEB=importlib.import_module('16_candidate_webpage')
TZ=ZoneInfo('America/Los_Angeles')


def field(name,unwise):
    path=OUT/'neighbour_cache'/f'{name}_screen.json'
    if not path.exists():return 'pending','Optical field not checked'
    d=json.loads(path.read_text());flags=list(d['flags']);uw=unwise.get(name,{})
    if uw.get('status')=='matched' and (uw.get('fracflux_w1',0)<.8 or uw.get('flags_unwise_w1',1)!=0):flags.append('unWISE deblending/quality inspection')
    if flags:return 'review','; '.join(flags)
    if d.get('sdss_query')!='available' or d.get('gaia_query')!='available' or uw.get('status')!='matched':return 'pending','Catalogue coverage incomplete'
    return 'clear','No catalogue flag; visual/slit inspection still required'


def science(target,sensitivity,known,neowise):
    name=target['name'];row=dict(name=name)
    ref=sensitivity.get(name,{})
    mjd=ref.get('reference_mjd',np.nan);row['reference_mjd']=mjd
    row['reference_date']=Time(mjd,format='mjd').strftime('%Y-%m-%d') if pd.notna(mjd) else ''
    ztf,info=WEB.review_ztf(name, bin_days=1);changes=[]
    row['ztf_status']=info.get('status');row['ztf_last_date']=info.get('last_date')
    for band in ['g','r']:
        a=np.asarray(ztf.get(band,[]),float)
        if len(a)<3:continue
        a=a[np.isfinite(a[:,:3]).all(axis=1)];latest=a[a[:,0]>=a[:,0].max()-180]
        row[f'ztf_{band}_latest180_mag']=float(np.median(latest[:,1]))
        row[f'ztf_{band}_span_days']=float(np.ptp(a[:,0]))
        if not np.isfinite(mjd):continue
        baseline=a[np.abs(a[:,0]-mjd)<=45]
        recent=latest[latest[:,0]>=mjd+90]
        if len(baseline)<3 or len(recent)<3:continue
        delta=float(np.median(recent[:,1])-np.median(baseline[:,1]))
        def uncertainty(x):return max(.03,1.4826*np.median(np.abs(x-np.median(x)))/np.sqrt(len(x)))
        err=float(np.hypot(uncertainty(baseline[:,1]),uncertainty(recent[:,1])))
        row[f'ztf_{band}_change_after_reference']=delta;row[f'ztf_{band}_change_uncertainty']=err
        if abs(delta)>=.3 and abs(delta)>=3*err:changes.append(delta)
    same_direction=bool(changes) and (len(changes)==1 or np.sign(changes[0])==np.sign(changes[1]))
    row['post_spectrum_optical_trigger']=same_direction and not info.get('association_warning',False)
    row['trigger_direction']='fading' if same_direction and np.median(changes)>0 else 'brightening' if same_direction else ''
    row['post_spectrum_ir_flag']=False
    neo=neowise.get(name)
    if neo is not None and len(neo)>=3:
        neo=neo.sort_values('mjd');tail=neo[neo.mjd>=neo.mjd.max()-3*365.25]
        row['neowise_last_mjd']=float(neo.mjd.max());row['neowise_last_date']=Time(neo.mjd.max(),format='mjd').strftime('%Y-%m-%d')
        row['neowise_amplitude_mag']=float(neo.w1.quantile(.95)-neo.w1.quantile(.05))
        row['neowise_recent_slope_mag_year']=float(np.polyfit((tail.mjd-tail.mjd.max())/365.25,tail.w1,1)[0]) if len(tail)>=4 else np.nan
        row['accepted_spectrum_after_wise_record']=bool(np.isfinite(mjd) and mjd>neo.mjd.max())
        if np.isfinite(mjd):
            before=neo[(neo.mjd<=mjd+180)&(neo.mjd>=mjd-365)]
            after=neo[(neo.mjd>=mjd+180)&(neo.mjd>=neo.mjd.max()-400)]
            if len(before)>=2 and len(after)>=2:
                delta=float(after.w1.median()-before.w1.median())
                row['neowise_change_after_reference_mag']=delta
                # IR echoes can lag an optical transition; this flags dated
                # evidence for inspection, not a new post-spectrum state.
                row['post_spectrum_ir_flag']=abs(delta)>=.3
    row['known_state_status']=known.get(name,'no match in checked catalogs')
    if target['pool_role']=='reserve':question='Comparison spectrum outside the selected manifold regions'
    elif row['post_spectrum_optical_trigger']:question='Test broad-line response to dated '+row['trigger_direction']+' after the accepted reference'
    elif row['post_spectrum_ir_flag']:question='Compare spectral epochs with dated IR evolution, allowing for a dust lag'
    elif row['known_state_status']=='catalog-confirmed CLAGN':question='Measure current state/persistence; recurrence is unverified'
    else:question='Test the manifold selection against the archival spectral baseline'
    row['science_question']=question
    row['state_classification']='Requires broad-line review; no state assigned automatically'
    return row


def dt(value):return datetime.fromisoformat(value).replace(tzinfo=timezone.utc)
def covers(value,start,end):
    # Saved Astropy ISO strings have millisecond precision, whereas datetime
    # retains microseconds. Permit serialization rounding, not extra sky time.
    epsilon=timedelta(milliseconds=1)
    return any(dt(r['start_utc'])<=start+epsilon and dt(r['end_utc'])>=end-epsilon for r in json.loads(value))


def main():
    targets=pd.read_csv(OUT/'compact_review_objects.csv');windows=pd.read_csv(OUT/'compact_review_object_nights.csv')
    etc=pd.read_csv(OUT/'three_night_2x600_sensitivity.csv');sens=etc.set_index('name').to_dict('index')
    known=pd.read_csv(OUT/'compact_spectral_audit.csv').set_index('name').known_state_status.to_dict()
    uw=pd.read_csv(OUT/'neighbour_unwise_manifest.csv').set_index('name').to_dict('index')
    frames=[]
    for tag in ['pool','poolall','zeltyn','three_night']+sorted(f.stem.removeprefix('neowise_visits_') for f in OUT.glob('neowise_visits_prepared*.csv')):
        path=(OUT if tag.startswith('prepared') else ROOT/'data')/f'neowise_visits_{tag}.csv'
        if path.exists():
            d=pd.read_csv(path);frames.append(d[d.name.isin(targets.name)])
    neo=pd.concat(frames,ignore_index=True).drop_duplicates(['name','mjd'],keep='last')
    neowise=dict(tuple(neo.groupby('name')))
    rows=[]
    for t in targets.to_dict('records'):
        row=science(t,sens,known,neowise);row['field_status'],row['field_notes']=field(t['name'],uw);rows.append(row)
    evidence=pd.DataFrame(rows);evidence=evidence.merge(etc.drop(columns=['reference_mjd'],errors='ignore'),on='name',how='left',validate='one_to_one')
    # Transparent review ordering: dated post-spectrum changes first, then
    # sensitivity within a tier. Spectral count is not a merit score.
    evidence['review_order_score']=100*evidence.post_spectrum_optical_trigger.astype(int)+25*evidence.post_spectrum_ir_flag.astype(int)+evidence.snr_sky18p5_X1p3.fillna(0).clip(upper=20)
    evidence.to_csv(DEST/'science_and_sensitivity.csv',index=False)
    allslots=[];bundles=[];shortlists=[]
    for night,(date,part) in OBS.NIGHTS.items():
        t0,t1,_,_=OBS.night_window(date,part);duration=(t1-t0).to_value(u.min)
        slots=[]
        for minute in np.arange(0,duration-30+1e-6,30):
            start=(t0+minute*u.min).to_datetime(timezone=timezone.utc);end=(t0+(minute+30)*u.min).to_datetime(timezone=timezone.utc)
            slots.append((start,end))
        prepared=set(targets.loc[targets.prepared_for_nights.fillna('').str.split(',').map(lambda ns:night in ns),'name']) if 'prepared_for_nights' in targets else set(targets.name)
        w=windows[windows.night.eq(night)&windows.name.isin(prepared)][['name','ranges_json','preferred_longest_minutes','min_airmass','minutes_airmass_le1p3']].merge(evidence,on='name')
        w=w.merge(targets[['name','ra','dec']],on='name')
        for start,end in slots:
            available=w[w.ranges_json.map(lambda s:covers(s,start,end))]
            m=available[available.pool_role.eq('manifold')];clear=m[m.field_status.eq('clear')]
            allslots.append(dict(night=night,start_pdt=start.astimezone(TZ).strftime('%Y-%m-%d %H:%M'),end_pdt=end.astimezone(TZ).strftime('%H:%M'),
                all_manifold=len(m),clear_manifold=len(clear),clear_nominal_snr10=int(clear.snr_sky18p5_X1p3.ge(10).sum()),clear_bright_sky_snr10=int(clear.snr_sky18_X1p5.ge(10).sum()),
                clear_reserves=int((available.pool_role.eq('reserve')&available.field_status.eq('clear')).sum())))
        main=w[w.pool_role.eq('manifold')&w.field_status.eq('clear')&w.status.eq('Scenario calculation')].copy().reset_index(drop=True)
        # Three distinct alternative objects per visit, with no object reused
        # among that night's bundles. Dummy assignments expose coverage gaps.
        costs=np.full((3*len(slots),len(main)+3*len(slots)),1e6)
        for i,(start,end) in enumerate(slots):
            valid=main.ranges_json.map(lambda value:covers(value,start,end)).to_numpy()
            for j in np.where(valid)[0]:costs[3*i:3*i+3,j]=-float(main.iloc[j].review_order_score)
        rr,cc=linear_sum_assignment(costs);selected=set()
        for r,c in zip(rr,cc):
            start,end=slots[r//3]
            name=main.iloc[c]['name'] if c<len(main) and costs[r,c]<1e6 else ''
            bundles.append(dict(night=night,start_pdt=start.astimezone(TZ).strftime('%Y-%m-%d %H:%M'),choice=r%3+1,name=name,status='alternative' if name else 'unfilled with distinct screened manifold targets'))
            if name:selected.add(name)
        # At least 45 unique September alternatives when the screened pool
        # permits it; longer nights require at least three per visit.
        desired=max(45,3*len(slots))
        for name in main.sort_values(['review_order_score','r_planning'],ascending=[False,True]).name:
            if len(selected)>=desired:break
            selected.add(name)
        short=main[main.name.isin(selected)].copy();short['night']=night;shortlists.append(short)
        plot=w.sort_values(['ranges_json','ra']);height=max(7,.19*len(plot)+2)
        fig,ax=plt.subplots(figsize=(13,height));colors={'manifold':'#497b9b','reserve':'#b47c42'}
        for i,r in enumerate(plot.itertuples()):
            for segment in json.loads(r.ranges_json):
                left,right=mdates.date2num(dt(segment['start_utc'])),mdates.date2num(dt(segment['end_utc']))
                ax.barh(i,right-left,left=left,height=.68,color=colors[r.pool_role])
        ax.set_yticks(range(len(plot)),[f'{r.name}{" [reserve]" if r.pool_role=="reserve" else ""}  r={r.r_planning:.1f}  S/N~{r.snr_sky18p5_X1p3:.0f}' for r in plot.itertuples()],fontsize=7)
        ax.invert_yaxis();ax.set_xlim(mdates.date2num(t0.to_datetime(timezone=timezone.utc)),mdates.date2num(t1.to_datetime(timezone=timezone.utc)))
        ax.xaxis.set_major_locator(mdates.HourLocator(tz=TZ));ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M',tz=TZ));ax.grid(axis='x',alpha=.25)
        ax.set_xlabel('Palomar local time (PDT)');ax.set_title(f'{date}: prepared alternatives (blue: manifold; amber: reserves)\nX ≤ 1.5, Moon ≥ 40°, r < 19; 30-minute minimum window')
        fig.text(.01,.012,'Availability, not assigned observations. S/N labels: archival continuum, sky V=18.5, X=1.3; not broad-line significance. Field/slit inspection remains required.',fontsize=8)
        fig.tight_layout(rect=(0,.035,1,1))
        for extension in ['png','pdf']:fig.savefig(DEST/f'{night}_visibility.{extension}',dpi=150)
        plt.close(fig)
    coverage=pd.DataFrame(allslots);coverage.to_csv(DEST/'slot_coverage.csv',index=False)
    pd.DataFrame(bundles).to_csv(DEST/'distinct_alternatives_by_slot.csv',index=False)
    pd.concat(shortlists,ignore_index=True).to_csv(DEST/'time_balanced_shortlist.csv',index=False)
    summary=coverage.groupby('night')[['clear_manifold','clear_nominal_snr10','clear_bright_sky_snr10','clear_reserves']].min()
    (DEST/'REVIEW.md').write_text('# Three-night review\n\nThe full parent is retained separately. These are provisional, time-balanced alternatives from the currently prepared manifold candidates, not a confirmed observing schedule. W1 projection coverage and acquisition manifests must be checked before calling the search complete.\n\n'+summary.to_markdown()+'\n\nTable values are the minimum available counts in any complete 30-minute block. The local bundle table assigns manifold-only alternatives without reusing an object across blocks; blank assignments expose shortages. The prepared_distinct_alternatives.csv file one directory above includes the separately labeled reserves. The visibility plots include both pools. S/N is continuum per binned pixel, not integrated broad-H-beta significance. Sensitivity is reported rather than imposed as an undocumented cut.\n\nReview ordering prioritizes a same-filter optical change of at least 0.3 mag and three robust standard errors after the latest accepted spectral reference; at least three nightly measurements near that reference and three at least 90 days later are required. Consistent g/r directions are required when both trigger, and association flags prevent promotion. This is a screening heuristic, not a CLAGN classification or an extrapolation to September/October 2026. Dated IR changes receive a smaller review weight, with dust lag explicitly unresolved; expected continuum S/N breaks ties. Spectral count never supplies a scientific reward or penalty. Manifold neighbor fractions are not transition probabilities.\n\nThe field-clear subset has no listed SDSS/Gaia/unWISE catalogue flags; this is not a contamination guarantee. Inspect the cutout and slit PA. The complete pool retains conditional fields and lower-S/N alternatives for review.\n')
    print(summary.to_string(),flush=True)


if __name__=='__main__':main()
