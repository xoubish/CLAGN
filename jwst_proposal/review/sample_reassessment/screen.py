"""Descriptive history screening, not a detection forecast or final selection.

Uses saved photometry only. Run with the project's scientific Python environment.
Outputs CSV/PNG review artifacts; does not edit sample membership or the proposal.
"""
from pathlib import Path
import os, json, hashlib
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/clagn-matplotlib')
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
INPUTS=ROOT/'jwst_proposal/inputs'
AUDIT=OUT.parent/'history_audit'
USED={}
def read(p, **kw):
    USED[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    return pd.read_csv(p, **kw)

def metrics(t,m):
    good=np.isfinite(t)&np.isfinite(m)
    t=np.asarray(t)[good];m=np.asarray(m)[good]
    ix=np.argsort(t);t=t[ix];m=m[ix]
    if len(t)<4:return {}
    early=np.median(m[t<=t[0]+400]);late=np.median(m[t>=t[-1]-400])
    delta=late-early;sign=np.sign(delta)
    sm=pd.Series(m).rolling(3,center=True,min_periods=2).median().to_numpy()
    span=np.ptp(sm)
    # Fraction of smoothed full excursion retained in endpoint change.
    retention=abs(delta)/span if span>0 else 0
    retreat=max(0,float(np.max(sign*sm)-sign*late))
    recent=t>=t[-1]-3*365.25
    slope=np.polyfit((t[recent]-t[recent][0])/365.25,m[recent],1)[0] if recent.sum()>=4 else np.nan
    return dict(n=len(t),first_year=float(Time(t[0],format='mjd').decimalyear),
                last_year=float(Time(t[-1],format='mjd').decimalyear),
                delta_mag=float(delta),smooth_span_mag=float(span),retention=float(retention),
                terminal_retreat_mag=retreat,late_slope_mag_yr=float(slope),
                p95_p05_mag=float(np.percentile(m,95)-np.percentile(m,5)))

def pass_wise(r, threshold=.30, retention=.65):
    return (r.get('W1_n',0)>=15 and r.get('W2_n',0)>=15
        and abs(r.get('W1_delta_mag',0))>=threshold
        and abs(r.get('W2_delta_mag',0))>=threshold
        and r['W1_delta_mag']*r['W2_delta_mag']>0
        and min(r['W1_retention'],r['W2_retention'])>=retention
        and max(r['W1_terminal_retreat_mag'],r['W2_terminal_retreat_mag'])<=.15
        and all(np.sign(r[b+'_delta_mag'])*r[b+'_late_slope_mag_yr']>=-.05 for b in ['W1','W2']))

def optical(name,ra,dec):
    paths=[ROOT/'data/ztf_cache'/tag/f'{name}.csv' for tag in
           ['review_dr24_20260920','calib','zeltyn','pool','v2']]
    result={};curves=[]
    for p in paths:
        if not p.exists() or p.stat().st_size<10:continue
        try:d=read(p,dtype={'oid':str,'expid':str},low_memory=False)
        except pd.errors.EmptyDataError:continue
        if not {'ra','dec','mjd','mag','magerr','catflags','oid','filtercode','expid'}.issubset(d):continue
        d=d[np.isfinite(d[['ra','dec','mjd','mag','magerr']]).all(axis=1)&(d.magerr>0)].copy()
        d=d[(d.catflags.fillna(32768).astype(int)&32768)==0]
        if d.empty:continue
        sep=SkyCoord(d.ra.to_numpy()*u.deg,d.dec.to_numpy()*u.deg).separation(SkyCoord(ra*u.deg,dec*u.deg)).arcsec
        d=d[sep<1.5].copy();d['bin']=np.floor(d.mjd/90)
        for fc,b in [('zg','g'),('zr','r')]:
            a=d[d.filtercode==fc]
            if a.empty:continue
            ranked=[]
            for oid,g in a.groupby('oid'):
                ranked.append(((g.groupby('bin').size()>=3).sum(),g.mjd.max()-g.mjd.min(),len(g),oid))
            oid=max(ranked)[-1]
            a=a[a.oid==oid].sort_values('magerr').drop_duplicates('expid')
            bins=a.groupby('bin').agg(mjd=('mjd','median'),mag=('mag','median'),n=('mag','size'))
            bins=bins[bins.n>=3].copy()
            bins['name']=name;bins['band']=b;bins['source']=str(p.relative_to(ROOT));bins['oid']=oid
            result.update({b+'_'+k:v for k,v in metrics(bins.mjd.to_numpy(),bins.mag.to_numpy()).items()})
            curves.append(bins.reset_index(drop=True))
        result['optical_source']=str(p.relative_to(ROOT))
        return result,curves
    return result,curves

def atlas(rows,curves,prefix):
    colors={'W1':'#782ba4','W2':'#bd7813','g':'#008562','r':'#c34054'}
    for start in range(0,len(rows),6):
        fig,axes=plt.subplots(3,2,figsize=(12,10),sharex=True)
        for ax,(_,r) in zip(axes.flat,rows.iloc[start:start+6].iterrows()):
            for b,d in curves[r['name']].groupby('band'):
                d=d.sort_values('mjd');ref=d[d.mjd>=d.mjd.max()-400].mag.median()
                ax.plot(Time(d.mjd.to_numpy(),format='mjd').decimalyear,10**(-.4*(d.mag-ref)),'.-',lw=1,ms=3,label=b,color=colors[b])
            ax.axhline(1,color='.6',ls=':',lw=.6);ax.set_yscale('log');ax.grid(alpha=.15)
            ax.set_title(f"{r['name']} | {r['direction']} | z={r['z']:.3f}\n{r.get('optical_status','')} / {r.get('ngps_status','')}",fontsize=9)
            ax.legend(ncol=4,fontsize=7);ax.set_xlim(2013.6,2026.8);ax.set_ylabel('Flux / late median')
        for ax in axes.flat[len(rows.iloc[start:start+6]):]:ax.set_visible(False)
        fig.suptitle('Replacement candidates: measured total-source histories',fontsize=14)
        fig.text(.5,.015,'Each band normalized separately; host included. Lines across gaps are not measurements. No warm-response prediction.',ha='center',fontsize=9)
        fig.tight_layout(rect=[0,.035,1,.95]);fig.savefig(OUT/f'{prefix}_{start//6+1}.png',dpi=130);plt.close(fig)

def main():
    sample=read(INPUTS/'jwst_sample_cycle6.csv')
    current=read(AUDIT/'binned_measurements.csv')
    rows=[]
    for _,r in sample.iterrows():
        row=dict(id=r.id,name=r.internal_id,target=r.target,z=r.z,selection_family=r.family)
        for b,d in current[current.id==r.id].groupby('band'):
            row.update({b+'_'+k:v for k,v in metrics(d.mjd.to_numpy(),d.mag.to_numpy()).items()})
        row['wise_screen_pass']=pass_wise(row)
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT/'current_metrics.csv',index=False)
    pool=read(INPUTS/'pool_forecast.csv').drop_duplicates('name')
    # This is an existing W3-matched preselected subset, not the entire parent catalogue.
    master=read(ROOT/'data/master_list_scored.csv',low_memory=False).drop_duplicates('name').set_index('name')
    visits=read(ROOT/'data/neowise_visits_poolall.csv')
    groups={k:v for k,v in visits.groupby('name')}
    ngps={p.stem for p in (ROOT/'sep23_data/reduction_20260924/products_p330e/spectra').glob('*.csv')}
    schedule=read(ROOT/'observing/oct26_27_2026/all_80_targets.csv')
    scheduled=set(schedule.NAME)
    rows=[];allcurves={}
    currentcoords=SkyCoord(sample.ra.to_numpy()*u.deg,sample.dec.to_numpy()*u.deg)
    for _,r in pool.iterrows():
        if not 0<r.z<.3 or r['name'] not in groups:continue
        d=groups[r['name']].sort_values('mjd')
        row={k:r[k] for k in ['name','ra','dec','z','r_mag','w3_mjy']}
        row['current_sample_match']=bool(SkyCoord(r.ra*u.deg,r.dec*u.deg).separation(currentcoords).arcsec.min()<3)
        for band,col in [('W1','w1'),('W2','w2')]:
            row.update({band+'_'+k:v for k,v in metrics(d.mjd.to_numpy(),d[col].to_numpy()).items()})
        if 'W1_delta_mag' not in row or 'W2_delta_mag' not in row:continue
        row['direction']='fade' if row['W1_delta_mag']>0 else 'rise'
        row['wise_screen_pass']=pass_wise(row)
        row['ngps_status']='local_reduced_spectrum' if r['name'] in ngps else 'proposed_Oct_queue' if r['name'] in scheduled else 'not_in_checked_followup'
        if r['name'] in master.index:
            m=master.loc[r['name']]
            for k in ['n_spec','blend_flag','blend_kind','r_mag','source_catalog']:row[k]=m.get(k,np.nan)
        rows.append(row)
        if row['wise_screen_pass'] and not row['current_sample_match']:
            allcurves[r['name']]=pd.concat([pd.DataFrame({'mjd':d.mjd,'mag':d[col],'band':band}) for band,col in [('W1','w1'),('W2','w2')]])
    screened=pd.DataFrame(rows)
    screened.to_csv(OUT/'parent_wise_screen.csv',index=False)
    eligible=screened[screened.wise_screen_pass&~screened.current_sample_match].copy()
    results=[]
    for _,row in eligible.iterrows():
        r=row.to_dict();stats,curves=optical(r['name'],r['ra'],r['dec']);r.update(stats)
        if curves:allcurves[r['name']]=pd.concat([allcurves[r['name']],*curves],ignore_index=True)
        direction=1 if r['direction']=='fade' else -1
        informative=[b for b in ['g','r'] if r.get(b+'_n',0)>=10 and r.get(b+'_last_year',0)>=2025]
        opposing=[b for b in informative if direction*r[b+'_delta_mag']<-.15 or direction*r[b+'_late_slope_mag_yr']<-.05]
        supporting=[b for b in informative if direction*r[b+'_delta_mag']>=.2 and r[b+'_terminal_retreat_mag']<=.2]
        r['optical_status']='conflict_or_reversal' if opposing else 'supports_direction' if supporting else 'weak_or_missing'
        r['optical_support_bands']=','.join(supporting)
        r['history_rank']=min(abs(r['W1_delta_mag']),abs(r['W2_delta_mag']))*min(r['W1_retention'],r['W2_retention'],1)
        results.append(r)
    result=pd.DataFrame(results).sort_values(['optical_status','history_rank'],ascending=[True,False])
    result.to_csv(OUT/'replacement_candidates.csv',index=False)
    # Review the best twelve supported candidates per direction; broader table is retained.
    top=result[result.optical_status=='supports_direction'].sort_values('history_rank',ascending=False).groupby('direction',group_keys=False).head(12).sort_values(['direction','history_rank'],ascending=[True,False])
    top.to_csv(OUT/'replacement_visual_review.csv',index=False)
    atlas(top,allcurves,'replacement_histories')
    pd.concat([d.assign(name=k) for k,d in allcurves.items()],ignore_index=True).to_csv(OUT/'replacement_binned_measurements.csv',index=False)
    sensitivity=[]
    for amp in [.25,.30,.40]:
        for retention in [.55,.65,.75]:
            ok=screened.apply(lambda r:pass_wise(r,amp,retention),axis=1)&~screened.current_sample_match
            sensitivity.append(dict(min_change_mag=amp,min_retention=retention,**screened[ok].direction.value_counts().to_dict()))
    pd.DataFrame(sensitivity).to_csv(OUT/'screen_threshold_sensitivity.csv',index=False)
    (OUT/'provenance.json').write_text(json.dumps(dict(inputs=USED,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),parent_rows=len(pool),screened_rows=len(screened),wise_candidates=len(eligible),optical_status_counts=result.groupby(['direction','optical_status']).size().to_dict().__str__(),notes='Descriptive screening only. Parent preselected for W1 variability and W3 availability. WISE visits end in 2024. No SPHEREx completeness or nuclear decomposition validation for replacements. Current sample untouched.'),indent=2)+'\n')
    print(result.groupby(['direction','optical_status']).size().to_string())
    print(top[['name','direction','z','W1_delta_mag','W2_delta_mag','optical_support_bands','ngps_status','n_spec']].to_string(index=False))

if __name__=='__main__':main()
