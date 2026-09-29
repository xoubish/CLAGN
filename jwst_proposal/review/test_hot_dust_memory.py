"""Conditional hot-dust echo test with last-year linear ZTF extrapolation.

Exploratory: no proposal edits; no claim of a detected lag or structural memory.
Prespecified lag grid, source-normalized temporal changes and held-out targets.
"""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/clagn-matplotlib')
import json,sys,hashlib,argparse
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from scipy.optimize import least_squares,nnls
from scipy.linalg import pinvh
from scipy.integrate import trapezoid
from scipy.stats import chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.cosmology import Planck18

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;OUT=HERE/'hot_dust_memory_test'
sys.path[:0]=[str(HERE),str(ROOT)]
from review_spherex_sample import read_target
from fit_spherex_hot_dust import blackbody,LINES
LAGS=np.array([0.,.1,.25,.5,1.,2.,3.])
USE_SAVED_ALERTS=False
RNG=np.random.default_rng(29092026)


def continuum_fit(df,z,template):
    """Constant host; separate disc amplitude, dust amplitude and T per visit."""
    d=df.copy();rest=d.wavelength_um/(1+z)
    valid_visits=[]
    for visit,g in d.groupby('visit'):
        w=g.wavelength_um/(1+z)
        if len(g)>=35 and ((w>=1)&(w<=1.4)).sum()>=3 and ((w>=3.4)&(w<=4.05)).sum()>=3:
            valid_visits.append(int(visit))
    d=d[d.visit.isin(valid_visits)].copy()
    if len(valid_visits)<2:return None
    lo=(d.wavelength_um-d.wavelength_half_width_um).to_numpy()/(1+z)
    hi=(d.wavelength_um+d.wavelength_half_width_um).to_numpy()/(1+z)
    keep=np.isfinite(lo)&np.isfinite(hi)&(lo>0)&(hi>=lo)
    for line in LINES:keep &= ~((lo<=line*1.02)&(hi>=line*.98))
    keep &= ~((lo<=3.35)&(hi>=3.20))
    d=d[keep].copy();lo=lo[keep];hi=hi[keep]
    if any((d.visit==v).sum()<25 for v in valid_visits):return None
    q,qw=np.polynomial.legendre.leggauss(12)
    waves=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*q/2
    tab=np.loadtxt(ROOT/f'inputs/host_templates/{template}_template_norm.sed')
    hw=tab[:,0]/1e4;hf=tab[:,1]*hw**2;hf/=np.interp(1.6,hw,hf)
    host=np.interp(waves,hw,hf)@(qw/2);disc=waves**(-1/3)@(qw/2)
    ep=np.array([valid_visits.index(int(v)) for v in d.visit]);n=len(valid_visits)
    scale=float(np.median(d.flux_mjy));y=d.flux_mjy.to_numpy()/scale;e=d.flux_err_mjy.to_numpy()/scale
    def model(p):
        dust=blackbody(waves,p[1+2*n:][ep,None]*1000)@(qw/2)
        return p[0]*host+p[1:1+n][ep]*disc+p[1+n:1+2*n][ep]*dust
    def residual(p):return (model(p)-y)/e
    fits=[]
    for t in [.85,1.3,1.8]:
        mat=np.zeros((len(y),1+2*n));mat[:,0]=host
        for j in range(n):
            take=ep==j;mat[take,1+j]=disc[take]
            mat[take,1+n+j]=blackbody(waves[take],t*1000)@(qw/2)
        amp,_=nnls(mat/e[:,None],y/e)
        p=np.r_[amp,np.full(n,t)]
        fit=least_squares(residual,p,bounds=(np.r_[np.zeros(1+2*n),np.full(n,.6)],
                                           np.r_[np.full(1+2*n,np.inf),np.full(n,2.4)]),
                          x_scale='jac',max_nfev=350,ftol=1e-8,xtol=1e-8,gtol=1e-8)
        fits.append(fit)
    fit=min(fits,key=lambda f:float(f.fun@f.fun));p=fit.x;dof=len(y)-len(p)
    red=float(fit.fun@fit.fun/dof)
    covariance=pinvh(fit.jac.T@fit.jac)*max(1,red)
    # Integral of a Planck Fnu shape over all frequency is proportional to
    # T^4 [exp(h nu_2/kT)-1]. Constants cancel in temporal ratios.
    temp=p[1+2*n:]*1000;amp=p[1+n:1+2*n]*scale
    bol=amp*(temp/1200)**4*np.expm1(14387.76877/(2*temp))/np.expm1(14387.76877/2400)
    jac=np.zeros((n,len(p)))
    for j in range(n):
        jac[j,1+n+j]=1/max(p[1+n+j],1e-12)
        x=14387.76877/(2*temp[j]);jac[j,1+2*n+j]=(4-x/(1-np.exp(-x)))/p[1+2*n+j]
    logcov=jac@covariance@jac.T
    # Prespecified sensitivity allowance, not an instrument calibration claim:
    # independent visit scale errors of 3%, coherent across wavelengths.
    logcov+=np.eye(n)*.03**2
    rows=[]
    for j,visit in enumerate(valid_visits):
        g=d[d.visit==visit]
        rows.append(dict(visit=visit,mjd=float(g.mjds.median()),first_mjd=float(g.mjds.min()),
                         last_mjd=float(g.mjds.max()),n_used=len(g),temperature_K=float(temp[j]),
                         dust_2um_mjy=float(amp[j]),bolometric_proxy=float(bol[j]),
                         log_bolometric_error=float(np.sqrt(max(0,logcov[j,j]))),
                         temperature_at_bound=bool(temp[j]<610 or temp[j]>2390)))
    result=dict(template=template,converged=bool(fit.success),reduced_chi2=red,dof=dof,
                supplied_error_adequacy_p=float(chi2.sf(red*dof,dof)),
                host_rest16_mjy=float(p[0]*scale),visits=rows,log_bolometric_covariance=logcov.tolist(),
                jacobian_rank=int(np.linalg.matrix_rank(fit.jac)),n_parameters=len(p))
    return result


def optical_driver(row,target):
    """Use the audited source and OID, convert to flux, bin into 14-day cells."""
    path=ROOT.parent/row.source
    raw=pd.read_csv(path,dtype={'oid':str,'expid':str},low_memory=False)
    raw=raw[raw.oid.eq(str(row.oid))&raw.filtercode.eq('z'+row.band)].copy()
    good=np.isfinite(raw[['mjd','mag','magerr','ra','dec']]).all(axis=1)&(raw.magerr>0)
    raw=raw[good&((raw.catflags.fillna(32768).astype(int)&32768)==0)]
    sep=SkyCoord(raw.ra.to_numpy()*u.deg,raw.dec.to_numpy()*u.deg).separation(SkyCoord(target.ra*u.deg,target.dec*u.deg)).arcsec
    raw=raw[sep<1.5].sort_values('magerr').drop_duplicates('expid').sort_values('mjd')
    if len(raw)<6:return None
    alert_meta={}
    apath=ROOT/f'inputs/lightcurves/{target.id}_alerce.csv'
    if USE_SAVED_ALERTS and apath.exists():
        al=pd.read_csv(apath).drop_duplicates('candid')
        al=al[al.fid.eq(1 if row.band=='g' else 2)&al.corrected.astype(bool)&~al.dubious.astype(bool)&np.isfinite(al.magpsf_corr)].copy()
        if len(al):
            distances=np.abs(raw.mjd.to_numpy()[None,:]-al.mjd.to_numpy()[:,None])
            near=distances.argmin(axis=1);matched=distances[np.arange(len(al)),near]<.02
            if matched.sum()>=5:
                off=float(np.median(al.magpsf_corr.to_numpy()[matched]-raw.mag.to_numpy()[near[matched]]))
                post=al[al.mjd>raw.mjd.max()].copy()
                error_name='sigmapsf_corr_ext' if 'sigmapsf_corr_ext' in post else 'sigmapsf_corr'
                good_error=np.isfinite(post[error_name])&(post[error_name]>0)
                post=post[good_error]
                extra=pd.DataFrame(dict(mjd=post.mjd,mag=post.magpsf_corr-off,magerr=post[error_name]))
                raw=pd.concat([raw,extra],ignore_index=True).sort_values('mjd')
                alert_meta=dict(file=str(apath.relative_to(ROOT.parent)),sha256=hashlib.sha256(apath.read_bytes()).hexdigest(),
                                matched=int(matched.sum()),offset_mag=off,added_points=len(extra))
    raw['flux']=3631000*10**(-.4*raw.mag);raw['err']=raw.flux*np.log(10)*.4*raw.magerr
    raw['cell']=np.floor(raw.mjd/14).astype(int)
    rows=[]
    for _,g in raw.groupby('cell'):
        w=1/g.err.to_numpy()**2;f=float(np.sum(w*g.flux)/w.sum())
        scatter=float(g.flux.std()/np.sqrt(len(g))) if len(g)>1 else 0
        rows.append(dict(mjd=float(g.mjd.median()),flux=f,error=max(float(1/np.sqrt(w.sum())),scatter),n=len(g)))
    d=pd.DataFrame(rows).sort_values('mjd');t=d.mjd.to_numpy();f=d.flux.to_numpy();e=d.error.to_numpy()
    # Exactly the last observed year, not the last year before the SPHEREx epoch.
    last=t>=t[-1]-365.25
    if last.sum()<3 or np.ptp(t[last])<90:return None
    x=(t[last]-t[-1])/365.25;design=np.column_stack([np.ones(last.sum()),x])
    # Equal-bin regression avoids letting a densely sampled week dominate the
    # annual trend; its residual scatter enters slope uncertainty.
    beta=np.linalg.lstsq(design,f[last],rcond=None)[0]
    residual=f[last]-design@beta;variance=max(float(residual@residual/(last.sum()-2)),float(np.mean(e[last]**2)))
    covariance=np.linalg.inv(design.T@design)*variance
    return dict(t=t,f=f,e=e,last=last,design=design,beta=beta,cov=covariance,
                source=str(path.relative_to(ROOT.parent)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),alerts=alert_meta)


def driver_samples(driver,ndraw=160):
    f=driver['f'];e=driver['e'];design=driver['design'];last=driver['last']
    draws=f[:,None]+e[:,None]*RNG.normal(size=(len(f),ndraw))
    fitdraw=np.linalg.lstsq(design,draws[last],rcond=None)[0]
    # Add annual-trend residual uncertainty beyond formal photometric noise.
    formal=np.linalg.inv(design.T@design)@design.T@np.diag(e[last]**2)@design@np.linalg.inv(design.T@design)
    extra=(driver['cov']-formal);eig,vec=np.linalg.eigh((extra+extra.T)/2)
    fitdraw+=vec@np.diag(np.sqrt(np.maximum(eig,0)))@RNG.normal(size=(2,ndraw))
    return np.column_stack([f,draws]),np.r_[driver['beta'][1],fitdraw[1]]


def response(driver,fluxes,slopes,when,z,tau,mode,host_fraction):
    """Linear interpolation; future line continuous at last measured bin.

    No positivity clipping. Invalid extrapolated illumination is rejected.
    Fixed host-fraction sensitivity is explicitly an assumption.
    """
    lag=tau*365.25*(1+z);lo=when-2*lag
    if lo<driver['t'][0]:return None
    knots=np.array([when]) if tau==0 else np.unique(np.r_[lo,driver['t'][(driver['t']>lo)&(driver['t']<when)],when])
    # Include the history/extrapolation boundary for exact trapezoid integration.
    if tau>0 and lo<driver['t'][-1]<when:knots=np.unique(np.r_[knots,driver['t'][-1]])
    values=np.column_stack([np.interp(knots,driver['t'],fluxes[:,j]) for j in range(fluxes.shape[1])])
    future=np.maximum(knots-driver['t'][-1],0)/365.25
    if mode=='linear':values+=future[:,None]*slopes[None,:]
    values-=host_fraction*np.min(driver['f'])
    acceptable=np.all(values>0,axis=0)
    result=values[0].copy() if tau==0 else trapezoid(values,knots,axis=0)/(2*lag)
    result[~acceptable]=np.nan
    return result


def score(table,lag_grid=LAGS,fit_responsivity=False):
    """Leave one target out for choosing a lag; score final/first changes."""
    ids=sorted(set.intersection(*(set(table[table.tau.eq(t)].id) for t in lag_grid)))
    t=table[table.id.isin(ids)].copy()
    if len(ids)<4:return dict(n=len(ids),reason='fewer than four common targets')
    obs=t[t.tau.eq(0)].set_index('id').loc[ids].observed_log_change.to_numpy()
    pred=np.array([t[t.tau.eq(lag)].set_index('id').loc[ids].predicted_log_change.to_numpy() for lag in lag_grid]).T
    optical_error=np.array([t[t.tau.eq(lag)].set_index('id').loc[ids].optical_prediction_error.to_numpy() for lag in lag_grid]).T
    se=t[t.tau.eq(0)].set_index('id').loc[ids].dust_change_error.to_numpy()
    # Equal target weighting; no preference for a lag merely because its
    # extrapolation uncertainty is larger. Errors shown separately.
    def gains(y,x):
        return np.clip(np.sum(y[:,None]*x,axis=0)/np.maximum(np.sum(x*x,axis=0),1e-15),0,3) if fit_responsivity else np.ones(x.shape[1])
    gain_all=gains(obs,pred)
    loss=(obs[:,None]-pred*gain_all)**2
    chosen=[];cv=[];base=[];cvpred=[];basepred=[];cvgain=[]
    for i in range(len(ids)):
        train_y=np.delete(obs,i);train_x=np.delete(pred,i,axis=0);gain=gains(train_y,train_x)
        j=int(np.argmin(np.mean((train_y[:,None]-train_x*gain)**2,axis=0)))
        chosen.append(float(lag_grid[j]));cvpred.append(float(pred[i,j]*gain[j]));basepred.append(float(pred[i,0]*gain[0]));cvgain.append(float(gain[j]))
        cv.append(float((obs[i]-cvpred[-1])**2));base.append(float((obs[i]-basepred[-1])**2))
    best=int(np.argmin(loss.mean(axis=0)))
    delta=np.array(base)-np.array(cv)
    boot=np.mean(delta[RNG.integers(0,len(ids),size=(4000,len(ids)))],axis=1)
    return dict(n=len(ids),ids=ids,fit_responsivity=fit_responsivity,rmse_by_tau={str(lag):float(np.sqrt(loss[:,j].mean())) for j,lag in enumerate(lag_grid)},
                best_grid_delay_parameter=float(lag_grid[best]),best_responsivity=float(gain_all[best]),rmse_instantaneous=float(np.sqrt(np.mean(base))),
                rmse_constant_dust=float(np.sqrt(np.mean(obs**2))),
                rmse_nested_heldout=float(np.sqrt(np.mean(cv))),heldout_selected_tau=chosen,
                fractional_rmse_reduction=float(1-np.sqrt(np.mean(cv)/np.mean(base))),
                mean_heldout_squared_error_improvement=float(delta.mean()),
                paired_bootstrap_improvement_16_84=np.quantile(boot,[.16,.84]).tolist(),
                median_dust_change_error=float(np.median(se)),
                per_target=[dict(id=tid,observed=float(obs[i]),instantaneous=basepred[i],
                                 heldout_prediction=cvpred[i],dust_error=float(se[i]),
                                 optical_error=float(optical_error[i,list(lag_grid).index(chosen[i])]*cvgain[i]),
                                 selected_tau=chosen[i],responsivity=cvgain[i]) for i,tid in enumerate(ids)])


def main():
    global OUT,USE_SAVED_ALERTS
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--use-saved-alerts',action='store_true')
    args=parser.parse_args();USE_SAVED_ALERTS=args.use_saved_alerts
    if USE_SAVED_ALERTS:OUT=OUT/'with_saved_alerts'
    OUT.mkdir(exist_ok=True)
    sample=pd.read_csv(ROOT/'inputs/jwst_sample_cycle6.csv')
    optical=pd.read_csv(HERE/'history_audit/optical_inventory.csv',dtype={'oid':str})
    optical=optical[optical.status.eq('available')]
    fitrows=[];excluded=[];models={};drivers={};driver_records=[];snapshots=[]
    for target in sample.itertuples():
        df,meta=read_target(target.id);snapshots.append(meta)
        if len(df):df.to_csv(OUT/f'{target.id}_spherex_snapshot.csv',index=False)
        for template in ['Ell5','S0']:
            cache=(OUT.parent if USE_SAVED_ALERTS else OUT)/f'{target.id}_{template}_fit.json'
            key=hashlib.sha256((meta['snapshot_sha256']+template+Path(__file__).read_text().split('def optical_driver')[0]).encode()).hexdigest()
            fit=None
            if cache.exists():
                previous=json.loads(cache.read_text())
                if previous.get('cache_key')==key:fit=previous
            if fit is None:
                fit=continuum_fit(df,target.z,template)
                if fit:
                    fit['cache_key']=key;cache.write_text(json.dumps(fit,indent=2)+'\n')
            if fit:
                models[(target.id,template)]=fit
                for v in fit['visits']:fitrows.append(dict(id=target.id,template=template,reduced_chi2=fit['reduced_chi2'],**v))
                print(target.id,template,'visits',len(fit['visits']),'reduced chi2',round(fit['reduced_chi2'],2),flush=True)
            else:excluded.append(dict(id=target.id,template=template,reason='fewer than two adequately covered spectral visits'))
        for r in optical[optical.id.eq(target.id)].itertuples():
            driver=optical_driver(r,target)
            if driver:
                drivers[(target.id,r.band)]=driver
                driver_records.append(dict(id=target.id,band=r.band,n_bins=len(driver['t']),last_year_bins=int(driver['last'].sum()),
                    last_mjd=float(driver['t'][-1]),last_flux_mjy=float(driver['f'][-1]),
                    slope_mjy_per_observed_year=float(driver['beta'][1]),slope_error=float(np.sqrt(driver['cov'][1,1])),
                    source=driver['source'],sha256=driver['sha256'],saved_alerts=json.dumps(driver['alerts'])))
                pd.DataFrame(dict(mjd=driver['t'],flux_mjy=driver['f'],error_mjy=driver['e'])).to_csv(OUT/f'{target.id}_{r.band}_optical_bins.csv',index=False)
            else:excluded.append(dict(id=target.id,band=r.band,reason='insufficient last-year optical bins/span'))
    pd.DataFrame(fitrows).to_csv(OUT/'spectral_fits.csv',index=False)
    pd.DataFrame(driver_records).to_csv(OUT/'optical_extrapolations.csv',index=False)
    pd.DataFrame(excluded).to_csv(OUT/'exclusions.csv',index=False)
    (OUT/'snapshot_manifest.json').write_text(json.dumps(snapshots,indent=2)+'\n')
    rows=[]
    for (tid,band),driver in drivers.items():
        z=float(sample.set_index('id').loc[tid,'z']);fluxes,slopes=driver_samples(driver)
        for template in ['Ell5','S0']:
            fit=models.get((tid,template))
            if not fit:continue
            a,b=fit['visits'][0],fit['visits'][-1]
            cov=np.array(fit['log_bolometric_covariance'])
            if min(a['bolometric_proxy'],b['bolometric_proxy'])<=0:continue
            obs=float(np.log(b['bolometric_proxy']/a['bolometric_proxy']))
            dy=float(np.sqrt(max(0,cov[0,0]+cov[-1,-1]-2*cov[0,-1])))
            # Report all fits; fixed reliability flags define a second cohort.
            reliable=bool(fit['converged'] and not a['temperature_at_bound'] and not b['temperature_at_bound']
                          and dy<.5 and fit['jacobian_rank']==fit['n_parameters'])
            for mode in ['linear','flat']:
                for host in [0.,.5]:
                    # Approximate Koshida V-band relation from a fixed median
                    # optical luminosity; host and optical SED assumptions are
                    # sensitivity cases. This is a lag prior, not a measurement.
                    fref=np.median(driver['f'])-host*np.min(driver['f'])
                    dl=Planck18.luminosity_distance(z).to_value(u.cm)
                    lnu=4*np.pi*dl**2*fref*1e-26/(1+z)*((.472 if band=='g' else .633)/((1+z)*.55))**(1/3)
                    absolute_v=-2.5*np.log10(lnu/(4*np.pi*(10*u.pc.to(u.cm))**2))-48.6
                    tau_ref=10**(-2.11-.2*absolute_v)/365.25
                    for family,grid in [('fixed_years',LAGS),('luminosity_scaled',np.array([0.,.5,1.,2.,4.]))]:
                      for parameter in grid:
                        tau=parameter if family=='fixed_years' else parameter*tau_ref
                        x=response(driver,fluxes,slopes,a['mjd'],z,tau,mode,host)
                        y=response(driver,fluxes,slopes,b['mjd'],z,tau,mode,host)
                        if x is None or y is None or not np.isfinite([x[0],y[0]]).all():continue
                        change=np.log(y/x);valid=np.isfinite(change[1:])
                        rows.append(dict(id=tid,band=band,template=template,continuation=mode,host_fraction=host,tau=float(parameter),
                            lag_family=family,tau_rest_years=float(tau),tau_reference_rest_years=float(tau_ref),
                            observed_log_change=obs,dust_change_error=dy,predicted_log_change=float(change[0]),
                            optical_prediction_error=float(np.nanstd(change[1:])) if valid.any() else np.nan,
                            valid_driver_draw_fraction=float(valid.mean()),screened_decomposition=reliable,
                            adequate_supplied_errors=bool(fit['supplied_error_adequacy_p']>.01),
                            reduced_chi2=fit['reduced_chi2'],first_mjd=a['mjd'],last_mjd=b['mjd'],
                            extrapolated_last_years=max(0,(b['mjd']-driver['t'][-1])/365.25),
                            history_index_last=float(np.log(y[0]/response(driver,fluxes[:,:1],slopes[:1],b['mjd'],z,0,mode,host)[0]))))
    tab=pd.DataFrame(rows);tab.to_csv(OUT/'conditional_predictions.csv',index=False)
    reports={}
    for keys,group in tab.groupby(['band','template','continuation','host_fraction','lag_family']):
        grids={'0to3yr':LAGS,'0to1yr':LAGS[LAGS<=1]} if keys[-1]=='fixed_years' else {'scaled':np.array([0.,.5,1.,2.,4.])}
        for grid_name,grid in grids.items():
          for cohort in ['all','screened','adequate']:
            g=group if cohort=='all' else group[group.screened_decomposition & group.valid_driver_draw_fraction.ge(.8)]
            if cohort=='adequate':g=g[g.adequate_supplied_errors]
            for free in [False,True]:
                reports['|'.join(map(str,keys))+f'|{grid_name}|{cohort}|'+('free_gain' if free else 'unit_gain')]=score(g,grid,free)
    summary=dict(created_utc=datetime.now(timezone.utc).isoformat(),use_saved_alerts=USE_SAVED_ALERTS,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),model='D_bolometric_proxy = constant_per_target * filtered_optical_flux',
        primary='g|Ell5|linear|0.0|fixed_years|0to1yr|screened|unit_gain',lag_grid_rest_years=LAGS.tolist(),
        spectrum_status='Exploratory single-blackbody decomposition; goodness of fit and covariance limitations retained.',
        extrapolation='Slope from final 365.25 observed days of 14-day flux bins; extend continuously from final bin. No negative-flux clipping.',
        uncertainty='Local spectral covariance scaled by max(1, reduced chi2), plus assumed 3% independent visit scale; optical Monte Carlo includes slope uncertainty. Conditional RMSE comparison does not marginalize all these uncertainties.',
        scoring='Equal-target log-change RMSE; first/last adequately covered visits; lag chosen on other targets for each held-out prediction. Same cohort across every lag.',
        limitations=['Optical driver includes host; 50%-of-minimum host subtraction is an assumed sensitivity case.',
          'No measured bolometric corrections or extinction corrections; blackbody extrapolation is model dependent.',
          'Poor spectral fits are not repaired by covariance inflation; these results are conditional exploratory tests.',
          'SPHEREx visits represented by median observation dates; intra-visit evolution is not modeled.',
          'Source normalization cancels in temporal changes; no naive correlation of shared-denominator ratios.',
          'No structural-radius evolution is inferred by this fixed-response experiment.',
          'Bootstrap intervals describe target resampling, not lag significance or a full posterior.'],reports=reports)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    # Predictive comparisons; the plotted delayed prediction is genuinely held out.
    fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
    for ax,band in zip(axes,['g','r']):
        rep=reports.get(f'{band}|Ell5|linear|0.0|fixed_years|0to1yr|screened|unit_gain',{})
        for r in rep.get('per_target',[]):
            ax.plot([r['instantaneous'],r['heldout_prediction']],[r['observed']]*2,color='#bbb',lw=.8)
            ax.scatter(r['instantaneous'],r['observed'],marker='x',color='#777',s=28)
            ax.scatter(r['heldout_prediction'],r['observed'],color='#c46a35',s=25)
            ax.errorbar(r['heldout_prediction'],r['observed'],xerr=r['optical_error'],yerr=r['dust_error'],fmt='none',ecolor='#c46a35',alpha=.3,lw=.7)
            ax.annotate(r['id'],(r['heldout_prediction'],r['observed']),xytext=(3,3),textcoords='offset points',fontsize=7)
        lim=max(.2,max([abs(r[k]) for r in rep.get('per_target',[]) for k in ['observed','instantaneous','heldout_prediction']],default=.2))*1.15
        ax.plot([-lim,lim],[-lim,lim],':',color='#aaa');ax.set(xlim=(-lim,lim),ylim=(-lim,lim),
            xlabel='Predicted ln(last / first dust flux)',ylabel='Fitted ln(last / first dust flux)',
            title=f'{band}-band estimated driver · n={rep.get("n",0)}')
        if 'rmse_instantaneous' in rep:
            ax.text(.03,.97,f'RMSE: current {rep["rmse_instantaneous"]:.3f}; delayed {rep["rmse_nested_heldout"]:.3f}',transform=ax.transAxes,va='top',fontsize=9)
    fig.suptitle('Conditional test: last-year linear optical extrapolation\nGrey ×: current state · orange: delay selected using other targets',fontsize=12)
    fig.savefig(OUT/'predictive_test.png',dpi=160);fig.savefig(OUT/'predictive_test.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5),constrained_layout=True)
    for band,color in [('g','#287d6d'),('r','#a34955')]:
        for mode,ls in [('linear','-'),('flat','--')]:
            rep=reports.get(f'{band}|Ell5|{mode}|0.0|fixed_years|0to1yr|screened|unit_gain',{})
            if 'rmse_by_tau' in rep:ax.plot(LAGS[LAGS<=1],list(rep['rmse_by_tau'].values()),ls+'o',color=color,label=f'{band}: {mode} continuation (n={rep["n"]})')
    ax.set(xlabel='Assumed mean response time (rest years)',ylabel='RMSE of ln dust-flux changes',
           title='Sensitivity to response time and unobserved optical continuation')
    ax.legend(fontsize=8);fig.savefig(OUT/'lag_sensitivity.png',dpi=160);fig.savefig(OUT/'lag_sensitivity.pdf');plt.close(fig)
    for key in [summary['primary'],'r|Ell5|linear|0.0|fixed_years|0to1yr|screened|unit_gain','g|Ell5|linear|0.0|luminosity_scaled|scaled|screened|free_gain']:
        print(key,json.dumps({k:v for k,v in reports.get(key,{}).items() if k!='per_target'}),flush=True)


if __name__=='__main__':main()
