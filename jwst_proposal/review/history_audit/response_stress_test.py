"""Conditional history contrasts, NOT calibrated warm-dust predictions.

Applies unit-gain linear smoothing to observed total-source flux histories.
Response times are free sensitivity parameters, not measured lags. Prehistory
and future are held constant. This isolates how much the proposed group test
depends on response time, host dilution, and unobserved illumination.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
LAGS=[.5,1,2,3,5,10,20]
DATES=['2027-07-01','2028-01-01','2028-06-30']

def exponential(t,f,now,tau):
    # Exact convolution of a piecewise-linear input with an exponential kernel;
    # the response is in equilibrium with the first input before the record.
    response=f[0]
    for a,b,x,y in zip(t[:-1],t[1:],f[:-1],f[1:]):
        dt=b-a;one_minus=-np.expm1(-dt/tau)
        response=response*(1-one_minus)+x*one_minus+(y-x)/dt*(dt-tau*one_minus)
    return response*np.exp(-(now-t[-1])/tau)+f[-1]*(-np.expm1(-(now-t[-1])/tau))

def uniform(t,f,now,tau):
    # A thin spherical shell has delays uniformly distributed from 0 to 2R/c;
    # this is used only as a contrasting linear-response kernel.
    lo=now-2*tau
    knots=np.r_[lo,t[(t>lo)&(t<now)],now]
    return np.trapz(np.interp(knots,t,f),knots)/(2*tau)

def weights(t,now,tau,kernel):
    if kernel=='exponential':
        pre=np.exp(-(now-t[0])/tau);future=1-np.exp(-(now-t[-1])/tau)
    else:
        pre=np.clip((t[0]-(now-2*tau))/(2*tau),0,1)
        future=np.clip((now-t[-1])/(2*tau),0,1)
    return float(pre),float(1-pre-future),float(future)

def main():
    curves=pd.read_csv(OUT/'binned_measurements.csv')
    sample=pd.read_csv(OUT/'history_summary.csv').set_index('id')
    rows=[]
    for (ident,band),d in curves.groupby(['id','band']):
        if band not in ['g','r','W1']:continue
        d=d.sort_values('mjd');z=d.z.iloc[0]
        t=(d.mjd.to_numpy()-d.mjd.min())/(365.25*(1+z))
        f=10**(-.4*(d.mag.to_numpy()-d.mag.median()))
        for anchor in ['last_bin','last400d_median']:
            a=f.copy()
            if anchor=='last400d_median':
                a[-1]=np.median(f[d.mjd>=d.mjd.max()-400])
            for host_min_fraction in [0,.5,.8]:
                # Host is a constant fraction of the minimum of the ORIGINAL
                # observed curve, not a measurement or an inferred prior.
                driver=a-host_min_fraction*f.min()
                for date in DATES:
                    now=(Time(date).mjd-d.mjd.min())/(365.25*(1+z))
                    for tau in LAGS:
                        for kernel,fn in [('exponential',exponential),('uniform',uniform)]:
                            response=fn(t,driver,now,tau)
                            pre,observed,future=weights(t,now,tau,kernel)
                            rows.append(dict(id=ident,family=sample.loc[ident,'family'],band=band,
                                anchor=anchor,host_min_fraction=host_min_fraction,date=date,
                                tau_rest_yr=tau,kernel=kernel,
                                history_contrast_dex=np.log10(response/driver[-1]),
                                prehistory_weight=pre,record_window_weight=observed,
                                future_weight=future,n_bins=len(d),
                                last_observed_year=float(Time(d.mjd.max(),format='mjd').decimalyear)))
    full=pd.DataFrame(rows);full.to_csv(OUT/'conditional_response_grid.csv',index=False)
    keys=['band','anchor','host_min_fraction','date','tau_rest_yr','kernel']
    summary=[]
    for vals,g in full.groupby(keys):
        fade=g[g.family=='fade'].history_contrast_dex;rise=g[g.family=='rise'].history_contrast_dex
        summary.append(dict(zip(keys,vals),mean_fade_minus_rise_dex=float(fade.mean()-rise.mean()),
            median_record_window_weight=float(g.record_window_weight.median()),
            median_prehistory_weight=float(g.prehistory_weight.median()),
            median_future_weight=float(g.future_weight.median()),
            n_fade_positive=int((fade>0).sum()),n_rise_negative=int((rise<0).sum())))
    result=pd.DataFrame(summary);result.to_csv(OUT/'conditional_group_contrasts.csv',index=False)
    base=result[(result.date=='2028-01-01')&(result.anchor=='last400d_median')&(result.host_min_fraction==0)]
    print(base[['band','kernel','tau_rest_yr','mean_fade_minus_rise_dex','median_record_window_weight']].round(3).to_string(index=False))
    fig,axes=plt.subplots(1,2,figsize=(11,4.7))
    colors={'g':'#008562','r':'#c34054','W1':'#782ba4'}
    for band,color in colors.items():
        for kernel,style in [('exponential','-'),('uniform','--')]:
            g=base[(base.band==band)&(base.kernel==kernel)].sort_values('tau_rest_yr')
            axes[0].plot(g.tau_rest_yr,g.mean_fade_minus_rise_dex,style+'o',color=color,label=f'{band}: {kernel}')
            axes[1].plot(g.tau_rest_yr,g.median_record_window_weight,style+'o',color=color)
    axes[0].axhline(.24,color='black',ls=':',label='Proposal planning threshold (0.24 dex)')
    axes[0].axhline(0,color='.6',lw=.6);axes[0].set_ylabel('Mean fade − rise history contrast (dex)')
    axes[1].set_ylabel('Median kernel weight inside measured time span')
    for ax in axes:ax.set_xscale('log');ax.set_xlabel('Assumed mean response time (rest years)');ax.grid(alpha=.15)
    axes[0].legend(fontsize=7)
    fig.suptitle('Sensitivity experiment at 2028 Jan 1 — not a MIRI forecast')
    fig.text(.5,.015,'Unit-gain response; host included; constant extrapolation; last-400-day endpoint. Within-span gaps are interpolated.',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.05,1,.95]);fig.savefig(OUT/'response_stress_test.pdf');fig.savefig(OUT/'response_stress_test.png',dpi=160)
    (OUT/'response_assumptions.json').write_text(json.dumps(dict(
        interpretation='Conditional linear-filter history contrasts, not fitted dust spectra or probabilities.',
        kernels=['exponential mean tau','uniform delays 0 to 2 tau'],
        mean_lags_rest_yr=LAGS,dates=DATES,
        host_min_fractions=[0,.5,.8],
        anchors=['last bin','last 400-day median substituted at final observed bin'],
        extrapolation='constant before first measurement and after final measurement',
        excluded='Unknown UV conversion, nonlinear dust response, actual warm geometry, covering factors, optical host measurements, SPHEREx for remaining 23, future variability and calibration covariance.',
        warning='The 0.24 dex line is only a scale comparison: these optical/W1 filters do not directly predict the warm-continuum observable.'),indent=2)+'\n')
    # Independent analytic checks for constant inputs and numerical integration.
    for fn in [exponential,uniform]:
        assert abs(fn(np.array([0.,1.,2.]),np.array([3.,3.,3.]),5.,2.)-3)<1e-12
    t=np.array([0.,1.,3.]);f=np.array([1.,4.,2.]);now=5.;tau=2.
    delay=np.linspace(0,100,100001)
    numeric=np.trapz(np.interp(now-delay,t,f)*np.exp(-delay/tau)/tau,delay)
    assert abs(exponential(t,f,now,tau)-numeric)<1e-6

if __name__=='__main__':main()
