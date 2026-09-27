"""Summarize adequacy and conditional forecasts; never label failed fits allowed."""
from pathlib import Path
import sys
import json
import hashlib
import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import minimize, nnls
from scipy.stats import chi2
from numpy.polynomial.legendre import leggauss
from astropy.cosmology import FlatLambdaCDM
from astropy.time import Time
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT))
from test_p9694_miri_predictivity import Audit, BOUNDS, Z, GRAINS


class Constrained:
    def __init__(self,a):
        self.a=a
        distance=FlatLambdaCDM(H0=70,Om0=.3).luminosity_distance(Z).to_value('m')
        self.bol=np.array([s['cooling']/s['norm'] for s in a.shells])*1e-29*4*np.pi*distance**2/(1+Z)
        g=next(g for g in GRAINS if g['name'].startswith('Gra'))
        heat=np.exp(np.interp(np.log(1500),g['log_t'],g['log_cool']))/g['q_heat']
        self.L0=16*np.pi**2*(9.46073047e15)**2*heat

    def matrices(self,p,eta):
        m,off=self.a.design(p,eta)
        aw=solve_triangular(self.a.chol5,m,lower=True,check_finite=False)
        yw=solve_triangular(self.a.chol5,self.a.y-off,lower=True,check_finite=False)
        w,wo=self.a.warm_bounds(p,eta)
        # Historical W3/W4 are upper constraints on nuclear flux, allowing 20%
        # calibration/colour-conversion slack, not 20% measured statistical errors.
        limits=np.array([23.52,74.81])*1.2-wo
        hs=self.a.driver(np.linspace(1900,2026.74,200),p)
        c=self.bol/(self.L0*10**(2*p[3]))
        inner=np.array([s['inner'] for s in self.a.shells])
        energy=np.array([c*np.where(inner,h**eta,1) for h in [hs.min(),1.,hs.max()]])
        inequality=np.vstack([w/limits[:,None],energy])
        # Remove shells exceeding illustrative sublimation limits in the fitted history.
        upper=[]
        for s in self.a.shells:
            g=s['grain']
            maximum=np.exp(np.interp(np.log(s['cooling']*hs.max()),g['log_cool'],g['log_t']))
            limit=1500. if g['name'].startswith('Sil') else 2000.
            upper.append(0. if maximum>limit else np.inf)
        return m,off,aw,yw,inequality,upper

    def fit(self,p,eta):
        m,off,aw,yw,ineq,upper=self.matrices(p,eta)
        if any(np.all(row<0) for row in ineq[:2]):
            return 1e12,np.zeros(12),False
        start,_=nnls(aw,yw,maxiter=300)
        start[np.array(upper)==0]=0
        start/=max(1.,float((ineq@start).max())*1.001)
        def fun(x):
            res=aw@x-yw
            return float(res@res)/len(yw)
        def jac(x):return 2*aw.T@(aw@x-yw)/len(yw)
        result=minimize(fun,start,jac=jac,bounds=list(zip(np.zeros(len(upper)),upper)),
                        constraints=[dict(type='ineq',fun=lambda x:1-ineq@x,jac=lambda x:-ineq)],
                        method='SLSQP',options={'ftol':1e-8,'maxiter':150})
        valid=bool(result.success and np.max(ineq@result.x)<=1.0001)
        return fun(result.x)*len(yw),result.x,valid

    def profile(self,p,eta,bestamps,objective,limit):
        m,off,aw,yw,ineq,upper=self.matrices(p,eta)
        constraints=[dict(type='ineq',fun=lambda x:1-ineq@x,jac=lambda x:-ineq),
                     dict(type='ineq',fun=lambda x:(limit-np.sum((aw@x-yw)**2))/len(yw),
                          jac=lambda x:-2*aw.T@(aw@x-yw)/len(yw))]
        found=[]
        objective=objective/max(np.max(np.abs(objective)),1e-30)
        for sign in [-1,1]:
            opt=minimize(lambda x:sign*objective@x,bestamps,jac=lambda x:sign*objective,
                         method='SLSQP',bounds=list(zip(np.zeros(len(upper)),upper)),
                         constraints=constraints,options={'ftol':1e-8,'maxiter':300})
            if opt.success and np.sum((aw@opt.x-yw)**2)<=limit+1e-3 and np.max(ineq@opt.x)<=1.0001:
                found.append(opt.x)
        return found

    def forecast_valid(self,p,eta,amplitudes,year,future):
        hs=self.a.driver(np.linspace(1900,year,300),p,future)
        c=self.bol/(self.L0*10**(2*p[3]))
        inner=np.array([s['inner'] for s in self.a.shells])
        for h in [hs.min(),1.,hs.max()]:
            if np.sum(c*amplitudes*np.where(inner,h**eta,1))>1.0001:return False
        for s,amp in zip(self.a.shells,amplitudes):
            if amp<1e-7:continue
            g=s['grain']
            maximum=np.exp(np.interp(np.log(s['cooling']*hs.max()),g['log_cool'],g['log_t']))
            if maximum>(1500. if g['name'].startswith('Sil') else 2000.):return False
        return True


def main():
    a=Audit()
    # Final adequacy/forecast analysis resolves the echo integral more finely
    # than the initial 8-node coarse parameter search.
    nodes,weights=leggauss(32);a.delay_nodes=(nodes+1)/2;a.delay_weights=weights/2
    constrained=Constrained(a)
    archive=json.loads((HERE/'profile_fits.json').read_text())
    n=len(a.y); dof=n-len(a.shells)-len(BOUNDS)
    threshold=float(chi2.ppf(.99,dof))
    results={'n_data':n,'nominal_dof':dof,'adequacy_chi2_threshold_99pct':threshold,
             'scope':'Exploratory optically thin shell families only; no clumpy/disc-wind posterior.',
             'primary':{},'sensitivity_5pct_extra_scatter':{}}
    refined=[]
    for eta in [0.,-.5,-1.]:
        candidates=sorted([r for r in archive['records'] if r['evolution']==eta],key=lambda r:r['chi2'])
        best=candidates[0]
        stat,amps,pred=a.fit(best['parameters'],eta)
        results['primary'][str(eta)]={**best,'chi2':stat,'reduced_chi2':stat/dof,'nominal_p_value':float(chi2.sf(stat,dof)),
            'passes_adequacy':bool(stat<=threshold),'residual_rms_fraction_by_dataset':
            {k:float(np.sqrt(np.mean(((pred[a.kind==k]-a.y[a.kind==k])/a.y[a.kind==k])**2))) for k in np.unique(a.kind)}}
        start=np.array(best['parameters']); start[3]=max(start[3],.05)
        opt=minimize(lambda p:constrained.fit(p,eta)[0],start,method='Powell',bounds=BOUNDS,
                     options={'maxfev':550,'ftol':.001,'xtol':.005})
        stat,amps,valid=constrained.fit(opt.x,eta)
        record=dict(evolution=eta,parameters=opt.x.tolist(),amplitudes=amps.tolist(),chi2=stat,
                    reduced_chi2=stat/dof,linear_constraints_satisfied=valid,
                    nonlinear_optimizer_success=bool(opt.success),
                    passes_conditional_adequacy=bool(valid and stat<=threshold))
        results['sensitivity_5pct_extra_scatter'][str(eta)]=record
        refined.append(record)
        print('CONSTRAINED +5% SCATTER',json.dumps(record),flush=True)
        (HERE/'adequacy_and_forecasts.json').write_text(json.dumps(results,indent=2)+'\n')
    wave=np.geomspace(4.9/(1+Z),27.9/(1+Z),180)
    observed=wave*(1+Z)
    # Continuous monochromatic forecast in the rest-frame MIRI wavelength range.
    waves=observed[:,None];weights=np.ones_like(waves)
    all_curves=[]; curve_labels=[]; examples=[]
    for fit in refined:
        if not fit['passes_conditional_adequacy']:
            continue
        eta=fit['evolution'];p=np.array(fit['parameters']);amp=np.array(fit['amplitudes'])
        m,off=a.design(p,eta,waves,weights,np.full(len(wave),2028.0))
        solutions=[amp]
        # Explore non-unique warm-shell weights as well as the nominal optimum.
        for rest in [6.,9.7,12.,18.]:
            j=np.argmin(abs(wave-rest))
            solutions.extend(constrained.profile(p,eta,amp,m[j],threshold))
        for future in [.5,1.,2.]:
            for t in [float(Time('2027-07-01').decimalyear),2028.,float(Time('2028-06-30').decimalyear)]:
                matrix,offset=a.design(p,eta,waves,weights,np.full(len(wave),t),future)
                for s,amplitude in enumerate(solutions):
                    if not constrained.forecast_valid(p,eta,amplitude,t,future):continue
                    all_curves.append(matrix@amplitude+offset)
                    curve_labels.append(dict(evolution=eta,future_driver_factor=future,year=t,profile_solution=s))
        examples.append((fit,solutions))
    curves=np.array(all_curves)
    np.savez_compressed(HERE/'conditional_curves.npz',rest_um=wave,flux_mjy=curves)
    results['conditional_curve_labels']=curve_labels
    results['forecast_date_range']=['2027-07-01','2028-06-30']
    results['future_illumination_factors']=[.5,1.,2.]
    results['extra_scatter_warning']=('5% added independently to every point is an explicit sensitivity assumption, '
        'not a measured calibration covariance or validation of the physical model. Envelopes are conditional '
        'scenario ranges, not posterior credible intervals.')
    results['historical_W3_W4_constraint']=('Nuclear band flux in 2010 <= 1.2 times total-source archival values '
        '(23.52,74.81 mJy); 20% is conservative calibration/colour slack, not a measured error.')
    results['energy_constraint']=('Dust luminosity does not exceed monochromatic heating luminosity at '
        'sampled fitted-history levels; necessary but not sufficient for radiative-transfer consistency.')
    results['curves_at_rest_um']={}
    if len(curves):
        for w in [6.,9.7,12.,18.]:
            j=np.argmin(abs(wave-w))
            results['curves_at_rest_um'][str(w)]={}
            for eta in [0.,-.5,-1.]:
                choose=np.array([l['evolution']==eta for l in curve_labels])
                if choose.any():
                    results['curves_at_rest_um'][str(w)][str(eta)]=[float(curves[choose,j].min()),float(curves[choose,j].max())]
    # Constructive whole-spectrum test: refit fixed shell masses to an evolving
    # prediction while retaining the existing-data and physical constraints.
    results['whole_spectrum_counterexamples']=[]
    fixed=results['sensitivity_5pct_extra_scatter']['0.0']
    if fixed['passes_conditional_adequacy']:
        fp=fixed['parameters'];fa=np.array(fixed['amplitudes'])
        _,_,aw,yw,ineq,upper=constrained.matrices(fp,0.)
        constraints=[dict(type='ineq',fun=lambda x:1-ineq@x,jac=lambda x:-ineq),
                     dict(type='ineq',fun=lambda x:(threshold-np.sum((aw@x-yw)**2))/len(yw),
                          jac=lambda x:-2*aw.T@(aw@x-yw)/len(yw))]
        for eta in [-.5,-1.]:
            example=results['sensitivity_5pct_extra_scatter'][str(eta)]
            if not example['passes_conditional_adequacy']:continue
            ep=example['parameters'];ea=np.array(example['amplitudes'])
            em,eo=a.design(ep,eta,waves,weights,np.full(len(wave),2028.),1.)
            target=em@ea+eo
            fm,fo=a.design(fp,0.,waves,weights,np.full(len(wave),2028.),1.)
            ar=fm/target[:,None];yr=(target-fo)/target
            opt=minimize(lambda x:np.mean((ar@x-yr)**2),fa,
                         jac=lambda x:2*ar.T@(ar@x-yr)/len(yr),
                         bounds=list(zip(np.zeros(len(upper)),upper)),constraints=constraints,
                         method='SLSQP',options={'ftol':1e-12,'maxiter':500})
            mismatch=(fm@opt.x+fo)/target-1
            datachi=float(np.sum((aw@opt.x-yw)**2))
            valid=bool(opt.success and datachi<=threshold+1e-3 and np.max(ineq@opt.x)<=1.0001
                       and constrained.forecast_valid(fp,0.,opt.x,2028.,1.)
                       and constrained.forecast_valid(ep,eta,ea,2028.,1.))
            item=dict(evolution=eta,year=2028.,future_factor_both=1.,valid=valid,
                      fixed_amplitudes=opt.x.tolist(),fixed_existing_data_chi2=datachi,
                      rms_fractional_mismatch=float(np.sqrt(np.mean(mismatch**2))),
                      maximum_fractional_mismatch=float(np.max(abs(mismatch))))
            results['whole_spectrum_counterexamples'].append(item)
            if eta==-.5 and valid:
                np.savetxt(HERE/'counterexample_spectra.csv',np.column_stack([wave,target,fm@opt.x+fo,mismatch]),
                           delimiter=',',header='rest_um,evolving_observed_mjy,fixed_observed_mjy,fixed_over_evolving_minus1',comments='')
                fig,axes=plt.subplots(2,1,figsize=(7,5),sharex=True,layout='constrained',
                                      gridspec_kw={'height_ratios':[3,1]})
                axes[0].plot(wave,target,color='#be513b',label='Evolving inner dust mass (η = −0.5)')
                axes[0].plot(wave,fm@opt.x+fo,'--',color='#24658a',label='Refitted fixed dust distribution')
                axes[0].set(yscale='log',ylabel='Predicted flux density (mJy)',
                            title='Different dust responses, nearly identical MIRI spectra')
                axes[0].legend(fontsize=9)
                axes[1].plot(wave,100*mismatch,color='#24658a');axes[1].axhline(0,color='.5',lw=.6)
                axes[1].set(xlabel='Rest wavelength (µm)',ylabel='Fixed / evolving\n− 1 (%)')
                for axis in axes:
                    for w in [9.7,18.]:axis.axvline(w,color='.7',ls=':',lw=.6)
                fig.suptitle('Same forecast date and future factor; an additional 5% historical scatter is assumed',fontsize=9)
                for ext in ['pdf','png']:fig.savefig(HERE/f'counterexample.{ext}',dpi=180)
                plt.close(fig)
    # Check the final echo integral against a doubled angular grid.
    numerical={}
    for record in refined:
        p=record['parameters'];eta=record['evolution'];amps=record['amplitudes']
        m,o=a.design(p,eta);p32=m@amps+o
        nodes,weights64=leggauss(64);a.delay_nodes=(nodes+1)/2;a.delay_weights=weights64/2
        m,o=a.design(p,eta);p64=m@amps+o
        numerical[str(eta)]=dict(max_fractional_32_to_64=float(np.max(abs(p32/p64-1))),
                                 rms_fractional_32_to_64=float(np.sqrt(np.mean((p32/p64-1)**2))))
        nodes,weights32=leggauss(32);a.delay_nodes=(nodes+1)/2;a.delay_weights=weights32/2
    results['numerical_checks']=numerical
    results['angular_quadrature_nodes']=32
    results['analysis_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/'adequacy_and_forecasts.json').write_text(json.dumps(results,indent=2)+'\n')
    # Evidence of adequacy failure, not an apparent evolutionary detection.
    fig,ax=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    for kind,marker in [('W1','o'),('W2','s')]:
        use=a.kind==kind
        ax[0].errorbar(a.times[use],a.y[use],yerr=a.err[use],fmt=marker,ms=3,label=f'{kind} measured')
    styles={0.:('#24658a','-'),-.5:('#be513b','--'),-1.:('#a6761d',':')}
    for eta in [0.,-.5,-1.]:
        b=results['primary'][str(eta)]; _,_,pred=a.fit(b['parameters'],eta)
        color,ls=styles[eta]
        for kind in ['W1','W2']:
            use=a.kind==kind;order=np.argsort(a.times[use])
            ax[0].plot(a.times[use][order],pred[use][order],ls,color=color,lw=1,
                       label=('Fixed' if eta==0 else f'Mass response η={eta}') if kind=='W1' else None)
    ax[0].set(xlabel='Year',ylabel='WISE flux density (mJy)',title='Fits using supplied errors + WISE calibration covariance')
    ax[0].legend(fontsize=7)
    for eta in [0.,-.5,-1.]:
        choose=np.array([l['evolution']==eta for l in curve_labels])
        if choose.any():
            color,ls=styles[eta]
            ax[1].fill_between(wave,curves[choose].min(axis=0),curves[choose].max(axis=0),
                               color=color,alpha=.15,label='Fixed' if eta==0 else f'Mass response η={eta}')
    ax[1].set(xlabel='Rest wavelength (µm)',ylabel='Predicted total flux density (mJy)',yscale='log',
              title='Conditional forecasts: additional 5% scatter assumed')
    if len(curves):ax[1].legend(fontsize=8)
    else:ax[1].text(.5,.5,'No family passes the adequacy test',ha='center',transform=ax[1].transAxes)
    for w in [9.7,18.]:ax[1].axvline(w,color='.6',ls=':',lw=.7)
    for ext in ['pdf','png']:fig.savefig(HERE/f'forecast_audit.{ext}',dpi=180)
    plt.close(fig)
    print('Wrote forecast audit',flush=True)


if __name__=='__main__':main()
