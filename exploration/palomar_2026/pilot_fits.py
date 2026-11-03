"""Exploratory, independently fitted H-beta decompositions for the Palomar pilot.

This is a model-sensitivity experiment, not a calibrated CLAGN classifier or a
joint multi-epoch inference. Fits use supplied diagonal pixel errors. Their
chi-square and the model envelope do NOT include flux calibration, aperture,
resampling covariance, or a posterior over stellar populations/Fe II shapes.
"""
import hashlib,json,sys,time
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares,nnls
from scipy.special import ndtr
from paper_protocol import OUT
from pilot_data import load
from analysis import integration_weights

C=299792.458
VARIANTS={
    'baseline':dict(young='19_ssp_1.4Gyr_z02.spec.gz',fe_sigma=1000.),
    'young_host':dict(young='28_ssp_290Myr_z02.spec.gz',fe_sigma=1000.),
    'wide_fe':dict(young='19_ssp_1.4Gyr_z02.spec.gz',fe_sigma=2000.),
}
NAMES=['agn','host_5gyr','host_young','feii','hb_narrow','oiii_core','oiii_wing',
       'hb_broad1','hb_broad2','hg_narrow','hg_broad','heii_broad','oiii4364']
PARAMS=['alpha_lambda','narrow_sigma_1000kms','narrow_v_1000kms',
        'broad1_sigma_1000kms','broad2_sigma_1000kms','broad1_v_1000kms','broad2_v_1000kms']
LOW=np.array([-3.5,.06,-.5,.5,1.,-2.,-3.])
HIGH=np.array([2.5,.6,.5,5.,8.,2.,3.])


class Basis:
    def __init__(self,wave,lsf,variant='baseline'):
        self.w=np.asarray(wave);self.lsf=np.asarray(lsf)
        self.edges=np.r_[self.w[0]-(self.w[1]-self.w[0])/2,(self.w[1:]+self.w[:-1])/2,self.w[-1]+(self.w[-1]-self.w[-2])/2]
        cfg=VARIANTS[variant];folder=OUT/'templates'
        # Template smoothing on a fine logarithmic grid. Fixed added host sigma
        # is a pilot assumption; template resolution is not deconvolved.
        logw=np.arange(np.log(3900),np.log(5800),30/C);grid=np.exp(logw)
        self.host=[]
        for name in ['13_ssp_5Gyr_z02.spec.gz',cfg['young']]:
            t=np.loadtxt(folder/name)
            f=np.interp(grid,t[:,0],t[:,1]);f=gaussian_filter1d(f,150/30)
            self.host.append(np.interp(self.w,grid,f)/np.interp(5100,grid,f))
        fe=np.loadtxt(folder/'fe_optical.txt')
        ff=np.interp(grid,10**fe[:,0],fe[:,1],left=0,right=0)
        ff=gaussian_filter1d(ff,cfg['fe_sigma']/30)
        ff/=np.mean(ff[(grid>4434)&(grid<4684)])
        self.fe=np.interp(self.w,grid,ff)

    def gaussian(self,center,sigma,v):
        mu=center*(1+v/C)
        inst=float(np.interp(mu,self.w,self.lsf))
        sig=np.hypot(center*sigma/C,inst)
        # Integrate across native pixel bounds; mask selection occurs later so
        # a masked gap cannot enlarge a neighboring pixel or bridge missing data.
        return np.diff(ndtr((self.edges-mu)/sig))/np.diff(self.edges)

    def matrix(self,p):
        alpha,ns,nv,s1,s2,v1,v2=p
        ns*=1000;nv*=1000;s1*=1000;s2*=1000;v1*=1000;v2*=1000
        g=self.gaussian
        core=10*(g(5008.24,ns,nv)+g(4960.30,ns,nv)/2.98)
        wing=10*(g(5008.24,800,nv-200)+g(4960.30,800,nv-200)/2.98)
        return np.column_stack([(self.w/5100)**alpha,*self.host,self.fe,
            10*g(4862.68,ns,nv),core,wing,100*g(4862.68,s1,v1),100*g(4862.68,s2,v2),
            10*g(4341.68,ns,nv),100*g(4341.68,s1,v1),100*g(4687.02,s1,v1),10*g(4364.44,ns,nv)])


def prepare(s,flux_key='flux'):
    w=s['wave']/(1+s['z']);m=(w>=4200)&(w<=5500)
    x=w[m];y=np.asarray(s[flux_key])[m];e=s['error'][m]
    lsf=s['lsf_sigma'][m]/(1+s['z'])
    lsf=np.where(np.isfinite(lsf)&(lsf>0),lsf,x/1650/2.35482)
    good=np.isfinite(y)&np.isfinite(e)&(e>0)
    observed=x*(1+s['z'])
    for a,b in [(5569,5585),(6293,6308),(6860,6960),(7580,7705),(9300,9650)]:
        good&=~((observed>=a)&(observed<=b))
    return x,y,e,lsf,good


def fit_arrays(x,y,e,lsf,good,variant='baseline',starts=2):
    basis=Basis(x,lsf,variant)
    yg=y[good];eg=e[good]
    def solve(p,return_parts=False):
        mat=basis.matrix(p);a=mat[good]/eg[:,None];b=yg/eg
        norms=np.linalg.norm(a,axis=0);norms=np.maximum(norms,1e-12)
        coef=nnls(a/norms,b,maxiter=500)[0]/norms
        residual=(mat[good]@coef-yg)/eg
        return (mat,coef,residual) if return_parts else residual
    seeds=[[-1.5,.18,0,1.3,3.,0,0],[-.5,.3,0,.8,4.5,-.2,.3]]
    attempts=[least_squares(solve,p,bounds=(LOW,HIGH),max_nfev=110,ftol=2e-6,xtol=2e-6,gtol=2e-6,x_scale='jac') for p in seeds[:starts]]
    opt=min(attempts,key=lambda r:np.sum(r.fun**2))
    mat,coef,res=solve(opt.x,True)
    return opt,mat,coef,res


def measurements(x,mat,coef,p,z):
    parts=mat*coef[None,:]
    continuum=parts[:,:3].sum(axis=1)
    broad=parts[:,7:9].sum(axis=1)
    # Model integrals use the continuous fitted profile; explicit coverage
    # diagnostics distinguish interpolation from directly supported data.
    ew=np.trapezoid(np.divide(broad,continuum,out=np.zeros_like(broad),where=continuum>0),x)
    hb=100*(coef[7]+coef[8])*(1+z)
    oiii=10*(coef[5]+coef[6])*(1+z)
    agn5100=coef[0];total5100=np.interp(5100,x,continuum)
    return dict(hb_flux=hb,oiii5008_flux=oiii,hb_oiii=hb/oiii if oiii>0 else np.nan,
        hb_ew_total_cont=ew,agn5100=agn5100,total5100=total5100,
        host_fraction5100=1-agn5100/total5100 if total5100>0 else np.nan,
        alpha_lambda=p[0],broad_component2_fraction=coef[8]/max(coef[7]+coef[8],1e-30))


def one(s,variant='baseline',save=True,tag=None):
    x,y,e,lsf,good=prepare(s)
    row={k:s[k] for k in ['name','epoch_id','source','night','mjd','z','error_kind','resolution_status']}
    row.update(variant=variant,n_pixels=int(good.sum()))
    if good.sum()<150:
        return row|dict(status='insufficient native error pixels')
    coverage=integration_weights(x,good,[(4750,4950)])[1]
    ocoverage=integration_weights(x,good,[(4985,5030)])[1]
    opt,mat,coef,res=fit_arrays(x,y,e,lsf,good,variant)
    near=(opt.x-LOW<.015*(HIGH-LOW))|(HIGH-opt.x<.015*(HIGH-LOW))
    dof=max(1,int(good.sum())-len(opt.x)-int((coef>1e-8).sum()))
    row.update(status='fit',optimizer_success=bool(opt.success),nfev=int(opt.nfev),
        reduced_chi2=float(np.sum(res**2)/dof),hb_window_coverage=coverage,oiii_window_coverage=ocoverage,
        parameter_bounds=';'.join(np.array(PARAMS)[near]),**measurements(x,mat,coef,opt.x,s['z']))
    row['hb_supported']=coverage>=.97
    row['oiii_supported']=ocoverage>=.97
    if save:
        folder=OUT/'fits';folder.mkdir(exist_ok=True)
        np.savez_compressed(folder/f"{s['epoch_id']}__{tag or variant}.npz",wave_rest=x,flux=y,error=e,good=good,
            lsf_sigma_rest=lsf,parameters=opt.x,coefficients=coef,components=mat*coef[None,:],component_names=NAMES)
    return row


def compare(results):
    rows=[]
    selected=pd.read_csv(OUT/'pilot_epochs.csv')
    missing=[]
    good=results[results.status.eq('fit')]
    for (name,variant),g in good.groupby(['name','variant']):
        old=g[g.source.ne('NGPS')].sort_values('mjd')
        new=g[g.source.eq('NGPS')].sort_values('mjd')
        for n in new.itertuples(index=False):
            before=selected[selected.name.eq(name)&selected.source.ne('NGPS')&selected.max_mjd.lt(n.mjd)].sort_values('max_mjd')
            refs=[]
            wanted=[]
            if len(before):wanted.append(('latest_archive',before.iloc[-1].epoch_id))
            if len(before)>1:wanted.append(('earliest_archive',before.iloc[0].epoch_id))
            for role,refid in wanted:
                fit=old[old.epoch_id.eq(refid)]
                if len(fit):refs.append((role,fit.iloc[0]))
                else:missing.append(dict(name=name,variant=variant,epoch_id=n.epoch_id,reference_epoch=refid,comparison=role,reason='date-selected reference lacks a usable fit; no older substitute'))
            repeat=new[new.mjd<n.mjd]
            if len(repeat):refs.append(('previous_ngps',repeat.iloc[-1]))
            for role,o in refs:
                r=dict(name=name,variant=variant,comparison=role,epoch_id=n.epoch_id,reference_epoch=o.epoch_id,
                    baseline_rest_years=(n.mjd-o.mjd)/365.25/(1+n.z),
                    both_hb_supported=bool(n.hb_supported and o.hb_supported),
                    both_oiii_supported=bool(n.oiii_supported and o.oiii_supported),
                    max_reduced_chi2=max(n.reduced_chi2,o.reduced_chi2),
                    any_parameter_bound=bool(n.parameter_bounds or o.parameter_bounds),
                    both_optimizers_success=bool(n.optimizer_success and o.optimizer_success))
                for key in ['hb_ew_total_cont','hb_oiii','hb_flux','agn5100','total5100']:
                    r[key+'_ratio']=getattr(n,key)/o[key] if o[key]>0 else np.nan
                rows.append(r)
    pairs=pd.DataFrame(rows);pairs.to_csv(OUT/'fit_comparisons.csv',index=False)
    pd.DataFrame(missing).to_csv(OUT/'unavailable_fit_comparisons.csv',index=False)
    env=[]
    for (name,epoch,reference,role),g in pairs.groupby(['name','epoch_id','reference_epoch','comparison']):
        r=dict(name=name,epoch_id=epoch,reference_epoch=reference,comparison=role,n_variants=len(g),
            both_hb_supported=bool(g.both_hb_supported.all()),both_oiii_supported=bool(g.both_oiii_supported.all()),
            any_parameter_bound=bool(g.any_parameter_bound.any()),max_reduced_chi2=g.max_reduced_chi2.max(),
            both_optimizers_success=bool(g.both_optimizers_success.all()))
        for key in ['hb_ew_total_cont','hb_oiii','hb_flux','agn5100','total5100']:
            r[key+'_ratio_min']=g[key+'_ratio'].min();r[key+'_ratio_max']=g[key+'_ratio'].max()
        env.append(r)
    pd.DataFrame(env).to_csv(OUT/'model_sensitivity_envelope.csv',index=False)
    return pairs


def main():
    spectra=load();path=OUT/'fit_results.csv'
    subset=set(sys.argv[1:]);spectra=[s for s in spectra if not subset or s['name'] in subset]
    results=pd.read_csv(path).to_dict('records') if path.exists() and subset else []
    for i,s in enumerate(spectra,1):
        for variant in VARIANTS:
            start=time.monotonic();r=one(s,variant)
            results=[a for a in results if (a['epoch_id'],a['variant'])!=(s['epoch_id'],variant)]+[r]
            pd.DataFrame(results).to_csv(path,index=False)
            print(i,len(spectra),s['epoch_id'],variant,r['status'],
                'EW',round(r.get('hb_ew_total_cont',np.nan),2),'chi2',round(r.get('reduced_chi2',np.nan),2),
                'seconds',round(time.monotonic()-start,1),flush=True)
    compare(pd.DataFrame(results).fillna({'parameter_bounds':''}))
    (OUT/'fit_method.json').write_text(json.dumps(dict(method=__doc__,variants=VARIANTS,parameter_names=PARAMS,
        lower_bounds=LOW.tolist(),upper_bounds=HIGH.tolist(),components=NAMES,
        fixed_oiii_wing='sigma=800 km/s, centroid=narrow centroid -200 km/s; free nonnegative amplitude',
        source_sha256=hashlib.sha256(open(__file__,'rb').read()).hexdigest(),
        limitations=['independent epochs, not a joint physical fit','no full LSF-matrix convolution',
        'host sigma fixed at 150 km/s with template resolution not deconvolved','Fe II added Gaussian sigma fixed per variant',
        'no calibration, aperture, telluric, covariance, or reddening nuisance posterior',
        'model envelope is not a confidence interval; no formal significance or type transitions inferred']),indent=2)+'\n')


if __name__=='__main__':main()
