"""Conditional R01 mock-MIRI refits and measurement precision.

This does not validate the shell physics or infer a dust mechanism. Historical
fits require an assumed 5% extra scatter. Forecast date, illumination and each
family's nonlinear parameters are held fixed; shell weights are refitted.
"""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
from pathlib import Path
import sys, json, hashlib
import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.optimize import minimize
from scipy.integrate import trapezoid
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE/'figure2_recovery'
sys.path[:0] = [str(ROOT), str(HERE/'p9694_predictivity')]
from analyse_forecasts import Audit, Constrained, Z


def band_weights(wave, lo=8., hi=13.):
    # Piecewise-linear integral Fnu dnu, including exact band edges.
    grid = np.r_[lo, wave[(wave>lo)&(wave<hi)], hi]
    mat = np.array([np.interp(grid, wave, row) for row in np.eye(len(wave))])
    return -trapezoid(mat, 299792458./(grid*1e-6), axis=1)


def diagnostics(wave, flux):
    """Total extracted spectral diagnostics; these do not separate a host."""
    flux = np.atleast_2d(flux)
    at = lambda x: np.array([np.interp(x, wave, f) for f in flux])
    baseline = lambda x, lo, hi: np.exp(np.log(at(lo))+(np.log(x/lo)/np.log(hi/lo))*np.log(at(hi)/at(lo)))
    warm = flux@band_weights(wave)*1e-26/(1+Z)
    return np.column_stack([warm, np.log(at(9.7)/baseline(9.7,6,14)),
                            np.log(at(18)/baseline(18,14,21))])


def main():
    OUT.mkdir(exist_ok=True)
    a = Audit(); nodes, ww = leggauss(64)
    a.delay_nodes = (nodes+1)/2; a.delay_weights = ww/2
    c = Constrained(a)
    archive = json.loads((HERE/'p9694_predictivity/adequacy_and_forecasts.json').read_text())
    limit = archive['adequacy_chi2_threshold_99pct']
    # Contiguous logarithmic top-hat bins with R approximately 100.
    edges = np.geomspace(4.9, 27.9, int(np.ceil(100*np.log(27.9/4.9)))+1)
    centers = np.sqrt(edges[:-1]*edges[1:]); rest = centers/(1+Z)
    q, qw = leggauss(8)
    waves = (edges[1:,None]+edges[:-1,None])/2+(edges[1:,None]-edges[:-1,None])*q/2
    weights = np.broadcast_to(qw/2, waves.shape)
    fits = {}
    for key, rec in archive['sensitivity_5pct_extra_scatter'].items():
        p = np.array(rec['parameters']); eta = rec['evolution']
        hist, offset, aw, yw, ineq, upper = c.matrices(p, eta)
        stat, amp, valid = c.fit(p, eta)
        m, off = a.design(p,eta,waves,weights,np.full(len(rest),2028.),future=1.)
        assert valid and stat < limit and c.forecast_valid(p,eta,amp,2028.,1.)
        fits[key] = dict(p=p, eta=eta, hist=hist, offset=offset, aw=aw, yw=yw,
                         ineq=ineq, upper=upper, amp=amp, m=m, off=off, chi2=stat)
    truth = fits['-0.5']['m']@fits['-0.5']['amp']+fits['-0.5']['off']
    # Precision at the REQUESTED S/N, not a new ETC prediction for this target.
    frac = np.interp(rest,[5.,12.,18.,23.],[.02,.02,1/30,1/30])
    staterr = frac*truth
    # Explicit sensitivity assumptions. These are coherent error modes; they
    # must not be averaged down by counting wavelength bins.
    scale_mode = .03*truth
    tilt_mode = .03*truth*np.log(rest/12.)/np.log(22./6.)
    # A smooth 200 K modified blackbody is an illustrative host-subtraction
    # residual shape, with uncertainty 10% of extracted flux at rest 12 um.
    hostshape = rest**(-4.5)/np.expm1(14387.76877/(rest*200.))
    hostshape /= np.interp(12.,rest,hostshape)
    host_mode = .10*np.interp(12.,rest,truth)*hostshape
    modes = np.column_stack([scale_mode,tilt_mode,host_mode])
    cov = np.diag(staterr**2)+modes@modes.T
    chol = np.linalg.cholesky(cov)
    rng = np.random.default_rng(96942026)
    noise = rng.normal(size=len(rest))*staterr+modes@rng.normal(size=3)
    mock = truth+noise
    recs = {}
    for key, fit in fits.items():
        m,off,aw,yw,ineq = [fit[k] for k in ('m','off','aw','yw','ineq')]
        mw = solve_triangular(chol,m,lower=True)
        target = solve_triangular(chol,mock-off,lower=True)
        cons = [dict(type='ineq',fun=lambda x:1-ineq@x,jac=lambda x:-ineq),
                dict(type='ineq',fun=lambda x:(limit-np.sum((aw@x-yw)**2))/len(yw),
                     jac=lambda x:-2*aw.T@(aw@x-yw)/len(yw))]
        opt = minimize(lambda x:np.sum((mw@x-target)**2)/len(rest),fit['amp'],
                       jac=lambda x:2*mw.T@(mw@x-target)/len(rest), method='SLSQP',
                       bounds=list(zip(np.zeros(12),fit['upper'])),constraints=cons,
                       options={'maxiter':500,'ftol':1e-10})
        valid = bool(opt.success and np.max(ineq@opt.x)<1.0001 and
                     np.sum((aw@opt.x-yw)**2)<=limit+.001 and
                     c.forecast_valid(fit['p'],fit['eta'],opt.x,2028.,1.))
        assert valid, (key,opt.message)
        model = m@opt.x+off
        nuisance = modes.T@np.linalg.solve(cov, mock-model)
        observed_model = model+modes@nuisance
        check_stat = np.sum(((observed_model-mock)/staterr)**2)+np.sum(nuisance**2)
        assert np.isclose(check_stat,np.sum((mw@opt.x-target)**2),rtol=1e-7)
        recs[key] = dict(eta=fit['eta'], chi2_history=float(np.sum((aw@opt.x-yw)**2)),
                         chi2_miri=float(np.sum((mw@opt.x-target)**2)), valid=valid,
                         amplitudes=opt.x.tolist(), parameters=fit['p'].tolist(),
                         nuisance_standard_deviations=nuisance.tolist(),
                         diagnostics=diagnostics(rest,model)[0].tolist())
        fit['model']=observed_model
        fit['intrinsic']=model
        print('refit',key,recs[key]['chi2_history'],recs[key]['chi2_miri'],flush=True)
    # Repeated noise realizations measure diagnostic precision only, not
    # posterior model probabilities or false-positive rates for mechanisms.
    draws = truth+rng.normal(size=(4000,len(rest)))*staterr+rng.normal(size=(4000,3))@modes.T
    diag = diagnostics(rest,draws)
    true_diag = diagnostics(rest,truth)[0]
    quantile = np.quantile(diag,[.16,.5,.84],axis=0)
    statdraws = truth+rng.normal(size=(4000,len(rest)))*staterr
    statdiag = diagnostics(rest,statdraws)
    # Numerical stability of forward spectra at twice the angular quadrature.
    n,w=leggauss(128);a.delay_nodes=(n+1)/2;a.delay_weights=w/2
    fit=fits['-0.5']
    m,o=a.design(fit['p'],fit['eta'],waves,weights,np.full(len(rest),2028.))
    quad=float(np.max(np.abs((m@fit['amp']+o)/truth-1)))
    data = dict(rest_um=rest, observed_um=centers, truth_mjy=truth, mock_mjy=mock,
                statistical_error_mjy=staterr, total_error_mjy=np.sqrt(np.diag(cov)))
    for key, fit in fits.items():
        data[f'fit_eta{key}_mjy']=fit['model']
        data[f'intrinsic_eta{key}_mjy']=fit['intrinsic']
    pd.DataFrame(data).to_csv(OUT/'mock_spectrum.csv',index=False)
    np.savez_compressed(OUT/'measurement_draws.npz',diagnostics=diag,statistical_diagnostics=statdiag,covariance=cov)
    report=dict(scope=__doc__,n_bins=len(rest),year=2028.,future_driver_factor=1.,seed=96942026,
                truth_eta=-.5,truth_diagnostics=true_diag.tolist(),measurement_quantiles_16_50_84=quantile.tolist(),
                measurement_std=np.std(diag,axis=0).tolist(),statistical_std=np.std(statdiag,axis=0).tolist(),
                fitted_families=recs,maximum_fractional_quadrature_change=quad,
                assumptions=dict(statistical='S/N 50 through rest 12 um, interpolating to 30 at 18 um; assumed throughout bins',
                                 calibration_scale_fraction=.03,calibration_tilt='3% coefficient times ln(lambda/12)/ln(22/6)',
                                 host_residual='200 K Fnu proportional to lambda^-4.5/(exp(hc/lambda kT)-1); 10% of flux at 12 um',
                                 historical_extra_scatter_fraction=.05,
                                 silicate_indices='Natural log of flux / local power law through 6,14 or 14,21 um; not full decomposition'),
                limitations=['No complete torus library; nonlinear parameters fixed to historical optima.',
                             'No model-selection significance or sample-wide recovery forecast follows.',
                             'Quoted errors are conditional measurement precision, excluding uncertainty in physical interpretation.',
                             'The assumed host shape does not span every compact host or PAH contribution.'],
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['truth_diagnostics','measurement_std','statistical_std','maximum_fractional_quadrature_change']},indent=2))

if __name__=='__main__':main()
