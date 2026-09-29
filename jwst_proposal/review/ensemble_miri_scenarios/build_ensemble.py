"""A uniform 24-target, WISE-conditioned MIRI sensitivity grid.

Reuses the established shell model, heating proxy and host/size assumptions.
These are forward scenarios, not new fits to all 24 optical/SPHEREx spectra.
Temperature-grid interpolation accelerates the unchanged shell integral.
"""
from pathlib import Path
import os,sys,json,hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/clagn-matplotlib')
HERE=Path(__file__).resolve().parent; BASE=HERE.parents[1];sys.path.insert(0,str(BASE))
import numpy as np,pandas as pd
from scipy.integrate import trapezoid
import forecast_warm_response as F
from make_fig2_diagnostics import GRAINS,HEATING,planck_nu,WAVE as QWAVE

WAVE=np.unique(np.r_[np.geomspace(2.,22.,180),4.,6.,8.,9.7,12.,13.,14.,18.,21.,22.])
LOGT=np.linspace(np.log(10.),np.log(6000.),2400)
PLANCK=planck_nu(WAVE,np.exp(LOGT))
HOST_OFFSETS=[-.15,0.,.15]
RADIUS_FACTORS=[.5,1.,2.]
HYPOTHESES=['delayed','evolving']
EPOCH=2028.


def spectra(t,lum,scale,z):
    radii=F.R_GRID
    delay=2*radii[:,None]*((np.arange(F.N_DELAY)+.5)/F.N_DELAY)*scale*(1+z)
    illumination=np.interp(EPOCH-delay,t,lum)
    out={h:np.zeros(len(WAVE)) for h in HYPOTHESES}
    for grain in GRAINS:
        power=HEATING*grain['q_heat']/radii[:,None]**2*illumination
        temperature=np.interp(np.log(power),grain['log_cool'],grain['log_t'])
        fractional=np.clip((temperature-LOGT[0])/(LOGT[1]-LOGT[0]),0,len(LOGT)-1)
        low=np.minimum(fractional.astype(int),len(LOGT)-2);frac=fractional-low
        for h in HYPOTHESES:
            mask=(radii[:,None]>=grain['inner_radius']*np.sqrt(illumination)) if h=='evolving' else (radii[:,None]>=grain['inner_radius'])
            weights=np.broadcast_to(mask*radii[:,None]*F.DLOGR*grain['mass_fraction']/grain['rho']/F.N_DELAY,illumination.shape)
            hist=np.bincount(low.ravel(),weights=(weights*(1-frac)).ravel(),minlength=len(LOGT))
            hist+=np.bincount((low+1).ravel(),weights=(weights*frac).ravel(),minlength=len(LOGT))
            q=np.interp(WAVE,QWAVE,grain['q'])
            out[h]+=hist@PLANCK*q
    return out


def diagnostic(flux,z):
    at=lambda w:float(np.interp(w,WAVE,flux))
    cont=lambda w,x:np.exp(np.interp(np.log(w),np.log(x),np.log([at(v) for v in x])))
    use=(WAVE>=8)&(WAVE<=13)
    # Dimensionless integrated warm energy / nu Fnu at the same model 4um epoch.
    warm=-trapezoid(flux[use],1/WAVE[use])/(at(4)/4)
    return dict(warm_over_hot_energy=float(warm),F14_over_F4=at(14)/at(4),
                F21_over_F14=at(21)/at(14),S9_7_local=float(np.log(at(9.7)/cont(9.7,[6,14]))),
                S18_local=float(np.log(at(18)/cont(18,[14,21]))),
                mrs_rest_lo=4.9/(1+z),mrs_rest_hi=27.9/(1+z))


def main():
    sample=pd.read_csv(BASE/'inputs/jwst_sample_cycle6.csv')
    history=pd.read_csv(BASE/'review/history_audit/binned_measurements.csv')
    allwise=pd.read_csv(BASE/'inputs/sample_allwise_psd.csv')
    meta=json.loads((BASE/'inputs/warm_response_forecast_summary.json').read_text())
    F.DRIVER_GAMMA=meta['driver_gamma']
    byname={g['name']:g for g in GRAINS}
    rin_ratio=byname['Sil_21.gz']['inner_radius']/byname['Gra_21.gz']['inner_radius']
    curves=[];rows=[];checks=[]
    for _,target in sample.iterrows():
        aw=allwise[allwise.id.eq(target.id)].iloc[0];host0=F.host_fraction(aw)
        tau=target.tau_dust_yr if np.isfinite(target.tau_dust_yr) else F.koshida_lag_years(target.r_mag,target.z)
        actual_hosts=np.unique([np.clip(host0+offset,0,.9) for offset in HOST_OFFSETS])
        for host in actual_hosts:
            t,lum,t0=F.history(target,history,allwise,float(host))
            # Record where the established heating-proxy floor was applied.
            neo=history[history.id.eq(target.id)&history.band.eq('W1')].sort_values('mjd')
            total_ratio=np.r_[1.,10**(-.4*(neo.mag.to_numpy()-aw.w1mpro))]
            clipped=int(np.sum((total_ratio-host)/(1-host)<.05))
            for radius_factor in RADIUS_FACTORS:
                scale=rin_ratio*tau*radius_factor
                model=spectra(t,lum,scale,target.z)
                for h,flux in model.items():
                    norm=float(np.interp(4.,WAVE,flux));assert norm>0
                    normalized=flux/norm
                    assert np.isfinite(normalized).all() and (normalized>0).all()
                    index=len(curves);curves.append(normalized)
                    rows.append(dict(curve_index=index,id=target.id,family=target.family,z=target.z,
                                     hypothesis=h,host_fraction=float(host),host_baseline=host0,
                                     radius_factor=radius_factor,rin_si_years=scale,
                                     heating_floor_count=clipped,n_heating_points=len(t),
                                     last_wise_year=float(t[-1]),**diagnostic(normalized,target.z)))
                # Independently compare the acceleration to the existing direct integrator.
                if np.isclose(host,host0) and radius_factor==1:
                    for h in HYPOTHESES:
                        for w in [4.,9.7,18.]:
                            direct=F.model_flux(np.array([EPOCH]),t,lum,w,scale,h,target.z)[0]
                            accelerated=float(np.interp(w,WAVE,model[h]))
                            error=float(abs(accelerated/direct-1))
                            assert error<.001,(target.id,h,w,error)
                            checks.append(dict(id=target.id,hypothesis=h,rest_um=w,fractional_error=error))
        print(target.id,'completed',flush=True)
    table=pd.DataFrame(rows);curves=np.array(curves)
    selected=[]
    for ident,group in table.groupby('id',sort=False):
        for label,index in [('weak',group.warm_over_hot_energy.idxmin()),('strong',group.warm_over_hot_energy.idxmax())]:
            record=table.loc[index].to_dict();record['response_case']=label;selected.append(record)
    endpoints=pd.DataFrame(selected)
    assert table.id.nunique()==24 and len(endpoints)==48
    table.to_csv(HERE/'scenario_grid.csv',index=False);endpoints.to_csv(HERE/'selected_endpoints.csv',index=False)
    np.savez_compressed(HERE/'spectra.npz',rest_um=WAVE,normalized_flux=curves)
    q=[]
    for case in ['weak','strong']:
        chosen=curves[endpoints.loc[endpoints.response_case.eq(case),'curve_index'].to_numpy()]
        quantiles=np.percentile(chosen,[0,10,25,50,75,90,100],axis=0)
        for i,w in enumerate(WAVE):
            q.append(dict(response_case=case,rest_um=w,**{f'p{p}':quantiles[j,i] for j,p in enumerate([0,10,25,50,75,90,100])}))
    pd.DataFrame(q).to_csv(HERE/'spectral_bands.csv',index=False)
    provenance=dict(targets=24,curves=len(table),epoch=EPOCH,
        input_scope='Measured AllWISE plus NEOWISE W1 histories for all 24; W1-W2 estimates host fraction. Not a joint optical/SPHEREx fit; extracted SPHEREx available for only three targets in these inputs.',
        model='Unchanged optically thin shell physics from forecast_warm_response.py. Fixed mass distribution vs immediate L^1/2 boundary response, both with light-travel delays.',
        nuisance_grid=dict(host_fraction_offsets=HOST_OFFSETS,host_clip=[0,.9],radius_factors=RADIUS_FACTORS,driver_gamma=F.DRIVER_GAMMA),
        extrapolation='Constant illumination before first and after last W1 datum.',
        normalization='Each nuclear model spectrum divided by its own 2028 rest-4um Fnu. No observed near-IR measurements rescaled into the bands.',
        endpoint_definition='For each target select the minimum and maximum integrated rest-8–13um energy / (nu Fnu at 4um), over all host/size/response scenarios.',
        spectral_shades='Pointwise interquartile (dark) and 10th–90th percentile (light) spreads across 24 object endpoints, each target weighted once in each case. Median curves. Not confidence intervals or realizable spectra.',
        diagnostic_shades='Convex hulls of the 24 endpoint positions for each response case. Geometric scenario spans, not likelihood or classification boundaries.',
        diagnostic_definitions='F14/F4 and F21/F14 are spectral ratios. S9.7=ln(F9.7/local power law through 6 and 14um); S18 uses 14 and 21um. Same model epoch. Full decomposition would replace these local indices.',
        common_mrs_rest_um=[float(4.9/(1+sample.z.min())),float(27.9/(1+sample.z.max()))],
        models_with_heating_floor=int(table.heating_floor_count.gt(0).sum()),
        numerical_check_max_error=max(r['fractional_error'] for r in checks),numerical_checks=checks,
        source_hashes={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [BASE/'forecast_warm_response.py',BASE/'make_fig2_diagnostics.py',BASE/'inputs/jwst_sample_cycle6.csv',BASE/'inputs/sample_allwise_psd.csv',BASE/'review/history_audit/binned_measurements.csv',BASE/'inputs/warm_response_forecast_summary.json']})
    (HERE/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps({k:provenance[k] for k in ['targets','curves','models_with_heating_floor','numerical_check_max_error']},indent=2))


if __name__=='__main__':main()
