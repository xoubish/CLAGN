"""Limited injection recovery and alternative-standard sensitivity, not completeness."""
import json
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from test_pilot_fits import synthetic,recovered
from pilot_fits import Basis,fit_arrays,measurements,one,VARIANTS
from pilot_data import load
from paper_protocol import OUT


def injections():
    rows=[]
    for seed in range(12):
        for truth in [.5,1.,2.]:
            a=synthetic(noise=1,seed=seed,resolution=2000)
            b=list(synthetic(noise=1,seed=100+seed,host=1.7,broad=truth,grey=.65))
            # An explicitly out-of-model continuum curvature and correlated
            # residual challenge is retained alongside the matched-model case.
            for challenge in ['matched_model','curvature_and_correlated_residual']:
                bb=list(b)
                if challenge!='matched_model':
                    w=bb[0];rng=np.random.default_rng(1000+seed)
                    residual=gaussian_filter1d(rng.normal(size=len(w)),1.2)
                    residual/=np.std(residual)
                    bb[1]=bb[1]*(1+.08*((w-4850)/650)**2)+residual*bb[2]
                ra=recovered(a);rb=recovered(bb)
                rows.append(dict(seed=seed,challenge=challenge,true_hb_oiii_ratio=truth,
                    recovered_hb_oiii_ratio=rb['hb_oiii']/ra['hb_oiii'],
                    fractional_error=(rb['hb_oiii']/ra['hb_oiii'])/truth-1))
    pd.DataFrame(rows).to_csv(OUT/'injection_recovery.csv',index=False)
    print(pd.DataFrame(rows).groupby(['challenge','true_hb_oiii_ratio']).fractional_error.agg(['median','std','min','max']).to_string())


def alternate_standards():
    spectra=load();base=pd.read_csv(OUT/'fit_results.csv');rows=[]
    for s in spectra:
        alt=s.get('alternate')
        if s['source']!='NGPS' or alt is None or not np.isfinite(alt).any():continue
        # These are two deterministic reductions of the same exposure, not
        # independent measurements. Scale the diagonal errors by the response
        # ratio solely to preserve fractional statistical weights in the refit.
        a=dict(s);ratio=np.divide(alt,s['flux'],out=np.full(len(alt),np.nan),where=s['flux']!=0)
        a['flux']=alt.copy();a['error']=s['error']*abs(ratio)
        for variant in VARIANTS:
            r=one(a,variant,tag=variant+'_alternate')
            old=base[base.epoch_id.eq(s['epoch_id'])&base.variant.eq(variant)].iloc[0]
            for key in ['hb_flux','hb_ew_total_cont','hb_oiii','agn5100','total5100']:
                r[key+'_alternate_over_primary']=r.get(key,np.nan)/old[key] if old[key]>0 else np.nan
            rows.append(r)
    pd.DataFrame(rows).to_csv(OUT/'alternate_standard_fits.csv',index=False)
    print('Alternative-standard fits:',len(rows))
    # Keep the same host/Fe family across epochs, allowing the standard choice
    # on each observing night to vary independently. This is a stress envelope,
    # not an assertion that all standards are equally plausible or a posterior.
    alt=pd.DataFrame(rows).set_index(['epoch_id','variant'])
    pairs=pd.read_csv(OUT/'fit_comparisons.csv');stress=[]
    keys=['hb_flux','hb_ew_total_cont','hb_oiii','agn5100','total5100']
    for p in pairs.to_dict('records'):
        nk=(p['epoch_id'],p['variant']);ok=(p['reference_epoch'],p['variant'])
        for nalt in ([False,True] if nk in alt.index else [False]):
            for oalt in ([False,True] if ok in alt.index else [False]):
                r=p.copy();r.update(new_standard='alternate' if nalt else 'primary',reference_standard='alternate' if oalt else 'primary')
                for key in keys:
                    n=alt.loc[nk,key+'_alternate_over_primary'] if nalt else 1.
                    o=alt.loc[ok,key+'_alternate_over_primary'] if oalt else 1.
                    r[key+'_ratio']=p[key+'_ratio']*n/o
                stress.append(r)
    stress=pd.DataFrame(stress);stress.to_csv(OUT/'calibration_model_comparisons.csv',index=False)
    envelope=[]
    for (name,epoch,ref,comparison),g in stress.groupby(['name','epoch_id','reference_epoch','comparison']):
        r=dict(name=name,epoch_id=epoch,reference_epoch=ref,comparison=comparison,n_fits=len(g),
               both_hb_supported=bool(g.both_hb_supported.all()),both_oiii_supported=bool(g.both_oiii_supported.all()))
        for key in keys:
            r[key+'_ratio_min']=g[key+'_ratio'].min();r[key+'_ratio_max']=g[key+'_ratio'].max()
        envelope.append(r)
    pd.DataFrame(envelope).to_csv(OUT/'calibration_model_envelope.csv',index=False)


if __name__=='__main__':
    import sys
    if '--calibration-only' not in sys.argv:injections()
    alternate_standards()
