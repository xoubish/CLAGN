"""Test wavelength-dependent sky residuals with localized, resolution-matched atlas fits."""
import importlib,json,sys
import shutil
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.optimize import least_squares
env=importlib.import_module('54_reduce_ngps_sep23')

def measure(ch):
    with fits.open(env.BASE/'software/PypeIt/pypeit/data/sky_spec/paranal_sky.fits') as h:
        aw=h[1].data.copy();af=h[0].data.copy()
    step=np.median(np.diff(aw));sigma=1.8 if ch=='r' else 2.2
    smooth=gaussian_filter1d(af,sigma/step)
    lo,hi=(6200,7870) if ch=='r' else (7700,10200)
    sel=(aw>lo)&(aw<hi);sw=aw[sel];sf=smooth[sel]
    peaks,props=find_peaks(sf,distance=int(23/step),prominence=np.percentile(sf,80)*.4)
    centers=sw[peaks];records=[]
    for path in sorted((env.NIGHT/'central_exposures'/ch).glob('spec1d*fits')):
        with fits.open(path) as h:
            d=h[1].data;vel=h[1].header['VEL_CORR'];w=d['OPT_WAVE']/vel;y=d['OPT_COUNTS_SKY']
            valid=d['OPT_MASK']&(w>0)
            for ref in centers:
                q=valid&(abs(w-ref)<9)
                if q.sum()<7:continue
                x=w[q];obs=y[q];height=np.ptp(obs)
                if height<15:continue
                template=np.interp(x,aw,smooth)
                scale=height/max(np.ptp(template),1e-8)
                def residual(v):
                    shift,amp,bg,slope=v
                    return (amp*np.interp(x-shift,aw,smooth)+bg+slope*(x-ref)-obs)/height
                fit=least_squares(residual,[0,scale,np.median(obs),0],
                    bounds=([-3,0,-np.inf,-np.inf],[3,scale*5,np.inf,np.inf]),loss='soft_l1',f_scale=.05)
                rms=np.std(residual(fit.x))
                if abs(fit.x[0])>2.8 or rms>.12:continue
                records.append(dict(file=path.name,target=h[0].header['TARGET'],channel=ch,
                    wave=float(ref),offset=float(fit.x[0]),fractional_fit_rms=float(rms)))
    Table(rows=records).write(env.BASE/f'audit/sky_template_local_{ch}.csv',overwrite=True)
    # Fit nightly wavelength dependence using many science frames. Calibration
    # standards are held out. Per-frame offsets are estimated after this fit.
    sci=[r for r in records if r['target'] not in ['P330E','BD284211']]
    x=np.array([r['wave'] for r in sci]);y=np.array([r['offset'] for r in sci]);pivot=(lo+hi)/2
    fun=lambda v:v[0]+v[1]*(x-pivot)/1000-y
    model=least_squares(fun,[np.median(y),0],loss='soft_l1',f_scale=.12)
    res=fun(model.x)
    report=dict(channel=ch,n=len(x),pivot=pivot,offset_at_pivot=float(model.x[0]),
        slope_angstrom_per_1000A=float(model.x[1]),
        residual_mad_angstrom=float(1.4826*np.median(abs(res-np.median(res)))),
        wave_min=float(min(x)),wave_max=float(max(x)))
    (env.BASE/f'audit/sky_template_model_{ch}.json').write_text(json.dumps(report,indent=2))
    print(report)

def apply_i():
    import pandas as pd
    ch='i';d=pd.read_csv(env.BASE/'audit/sky_template_local_i.csv')
    grouped=d[~d.target.isin(['P330E','BD284211'])].groupby('wave').offset.agg(['count','median','std'])
    grouped=grouped[grouped['count']>=15]
    x=(grouped.index.to_numpy()-8950)/1000;y=grouped['median'].to_numpy()
    trials=[]
    for degree in [1,2,3]:
        f=least_squares(lambda v:np.polynomial.polynomial.polyval(x,v)-y,
            np.zeros(degree+1),loss='soft_l1',f_scale=.12)
        errors=[]
        for k in range(len(x)):
            q=np.arange(len(x))!=k
            ff=least_squares(lambda v:np.polynomial.polynomial.polyval(x[q],v)-y[q],
                f.x,loss='soft_l1',f_scale=.12)
            errors.append(float(np.polynomial.polynomial.polyval(x[k],ff.x)-y[k]))
        trials.append(dict(degree=degree,coefficients=f.x.tolist(),
            leave_one_line_out_median_abs=float(np.median(abs(np.array(errors)))),
            leave_one_line_out_p84_abs=float(np.percentile(abs(np.array(errors)),84))))
    choice=min(trials,key=lambda r:r['leave_one_line_out_median_abs'])
    assert choice['leave_one_line_out_median_abs']<.2
    coeff=np.array(choice['coefficients']);lo=float(grouped.index.min());hi=float(grouped.index.max())
    indir=env.NIGHT/'central_exposures/i';backup=env.BASE/'audit/i_before_sky_refinement'
    assert not backup.exists(),'Refinement already applied; do not apply twice.'
    shutil.copytree(indir,backup)
    report=dict(pivot=8950,valid_min=lo,valid_max=hi,trials=trials,selected=choice,frames=[])
    for p in sorted(indir.glob('*fits')):
        q=d[d.file==p.name];res=q.offset-np.polynomial.polynomial.polyval((q.wave-8950)/1000,coeff)
        delta=float(np.median(res)) if len(res)>=5 else 0.
        assert abs(delta)<1.5,(p,delta)
        with fits.open(p,mode='update') as h:
            vel=h[1].header['VEL_CORR']
            for ext in ['OPT','BOX']:
                w=h[1].data[ext+'_WAVE'];valid=w>0;obs=w[valid]/vel
                correction=np.polynomial.polynomial.polyval((np.clip(obs,lo,hi)-8950)/1000,coeff)+delta
                w[valid]-=correction*vel
            h[0].header['SKYREFIN']='Cubic atlas residual + frame median; clipped at fitted endpoints'
            h[0].header['SKYEXTRA']=delta
        report['frames'].append(dict(file=p.name,n_lines=len(res),offset=delta,
            residual_mad=float(1.4826*np.median(abs(res-delta)))))
    (env.BASE/'audit/i_sky_refinement_applied.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({k:v for k,v in report.items() if k!='frames'},indent=2))

if __name__=='__main__':
    if '--apply-i' in sys.argv:apply_i()
    else:
        for ch in sys.argv[1:] or ['r','i']:measure(ch)
