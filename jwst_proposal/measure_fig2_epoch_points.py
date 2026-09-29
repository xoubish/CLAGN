"""Measured SPHEREx synthetic W1 and optical continuum epochs for Figure 2."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from astropy.time import Time
from spherex_data import load_spherex, CSV_PATHS

HERE=Path(__file__).resolve().parent
Z=.2377224

def spherex_w1():
    data,groups,_=load_spherex('P9694')
    raw=pd.read_csv(CSV_PATHS['P9694'])
    response=np.loadtxt(HERE/'inputs/wise_response/W1.txt')
    grid=np.linspace(response[0,0],response[-1,0],8001)
    throughput=np.interp(grid,response[:,0],response[:,1])/grid
    full=trapezoid(throughput,grid)
    rows=[]
    for visit,g in enumerate(groups,1):
        idx=g['indices'];wave=data['wavelength_um'][idx]
        assert np.all(np.diff(wave)>0), 'Duplicate wavelengths require combination before interpolation'
        use=(grid>=wave.min())&(grid<=wave.max())
        x,t=grid[use],throughput[use]
        coverage=trapezoid(t,x)/full
        assert coverage>.999, 'Incomplete W1 spectral coverage'
        # Integrate a piecewise-linear Fnu spectrum with the photon response.
        # Flat-Fnu band-average convention matches the catalogue conversion
        # used for the W1 light curve; do not apply an additional color correction.
        basis=np.array([np.interp(x,wave,e) for e in np.eye(len(wave))])
        coeff=trapezoid(basis*t,x,axis=1)/trapezoid(t,x)
        assert np.isclose(coeff.sum(),1.)
        flux=float(coeff@data['flux_mjy'][idx])
        stat=float(np.sqrt(coeff**2@raw.variance_measurement_mjy2.to_numpy()[idx]))
        # The supplied calibration variance is propagated coherently within
        # each visit, rather than averaged down as independent spectral noise.
        cal=float(coeff@np.sqrt(raw.variance_calibration_mjy2.to_numpy()[idx]))
        mjd=float(coeff@data['mjds'][idx])
        rows.append(dict(kind='SPHEREx synthetic W1',visit=visit,label=g['label'],
                         mjd=mjd,year=float(Time(mjd,format='mjd').decimalyear),
                         flux_mjy=flux,error_mjy=float(np.hypot(stat,cal)),
                         statistical_error_mjy=stat,calibration_error_mjy=cal,
                         response_coverage=coverage,n_measurements=len(idx)))
    assert len(rows)==3
    return rows

def continuum(wave,flux,error,good):
    rest=wave/(1+Z);conversion=wave**2/2.99792458e18*1e9
    f=flux*conversion;e=error*conversion
    values=[];variances=[];centers=[];counts=[]
    for lo,hi in [(6400.,6450.),(6750.,6800.)]:
        keep=good&(rest>=lo)&(rest<=hi)&np.isfinite(f)&np.isfinite(e)&(e>0)
        assert keep.sum()>=8, 'Insufficient continuum pixels'
        values.append(np.mean(f[keep]));variances.append(np.sum(e[keep]**2)/keep.sum()**2)
        centers.append(np.mean(rest[keep]));counts.append(int(keep.sum()))
    fraction=(6564.61-centers[0])/(centers[1]-centers[0])
    coeff=np.array([1-fraction,fraction])
    return float(coeff@values),float(np.sqrt(coeff**2@variances)),counts

def optical_points():
    archival=pd.read_csv(HERE/'inputs/archival_spectra/R01/R01_mjd58424_eboss.csv')
    ngps=pd.read_csv(HERE.parent/'sep23_data/reduction_20260924/products_p330e/spectra/P9694.csv')
    error=np.full(len(archival),np.nan);good=archival.ivar.to_numpy()>0
    error[good]=1/np.sqrt(archival.ivar.to_numpy()[good])
    base,baseerr,nbase=continuum(archival.wave_A.to_numpy(),archival.flux.to_numpy(),error,good)
    flux,err,n=continuum(ngps.WAVE_VAC_HELIO_A.to_numpy(),ngps.FLUX.to_numpy(),
                         ngps.FLUX_ERR.to_numpy(),ngps.MASK.to_numpy().astype(bool))
    target=json.loads((HERE/'inputs/P9694_figure_target.json').read_text())
    epoch=next(e for e in target['spec']['epochs'] if e.get('instrument')=='NGPS')
    ratio=flux/base;ratioerr=ratio*np.hypot(err/flux,baseerr/base)
    return dict(kind='NGPS continuum near Halpha',mjd=epoch['mjd'],
                year=float(Time(epoch['mjd'],format='mjd').decimalyear),flux_mjy=flux,error_mjy=err,
                baseline_mjd=58424.,baseline_flux_mjy=base,baseline_error_mjy=baseerr,
                relative_flux=ratio,relative_error=ratioerr,pixels_per_window=n,baseline_pixels_per_window=nbase,
                rest_windows_A=[[6400,6450],[6750,6800]],
                uncertainty_scope='Formal propagated pixel errors only; inter-pixel covariance, absolute calibration and aperture/slit-loss differences are not quantified.')

def measurements():
    spx=spherex_w1()
    result=dict(spherex=spx,
                spherex_method='Piecewise-linear native Fnu integrated against W1 photon response / wavelength; flat-Fnu band average. Supplied measurement errors plus fully correlated within-visit calibration term.',
                reference='https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/sec4_4h.html')
    (HERE/'inputs/fig2_epoch_points.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':print(json.dumps(measurements(),indent=2))
