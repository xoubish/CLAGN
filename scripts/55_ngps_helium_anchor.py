"""Validate and add the observed He I 3188.67 vacuum-A anchor to the U central slice.

Preserves the original ThAr identification/rejection mask and polynomial order.
Writes a proposed WaveCalib under audit/helium_refit; does not alter a running reduction.
"""
import importlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit
from numpy.polynomial.legendre import legfit, legval

env = importlib.import_module('54_reduce_ngps_sep23')
from pypeit.wavecalib import WaveCalib


def main():
    path=env.pypeit_file('u').parent/'Calibrations/WaveCalib_A_0_DET01.fits'
    wavecal=WaveCalib.from_file(str(path))
    central=int(np.argsort(wavecal.spat_ids)[len(wavecal.spat_ids)//2])
    fit=wavecal.wv_fits[central]
    # The isolated blue peak was detected by the automatic fit but lay ~3 pixels
    # from its extrapolated wavelength, outside the automatic 0.3-pixel ID tolerance.
    centers=np.asarray(fit.tcent)
    candidates=centers[(centers>150)&(centers<195)]
    assert len(candidates)==1, 'The blue-line identification requires inspection.'
    mu=float(candidates[0]);pixels=np.arange(int(mu)-9,int(mu)+10)
    def profile(x,amp,center,sigma,bg,slope):
        return amp*np.exp(-.5*((x-center)/sigma)**2)+bg+slope*(x-mu)
    parameters,cov=curve_fit(profile,pixels,fit.spec[pixels],p0=[10,mu,1,2,0],
        bounds=([0,mu-3,.5,-20,-2],[100,mu+3,6,20,2]))
    centroid=float(parameters[1]);error=float(np.sqrt(cov[1,1]))
    assert error < .2 and .5 < parameters[2] < 3
    old_mask=fit.pypeitfit.gpm.astype(bool)
    x=np.r_[centroid/fit.xnorm,fit.pypeitfit.xval]
    y=np.r_[3188.67,fit.pypeitfit.yval]
    mask=np.r_[True,old_mask]
    coeff=legfit(2*x[mask]-1,y[mask],4)
    residual=legval(2*x-1,coeff)-y
    assert abs(residual[0])<.1
    assert np.std(residual[1:][old_mask])<.15
    old_prediction=float(fit.pypeitfit.eval(centroid/fit.xnorm))
    fit.pypeitfit.xval=x
    fit.pypeitfit.yval=y
    fit.pypeitfit.gpm=mask.astype(int)
    fit.pypeitfit.weights=np.ones_like(x)
    fit.pypeitfit.fitc=coeff
    fit.pixel_fit=np.r_[centroid,fit.pixel_fit]
    fit.wave_fit=y
    fit.ion_bits=np.r_[fit.bitmask.turn_on(np.int64(0),'HeI'),fit.ion_bits].astype(fit.ion_bits.dtype)
    fit.wave_soln=fit.pypeitfit.eval(np.arange(len(fit.spec))/fit.xnorm)
    fit.rms=float(np.sqrt(np.mean(residual[mask]**2))/fit.cen_disp)
    out=env.BASE/'audit/helium_refit';out.mkdir(parents=True,exist_ok=True)
    wavecal.to_file(str(out/path.name),overwrite=True)
    report=dict(spat_id=int(fit.spat_id),line_vacuum_angstrom=3188.67,
        centroid_pixel=centroid,centroid_error_pixel=error,line_fwhm_pixel=float(2.35482*parameters[2]),
        initial_predicted_angstrom=old_prediction,initial_error_angstrom=old_prediction-3188.67,
        revised_anchor_residual_angstrom=float(residual[0]),
        revised_thar_rms_angstrom=float(np.std(residual[1:][old_mask])),
        revised_rms_pixel=fit.rms,n_original_accepted_thar=int(old_mask.sum()),
        model='Quartic Legendre; equal weight for new He anchor and previously accepted ThAr lines; original ThAr rejection mask preserved.',
        scope='Central U slice only. Side slices retain the original unanchored blue extrapolation.')
    (out/'anchor_report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
