"""Physical-invariance and recovery checks for the conditional pilot fitter."""
import unittest
import numpy as np
from pilot_fits import Basis,fit_arrays,measurements,prepare


def synthetic(host=1.,broad=1.,grey=1.,resolution=1650,noise=0,seed=5):
    w=np.arange(4200.,5500.,1.2);lsf=w/resolution/2.35482
    p=np.array([-1.4,.2,.02,1.4,3.4,-.12,.23])
    c=np.array([15.,8.*host,3.*host,1.6,2.,10.,2.,7.*broad,3.*broad,1.,3.,.6,.4])
    mat=Basis(w,lsf).matrix(p);e=np.full(len(w),.8)
    y=mat@c+noise*np.random.default_rng(seed).normal(size=len(w))*e
    return w,y*grey,e*grey,lsf,np.ones(len(w),bool),p,c*grey


def recovered(s):
    x,y,e,lsf,g,_,_=s
    opt,mat,c,_=fit_arrays(x,y,e,lsf,g)
    return measurements(x,mat,c,opt.x,.2)


class PilotFitTests(unittest.TestCase):
    def test_grey_calibration_invariance(self):
        a=recovered(synthetic(noise=1));b=recovered(synthetic(noise=1,grey=.37))
        for key in ['hb_ew_total_cont','hb_oiii']:
            self.assertAlmostEqual(a[key]/b[key],1.,places=3)
        self.assertAlmostEqual(b['hb_flux']/a['hb_flux'],.37,places=3)

    def test_host_dilution_does_not_create_line_flux_change(self):
        a=recovered(synthetic(host=.5));b=recovered(synthetic(host=2.))
        self.assertLess(b['hb_ew_total_cont'],a['hb_ew_total_cont'])
        self.assertLess(abs(b['hb_oiii']/a['hb_oiii']-1),.015)

    def test_known_line_change_across_resolution_and_calibration(self):
        a=recovered(synthetic(resolution=2000));b=recovered(synthetic(broad=.5,grey=1.8,resolution=1650))
        self.assertLess(abs(b['hb_oiii']/a['hb_oiii']-.5),.015)
        self.assertLess(abs(b['hb_ew_total_cont']/a['hb_ew_total_cont']-.5),.015)

    def test_observed_flux_redshift_factor(self):
        x,_,_,lsf,_,p,c=synthetic();mat=Basis(x,lsf).matrix(p)
        a=measurements(x,mat,c,p,0);b=measurements(x,mat,c,p,1)
        self.assertAlmostEqual(b['hb_flux']/a['hb_flux'],2.)
        self.assertAlmostEqual(b['hb_ew_total_cont']/a['hb_ew_total_cont'],1.)

    def test_missing_errors_and_masked_samples_are_not_imputed(self):
        x,y,e,lsf,_,_,_=synthetic();e[:20]=np.nan;y[400:410]=np.nan
        s=dict(wave=x,flux=y,error=e,z=0,lsf_sigma=lsf)
        _,_,_,_,good=prepare(s)
        self.assertFalse(good[:20].any());self.assertFalse(good[400:410].any())


if __name__=='__main__':unittest.main()
