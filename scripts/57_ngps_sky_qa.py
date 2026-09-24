"""Measure observed-frame sky-line offsets without modifying spectra."""
import importlib
import json
import sys
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import fits
from astropy.table import Table
env=importlib.import_module('54_reduce_ngps_sep23')
from pypeit.core.wave import airtovac
from astropy import units as u

def main():
    lines=np.loadtxt(env.BASE/'software/PypeIt/pypeit/data/sky_spec/sky_single_mg.dat',skiprows=1)
    rows=[]
    corrected='--corrected' in sys.argv
    for ch in 'gri':
        refs=np.array([airtovac(5577.338*u.AA).value]) if ch=='g' else np.r_[lines,[6302.046,6365.536]]
        directory=env.NIGHT/'central_exposures'/ch if corrected else env.pypeit_file(ch).parent/'Science'
        for p in sorted(directory.glob('spec1d*fits')):
            with fits.open(p) as h:
                a=[x for x in h[1:] if x.header.get('SLITID')=={'g':240,'r':243,'i':236}[ch]]
                if len(a)!=1:continue
                d=a[0].data;w=d['OPT_WAVE']/a[0].header.get('VEL_CORR',1.)
                sky=d['OPT_COUNTS_SKY'];mask=d['OPT_MASK']&(w>0)
                for ref in refs:
                    sel=mask&(abs(w-ref)<7)
                    if sel.sum()<7:continue
                    x=w[sel]-ref;y=sky[sel]
                    amp=np.max(y)-np.median(y)
                    if amp<10:continue
                    def fun(x,a,mu,s,b,m):return a*np.exp(-.5*((x-mu)/s)**2)+b+m*x
                    try:
                        par,cov=curve_fit(fun,x,y,p0=[amp,0,1.5,np.median(y),0],
                            bounds=([0,-3,.4,-np.inf,-np.inf],[np.inf,3,4,np.inf,np.inf]))
                        err=np.sqrt(cov[1,1]);rms=np.std(y-fun(x,*par))
                    except (ValueError,RuntimeError):continue
                    if err>.3 or par[0]<10*rms or abs(par[1])>2.9:continue
                    rows.append(dict(file=p.name,target=h[0].header['TARGET'],channel=ch,
                        mjd=h[0].header['MJD'],line=ref,offset=par[1],error=err,
                        fwhm=2.35482*par[2],amplitude=par[0],fit_rms=rms))
    suffix='_corrected' if corrected else ''
    Table(rows=rows).write(env.BASE/f'audit/sky_line_measurements{suffix}.csv',overwrite=True)
    summary={}
    for ch in 'gri':
        q=[r for r in rows if r['channel']==ch and r['target'] not in ['P330E','BD284211']]
        z=np.array([r['offset'] for r in q])
        summary[ch]=dict(n=len(z),median=float(np.median(z)),std=float(np.std(z)),
            p16=float(np.percentile(z,16)),p84=float(np.percentile(z,84)))
    (env.BASE/f'audit/sky_line_summary{suffix}.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
