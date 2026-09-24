"""Central-slice selection, measured sky flexure, fluxing, and native-channel coadds."""
import importlib
import sys
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
env=importlib.import_module('54_reduce_ngps_sep23')
from pypeit import specobjs, coadd1d, fluxcalibrate
from pypeit.onespec import OneSpec
from pypeit.core.flexure import spec_flex_shift, flexure_interp
from pypeit.core.wave import airtovac
from pypeit.spectrographs.util import load_spectrograph

SLITS=dict(u=550,g=240,r=243,i=236)
SPATS=dict(u=545,g=235,r=236,i=229)

def prepare(ch):
    out=env.NIGHT/'central_exposures'/ch;out.mkdir(parents=True,exist_ok=True)
    records=[];objects=[]
    for path in sorted((env.pypeit_file(ch).parent/'Science').glob('spec1d*fits')):
        sobjs=specobjs.SpecObjs.from_fitsfile(str(path))
        selected=[o for o in sobjs if o.SLITID==SLITS[ch]]
        assert len(selected)==1,(path,[o.NAME for o in selected])
        obj=selected[0]
        assert abs(obj.SPAT_PIXPOS-SPATS[ch])<10,(path,obj.SPAT_PIXPOS)
        w=obj.OPT_WAVE/obj.VEL_CORR;sky=obj.OPT_COUNTS_SKY
        mask=obj.OPT_MASK&(w>0)&np.isfinite(sky)
        shift=0.;method='none: U has no reliable sky-line constraint'
        if ch=='g':
            ref=airtovac(5577.338*u.AA).value
            q=mask&(abs(w-ref)<10);x=w[q]-ref;y=sky[q]
            def fn(x,a,mu,s,b,m):return a*np.exp(-.5*((x-mu)/s)**2)+b+m*x
            try:
                pars,cov=curve_fit(fn,x,y,p0=[np.ptp(y),-1.5,1.3,np.median(y),0],
                    bounds=([0,-5,.4,-np.inf,-np.inf],[np.inf,5,4,np.inf,np.inf]))
                assert np.sqrt(cov[1,1])<.25 and pars[0]>10*np.std(y-fn(x,*pars))
                shift=-pars[1]/np.median(np.diff(w[q]));method='OI 5577.338 air; Gaussian + linear background'
            except (AssertionError,ValueError,RuntimeError):shift=np.nan;method='insufficient sky S/N'
        elif ch in 'ri':
            spec=OneSpec(wave=w[mask],wave_grid_mid=w[mask],flux=sky[mask],
                PYP_SPEC=f'p200_ngps_{ch}',ext_mode='OPT',fluxed=False,mask=np.ones(mask.sum(),dtype=int))
            try:
                result=spec_flex_shift(spec,sky_file='paranal_sky.fits',spec_fwhm_pix=2.5,
                    mxshft=5,minwave=6400 if ch=='r' else 8000,maxwave=7900 if ch=='r' else 9500)
                shift=float(result['shift']);assert abs(shift)<4
                method='PypeIt sky cross-correlation against paranal_sky.fits'
            except (Exception,) as ex:
                shift=np.nan;method='insufficient sky constraint: '+str(ex)
        records.append(dict(file=path.name,target=sobjs.header['TARGET'],mjd=float(sobjs.header['MJD']),
            channel=ch,spat=float(obj.SPAT_PIXPOS),slit=int(obj.SLITID),
            snr=float(np.median((obj.OPT_COUNTS*np.sqrt(obj.OPT_COUNTS_IVAR))[mask])),
            shift_pixel=shift,method=method,heliocentric_factor=float(obj.VEL_CORR)))
        objects.append((sobjs,obj))
    finite=[r['shift_pixel'] for r in records if np.isfinite(r['shift_pixel']) and r['target'] not in ['P330E','BD284211']]
    assert len(finite)>30
    for r,(sobjs,obj) in zip(records,objects):
        if not np.isfinite(r['shift_pixel']):
            r['shift_pixel']=float(np.median(finite));r['method']+='; used science night median'
        shift=r['shift_pixel']
        if shift:
            for ext in ['OPT','BOX']:
                wave=obj[ext+'_WAVE']
                # Invalid wavelength zeros at the edges must not enter interpolation.
                good=wave>0
                fixed=wave.copy();fixed[good]=flexure_interp(shift,wave[good])
                obj[ext+'_WAVE']=fixed
            obj.FLEX_SHIFT_LOCAL=shift
            obj.FLEX_SHIFT_TOTAL=shift
        header=sobjs.header.copy()
        header['SKYSHFT']=(shift,'Applied extra spectral flexure, pixels')
        header['SLICE']='CENTRAL'
        header['HISTORY']='Sky reference: observed-frame wavelengths; heliocentric factor retained.'
        one=specobjs.SpecObjs(specobjs=np.array([obj]),header=header)
        one.write_to_fits(header,str(out/r['file']),overwrite=True)
    (env.BASE/f'audit/central_selection_{ch}.json').write_text(json.dumps(records,indent=2))
    print('PREPARED',ch,len(records),flush=True)

def flux_and_coadd(ch, standard=None, variant=None):
    suffix='' if variant is None else '_'+variant
    indir=env.NIGHT/'central_exposures'/ch
    out=env.NIGHT/('fluxed_exposures'+suffix)/ch;out.mkdir(parents=True,exist_ok=True)
    standards=env.BASE/'audit/standards'/ch
    paths=sorted(indir.glob('spec1d*fits'));groups=defaultdict(list)
    for p in paths:groups[fits.getheader(p)['TARGET'].upper()].append(p)
    sens=[]
    for p in sorted(standards.glob('sens_spec1d*fits')):
        h=fits.getheader(standards/p.name.removeprefix('sens_'))
        sens.append((float(h['MJD']),p,h['TARGET']))
    assert len(sens)==4
    if standard is not None:
        sens=[row for row in sens if row[2]==standard]
        assert len(sens)==2, (standard,sens)
    assignments=[]
    for target,files in groups.items():
        # Select once per target so a pair straddling the midpoint never uses
        # incompatible flux scales. Both standards' responses are delivered for QA.
        mjd=np.mean([fits.getheader(p)['MJD'] for p in files])
        _,sf,star=min(sens,key=lambda row:abs(row[0]-mjd))
        fluxcalibrate.flux_calibrate([str(p) for p in files],[str(sf)]*len(files),
            outfiles=[str(out/p.name) for p in files])
        rows=[]
        for p in files:
            fp=out/p.name;h=fits.getheader(fp,1);rows.append((str(fp),h['NAME']))
            assignments.append(dict(target=target,channel=ch,file=p.name,standard=star,sensfile=str(sf)))
        spec=load_spectrograph(f'p200_ngps_{ch}')
        par=spec.default_pypeit_par()['coadd1d']
        par['scale_method']='none';par['ex_value']='OPT';par['flux_value']=True
        c=coadd1d.CoAdd1D.get_instance([r[0] for r in rows],[r[1] for r in rows],spectrograph=spec,par=par)
        c.run()
        coadddir=env.NIGHT/('coadded_channels'+suffix);coadddir.mkdir(exist_ok=True)
        cp=coadddir/f'{target}_{ch.upper()}.fits';c.save(str(cp))
        with fits.open(cp,mode='update') as h:
            h[0].header['FLUXSTD']=star;h[0].header['NEXP']=len(files)
            h[0].header['TEXPTIME']=sum(fits.getheader(p)['EXPTIME'] for p in files)
            h[0].header['AIRMASS']=float(np.mean([fits.getheader(p)['AIRMASS'] for p in files]))
            h[0].header['WAVEREF']='HELIOCENTRIC VACUUM'
            h[0].header['SLICE']='CENTRAL'
            h[0].header['MEANVEL']=float(np.mean([fits.getheader(p,1)['VEL_CORR'] for p in files]))
        print('COADDED',target,ch,star,flush=True)
    (env.BASE/f'audit/flux_assignments_{ch}{suffix}.json').write_text(json.dumps(assignments,indent=2))

if __name__=='__main__':
    fn={'prepare':prepare,'flux':flux_and_coadd}[sys.argv[1]]
    for ch in sys.argv[2:]:fn(ch)
