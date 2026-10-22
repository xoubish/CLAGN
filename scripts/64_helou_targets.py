"""Optical aperture photometry and October feasibility for the five Helou targets.

Consumes cached public SIMBAD, IRAS FSC, ZTF and Legacy Survey responses in
data/additional_targets. Does not treat total-galaxy or mid-IR flux as slit flux.
"""
from pathlib import Path
import base64
import importlib
import json
import sys

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits, ascii
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'data/additional_targets'
REVIEW = BASE/'helou_review'
H = importlib.import_module('63_5muses327')


def photometry(row):
    """Subpixel circular apertures, no aperture/dereddening correction.

    Coadds are sky-subtracted nanomaggies/pixel. These small apertures include
    central host light; they are NOT an isolated AGN or nuclear point source.
    Resampling, seeing and centering dominate precision; quote 0.1 mag only.
    """
    path = BASE/row['name']/'legacy.fits'
    with fits.open(path) as hdus:
        data, header = hdus[0].data.copy(), hdus[0].header.copy()
    assert header['BANDS'].strip() == 'grz' and np.isfinite(data).all()
    w = WCS(header, naxis=2)
    cx, cy = w.world_to_pixel(SkyCoord(row['ra']*u.deg, row['dec']*u.deg))
    scale = np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))*3600
    y, x = np.indices(data.shape[1:])
    result = dict(survey='Legacy Surveys '+header['VERSION'].strip(), apertures={})
    for diameter in [1.5, 3.0]:
        mask = np.zeros(x.shape)
        for dx in (np.arange(20)+.5)/20-.5:
            for dy in (np.arange(20)+.5)/20-.5:
                mask += ((x+dx-cx)**2+(y+dy-cy)**2)*scale**2 < (diameter/2)**2
        mask /= 400
        flux = (data*mask).sum(axis=(1,2))
        result['apertures'][str(diameter)] = {b:float(22.5-2.5*np.log10(f)) if f>0 else None for b,f in zip('grz',flux)}
    # A 1.5 x 3 arcsec central extraction is a better match to the planned
    # slit than a circular aperture. Use the faintest orientation, not the
    # brightest knot. Keep the circular values for comparison.
    fluxes=[]
    for angle in np.arange(0,180,5):
        pa=np.deg2rad(angle); mask=np.zeros(x.shape)
        for dx in (np.arange(8)+.5)/8-.5:
            for dy in (np.arange(8)+.5)/8-.5:
                xx=(x+dx-cx)*scale; yy=(y+dy-cy)*scale
                mask += (abs(xx*np.cos(pa)+yy*np.sin(pa))<=.75)&(abs(-xx*np.sin(pa)+yy*np.cos(pa))<=1.5)
        fluxes.append((data*mask/64).sum(axis=(1,2)))
    result['slit_1p5_by_3_min_pa']={b:float(22.5-2.5*np.log10(f)) if f>0 else None for b,f in zip('grz',np.min(fluxes,axis=0))}
    result['method'] = 'Sky-subtracted coadd; circular diameter in arcsec; includes host; no aperture or extinction correction; systematic uncertainty not quantified.'
    return result


def etc(mag, wave, airmass, seeing=1.3, sky=18.5):
    sys.path.insert(0,str(ROOT/'data/reselection_2026-09-20/ngps_etc'))
    E = importlib.import_module('ETC_main')
    C = importlib.import_module('ETC_config')
    # Compatibility with SciPy's relocated low-level spline wrapper. The
    # official ETC's interpolation and throughput model remain unchanged.
    from scipy.interpolate import dfitpack, _fitpack
    if not hasattr(dfitpack,'bispeu'):
        dfitpack.bispeu = lambda tx,ty,c,kx,ky,x,y: _fitpack.bispeu(tx,ty,c,kx,ky,np.atleast_1d(x),np.atleast_1d(y))
    for ch,rn,scale in zip(C.channels,[2.8,7.8,3.7,4.6],[.193,.189,.186,.186]):
        C.readnoise[ch]=rn*u.count/u.pix
        C.platescale[ch]=scale*u.arcsec/u.pix
    ch=next(c for c in ['G','R','I','U'] if C.channelRange[c][0].to_value(u.nm)<wave-4 and C.channelRange[c][1].to_value(u.nm)>wave+4)
    args=E.parser.parse_args([ch,str(wave-4),str(wave+4),'EXPTIME','300',
        '-slit','SET','1.5','-binspect','3','-binspat','2','-seeing',str(seeing),'500',
        '-airmass',str(airmass),'-skymag',str(sky),'-mag',str(mag),'-magsystem','AB','-magfilter','match','-noslicer'])
    E.check_inputs_add_units(args)
    r=E.main(args,quiet=True)
    return float(r['SNR'].value*np.sqrt(2)/np.sqrt((C.dLambda[ch]*3).to_value(u.AA)))


def main():
    rows=json.loads((REVIEW/'targets.json').read_text())
    output=[]
    for row in rows:
        name=row['name']; dest=BASE/name
        photo=photometry(row)
        H.save_json(dest/'aperture_photometry.json',photo)
        ap=photo['slit_1p5_by_3_min_pa']
        nights=H.geometry(row['ra'],row['dec'],dest,26 if name=='HIRAS-F23290-1252' else 14)
        # Faint central-aperture flux as a point source, with an additional
        # one-mag-loss stress case. A scenario, not an extended-source ETC.
        deeper=name=='HIRAS-F23290-1252'
        xlimit=1.5 if deeper else 1.8
        nexp=4 if deeper else 2
        if deeper:
            for n in nights:
                allowed=n['tier_ranges']['preferred']
                n['tier_ranges']['extended']=allowed
                n['tier_ranges']['fallback']=allowed
                n['ranges']=allowed
                n['longest_minutes']=max((r['minutes'] for r in allowed),default=0)
        snr=etc(ap['g'],486.27*(1+row['z']),xlimit)*np.sqrt(nexp/2) if ap['g'] is not None else None
        stress=etc(ap['g']+1,486.27*(1+row['z']),xlimit,1.8,18.)*np.sqrt(nexp/2) if ap['g'] is not None else None
        feasible=[n['night'] for n in nights if n['ranges']]
        note=(f"Central 1.5 x 3 arcsec extraction g/r = {ap['g']:.1f}/{ap['r']:.1f} AB ({photo['survey']}); faintest sampled orientation, includes host light. "
              f"{nexp}x300 s H-beta continuum ETC scenarios: {snr:.1f} per Angstrom at X={xlimit}, sky V=18.5, zenith seeing 1.3 arcsec; "
              f"{stress:.1f} with one extra magnitude of loss, sky V=18 and seeing 1.8 arcsec. "
              'Broadband-normalized point-source approximation, not a line-detection guarantee. Center on the nucleus; inspect first exposure and adjust for host dilution and counts.') if snr is not None else 'Legacy DR9 cutout has no g/r coverage at this position; nuclear g/r aperture brightness and ETC not established. ZTF object catalogue is cached separately.'
        identity=''
        if name=='HUGC-03097':
            identity='Supplied z=0.01917 conflicts with SIMBAD z=0.012014 and LASr z=0.01197; use SIMBAD here. No Moon-safe window on either night.'
        elif name=='HIRAS-F23290-1252':
            identity='Abbreviated name resolved by IRAS FSC prefix search (11 matches): F23290-1252 has 25um=0.3778 Jy and SIMBAD z=0.0212045, matching the supplied table. Optical nucleus coordinates used instead of the IRAS centroid.'
        c=SkyCoord(row['ra']*u.deg,row['dec']*u.deg)
        alias=name[1:].replace('NGC-','NGC ').replace('IC-','IC ').replace('UGC-','UGC ').replace('IRAS-','IRAS ')
        source='https://simbad.cds.unistra.fr/simbad/sim-id?Ident='+alias.replace(' ','%20')
        target=dict(name=name,jname=alias,aliases=[alias],ra=row['ra'],dec=row['dec'],z=row['z'],
            sex=c.to_string('hmsdms',sep=':',precision=2),rmag=ap['r'],r_source=photo['survey']+' central 1.5 x 3 arcsec extraction',
            ux=None,uy=None,region='not projected',on_fraction=None,off_fraction=None,
            status='unmatched',known=[],n_spec=0,n_spec_available_dates=0,epochs=[],nights=nights,
            ztf={},ztf_status={'status':'object catalogue queried; no time-series variability assessment'},wise={},neo=[],
            pool_role='requested',prepared_nights=feasible,field_status='review',
            neighbour_screen=dict(status='review',flags=['Resolved galaxy: confirm nuclear acquisition and extraction']),
            spec=None,lines=[],exposure_plans={},
            science=dict(reference_private=False,state_classification='Radio-deficient proto-starburst / obscured-AGN candidate; optical core classification not established',
                science_question='Measure the central optical continuum and emission/absorption lines of the radio-deficient galaxy; deeply buried activity may be invisible in the optical.',
                field_notes=note+' '+identity),
            archive_review=dict(summary=note,note=identity+' The source description reports no visible variability; this is not independently established by the ZTF object-catalogue scatter. No W1 manifold projection has been made for this galaxy.',
                sources=[dict(label='SIMBAD identity',url=source),dict(label='Legacy Survey field',url=f"https://www.legacysurvey.org/viewer?ra={row['ra']}&dec={row['dec']}&layer=ls-dr9&zoom=15")]))
        if ap['r'] is None:
            ztf=ascii.read(dest/'ztf_objects.tbl',format='ipac').to_pandas()
            matches=ztf[(ztf.filtercode=='zr')&(ztf.dist<1)&(ztf.ngoodobs>20)].sort_values('ngoodobs',ascending=False)
            if len(matches):
                target['rmag']=float(matches.iloc[0].medianmag)
                target['r_source']='ZTF DR24 central source PSF median; host-dependent, not matched slit photometry'
                target['archive_review']['summary']+=f" Closest well-sampled central ZTF r source has median {target['rmag']:.2f}; PSF fitting to a resolved galaxy is seeing-dependent."
        wide=dest/'legacy.wide.jpg'; small=dest/'legacy.jpg'
        if wide.exists():
            from PIL import Image
            with Image.open(wide) as im: im.verify()
            meta=json.loads(wide.with_suffix('.jpg.json').read_text())
            survey='Legacy Surveys'  # exact DR9/DR10 products are in sidecar provenance
            target.update(image_status='available',cut=dict(source=survey,field_arcsec=240,image='data:image/jpeg;base64,'+base64.b64encode(wide.read_bytes()).decode()))
        if small.exists():
            target['cut40']='data:image/jpeg;base64,'+base64.b64encode(small.read_bytes()).decode()
        H.save_json(dest/'target_payload.json',target)
        output.append(dict(**row,central_g=ap['g'],central_r=ap['r'],snr_scenario=snr,snr_stress=stress,exposures=nexp,airmass_limit=xlimit,
            viable_nights=','.join(feasible),identity_note=identity))
        print(name,ap,'SNR',snr,stress,feasible,flush=True)
    pd.DataFrame(output).to_csv(REVIEW/'feasibility.csv',index=False)
    hs=json.loads((BASE/'H5MUSES-327/summary.json').read_text())
    result=dict(g=hs['ztf']['g']['median_mag'],r=hs['ztf']['r']['median_mag'])
    result['snr486_scenario']=etc(result['g'],486.27,1.5)
    result['snr486_stress']=etc(result['g']+1,486.27,1.8,1.8,18.)
    H.save_json(BASE/'H5MUSES-327/optical_etc.json',result)


if __name__=='__main__':main()
