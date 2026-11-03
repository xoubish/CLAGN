"""Native spectra for the Palomar paper pilot; never modify delivered products."""
import json
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import analysis as A
from paper_protocol import OUT


def load():
    selected=pd.read_csv(OUT/'pilot_epochs.csv')
    objects,visits,_,_=A.census()
    objects=objects[objects.name.isin(selected.name)]
    spectra,_=A.load_spectra(objects,visits)
    spectra=[s for s in spectra if s['epoch_id'] in set(selected.epoch_id)]
    audit=[]
    for s in spectra:
        status='original native spectrum'
        if s['source']=='DESI':
            obj=objects.set_index('name').loc[s['name']]
            matches=[]
            for p in (A.ROOT/'data/reselection_2026-09-20/desi_direct').glob(s['name']+'_*.fits'):
                with fits.open(p,memmap=False) as h:
                    hdr=h[0].header;fm=h['FIBERMAP'].data
                    sep=SkyCoord(obj.ra*u.deg,obj.dec*u.deg).separation(SkyCoord(float(fm['TARGET_RA'][0])*u.deg,float(fm['TARGET_DEC'][0])*u.deg)).arcsec
                    if sep<2 and abs(hdr.get('MINMJD',-1)-s['min_mjd'])<.001 and abs(hdr.get('MAXMJD',-1)-s['max_mjd'])<.001:
                        matches.append(p)
            if len(matches)==1:
                p=matches[0];arrays=[]
                with fits.open(A.track(p),memmap=False) as h:
                    for band,lo,hi in [('B',0,5750),('R',5750,7550),('Z',7550,np.inf)]:
                        w=h[band+'_WAVELENGTH'].data.astype(float)
                        f=h[band+'_FLUX'].data[0].astype(float)
                        iv=h[band+'_IVAR'].data[0].astype(float)
                        good=(iv>0)&(h[band+'_MASK'].data[0]==0)
                        err=np.where(good,1/np.sqrt(np.maximum(iv,1e-300)),np.nan)
                        f[~good]=np.nan
                        # Resolution matrix diagonal second moment is approximate,
                        # not a full forward convolution through the DESI matrix.
                        r=np.maximum(h[band+'_RESOLUTION'].data[0].astype(float),0)
                        off=np.arange(r.shape[0])-r.shape[0]//2
                        norm=np.maximum(r.sum(axis=0),1e-20)
                        mean=(r*off[:,None]).sum(axis=0)/norm
                        sig=np.sqrt(np.maximum((r*off[:,None]**2).sum(axis=0)/norm-mean**2,0))*np.median(np.diff(w))
                        m=(w>=lo)&(w<hi)
                        arrays.append((w[m],f[m],err[m],sig[m],np.full(m.sum(),band)))
                s.update(wave=np.concatenate([a[0] for a in arrays]),flux=np.concatenate([a[1] for a in arrays]),
                    error=np.concatenate([a[2] for a in arrays]),lsf_sigma=np.concatenate([a[3] for a in arrays]),
                    channel=np.concatenate([a[4] for a in arrays]),path=str(p.relative_to(A.ROOT)),error_kind='native diagonal',
                    alternate=None,uncorrected=None)
                s['quality']=np.zeros(len(s['wave']),int)
                status='DESI native recovered: exact interval and coordinate match'
                s['resolution_status']='DESI resolution-matrix second-moment approximation'
            else:
                status=f'DESI native not uniquely identified: {len(matches)} matches'
        if s['source']=='SDSS' and s['path'].endswith('.fits'):
            with fits.open(A.ROOT/s['path'],memmap=False) as h:
                d=h[1].data
                if 'wdisp' in [c.lower() for c in d.names]:
                    s['lsf_sigma']=np.asarray(d['wdisp'],float)*np.median(np.diff(d['loglam']))*np.log(10)*s['wave']
                    s['resolution_status']='SDSS native wdisp'
        if 'lsf_sigma' not in s:
            resolution=1650 if s['source']=='NGPS' else 2000
            s['lsf_sigma']=s['wave']/resolution/2.35482
            s['resolution_status']=f'assumed constant R={resolution}; needs empirical verification'
        audit.append({k:s[k] for k in ['name','epoch_id','source','path','error_kind','resolution_status']}|dict(status=status))
    pd.DataFrame(audit).to_csv(OUT/'native_data_audit.csv',index=False)
    (OUT/'spectrum_inputs.json').write_text(json.dumps(A.SOURCES,indent=2)+'\n')
    return spectra


if __name__=='__main__':
    spectra=load()
    print(pd.read_csv(OUT/'native_data_audit.csv').groupby(['source','status']).size().to_string())
