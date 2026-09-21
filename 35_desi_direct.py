"""Retrieve individual DESI DR1 coadds directly when SPARCL is unavailable.

HTTP range access reads only the target's rows. Preserve native band spectra,
inverse variances, masks and resolution matrices in a standalone FITS file.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,importlib,json,fcntl,time
import numpy as np
import pandas as pd
import hpgeom,aiohttp
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'data/reselection_2026-09-20'
DEST=OUT/'desi_direct';DEST.mkdir(exist_ok=True)
API=importlib.import_module('29_three_night_desi');PARSE=importlib.import_module('03d_fetch_spectra');NATIVE=importlib.import_module('17_fetch_review_spectra').finite


def one(item,attempt=0):
    name,epochs=item;records=[];failed=[]
    for e in epochs.rename(columns={'class':'spectral_class'}).itertuples():
        stem=f'{name}_{e.targetid}_{e.survey}_{e.program}'
        saved=API.CACHE/f'{stem}.v2.json';local=DEST/f'{stem}.json'
        if local.exists():records.extend(json.loads(local.read_text()));continue
        if saved.exists():records.extend(json.loads(saved.read_text()));continue
        value=getattr(e,'healpix',np.nan)
        pix=int(value) if pd.notna(value) else int(hpgeom.angle_to_pixel(64,e.ra,e.dec,nest=True,lonlat=True))
        url=f'https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/healpix/{e.survey}/{e.program}/{pix//100}/{pix}/coadd-{e.survey}-{e.program}-{pix}.fits'
        try:
            path=DEST/f'{stem}.fits'
            if not path.exists():
                with fits.open(url,use_fsspec=True,fsspec_kwargs={'block_size':262144,'cache_type':'bytes','client_kwargs':{'timeout':aiohttp.ClientTimeout(total=75)}}) as hdus:
                    fmap=hdus['FIBERMAP'].data;idx=np.where(fmap['TARGETID']==int(e.targetid))[0]
                    if len(idx)!=1:raise ValueError('Target ID does not identify exactly one coadd row')
                    i=int(idx[0]);r=fmap[i]
                    sep=SkyCoord(float(r['TARGET_RA'])*u.deg,float(r['TARGET_DEC'])*u.deg).separation(SkyCoord(e.ra*u.deg,e.dec*u.deg)).arcsec
                    if sep>2:raise ValueError('Coadd target position does not match the catalog')
                    primary=fits.PrimaryHDU();primary.header['ORIGURL']=url;primary.header['MJD']=float(e.mjd)
                    primary.header['MINMJD']=float(e.min_mjd);primary.header['MAXMJD']=float(e.max_mjd)
                    primary.header['TARGETID']=int(e.targetid);primary.header['SURVEY']=str(e.survey);primary.header['PROGRAM']=str(e.program)
                    copied=[primary,fits.BinTableHDU(fmap[i:i+1],name='FIBERMAP')]
                    for band in ['B','R','Z']:
                        copied.append(fits.ImageHDU(np.asarray(hdus[band+'_WAVELENGTH'].data),name=band+'_WAVELENGTH'))
                        for kind in ['FLUX','IVAR','MASK','RESOLUTION']:
                            ext=band+'_'+kind
                            data=np.asarray(hdus[ext].section[i,...])[None,...]
                            copied.append(fits.ImageHDU(data,name=ext))
                    tmp=path.with_suffix('.tmp');fits.HDUList(copied).writeto(tmp,overwrite=True);tmp.replace(path)
            waves=[];fluxes=[];weights=[]
            with fits.open(path) as h:
                for band in ['B','R','Z']:
                    wave=np.asarray(h[band+'_WAVELENGTH'].data,float);flux=np.asarray(h[band+'_FLUX'].data[0],float)
                    ivar=np.asarray(h[band+'_IVAR'].data[0],float);mask=np.asarray(h[band+'_MASK'].data[0])
                    waves.append(wave);fluxes.append(flux);weights.append(np.where((mask==0)&np.isfinite(flux),ivar,0))
            wave=np.concatenate(waves);flux=np.concatenate(fluxes);weight=np.concatenate(weights)
            idx=np.searchsorted(PARSE.GRID,wave)-1;rebinned=np.full(len(PARSE.GRID),np.nan)
            for i in np.unique(idx[(idx>=0)&(idx<len(rebinned))]):
                valid=(idx==i)&(weight>0)
                if valid.any():rebinned[i]=np.average(flux[valid],weights=weight[valid])
            record=NATIVE(dict(wave=PARSE.GRID.tolist(),flux=rebinned.tolist(),source='DESI',coadd=True,proprietary=False,
                mjd=float(e.mjd),min_mjd=float(e.min_mjd),max_mjd=float(e.max_mjd),coadd_numnight=int(e.coadd_numnight),survey=str(e.survey),program=str(e.program),
                url=url,archive_file=url,date_verified=True,metadata_quality_ok=bool(e.zwarning==0),lines={},ew={},
                meta={'class':str(e.spectral_class),'z':float(e.z),'specid':str(e.targetid),'zwarning':int(e.zwarning)}))
            local.write_text(json.dumps([record],separators=(',',':'),allow_nan=False));records.append(record)
        except Exception as exc:failed.append(dict(targetid=str(e.targetid),survey=str(e.survey),program=str(e.program),reason=type(exc).__name__,detail=str(exc)[:240]))
    if records:
        path=ROOT/'data/spectra_dl'/f'{name}.json'
        with path.with_suffix('.lock').open('w') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            old=json.loads(path.read_text()) if path.exists() else []
            ids={r['meta']['specid'] for r in records}
            # Replace ambiguous old DESI date assignments only after every
            # matched coadd for this target has been verified and cached.
            if not failed:old=[r for r in old if not (r.get('source')=='DESI' and str(r.get('meta',{}).get('specid')) in ids)]
            urls={r.get('url') for r in records};old=[r for r in old if r.get('url') not in urls]
            tmp=path.with_suffix('.desi.tmp');tmp.write_text(json.dumps(NATIVE(old+records),separators=(',',':'),allow_nan=False));tmp.replace(path)
    result=dict(name=name,records=len(records),failed=failed);(DEST/f'{name}_status.json').write_text(json.dumps(result,indent=2))
    if failed and attempt<2:
        time.sleep(3*(attempt+1))
        return one(item,attempt+1)
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--targets');ap.add_argument('--inventory-only',action='store_true');args=ap.parse_args()
    targets=pd.read_csv(args.targets or OUT/'compact_review_objects.csv');frames=[]
    with ThreadPoolExecutor(max_workers=3) as ex:
        futures=[ex.submit(API.query,targets.iloc[i:i+20]) for i in range(0,len(targets),20)]
        for f in as_completed(futures):frames.append(f.result())
    good=[d for d in frames if len(d)]
    if not good:print('No catalog matches; inspect query manifests for failures.');return
    epochs=pd.concat(good,ignore_index=True).drop_duplicates(['name','targetid','survey','program'])
    dest=ROOT/'data/spectra_epochs_three_night_desi.csv'
    if dest.exists():
        old=pd.read_csv(dest,dtype={'targetid':str});merged=pd.concat([old[~old.name.isin(targets.name)],epochs],ignore_index=True)
    else:merged=epochs
    merged.to_csv(dest,index=False)
    if args.inventory_only:return
    results=[]
    with ThreadPoolExecutor(max_workers=4) as ex:
        futures={ex.submit(one,item):item[0] for item in epochs.groupby('name')}
        for i,f in enumerate(as_completed(futures),1):
            try:results.append(f.result())
            except Exception as exc:results.append(dict(name=futures[f],status='worker failed',reason=type(exc).__name__))
            if i%10==0 or i==len(futures):
                (OUT/'desi_direct_status.json').write_text(json.dumps(results,indent=2));print('Direct DESI',i,'/',len(futures),flush=True)


if __name__=='__main__':main()
