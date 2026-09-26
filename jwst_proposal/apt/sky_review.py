"""Select and document nearby MRS sky fields from IRSA images and catalogs.

Survey-level screening, not a claim of zero emission at JWST depth. Downloads
are cached with request URLs; failures are never interpreted as blank sky.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import io
import json
import os
from pathlib import Path
import threading

os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import map_coordinates
from PIL import Image

HERE = Path(__file__).resolve().parent
OUT = HERE / 'sky_fields'
CACHE = OUT / 'cache'
LOCAL = threading.local()
IRSA = 'https://irsa.ipac.caltech.edu'


def fetch(path, url, params=None):
    if path.exists():
        return path.read_bytes()
    if not hasattr(LOCAL, 'session'):
        LOCAL.session = requests.Session()
    error = None
    for _ in range(2):
        try:
            r = LOCAL.session.get(url, params=params, timeout=(15, 65))
            r.raise_for_status()
            if len(r.content) < 50:
                raise ValueError('Unexpected empty response')
            if path.suffix == '.fits':
                with fits.open(io.BytesIO(r.content)) as hd:
                    assert hd[0].data is not None
            elif path.suffix == '.tbl':
                Table.read(r.text, format='ascii.ipac')
            elif path.suffix == '.jpg':
                im = Image.open(io.BytesIO(r.content)); im.load()
                assert np.asarray(im).std() > 1
            path.write_bytes(r.content)
            path.with_suffix(path.suffix+'.json').write_text(json.dumps({'url':r.url},indent=2)+'\n')
            return r.content
        except Exception as exc:
            error = exc
    raise RuntimeError(f'{path.name}: {error}')


def coord(row):
    # Use the actual sexagesimal science coordinates entered in APT.
    return SkyCoord(row['ra_hms'], row['dec_dms'], unit=(u.hourangle,u.deg))


def catalog(row, name):
    c=coord(row); path=CACHE/f"{row['id']}_{name}.tbl"
    catalog_name, columns = ('allwise_p3as_psd','designation,ra,dec,w1mpro,w3mpro,w4mpro,cc_flags,ext_flg') if name=='allwise' else ('fp_psc','designation,ra,dec,j_m,k_m')
    data=fetch(path,IRSA+'/cgi-bin/Gator/nph-query',dict(catalog=catalog_name,spatial='Cone',objstr=f'{c.ra.deg} {c.dec.deg}',radius=360 if row['id']=='R12' else 180,radunits='arcsec',outfmt=1,selcols=columns))
    t=Table.read(data.decode(),format='ascii.ipac')
    return t,SkyCoord(np.array(t['ra'])*u.deg,np.array(t['dec'])*u.deg)


def image_at(path):
    with fits.open(path) as h:
        return np.asarray(h[0].data,float),WCS(h[0].header,naxis=2)


def wise_images(row):
    c=coord(row); ident=row['id']; params={'POS':f'{c.ra.deg},{c.dec.deg}','ct':'csv'}
    raw=fetch(CACHE/f'{ident}_wise_index.csv',IRSA+'/ibe/search/wise/allwise/p3am_cdd',params)
    index=list(csv.DictReader(io.StringIO(raw.decode())))
    for band in [1,3,4]:
        path=CACHE/f'{ident}_W{band}.fits'
        if path.exists():continue
        entries=[r for r in index if int(r['band'])==band]
        entries.sort(key=lambda r:(float(r['crval1'])-c.ra.deg)**2*np.cos(c.dec.rad)**2+(float(r['crval2'])-c.dec.deg)**2)
        error=None
        for r in entries:
            tile=r['coadd_id']
            url=f'{IRSA}/ibe/data/wise/allwise/p3am_cdd/{tile[:2]}/{tile[:4]}/{tile}/{tile}-w{band}-int-3.fits'
            try:
                fetch(path,url,dict(center=f'{c.ra.deg},{c.dec.deg}',size='600arcsec' if ident=='R12' else '330arcsec',gzip='false'))
                break
            except Exception as exc:error=exc
        else:raise RuntimeError(f'{ident} W{band}: {error}')


def patch_metric(data,wcs,position):
    x,y=wcs.world_to_pixel(position)
    scale=float(np.mean(proj_plane_pixel_scales(wcs)))*3600
    yy,xx=np.indices(data.shape); dist=np.hypot(xx-x,yy-y)*scale
    core=data[dist<=20]; ring=data[(dist>=30)&(dist<=60)]
    if len(core)<30 or np.isfinite(core).mean()<1 or np.isfinite(ring).mean()<.95:return 999.
    med=np.nanmedian(ring); sigma=1.4826*np.nanmedian(abs(ring-med))
    return float((np.nanpercentile(core,98)-med)/max(sigma,1e-8))


def select(row):
    c=coord(row); ident=row['id']
    wise_images(row)
    wt,wc=catalog(row,'allwise');mt,mc=catalog(row,'2mass')
    images=[image_at(CACHE/f'{ident}_W{b}.fits') for b in [1,3,4]]
    candidates=[]
    for radius in (list(range(60,241,15)) if ident=='R12' else [60,75,90,105,120]):
        for pa in range(0,360,5):
            s=c.directional_offset_by(pa*u.deg,radius*u.arcsec)
            dw=s.separation(wc).arcsec; dm=s.separation(mc).arcsec
            nearw=float(min(dw,default=180)); nearm=float(min(dm,default=180))
            bright=np.array(np.ma.filled(wt['w1mpro'],99),float)<10
            bright_clear=float(min(dw[bright],default=180))
            metrics=[patch_metric(d,w,s) for d,w in images]
            # A 25-arcsec catalog clearance encloses the 20-arcsec inspection
            # region, itself larger than the dithered MRS field of view.
            bright_limit=60 if ident=='R12' else 90
            passed=nearw>=25 and nearm>=25 and bright_clear>=bright_limit and max(metrics)<5
            candidates.append(dict(radius_arcsec=radius,pa_deg=pa,ra_deg=float(s.ra.deg),dec_deg=float(s.dec.deg),
                                   nearest_allwise_arcsec=nearw,nearest_2mass_arcsec=nearm,
                                   nearest_bright_w1_arcsec=bright_clear,wise_patch_z=metrics,passed=bool(passed),
                                   score=min(nearw,nearm,75)-.04*radius-max(0,max(metrics)-2)*3))
    passing=[x for x in candidates if x['passed']]
    if not passing:
        (OUT/f'{ident}_failed_candidates.json').write_text(json.dumps(candidates,indent=2)+'\n')
        raise RuntimeError(f'{ident}: no candidate meets catalog/image screening')
    original=candidates[0]
    best=original if original['passed'] else max(passing,key=lambda x:x['score'])
    result=dict(id=ident,target=row['target'],science_ra_deg=float(c.ra.deg),science_dec_deg=float(c.dec.deg),
                selected=best,original=original,moved=(best is not original),catalog_counts=dict(allwise=len(wt),twomass=len(mt)),
                visual_review='pending')
    (OUT/f'{ident}_candidates.json').write_text(json.dumps(dict(result=result,candidates=candidates),indent=2)+'\n')
    s=SkyCoord(best['ra_deg']*u.deg,best['dec_deg']*u.deg)
    # Inspect 2MASS imaging at the selected sky location, selecting a frame
    # that fully contains a 60-arcsec cutout around the actual coordinate.
    raw=fetch(CACHE/f'{ident}_2mass_index.csv',IRSA+'/ibe/search/twomass/allsky/allsky',dict(POS=f'{s.ra.deg},{s.dec.deg}',ct='csv'))
    entries=[r for r in csv.DictReader(io.StringIO(raw.decode())) if r['filter']=='j']
    entries.sort(key=lambda r:(float(r['crval1'])-s.ra.deg)**2*np.cos(s.dec.rad)**2+(float(r['crval2'])-s.dec.deg)**2)
    for k,r in enumerate(entries):
        path=CACHE/f'{ident}_2massJ_{k}.fits'
        url=f"{IRSA}/ibe/data/twomass/allsky/allsky/{int(r['ordate']):06d}{r['hemisphere']}/s{int(r['scanno']):03d}/image/{r['fname']}"
        try:
            fetch(path,url,dict(center=f'{s.ra.deg},{s.dec.deg}',size='70arcsec',gzip='false'))
            d,w=image_at(path);x,y=w.world_to_pixel(s)
            if min(x,y,d.shape[1]-1-x,d.shape[0]-1-y)<29:continue
            result['twomass_image']=str(path.relative_to(OUT));break
        except Exception:continue
    else:raise RuntimeError(f'{ident}: no complete 2MASS J cutout')
    # Deeper optical image is an additional visual check, not the IR criterion.
    try:
        fetch(CACHE/f'{ident}_sdss.jpg','https://skyserver.sdss.org/dr18/SkyServerWS/ImgCutout/getjpeg',dict(ra=s.ra.deg,dec=s.dec.deg,scale=.25,width=280,height=280))
        result['sdss_image']=f'cache/{ident}_sdss.jpg'
    except Exception as exc:result['sdss_note']=str(exc)
    (OUT/f'{ident}_selection.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def northup(data,wcs,center,size=300,scale=1.):
    nn=round(size/scale);ww=WCS(naxis=2);ww.wcs.ctype=['RA---TAN','DEC--TAN'];ww.wcs.crval=[center.ra.deg,center.dec.deg];ww.wcs.crpix=[(nn+1)/2]*2;ww.wcs.cdelt=[-scale/3600,scale/3600]
    yy,xx=np.indices((nn,nn));px,py=wcs.world_to_pixel(ww.pixel_to_world(xx,yy))
    img=map_coordinates(data,[py,px],order=1,mode='constant',cval=np.nan)
    return img,ww


def plot(results):
    with PdfPages(OUT/'sky_field_review.pdf') as pdf:
        for start in range(0,len(results),4):
            part=results[start:start+4];fig,axes=plt.subplots(len(part),5,figsize=(17,3.4*len(part)),squeeze=False)
            for i,r in enumerate(part):
                ident=r['id'];c=SkyCoord(r['science_ra_deg']*u.deg,r['science_dec_deg']*u.deg);s=SkyCoord(r['selected']['ra_deg']*u.deg,r['selected']['dec_deg']*u.deg);old=c.directional_offset_by(0*u.deg,60*u.arcsec)
                for j,b in enumerate([1,3,4]):
                    size=600 if ident=='R12' else 300
                    d,w=image_at(CACHE/f'{ident}_W{b}.fits');img,ww=northup(d,w,c,size)
                    lo,hi=np.nanpercentile(img,[5,99.6]);ax=axes[i,j]
                    ax.imshow(np.arcsinh(np.clip((img-lo)/max(hi-lo,1e-5),0,None)*8),origin='lower',cmap='gray',vmin=0,vmax=np.arcsinh(8))
                    x,y=ww.world_to_pixel(c);ax.plot(x,y,'+',color='#4cc9ff',ms=10)
                    x,y=ww.world_to_pixel(old);ax.add_patch(plt.Circle((x,y),20,fill=False,color='#ffae42',lw=1.2,ls='--'))
                    x,y=ww.world_to_pixel(s);ax.add_patch(plt.Circle((x,y),20,fill=False,color='#20ef80',lw=1.8))
                    ax.set_title(f'{ident} W{b} | {size} arcsec',fontsize=10)
                d,w=image_at(OUT/r['twomass_image']);img,ww=northup(d,w,s,60,.5)
                lo,hi=np.nanpercentile(img,[2,99.5]);axes[i,3].imshow(img,origin='lower',cmap='gray',vmin=lo,vmax=hi)
                axes[i,3].add_patch(plt.Circle((59.5,59.5),40,fill=False,color='#20ef80',lw=1.5));axes[i,3].set_title(f'{ident} 2MASS J | 60 arcsec',fontsize=10)
                if 'sdss_image' in r:
                    axes[i,4].imshow(Image.open(OUT/r['sdss_image']))
                    axes[i,4].add_patch(plt.Circle((139.5,139.5),80,fill=False,color='#20ef80',lw=1.5))
                else:axes[i,4].text(.5,.5,'No validated SDSS cutout',ha='center',transform=axes[i,4].transAxes)
                axes[i,4].set_title(f"{ident} SDSS | 70 arcsec",fontsize=10)
                sel=r['selected'];axes[i,0].set_ylabel(f"{r['target']}\n{sel['radius_arcsec']} arcsec, PA {sel['pa_deg']} deg\nWISE/2MASS clearance {sel['nearest_allwise_arcsec']:.0f}/{sel['nearest_2mass_arcsec']:.0f} arcsec",fontsize=9)
                for ax in axes[i]:ax.set_xticks([]);ax.set_yticks([])
            fig.suptitle('Sky-field review: green = selected 20-arcsec safety region; orange = old north offset; blue = AGN\nWISE/2MASS panels north up, east left. Catalog screening plus visual inspection; survey resolution limits apply.',fontsize=12)
            fig.tight_layout(rect=(0,0,1,.95));pdf.savefig(fig);fig.savefig(OUT/f'review_{start//4+1}.png',dpi=125);plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plot-only',action='store_true');ap.add_argument('--ids',nargs='*');a=ap.parse_args()
    CACHE.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader((HERE.parent/'inputs/jwst_sample_cycle6.csv').open()))
    if a.ids:rows=[r for r in rows if r['id'] in a.ids]
    failures=[];results=[]
    if not a.plot_only:
        with ThreadPoolExecutor(max_workers=4) as executor:
            tasks={executor.submit(select,r):r['id'] for r in rows}
            for f in as_completed(tasks):
                try:
                    r=f.result();results.append(r);print(r['id'],'moved' if r['moved'] else 'kept','offset',r['selected']['radius_arcsec'],flush=True)
                except Exception as exc:failures.append(str(exc));print('FAILED',exc,flush=True)
        (OUT/'fetch_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    else:
        results=[json.loads((OUT/f"{r['id']}_selection.json").read_text()) for r in rows]
    results.sort(key=lambda r:r['id']);(OUT/'selections.json').write_text(json.dumps(results,indent=2)+'\n')
    if results:plot(results)
    print(f'{len(results)} selected, {len(failures)} failures',flush=True)


if __name__=='__main__':main()
