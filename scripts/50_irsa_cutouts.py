"""IRSA imaging cutouts for a set of targets: everything the archive holds at each position.

Discovers coverage with the IRSA SIA2 service (all collections) plus the ZTF and PTF
reference-image indexes, fetches 90-arcsecond cutouts of the science images (IBE cutout
parameters where the service supports them, whole-file download otherwise), resamples every
cutout onto a common north-up grid through its WCS, and saves a mosaic plus the FITS cutouts
and a coverage table. Collections with degree-scale pixels (IRAS, AKARI) and the thousands of
single NEOWISE frames are counted but not drawn.

Usage: python scripts/50_irsa_cutouts.py P05:P2921 P06:P3024 P10:P1793 P11:P1823 [--size 90]
"""
import argparse, gzip, io, json, re, warnings
from pathlib import Path
import numpy as np, pandas as pd, requests
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import AsinhStretch, PercentileInterval, ImageNormalize
from scipy.ndimage import map_coordinates
import pyvo
warnings.filterwarnings('ignore')
ROOT=Path(__file__).resolve().parents[1]; DEST=ROOT/'observing/sep23/figures'; CUT=DEST/'irsa_cutouts'; CUT.mkdir(parents=True,exist_ok=True)
S=requests.Session(); S.headers['User-Agent']='CLAGN-cutouts/1.0'
COLUMNS=[('ZTF g','ztf','zg'),('ZTF r','ztf','zr'),('ZTF i','ztf','zi'),('2MASS J','twomass_allsky','J'),('2MASS H','twomass_allsky','H'),('2MASS K','twomass_allsky','K'),
         ('IRAC 3.6','spitzer_sha','IRAC1'),('IRAC 4.5','spitzer_sha','IRAC2'),('unWISE W1','wise_unwise','W1'),('unWISE W2','wise_unwise','W2'),('AllWISE W3','wise_allwise','W3'),('AllWISE W4','wise_allwise','W4'),
         ('SPHEREx D1 0.75-1.1µm','spherex_qr3','SPHEREx-D1'),('SPHEREx D4 2.4-3.8µm','spherex_qr3','SPHEREx-D4')]
NOT_DRAWN={'neowiser':'single NEOWISE frames (time series lives in unWISE/NEOWISE coadds and the catalogue)','wise_allsky':'WISE all-sky single frames/atlas (superseded by AllWISE/unWISE)','wise_prelim':'WISE preliminary release','wise_2bandcryo':'WISE 2-band cryo','wise_3bandcryo':'WISE 3-band cryo','wise_prelim_2bandcryo':'WISE prelim 2-band','iras_iris':'IRAS IRIS (arcmin pixels)','iras_issa':'IRAS ISSA','akari_allskymaps':'AKARI far-IR maps','twomass_full':'2MASS full (duplicate of allsky)','twomass_sixdeg':'2MASS 6x calibration scans','spherex_qr2':'SPHEREx QR2 spectral images (drawn from QR3 instead)'}

def fetch(url,ra,dec,size,cut_ok=True):
    """Return (data, wcs) of a cutout; try the IBE cutout syntax first, then the whole file."""
    tries=[]
    if cut_ok and '/ibe/data/' in url: tries.append(f"{url}?center={ra},{dec}&size={size}arcsec&gzip=false")
    tries.append(url)
    for u in tries:
        try:
            r=S.get(u,timeout=180)
            if not r.ok or len(r.content)<2880: continue
            raw=r.content
            if raw[:2]==b'\x1f\x8b': raw=gzip.decompress(raw)
            with fits.open(io.BytesIO(raw),memmap=False) as h:
                hdu=next(x for x in h if x.data is not None and x.data.ndim>=2)
                data=np.asarray(hdu.data,float); hdr=hdu.header
                if data.ndim==3: data=data[0]
                return data,WCS(hdr,naxis=2)
        except Exception as e:
            continue
    return None,None

def northup(data,wcs,ra,dec,size,scale=0.75):
    """Resample onto a north-up, east-left TAN grid centred on the target."""
    n=int(round(size/scale)); w=WCS(naxis=2); w.wcs.ctype=['RA---TAN','DEC--TAN']; w.wcs.crval=[ra,dec]; w.wcs.crpix=[(n+1)/2,(n+1)/2]; w.wcs.cdelt=[-scale/3600,scale/3600]
    yy,xx=np.mgrid[0:n,0:n]; sky=w.pixel_to_world(xx,yy)
    px,py=wcs.world_to_pixel(sky)
    good=np.isfinite(data); filled=np.where(good,data,0.)
    out=map_coordinates(filled,[py,px],order=1,mode='constant',cval=np.nan)
    valid=map_coordinates(good.astype(float),[py,px],order=1,mode='constant',cval=0.)
    out[valid<0.99]=np.nan   # never paint padding or missing pixels with interpolated values
    return out,w

def discover(code,ra,dec):
    svc=pyvo.dal.SIA2Service('https://irsa.ipac.caltech.edu/SIA')
    sia=svc.search(pos=(ra,dec,0.002)).to_table().to_pandas()
    ztf=pd.DataFrame(); ptf=pd.DataFrame()
    for label,url in [('ztf','https://irsa.ipac.caltech.edu/ibe/search/ztf/products/ref'),('ptf','https://irsa.ipac.caltech.edu/ibe/search/ptf/images/level2')]:
        r=S.get(url,params=dict(POS=f'{ra},{dec}',ct='csv'),timeout=90)
        df=pd.read_csv(io.StringIO(r.text)) if r.ok and r.text.strip() and not r.text.startswith('<') else pd.DataFrame()
        if label=='ztf': ztf=df
        else: ptf=df
    return sia,ztf,ptf

def pick(sia,ztf,ptf,coll,band):
    """Ranked candidate URLs for one panel; the caller keeps the first with full coverage."""
    if coll=='ztf':
        s=ztf[ztf.filtercode==band].sort_values('nframes',ascending=False)
        return [f"https://irsa.ipac.caltech.edu/ibe/data/ztf/products/ref/{int(r.field)//1000:03d}/field{int(r.field):06d}/{r.filtercode}/ccd{int(r.ccdid):02d}/q{int(r.qid)}/ztf_{int(r.field):06d}_{r.filtercode}_c{int(r.ccdid):02d}_q{int(r.qid)}_refimg.fits" for _,r in s.iterrows()]
    if coll=='ptf':
        out=[]
        for _,r in ptf.iterrows():
            m=re.match(r'PTF_(d\d+)_(f\d+)_(c\d+)_(u\d+)_(p\d+)_refimg\.fits',r.filename)
            if m: out.append(f"https://irsa.ipac.caltech.edu/ibe/data/ptf/images/level2/{'/'.join(m.groups())}/{r.filename}")
        return out
    s=sia[(sia.obs_collection==coll)&(sia.energy_bandpassname==band)&(sia.dataproduct_subtype=='science')]
    if coll=='wise_unwise':
        m=s[s.access_url.str.contains('img-m')]; s=m if not m.empty else s
    if coll.startswith('spherex'): s=s.sort_values('t_min',ascending=False)
    return list(dict.fromkeys(s.access_url))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('targets',nargs='+',help='CODE:NAME'); ap.add_argument('--size',type=float,default=90.); a=ap.parse_args()
    cat=pd.read_csv(ROOT/'data/reselection_2026-09-20/compact_review_objects.csv').set_index('name')
    rows=[]; coverage=[]; panels={}
    targets=[(c.split(':')[0],c.split(':')[1]) for c in a.targets]
    for code,name in targets:
        ra,dec,z=float(cat.loc[name].ra),float(cat.loc[name].dec),float(cat.loc[name].z)
        sia,ztf,ptf=discover(code,ra,dec)
        counts=sia[sia.dataproduct_subtype=='science'].groupby('obs_collection').size().to_dict()
        for coll,n in counts.items(): coverage.append(dict(code=code,name=name,collection=coll,science_images=n,drawn=coll in {c for _,c,_ in COLUMNS},note=NOT_DRAWN.get(coll,'')))
        coverage.append(dict(code=code,name=name,collection='ztf_ref',science_images=len(ztf),drawn=True,note=', '.join(sorted(ztf.filtercode.unique())) if len(ztf) else 'none'))
        coverage.append(dict(code=code,name=name,collection='ptf_level2',science_images=len(ptf),drawn=False,note=('R-band reference exists but PTF Level 2 images are proprietary at IRSA' if len(ptf) else 'none')))
        for coll in ['euclid_DpdMerBksMosaic','euclid_q2','euclid_ero','spitzer_seip']:
            if coll not in counts: coverage.append(dict(code=code,name=name,collection=coll,science_images=0,drawn=False,note='no coverage at this position'))
        for label,coll,band in COLUMNS:
            urls=pick(sia,ztf,ptf,coll,band)
            if not urls: panels[(code,label)]=None; rows.append(dict(code=code,name=name,panel=label,status='no coverage',url='')); continue
            outfits=CUT/f"{name}_{label.split()[0]}_{label.split()[1] if len(label.split())>1 else ''}.fits".replace('__','_')
            if outfits.exists():
                panels[(code,label)]=np.asarray(fits.getdata(outfits),float); rows.append(dict(code=code,name=name,panel=label,status='ok (cached)',url=fits.getheader(outfits).get('SRCURL',''))); continue
            best=None
            for url in urls[:6]:
                data,wcs=fetch(url,ra,dec,a.size)
                if data is None: continue
                try:
                    img,w=northup(data,wcs,ra,dec,a.size)
                except Exception:
                    continue
                n=img.shape[0]; c=n//2; centre_ok=np.isfinite(img[c-3:c+4,c-3:c+4]).all()
                frac=float(np.isfinite(img).mean())
                if best is None or (centre_ok,frac)>(best[0],best[1]): best=(centre_ok,frac,img,w,url)
                if centre_ok and frac>0.98: break
            if best is None: panels[(code,label)]=None; rows.append(dict(code=code,name=name,panel=label,status='download failed',url=urls[0])); print(code,label,'FAILED',urls[0],flush=True); continue
            centre_ok,frac,img,w,url=best
            if not centre_ok:
                panels[(code,label)]='footprint'; rows.append(dict(code=code,name=name,panel=label,status='in archive but target outside the image footprint',url=url,coverage=frac)); print(code,label,'OUTSIDE FOOTPRINT',f'{frac:.2f}',flush=True); continue
            hdr=w.to_header(); hdr['SRCURL']=url[:68]; hdr['COVERAGE']=frac
            fits.PrimaryHDU(img.astype(np.float32),header=hdr).writeto(outfits,overwrite=True)
            panels[(code,label)]=img; rows.append(dict(code=code,name=name,panel=label,status='ok' if centre_ok else 'target at image edge',url=url,coverage=frac)); print(code,label,'ok' if centre_ok else 'EDGE',f'{frac:.2f}',flush=True)
    pd.DataFrame(rows).to_csv(DEST/'irsa_cutouts_status.csv',index=False); pd.DataFrame(coverage).to_csv(DEST/'irsa_coverage.csv',index=False)
    # Mosaic
    nr,nc=len(targets),len(COLUMNS)
    fig,axes=plt.subplots(nr,nc,figsize=(1.55*nc,1.75*nr+0.6),squeeze=False)
    for i,(code,name) in enumerate(targets):
        z=float(cat.loc[name].z)
        for j,(label,coll,band) in enumerate(COLUMNS):
            ax=axes[i,j]; img=panels.get((code,label)); ax.set_xticks([]); ax.set_yticks([])
            if isinstance(img,str):
                ax.set_facecolor('#f3e9dc'); ax.text(.5,.5,'in archive,\ntarget outside\nmosaic footprint',ha='center',va='center',fontsize=6.5,color='0.35',transform=ax.transAxes)
            elif img is None or not np.isfinite(img).any():
                ax.set_facecolor('#e8e8e8'); ax.text(.5,.5,'no\ncoverage',ha='center',va='center',fontsize=7,color='0.4',transform=ax.transAxes)
            else:
                norm=ImageNormalize(img,interval=PercentileInterval(99.3),stretch=AsinhStretch(0.1))
                ax.imshow(img,origin='lower',cmap='gray_r',norm=norm)
                n=img.shape[0]; c=(n-1)/2; r=5/0.75
                ax.add_patch(plt.Circle((c,c),r,fill=False,color='#e6550d',lw=1))
                if i==nr-1: ax.plot([n*.08,n*.08+30/0.75],[n*.07,n*.07],color='#e6550d',lw=2); ax.text(n*.08,n*.11,'30″',color='#e6550d',fontsize=6)
            if i==0: ax.set_title(label,fontsize=7.5)
            if j==0: ax.set_ylabel(f'{code} · {name}\nz={z:.3f}',fontsize=8)
    fig.suptitle(f'IRSA imaging at the four positions · {a.size:.0f}″ fields, north up east left, circle 5″ · asinh stretch',fontsize=10,y=0.995)
    fig.tight_layout(rect=(0,0,1,0.97)); out=DEST/f"irsa_cutouts_{'_'.join(c for c,_ in targets)}.png"
    fig.savefig(out,dpi=150); print('wrote',out)
    print(pd.DataFrame(rows).pivot(index='code',columns='panel',values='status').reindex(columns=[l for l,_,_ in COLUMNS]).to_string())

if __name__=='__main__': main()
