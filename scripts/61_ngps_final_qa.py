"""Extraction diagnostics, exposure consistency, and raw-file preservation check."""
import importlib,json,hashlib,csv
import numpy as np
from astropy.io import fits
from astropy.table import Table
env=importlib.import_module('54_reduce_ngps_sep23')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main():
    out=env.BASE/'products';(out/'plots').mkdir(parents=True,exist_ok=True)
    (out/'calibration').mkdir(exist_ok=True)
    (out/'plots/extraction_checks').mkdir(exist_ok=True)
    rows=[];targets={}
    for ch in 'ugri':
        for p in sorted((env.NIGHT/'central_exposures'/ch).glob('*fits')):
            with fits.open(p) as h:
                d=h[1].data;head=h[0].header;name=head['TARGET'].upper();w=d['OPT_WAVE']
                lo,hi={'u':(3600,4000),'g':(4800,5200),'r':(6300,6700),'i':(8400,8800)}[ch]
                q=(w>lo)&(w<hi)&d['OPT_MASK']&d['BOX_MASK']
                rate=np.median(d['OPT_COUNTS'][q])/head['EXPTIME']
                ratio=np.median(d['OPT_COUNTS'][q])/np.median(d['BOX_COUNTS'][q])
                spat=float(h[1].header['SPAT_PIXPOS'])
                rec=dict(NAME=name,CHANNEL=ch.upper(),FILE=p.name,COUNT_RATE=float(rate),
                    OPT_OVER_BOX=float(ratio),SPAT=spat)
                raw2=env.pypeit_file(ch).parent/'Science'/p.name.replace('spec1d_','spec2d_')
                with fits.open(raw2) as hh:
                    sci=hh['DET01-SCIIMG'].data;sky=hh['DET01-SKYMODEL'].data
                    obj=hh['DET01-OBJMODEL'].data;iv=hh['DET01-IVARMODEL'].data
                    valid=(hh['DET01-BPMMASK'].data==0)&(iv>0)
                    xx=np.arange(sci.shape[1])[None,:];ss=(abs(xx-spat)>15)&(abs(xx-spat)<40)
                    z=((sci-sky-obj)*np.sqrt(iv))[valid&ss]
                    rec.update(SKY_RESID_MEDIAN=float(np.median(z)),
                        SKY_RESID_ROBUST_SIGMA=float(1.4826*np.median(abs(z-np.median(z)))),
                        SKY_RESID_FRAC_GT5=float(np.mean(abs(z)>5)))
                    if name not in ['P330E','BD284211'] and ch not in targets.setdefault(name,{}):
                        sl=slice(int(spat)-15,int(spat)+16)
                        signal=((sci-sky)*np.sqrt(iv))[:,sl].T
                        targets[name][ch]=(signal,raw2.name)
                rows.append(rec)
    Table(rows=rows).write(out/'calibration/extraction_qa.csv',overwrite=True)
    with PdfPages(out/'plots/extraction_checks.pdf') as pdf:
        for target,chans in sorted(targets.items()):
            fig,axes=plt.subplots(4,1,figsize=(12,7))
            for ch,ax in zip('ugri',axes):
                im,name=chans[ch];ax.imshow(im,origin='lower',aspect='auto',vmin=-2,vmax=10,cmap='gray_r')
                ax.set_ylabel(ch.upper());ax.set_yticks([5,15,25],[-10,0,10])
                ax.set_title(name,fontsize=7)
            axes[-1].set_xlabel('Spectral pixel (blue → red); spatial pixels relative to trace')
            fig.suptitle(target+' — sky-subtracted first exposure, scaled by pixel noise')
            fig.tight_layout();pdf.savefig(fig);fig.savefig(out/f'plots/extraction_checks/{target}.png',dpi=110);plt.close(fig)
    manifest=env.BASE/'audit/raw_sha256_manifest.csv';fail=[];count=0
    with manifest.open() as f:
        for r in csv.DictReader(f):
            p=env.ROOT/'sep23_data/shemmati'/r['path'];d=hashlib.sha256()
            with p.open('rb') as fp:
                for chunk in iter(lambda:fp.read(1024*1024),b''):d.update(chunk)
            if d.hexdigest()!=r['sha256']:fail.append(r['path'])
            count+=1
    report=dict(files_checked=count,changed_raw_files=fail)
    (env.BASE/'audit/final_raw_integrity.json').write_text(json.dumps(report,indent=2))
    assert not fail,fail
    print(report)

if __name__=='__main__':main()
