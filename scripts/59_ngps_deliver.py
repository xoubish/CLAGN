"""Export calibrated NGPS spectra, empirical tellurics, plots, and explicit quality flags."""
import importlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy import units as units
from astropy.table import Table, vstack
from scipy.ndimage import median_filter
env=importlib.import_module('54_reduce_ngps_sep23')
from ngps_pipeline.finalize import empirical_T_per_channel,TELLURIC_WINDOWS,load_1d
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

OUT=env.BASE/'products'
COADDS=env.NIGHT/'coadded_channels'
TELLURIC_COADDS=env.NIGHT/'coadded_channels'
ASSIGNMENTS_SUFFIX=''
COLORS=dict(u='#7851a9',g='#16843b',r='#d94c36',i='#8a451e')
WINDOWS=dict(u=(3500,4200),g=(4500,5500),r=(6000,6800),i=(8000,8800))
LIMITS=dict(u=(3241,4340),g=(4340,5790),r=(5790,7740),i=(7740,10234))
UNIT='1e-17 erg s-1 cm-2 Angstrom-1'

def read(p):
    with fits.open(p) as h:
        d={k:np.asarray(h[1].data[k]).ravel().copy() for k in h[1].data.names}
        keep=np.flatnonzero(np.isfinite(d['wave'])&(d['wave']>100))
        order=keep[np.argsort(d['wave'][keep])]
        return h[0].header.copy(),{k:v[order] for k,v in d.items()}

def responses(ch,w):
    result={}
    for star in ['P330E','BD284211']:
        arr=[]
        for p in sorted((env.BASE/f'audit/standards/{ch}').glob(f'sens*{star}*fits')):
            with fits.open(p) as h:
                d=h['SENS'].data[0]
                arr.append(np.interp(w,d['SENS_WAVE'],d['SENS_ZEROPOINT_FIT']))
        assert len(arr)==2
        result[star]=np.array(arr)
    return result

def build_transmission(ch):
    p=TELLURIC_COADDS/f'BD284211_{ch.upper()}.fits'
    h,_=read(p);d=load_1d(p)
    sens=sorted((env.BASE/f'audit/standards/{ch}').glob('sens*BD284211*fits'))[0]
    w,t=empirical_T_per_channel(p,sens,ch)
    # Telluric lines are terrestrial: remove the standard's heliocentric shift
    # before evaluating on each science spectrum's observed-frame grid.
    observed=w/h['MEANVEL']
    frac=np.divide(d['sigma'],np.abs(d['flux']),out=np.full(len(w),np.inf),where=abs(d['flux'])>0)
    terr=np.where(t<1,t*frac,0)
    good=d['mask'].astype(bool)&np.isfinite(t)&np.isfinite(terr)
    Table(dict(WAVE_VAC_OBS_A=observed,TRANSMISSION=t,STAT_ERR=terr,MASK=good)).write(
        OUT/f'calibration/telluric_{ch.upper()}.csv',overwrite=True)
    return observed,t,terr,good

def main():
    for directory in ['spectra','channels','plots','calibration']:(OUT/directory).mkdir(parents=True,exist_ok=True)
    allpaths=sorted(COADDS.glob('*_G.fits'))
    targets=[p.name.removesuffix('_G.fits') for p in allpaths if not p.name.startswith(('P330E_','BD284211_'))]
    assert len(targets)==18
    trans={ch:build_transmission(ch) for ch in 'ugri'}
    standards={ch:read(TELLURIC_COADDS/f'BD284211_{ch.upper()}.fits')[0] for ch in 'ugri'}
    summaries=[];overlap=[];calibration=[]
    fig,axes=plt.subplots(2,2,figsize=(11,7))
    for ch,ax in zip('ugri',axes.flat):
        grid=np.linspace(*WINDOWS[ch],300);z=responses(ch,grid)
        ratio=10**(.4*(z['BD284211'].mean(axis=0)-z['P330E'].mean(axis=0)))
        repeat={s:float(np.median(abs(10**(.4*(a[0]-a[1]))-1))) for s,a in z.items()}
        calibration.append(dict(CHANNEL=ch.upper(),P330E_OVER_BD_FLUX=float(np.median(ratio)),
            P330E_REPEAT_FRAC=repeat['P330E'],BD_REPEAT_FRAC=repeat['BD284211']))
        ax.plot(grid,ratio);ax.axhline(1,c='gray',ls=':');ax.set_title(ch.upper())
        ax.set(xlabel='Vacuum wavelength (Å)',ylabel='P330E / BD-calibrated flux')
    fig.tight_layout();fig.savefig(OUT/'calibration/standard_comparison.png',dpi=160);plt.close(fig)
    Table(rows=calibration).write(OUT/'calibration/standard_comparison.csv',overwrite=True)
    with PdfPages(OUT/'plots/all_targets.pdf') as pdf:
        for target in targets:
            channels={};snrs={};header=None
            fig,axes=plt.subplots(2,1,figsize=(13,7),gridspec_kw={'height_ratios':[3,1]},sharex=True)
            for ch in 'ugri':
                h,d=read(COADDS/f'{target}_{ch.upper()}.fits');header=h
                w=d['wave'];f=d['flux'];err=d['sigma'];valid=d['mask'].astype(bool)&np.isfinite(f)&(d['ivar']>0)
                err[~valid]=np.inf
                wt,t,et,mt=trans[ch];wobs=w/h['MEANVEL']
                ratio=h['AIRMASS']/standards[ch]['AIRMASS']
                t0=np.interp(wobs,wt,t,left=1,right=1);te0=np.interp(wobs,wt,et,left=0,right=0)
                tt=t0**ratio;te=ratio*t0**(ratio-1)*te0
                tellvalid=np.interp(wobs,wt,mt.astype(float),left=1,right=1)>.99
                fc=f/tt;ec=np.sqrt((err/tt)**2+(f*te/tt**2)**2)
                quality=np.zeros(len(w),dtype=np.uint16)
                quality[~valid]|=1
                if ch=='u':quality[w<3400]|=2
                twindow=np.zeros(len(w),dtype=bool)
                for lo,hi in TELLURIC_WINDOWS[ch]:twindow|=(wobs>=lo)&(wobs<=hi)
                quality[twindow]|=4
                quality[(tt<.3)|(~tellvalid&twindow)]|=8
                if ch=='i':quality[(wobs<7717.6501)|(wobs>9975.2001)]|=16
                # The alternate standard measures a calibration choice, not a 1-sigma error.
                zp=responses(ch,w);main=h['FLUXSTD'];other='BD284211' if main=='P330E' else 'P330E'
                assignments=json.loads((env.BASE/f'audit/flux_assignments_{ch}{ASSIGNMENTS_SUFFIX}.json').read_text())
                used=next(a['sensfile'] for a in assignments if a['target']==target)
                with fits.open(used) as sh:
                    sd=sh['SENS'].data[0];zmain=np.interp(w,sd['SENS_WAVE'],sd['SENS_ZEROPOINT_FIT'])
                altfactor=10**(-.4*(zp[other].mean(axis=0)-zmain))
                mask=valid&((quality&10)==0)
                tab=Table(dict(WAVE_VAC_HELIO_A=w,FLUX=fc,FLUX_ERR=ec,FLUX_UNCORR=f,
                    ERR_UNCORR=err,TELLURIC=tt,TELLURIC_ERR=te,FLUX_ALTSTD=fc*altfactor,
                    MASK=mask,QUALITY=quality,CHANNEL=np.full(len(w),ch.upper())))
                tab.meta=dict(TARGET=target,FLUXSTD=main,ALTSTD=other,BUNIT=UNIT,WAVEREF='HELIOCENTRIC VACUUM',
                    NEXP=h['NEXP'],TEXPTIME=h['TEXPTIME'],SLICE='CENTRAL',RA=float(h['RA']),DEC=float(h['DEC']),
                    MJD=float(h['MJD']),AIRMASS=float(h['AIRMASS']),MEANVEL=float(h['MEANVEL']),
                    SLITASEC=1.5,EXTRACT='OPT',BINSPAT=2,BINSPEC=3,TELLSTD='BD284211')
                tab['WAVE_VAC_HELIO_A'].unit=units.AA
                for column in ['FLUX','FLUX_ERR','FLUX_UNCORR','ERR_UNCORR','FLUX_ALTSTD']:
                    tab[column].unit=units.Unit(1e-17*units.erg/units.s/units.cm**2/units.AA)
                tab.write(OUT/f'channels/{target}_{ch.upper()}.fits',overwrite=True)
                channels[ch]=tab
                clean=valid&(w>WINDOWS[ch][0])&(w<WINDOWS[ch][1])&~twindow
                snrs[ch]=float(np.median(f[clean]/err[clean]))
                q=(w>LIMITS[ch][0])&(w<LIMITS[ch][1])&mask
                axes[0].plot(w[q],fc[q],color=COLORS[ch],lw=.55,label=f'{ch.upper()} (S/N {snrs[ch]:.1f})')
                axes[0].fill_between(w[q],fc[q]-ec[q],fc[q]+ec[q],color=COLORS[ch],alpha=.13,lw=0)
                axes[1].plot(w[q],fc[q]/ec[q],color=COLORS[ch],lw=.6)
            for blue,red,lo,hi in [('u','g',4270,4400),('g','r',5680,5900),('r','i',7740,7880)]:
                a=channels[blue];b=channels[red];wb=b['WAVE_VAC_HELIO_A'];q=b['MASK']&(wb>lo)&(wb<hi)
                fa=np.interp(wb[q],a['WAVE_VAC_HELIO_A'],a['FLUX_UNCORR']);fb=np.asarray(b['FLUX_UNCORR'][q]);
                ok=np.isfinite(fa)&np.isfinite(fb)&(fb>0)&(fa>0)
                overlap.append(dict(TARGET=target,CHANNELS=blue.upper()+'/'+red.upper(),
                    BLUE_OVER_RED=float(np.median(fa[ok])/np.median(fb[ok])) if ok.sum()>10 else np.nan,NPIX=int(ok.sum())))
            segments=[]
            for ch,tab in channels.items():
                sel=(tab['WAVE_VAC_HELIO_A']>=LIMITS[ch][0])&(tab['WAVE_VAC_HELIO_A']<LIMITS[ch][1])
                segments.append(tab[sel])
            merged=vstack(segments,metadata_conflicts='silent');merged.sort('WAVE_VAC_HELIO_A')
            merged.meta.update(TARGET=target,MERGEMTH='Native channel segments; no relative scaling',NCHAN=4)
            merged.write(OUT/f'spectra/{target}.fits',overwrite=True)
            merged.write(OUT/f'spectra/{target}.csv',format='ascii.csv',overwrite=True)
            ff=np.asarray(merged['FLUX']);mm=merged['MASK'];yy=ff[mm]
            low,high=np.percentile(yy,[1,99]);axes[0].set_ylim(min(0,low)-.05*(high-low),high*1.13)
            axes[0].set_title(f'{target} — {header["NEXP"]} exposures, {header["TEXPTIME"]:.0f} s; flux standard {header["FLUXSTD"]}')
            axes[0].set_ylabel('Flux density (10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹)');axes[0].legend(ncol=4,fontsize=8)
            for lo,hi in [(6850,6940),(7160,7330),(7570,7700),(8120,8350),(8900,9700)]:
                for ax in axes:ax.axvspan(lo,hi,color='gray',alpha=.08)
            axes[1].set(xlabel='Vacuum heliocentric wavelength (Å)',ylabel='S/N per pixel',xlim=(3240,10235))
            axes[1].axhline(0,c='gray',lw=.5)
            fig.text(.07,.01,'Shading: statistical error. Gray bands: telluric regions. Absolute calibration/slit-loss uncertainty is not in the error bars.',fontsize=8)
            fig.tight_layout(rect=(0,.03,1,1));pdf.savefig(fig);fig.savefig(OUT/f'plots/{target}.png',dpi=150);plt.close(fig)
            summaries.append(dict(NAME=target,NEXP=int(header['NEXP']),EXPTIME_S=float(header['TEXPTIME']),
                RA_DEG=float(header['RA']),DEC_DEG=float(header['DEC']),FLUX_STANDARD=header['FLUXSTD'],
                **{f'SNR_{ch.upper()}':snrs[ch] for ch in 'ugri'}))
    Table(rows=summaries).write(OUT/'target_summary.csv',overwrite=True)
    Table(rows=overlap).write(OUT/'calibration/channel_overlap.csv',overwrite=True)
    print('EXPORTED',len(targets),'targets to',OUT)

if __name__=='__main__':main()
