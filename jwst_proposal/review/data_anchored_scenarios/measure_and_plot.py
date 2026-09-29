"""Measure saved sample data and illustrate conditional MIRI outcomes for R01.

No proposal file is edited. No new physical fits or synthetic detections.
"""
from pathlib import Path
import sys, os, json, hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/clagn-matplotlib')
HERE=Path(__file__).resolve().parent; BASE=HERE.parents[1]; sys.path.insert(0,str(BASE))
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from astropy.time import Time
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

ZERO={'g':3631000.,'r':3631000.,'W1':309540.,'W2':171787.,'W3':31674.,'W4':8363.}
BLUE,RED,INK,GREEN,GOLD='#24658a','#bc5039','#283542','#42806c','#a37a25'


def measurements():
    sample=pd.read_csv(BASE/'inputs/jwst_sample_cycle6.csv')
    history=pd.read_csv(BASE/'review/history_audit/binned_measurements.csv')
    aw=pd.read_csv(BASE/'inputs/sample_allwise_psd.csv').set_index('id')
    arch=pd.read_csv(BASE/'inputs/archival_spectra/index.csv')
    rows=[]; optical=[]; spx=[]
    for target in sample.itertuples():
        row=dict(id=target.id,target=target.target,family=target.family,z=target.z)
        for band in ['g','r','W1','W2']:
            d=history[history.id.eq(target.id)&history.band.eq(band)].sort_values('mjd')
            f=ZERO[band]*10**(-.4*d.mag.to_numpy())
            row[band+'_n_epochs']=len(d)
            if len(d):
                row[band+'_first_year']=float(Time(d.mjd.iloc[0],format='mjd').decimalyear)
                row[band+'_last_year']=float(Time(d.mjd.iloc[-1],format='mjd').decimalyear)
                row[band+'_p95_over_p05']=float(np.percentile(f,95)/np.percentile(f,5))
                row[band+'_last3_over_first3']=float(np.median(f[-3:])/np.median(f[:3]))
                row[band+'_last3_mjy']=float(np.median(f[-3:]))
        for band in ['W1','W2','W3','W4']:
            mag=aw.loc[target.id,band.lower()+'mpro']; em=aw.loc[target.id,band.lower()+'sigmpro']
            flux=ZERO[band]*10**(-.4*mag)
            row[band+'_2010_catalog_mjy']=flux
            row[band+'_2010_catalog_error_mjy']=flux*np.log(10)/2.5*em
        for band in ['W1','W2']:
            row[band+'_last3_over_2010']=row[band+'_last3_mjy']/row[band+'_2010_catalog_mjy']
        path=BASE/f'inputs/spherex_extracted/{target.id}.csv'
        if target.id=='R01':path=BASE/'inputs/spaxel_scryer_P9694_ra11p5333_dec9p1225_cleaned.csv'
        row['spherex_extracted_available']=path.exists()
        if path.exists():
            d=pd.read_csv(path)
            if 'spectrum_type' in d:d=d[d.spectrum_type.eq('cleaned')]
            d=d[np.isfinite(d.flux_mjy)&(d.flux_err_mjy>0)].sort_values('mjds')
            groups=np.split(np.arange(len(d)),np.flatnonzero(np.diff(d.mjds)>45)+1)
            for group in groups:
                e=d.iloc[group];rest=e.wavelength_um/(1+target.z)
                out=dict(id=target.id,first_mjd=float(e.mjds.min()),last_mjd=float(e.mjds.max()),
                         n_native_points=len(e),source=str(path.relative_to(BASE)))
                for label,lo,hi in [('short',2.23,2.45),('long',3.55,3.90)]:
                    b=e[(rest>=lo)&(rest<=hi)]
                    out[label+'_n']=len(b)
                    out[label+'_median_mjy']=float(b.flux_mjy.median()) if len(b) else None
                    out[label+'_median_point_error_mjy']=float(b.flux_err_mjy.median()) if len(b) else None
                out['long_over_short_total_flux']=out['long_median_mjy']/out['short_median_mjy'] if out['long_n'] and out['short_n'] else None
                spx.append(out)
        paths=[]
        for a in arch[arch.id.eq(target.id)].itertuples():
            d=pd.read_csv(BASE/a.file)
            paths.append((a.file,float(a.mjd),d.wave_A.to_numpy(),d.flux.to_numpy(),d.ivar.to_numpy()>0))
        npath=BASE.parent/f'sep23_data/reduction_20260924/products_p330e/spectra/{target.internal_id}.csv'
        row['ngps_reduced_available']=npath.exists()
        if npath.exists():
            d=pd.read_csv(npath); t=float(fits.getheader(npath.with_suffix('.fits'),1).get('MJD',61307.))
            paths.append((str(npath.relative_to(BASE.parent)),t,d.WAVE_VAC_HELIO_A.to_numpy(),d.FLUX.to_numpy(),d.MASK.to_numpy().astype(bool)))
        row['optical_spectra_available']=len(paths)
        for source,mjd,w,f,mask in paths:
            use=mask&np.isfinite(f)&(w/(1+target.z)>=5050)&(w/(1+target.z)<=5150)
            flux=f*1e-17*w*w/2.99792458e18/1e-26
            if use.sum()>=5:
                optical.append(dict(id=target.id,mjd=mjd,year=float(Time(mjd,format='mjd').decimalyear),
                                    rest_window_A='5050-5150',n_pixels=int(use.sum()),
                                    median_total_fnu_mjy=float(np.median(flux[use])),source=source))
        rows.append(row)
    table=pd.DataFrame(rows);table.to_csv(HERE/'sample_measurements.csv',index=False)
    pd.DataFrame(spx).to_csv(HERE/'spherex_epoch_measurements.csv',index=False)
    pd.DataFrame(optical).sort_values(['id','mjd']).to_csv(HERE/'optical_continuum_measurements.csv',index=False)
    return table,history,aw


def scenarios():
    root=BASE/'review/p9694_predictivity'
    meta=json.loads((root/'adequacy_and_forecasts.json').read_text())
    arrays=np.load(root/'conditional_curves.npz');wave=arrays['rest_um'];flux=arrays['flux_mjy']
    labels=meta['conditional_curve_labels']
    choose=np.array([r['evolution']==0 and r['year']==2028. and r['future_driver_factor']==1. for r in labels])
    idx=np.flatnonzero(choose)
    # Integrate Fnu over frequency with exact 8 and 13 micron endpoints.
    band=np.r_[8.,wave[(wave>8)&(wave<13)],13.]
    nuclear_plus_host=np.array([np.interp(band,wave,f) for f in flux])
    integrated=-trapezoid(nuclear_plus_host,299792458./(band*1e-6),axis=1)
    select=idx[[np.argmin(integrated[idx]),np.argmax(integrated[idx])]]
    records=[]
    for i,name in zip(select,['weak_warm','strong_warm']):
        records.append(dict(label=name,archive_curve_index=int(i),**labels[i],
                            Fnu_rest12_mjy=float(np.interp(12,wave,flux[i])),
                            Fnu_rest18_mjy=float(np.interp(18,wave,flux[i])),
                            integrated_8_13_observed_erg_s_cm2=float(integrated[i]*1e-26/(1+.2377224))))
    pd.DataFrame(dict(rest_um=wave,weak_warm_observed_mjy=flux[select[0]],strong_warm_observed_mjy=flux[select[1]])).to_csv(HERE/'r01_miri_scenarios.csv',index=False)
    return wave,flux[select],records,meta


def main():
    table,history,aw=measurements()
    wave,scenario,records,meta=scenarios()
    row=table.set_index('id').loc['R01'];z=float(row.z)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(8,6.2))
    ax=fig.add_axes([.09,.60,.86,.29]);sp=fig.add_axes([.09,.13,.86,.34])
    ax.set_title('R01 / P9694: measured history anchors the experiment',loc='left',fontsize=12,weight='bold')
    for band,color,marker in [('W1',BLUE,'o'),('W2',RED,'s')]:
        d=history[history.id.eq('R01')&history.band.eq(band)].sort_values('mjd')
        f=ZERO[band]*10**(-.4*d.mag.to_numpy());ref=ZERO[band]*10**(-.4*aw.loc['R01',band.lower()+'mpro'])
        y=Time(d.mjd.to_numpy(),format='mjd').decimalyear
        err=f*np.log(10)/2.5*d.mag_se.to_numpy()/ref
        ax.errorbar(y,f/ref,yerr=err,fmt=marker,ms=3,color=color,label=band+' measured',lw=.7)
        ax.scatter([2010],[1],s=22,facecolors='white',edgecolors=color,marker=marker)
    ax.axhline(1,color='.7',lw=.7)
    ax.set(xlim=(2009.5,2027.5),ylabel='Total flux / 2010 catalogue',xlabel='Year')
    ax.legend(loc='upper left',frameon=False,ncol=2,fontsize=9)
    ax.text(.98,.08,f"Latest three W1/W2 visits: {row.W1_last3_over_2010:.2f} / {row.W2_last3_over_2010:.2f} × 2010",transform=ax.transAxes,ha='right',fontsize=9)
    sp.set_title('Same existing data; contrasting possible MIRI outcomes',loc='left',fontsize=12,weight='bold')
    d=pd.read_csv(BASE/'inputs/spaxel_scryer_P9694_ra11p5333_dec9p1225_cleaned.csv')
    sp.errorbar(d.wavelength_um/(1+z),d.flux_mjy,yerr=d.flux_err_mjy,fmt='.',ms=2,color=INK,alpha=.5,lw=.4,
                label='SPHEREx 2025–26: measured')
    optical=pd.read_csv(HERE/'optical_continuum_measurements.csv')
    optical=optical[optical.id.eq('R01')].sort_values('mjd')
    for i, point in enumerate(optical.itertuples()):
        sp.plot(.51,point.median_total_fnu_mjy,'D',ms=4,color='#805298',
                label='Optical 2018 / 2026: measured' if i==0 else None)
    sp.axvspan(wave.min(),wave.max(),color='#eef3f6',zorder=0)
    for f,color,name in zip(scenario,[BLUE,RED],['Weak warm response','Strong warm response']):
        sp.plot(wave,f,color=color,lw=2.1,label=name+' (scenario)')
    for band,w in [('W3',11.5608),('W4',22.0883)]:
        f=row[band+'_2010_catalog_mjy'];e=row[band+'_2010_catalog_error_mjy']
        sp.errorbar(w/(1+z),f,yerr=e,fmt='s',mfc='white',mec=GOLD,color=GOLD,ms=5,capsize=2)
        sp.text(w/(1+z),f*1.18,'2010 '+band,ha='center',color=GOLD,fontsize=8)
    for w in [9.7,18.]:sp.axvline(w,color=GOLD,lw=.7,ls=':',alpha=.6)
    sp.set(xscale='log',yscale='log',xlim=(.45,23),ylim=(.45,115),xlabel='Rest wavelength (µm)',ylabel='Observed flux density (mJy)')
    sp.set_xticks([1,2,3,5,10,20]);sp.xaxis.set_major_formatter(ScalarFormatter())
    sp.set_yticks([1,3,10,30,100]);sp.yaxis.set_major_formatter(ScalarFormatter())
    sp.legend(loc='upper left',fontsize=8,frameon=False)
    fig.text(.09,.040,'Conditional fixed-dust scenarios: 2028, same future illumination, assumed 5% extra historical scatter.',fontsize=8.5)
    fig.savefig(HERE/'data_anchored_preview.pdf');fig.savefig(HERE/'data_anchored_preview.png',dpi=160);plt.close(fig)
    summary=dict(targets_measured=len(table),available_spherex_extractions=table.loc[table.spherex_extracted_available,'id'].tolist(),
                 available_ngps_reductions=table.loc[table.ngps_reduced_available,'id'].tolist(),
                 target_optical_spectra_available=int(table.optical_spectra_available.sum()),
                 W1_p95_over_p05_range=[float(table.W1_p95_over_p05.min()),float(table.W1_p95_over_p05.max())],
                 r01_latest3_over_2010=dict(W1=float(row.W1_last3_over_2010),W2=float(row.W2_last3_over_2010)),
                 scenario_scope='Two endpoints of archived fixed-dust shell profiles at 2028.0 with identical unit future driver. Conditional on an assumed additional 5% independent scatter; neither is an evolution-exclusive outcome.',
                 scenarios=records,contrast_8_13=records[1]['integrated_8_13_observed_erg_s_cm2']/records[0]['integrated_8_13_observed_erg_s_cm2'],
                 historical_fit_scope='Existing optical/WISE/SPHEREx shell fits; not newly refitted here. Primary fits fail supplied-error adequacy; scenario families use explicit extra scatter. See existing predictivity audit.',
                 source_hashes={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [BASE/'inputs/jwst_sample_cycle6.csv',BASE/'review/history_audit/binned_measurements.csv',BASE/'inputs/sample_allwise_psd.csv',BASE/'review/p9694_predictivity/conditional_curves.npz',BASE/'review/p9694_predictivity/adequacy_and_forecasts.json']})
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['source_hashes','historical_fit_scope']},indent=2))


if __name__=='__main__':main()
