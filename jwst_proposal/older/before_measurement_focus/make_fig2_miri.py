"""Figure 2: measured history, conditional MIRI experiment, spectral constraints.

Run review/validate_fig2_recovery.py first. Does not compile the proposal.
"""
from pathlib import Path
import os,json,sys,hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/clagn-matplotlib')
import numpy as np,pandas as pd
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'review'))
from validate_fig2_recovery import diagnostics
from measure_fig2_epoch_points import measurements
BLUE,ORANGE,INK,GOLD,GREEN='#24658a','#bc5039','#283542','#a37a25','#34836b'

def main():
    out=BASE/'review/figure2_recovery'
    report=json.loads((out/'summary.json').read_text())
    mock=pd.read_csv(out/'mock_spectrum.csv')
    archive=BASE/'review/p9694_predictivity'
    arrays=np.load(archive/'conditional_curves.npz')
    meta=json.loads((archive/'adequacy_and_forecasts.json').read_text())
    keep=[i for i,r in enumerate(meta['conditional_curve_labels']) if r['year']==2028. and r['future_driver_factor']==1.]
    w,curves=arrays['rest_um'],arrays['flux_mjy'][keep]
    explored=diagnostics(w,curves)
    draws=np.load(out/'measurement_draws.npz')['diagnostics']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(8,7.0))
    a=fig.add_axes([.10,.78,.86,.145]);b=fig.add_axes([.10,.385,.86,.28])
    residual=fig.add_axes([.10,.275,.86,.083])
    a.set_title('A  A measured history gives the spectrum its time context',loc='left',weight='bold',fontsize=11,pad=9)
    history=pd.read_csv(BASE/'review/history_audit/binned_measurements.csv')
    aw=pd.read_csv(BASE/'inputs/sample_allwise_psd.csv').set_index('id').loc['R01']
    normals={}
    for band,color,zero in [('g',GREEN,3631e3),('W1',BLUE,309540.),('W2',ORANGE,171787.)]:
        d=history[history.id.eq('R01')&history.band.eq(band)].sort_values('mjd')
        t=Time(d.mjd.values,format='mjd').decimalyear;f=zero*10**(-.4*d.mag.values)
        normal=np.median(f[(t>=2018)&(t<=2020)])
        normals[band]=normal
        a.plot(t,f/normal,'o-',color=color,ms=2.2,lw=.7,label=band,alpha=.85)
        if band.startswith('W'):
            cat=zero*10**(-.4*aw[band.lower()+'mpro'])
            a.plot(Time(aw['w1mjdmean'],format='mjd').decimalyear,cat/normal,
                   'o',mfc='white',mec=color,ms=4)
    epochs=measurements()
    sp=epochs['spherex']
    a.errorbar([r['year'] for r in sp],[r['flux_mjy']/normals['W1'] for r in sp],
               yerr=[r['error_mjy']/normals['W1'] for r in sp],fmt='s',ms=4.5,
               color=GOLD,mec='white',mew=.6,capsize=2,lw=.9,zorder=6,label='SPHEREx → W1')
    a.axvline(2028.,color=INK,ls='--',lw=1)
    a.text(2028.,1.03,'MIRI*',transform=a.get_xaxis_transform(),ha='center',fontsize=8,color=INK)
    a.set(xlim=(2009.5,2029.),ylabel='Relative flux',xlabel='Observed year')
    a.set_xticks([2010,2014,2018,2022,2026,2028]);a.margins(y=.2)
    a.legend(ncol=2,loc='upper left',frameon=False,fontsize=7.5,handlelength=1.6,columnspacing=1.3,labelspacing=.2)
    a.text(.01,.03,'R01 · SPHEREx shares the W1 normalization',transform=a.transAxes,fontsize=7.5,color=INK)
    b.set_title('B  MIRI measures the warm continuum and both silicate features',loc='left',weight='bold',fontsize=11,pad=12)
    spx=pd.read_csv(BASE/'inputs/spaxel_scryer_P9694_ra11p5333_dec9p1225_cleaned.csv')
    b.errorbar(spx.wavelength_um/1.2377224,spx.flux_mjy,yerr=spx.flux_err_mjy,
               fmt='.',ms=2,color=GREEN,alpha=.55,lw=.35,label='Measured SPHEREx')
    b.fill_between(w,curves.min(axis=0),curves.max(axis=0),color='#dce2e6',alpha=.8,
                   label='Explored predictions before MIRI')
    b.errorbar(mock.rest_um,mock.mock_mjy,yerr=mock.statistical_error_mjy,
               fmt='.',ms=2.5,lw=.45,color=INK,alpha=.65,label='Simulated MIRI*',zorder=4)
    b.plot(mock.rest_um,mock['fit_eta0.0_mjy'],color=BLUE,lw=1.6,label='Refitted fixed dust',zorder=5)
    b.plot(mock.rest_um,mock['fit_eta-0.5_mjy'],color=ORANGE,lw=1.6,ls=(0,(4,2)),label='Refitted evolving dust',zorder=6)
    for center in [9.7,18.]:
        b.axvline(center,color=GOLD,lw=.7,ls=':',alpha=.8)
        b.text(center,.97,f'{center:g} µm',ha='center',va='top',transform=b.get_xaxis_transform(),fontsize=8,color=GOLD)
    pah_bands=[6.2,7.7,8.6,11.3,12.7,17.]
    for center in pah_bands:
        b.plot([center,center],[.79,.84],transform=b.get_xaxis_transform(),
               color='#82549a',lw=.9)
        label_y=.90 if center==12.7 else .86
        b.text(center,label_y,f'{center:g}',ha='center',va='bottom',
               transform=b.get_xaxis_transform(),fontsize=7.5,color='#82549a')
    b.text(5.65,.815,'PAH',ha='right',va='center',transform=b.get_xaxis_transform(),
           fontsize=8,color='#82549a')
    b.set(xscale='log',yscale='log',xlim=(.65,23),ylim=(.7,100),ylabel='Flux density (mJy)')
    b.set_xticks([1,2,3,5,10,20]);b.xaxis.set_major_formatter(ScalarFormatter())
    b.set_yticks([1,3,10,30,100]);b.yaxis.set_major_formatter(ScalarFormatter())
    b.tick_params(labelbottom=False)
    b.legend(loc='upper left',fontsize=7.5,frameon=False,labelspacing=.25)
    b.text(.98,.06,'Mock excludes PAH emission\nBoth response models fit this mock spectrum',
           transform=b.transAxes,ha='right',color=INK,fontsize=8,linespacing=1.4)
    error=100*mock.total_error_mjy/mock.truth_mjy
    residual.fill_between(mock.rest_um,-error,error,color='#e1e6e9')
    for key,color,ls in [('0.0',BLUE,'-'),('-0.5',ORANGE,'--')]:
        residual.plot(mock.rest_um,100*(mock[f'fit_eta{key}_mjy']/mock.mock_mjy-1),color=color,lw=.85,ls=ls)
    residual.axhline(0,color=INK,lw=.5)
    residual.set(xscale='log',xlim=(.65,23),ylim=(-23,23),yticks=[-15,0,15],
                 xlabel='Rest wavelength (µm)',ylabel='Residual (%)')
    residual.set_xticks([1,2,3,5,10,20]);residual.xaxis.set_major_formatter(ScalarFormatter())
    residual.text(.02,.78,'Shading: statistical + correlated calibration / host errors',fontsize=7.2,transform=residual.transAxes,color=INK)
    fig.text(.10,.185,'C  Precision on three observables, even when response models overlap',fontsize=11,weight='bold')
    titles=[r'8–13 µm flux ($10^{-12}$ erg s$^{-1}$ cm$^{-2}$)',r'$S_{9.7}$ (local continuum)',r'$S_{18}$ (local continuum)']
    for j in range(3):
        axis=fig.add_axes([.10+j*.305,.093,.25,.058]);scale=1e12 if j==0 else 1
        lo,hi=explored[:,j].min()*scale,explored[:,j].max()*scale
        q=np.quantile(draws[:,j],[.16,.5,.84])*scale
        axis.plot([lo,hi],[1,1],color='#aeb8be',lw=4,solid_capstyle='round')
        axis.errorbar(q[1],0,xerr=[[q[1]-q[0]],[q[2]-q[1]]],fmt='o',color=ORANGE,ms=4,lw=2,capsize=3)
        precision=['≈5% flux uncertainty','σ(S₉.₇) ≈ 0.08','σ(S₁₈) ≈ 0.04'][j]
        axis.text(.5,1.10,precision,transform=axis.transAxes,ha='center',color=ORANGE,fontsize=8)
        axis.set(ylim=(-.7,1.7),yticks=[],xlabel=titles[j])
        axis.tick_params(axis='x',labelsize=8,pad=1);axis.xaxis.label.set_size(8)
        axis.spines['left'].set_visible(False);axis.set_xlim(lo-(hi-lo)*.07,hi+(hi-lo)*.07)
    fig.text(.10,.025,'Grey: explored model span   Orange: simulated measurement, 16–84% range   *Illustrative 2028 observation',fontsize=7.6,color=INK)
    fig.savefig(BASE/'fig2_miri.pdf');fig.savefig(BASE/'review/fig2_miri_preview.png',dpi=180);plt.close(fig)
    prov=dict(generator=Path(__file__).name,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              simulation='review/figure2_recovery/summary.json',conditional_models=len(keep),
              measured_epoch_points='inputs/fig2_epoch_points.json',
              pah_markers_rest_um=pah_bands,pah_emission_in_mock=False,
              interval_meaning='Grey finite model range, not a posterior; orange repeated-measurement quantiles, not a model posterior.',
              diagnostics=report['truth_diagnostics'],measurement_std=report['measurement_std'],
              scope='One conditional R01 experiment; no sample-wide mechanism recovery claimed.')
    (BASE/'inputs/fig2_miri_provenance.json').write_text(json.dumps(prov,indent=2)+'\n')
    print(json.dumps(prov,indent=2))

if __name__=='__main__':main()
