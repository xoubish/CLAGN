"""Measured SED and the spectral diagnostic planes for the same MIRI scenarios.

All points in the diagnostic planes are archived conditional models, not data.
Silicate indices are explicit local-continuum contrasts, not full decomposition.
"""
from pathlib import Path
import os,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/clagn-matplotlib')
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
BASE=Path(__file__).resolve().parent
HERE=BASE/'review/data_anchored_scenarios'
BLUE,RED,INK,GOLD='#24658a','#bc5039','#283542','#a37a25'


def indices(w,f):
    at=lambda q:float(np.interp(q,w,f))
    baseline=lambda q,x:float(np.exp(np.interp(np.log(q),np.log(x),np.log([at(k) for k in x]))))
    return dict(F14_over_F4=at(14)/at(4),F22_over_F14=at(22)/at(14),
                S9_7_local=float(np.log(at(9.7)/baseline(9.7,[6,14]))),
                S18_local=float(np.log(at(18)/baseline(18,[14,22]))))


def main():
    archive=BASE/'review/p9694_predictivity';a=np.load(archive/'conditional_curves.npz');w=a['rest_um'];f=a['flux_mjy']
    meta=json.loads((archive/'adequacy_and_forecasts.json').read_text())
    summary=json.loads((HERE/'summary.json').read_text()); selected=[r['archive_curve_index'] for r in summary['scenarios']]
    rows=[]
    for i,r in enumerate(meta['conditional_curve_labels']):
        if r['year']==2028 and r['future_driver_factor']==1:
            rows.append(dict(archive_curve_index=i,**r,**indices(w,f[i])))
    tab=pd.DataFrame(rows);tab.to_csv(BASE/'inputs/fig2_miri_diagnostics.csv',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(8,6.35))
    ax=fig.add_axes([.09,.54,.86,.36]);b=fig.add_axes([.09,.18,.35,.21]);c=fig.add_axes([.60,.18,.35,.21])
    ax.set_title('A  What MIRI adds: the warm-dust spectrum of R01',loc='left',fontsize=12,weight='bold',pad=10)
    z=.2377224
    spx=pd.read_csv(BASE/'inputs/spaxel_scryer_P9694_ra11p5333_dec9p1225_cleaned.csv')
    ax.errorbar(spx.wavelength_um/(1+z),spx.flux_mjy,yerr=spx.flux_err_mjy,
                fmt='.',ms=2,color=INK,alpha=.55,lw=.4,label='SPHEREx 2025–26')
    opt=pd.read_csv(HERE/'optical_continuum_measurements.csv');opt=opt[opt.id.eq('R01')].sort_values('mjd')
    ax.plot([.51]*len(opt),opt.median_total_fnu_mjy,'D',ms=4,color='#805298',label='Optical 2018 / 2026')
    ax.axvspan(w.min(),w.max(),color='#eef3f6',zorder=0)
    for i,color,label in zip(selected,[BLUE,RED],['Weak warm emission','Strong warm emission']):
        ax.plot(w,f[i],lw=2.2,color=color,label=label)
    row=pd.read_csv(HERE/'sample_measurements.csv').set_index('id').loc['R01']
    for band,lam in [('W3',11.5608),('W4',22.0883)]:
        flux=row[band+'_2010_catalog_mjy'];err=row[band+'_2010_catalog_error_mjy']
        ax.errorbar(lam/(1+z),flux,yerr=err,fmt='s',mfc='white',mec=GOLD,color=GOLD,ms=5,capsize=2)
        ax.text(lam/(1+z),flux*1.18,'2010 '+band,ha='center',color=GOLD,fontsize=8)
    ax.text(.98,.04,'MIRI/MRS coverage',transform=ax.transAxes,ha='right',color=BLUE,fontsize=10,weight='bold')
    for lam in [9.7,18.]:ax.axvline(lam,color=GOLD,lw=.7,ls=':',alpha=.55)
    ax.set(xscale='log',yscale='log',xlim=(.45,23),ylim=(.45,115),xlabel='Rest wavelength (µm)',ylabel='Observed flux density (mJy)')
    ax.set_xticks([1,2,3,5,10,20]);ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_yticks([1,3,10,30,100]);ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.legend(loc='upper left',fontsize=8,frameon=False)
    b.set_title('B  Warm continuum shape',loc='left',fontsize=11,weight='bold',pad=10)
    c.set_title('C  Paired silicate contrasts',loc='left',fontsize=11,weight='bold',pad=10)
    for eta in [0,-.5,-1]:
        t=tab[tab.evolution.eq(eta)]
        for axis,x,y in [(b,'F14_over_F4','F22_over_F14'),(c,'S9_7_local','S18_local')]:
            axis.scatter(t[x],t[y],s=24,marker='o' if eta==0 else '^',facecolor='none',edgecolor='#a5afb5',lw=.8,zorder=2)
    for i,color,name in zip(selected,[BLUE,RED],['Weak','Strong']):
        t=tab.set_index('archive_curve_index').loc[i]
        for axis,x,y in [(b,'F14_over_F4','F22_over_F14'),(c,'S9_7_local','S18_local')]:
            axis.scatter(t[x],t[y],s=48,color=color,zorder=4)
            axis.annotate(name,(t[x],t[y]),xytext=(7,7 if name=='Strong' else -14),textcoords='offset points',fontsize=9,color=color,weight='bold')
    b.set(xlabel=r'$F_\nu(14\,\mu\mathrm{m})\,/\,F_\nu(4\,\mu\mathrm{m})$',ylabel=r'$F_\nu(22\,\mu\mathrm{m})\,/\,F_\nu(14\,\mu\mathrm{m})$')
    c.set(xlabel=r'$S_{9.7}$  (local continuum)',ylabel=r'$S_{18}$  (local continuum)')
    for axis in [b,c]:axis.margins(x=.2,y=.25);axis.tick_params(labelsize=9)
    fig.text(.09,.085,'Lower panels: conditional model locations; circles = fixed dust, triangles = evolving dust.',fontsize=8.5)
    fig.text(.09,.047,'One spectrum measures the warm continuum and both silicate features together.',fontsize=9,weight='bold',color=INK)
    fig.savefig(BASE/'fig2_miri.pdf');fig.savefig(BASE/'review/fig2_miri_preview.png',dpi=170);plt.close(fig)
    result=dict(selected_indices=selected,conditional_models=len(tab),observables=['Fnu14/Fnu4','Fnu22/Fnu14','S9.7_local','S18_local'],
                feature_definition='S=ln(Fnu/Fnu_local). Log-linear baselines through rest 6 and 14 microns for 9.7; 14 and 22 microns for 18. Illustrative contrast indices, not full nuclear/host/feature decomposition.',
                historical_scatter_fraction=0.05, source_archive='review/p9694_predictivity/conditional_curves.npz', source_labels='review/p9694_predictivity/adequacy_and_forecasts.json', scope='Same saved spectra as previous preview. Diagnostic points are predictions, not measurements. All models use the same forecast date and future driver factor. Neither color denotes dust evolution.',
                selected=[tab.set_index('archive_curve_index').loc[i].to_dict() for i in selected])
    (BASE/'inputs/fig2_miri_provenance.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
