"""Nested spectral bands and diagnostic hulls for the 24-target scenario grid."""
from pathlib import Path
import os,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/clagn-matplotlib')
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.spatial import ConvexHull
HERE=Path(__file__).resolve().parent
BLUE,ORANGE,INK,GOLD='#2873a0','#d1683b','#283542','#a37a25'
COLORS={'weak':BLUE,'strong':ORANGE}


def hull(ax,points,color):
    points=np.unique(points,axis=0)
    if len(points)>=3 and np.linalg.matrix_rank(points-points.mean(axis=0))==2:
        h=ConvexHull(points); p=points[h.vertices]
        ax.fill(p[:,0],p[:,1],facecolor=color,edgecolor=color,alpha=.14,lw=1.1,zorder=1)


def main():
    end=pd.read_csv(HERE/'selected_endpoints.csv');grid=pd.read_csv(HERE/'scenario_grid.csv')
    bands=pd.read_csv(HERE/'spectral_bands.csv');data=np.load(HERE/'spectra.npz');w=data['rest_um'];curves=data['normalized_flux']
    meta=json.loads((HERE/'provenance.json').read_text());common=meta['common_mrs_rest_um']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':INK,'text.color':INK})
    fig=plt.figure(figsize=(8.2,6.8))
    a=fig.add_axes([.115,.55,.84,.32]);b=fig.add_axes([.115,.17,.325,.235]);c=fig.add_axes([.61,.17,.345,.235])
    fig.text(.095,.974,'24 measured histories → a range of MIRI outcomes',weight='bold',fontsize=13,va='top')
    fig.text(.095,.930,'Common WISE-driven model grid • each spectrum normalized at rest 4 µm',fontsize=10)
    for case in ['weak','strong']:
        q=bands[bands.response_case.eq(case)];color=COLORS[case]
        a.fill_between(q.rest_um,q.p10,q.p90,color=color,alpha=.13,lw=0)
        a.fill_between(q.rest_um,q.p25,q.p75,color=color,alpha=.30,lw=0)
        a.plot(q.rest_um,q.p50,color=color,lw=2.3,label=('Weaker' if case=='weak' else 'Stronger')+' warm endpoint')
    a.axhline(1,color='#c1c7ca',lw=.7,ls=':')
    for wave in [9.7,18.]:
        a.axvline(wave,color=GOLD,ls=':',lw=.8,zorder=0)
        a.text(wave,.975,f'{wave:g} µm silicate',transform=a.get_xaxis_transform(),ha='center',va='top',fontsize=9,color=GOLD)
    a.set(xlim=(4,common[1]),ylim=(0,2.7),xlabel='Rest wavelength (µm)',ylabel='Model flux / rest-4 µm flux')
    a.set_xticks([4,6,8,10,13,16,18,21])
    a.legend(loc='lower left',frameon=False,fontsize=9)
    a.text(.98,.06,'Dark: middle 50%\nLight: middle 80% of targets',transform=a.transAxes,ha='right',fontsize=8.5,color=INK)
    b.set_title('Warm spectral shape',loc='left',fontsize=11.5,weight='bold',pad=10)
    c.set_title('Paired silicate contrasts',loc='left',fontsize=11.5,weight='bold',pad=10)
    for axis,x,y in [(b,'F14_over_F4','F21_over_F14'),(c,'S9_7_local','S18_local')]:
        # All 24 endpoint locations contribute; no thinning or omitted outliers.
        for case in ['weak','strong']:
            d=end[end.response_case.eq(case)];color=COLORS[case]
            hull(axis,d[[x,y]].to_numpy(),color)
            for hypothesis,marker in [('delayed','o'),('evolving','^')]:
                chosen=d[d.hypothesis.eq(hypothesis)]
                axis.scatter(chosen[x],chosen[y],s=22,marker=marker,color=color,edgecolor='white',lw=.35,zorder=3,alpha=.85)
        axis.margins(x=.13,y=.18);axis.tick_params(labelsize=9)
    b.set(xlabel=r'$F_\nu(14\,\mu\mathrm{m})\,/\,F_\nu(4\,\mu\mathrm{m})$',ylabel=r'$F_\nu(21\,\mu\mathrm{m})\,/\,F_\nu(14\,\mu\mathrm{m})$')
    c.set(xlabel=r'$S_{9.7}$  (local continuum)',ylabel=r'$S_{18}$  (local continuum)')
    handles=[Line2D([],[],marker='o',ls='',color=INK,label='Fixed dust'),Line2D([],[],marker='^',ls='',color=INK,label='Evolving boundary'),Patch(facecolor='#a9b3ba',alpha=.22,label='Span of 24 target endpoints')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.525,.060),frameon=False,ncol=3,fontsize=8.5)
    fig.text(.095,.027,'Colors identify weaker/stronger warm outcomes; the shaded spans are not confidence or classification regions.',fontsize=8.4)
    fig.savefig(HERE/'ensemble_scenarios.pdf');fig.savefig(HERE/'ensemble_scenarios.png',dpi=170);plt.close(fig)

    # A complete atlas makes every object and its individual endpoints inspectable.
    fig,axes=plt.subplots(4,6,figsize=(13,8.1),sharex=True,sharey=True,layout='constrained')
    for ax,ident in zip(axes.flat,grid.id.drop_duplicates()):
        selected=end[end.id.eq(ident)]
        all_indices=grid.loc[grid.id.eq(ident),'curve_index'].to_numpy()
        ax.fill_between(w,curves[all_indices].min(axis=0),curves[all_indices].max(axis=0),color='#d6dce0',alpha=.45,lw=0)
        for row in selected.itertuples():
            ax.plot(w,curves[row.curve_index],color=COLORS[row.response_case],lw=1.4)
        ax.set(xlim=(4,common[1]),yscale='log',xticks=[5,10,20])
        ax.axvline(9.7,color=GOLD,ls=':',lw=.5);ax.axvline(18,color=GOLD,ls=':',lw=.5)
        ax.set_title(ident+'  '+('fading' if selected.family.iloc[0]=='fade' else 'rising')+' selection',fontsize=9)
        ax.tick_params(labelsize=8)
    fig.supxlabel('Rest wavelength (µm)',fontsize=11)
    fig.supylabel('Model Fν / Fν(rest 4 µm)',fontsize=11)
    fig.suptitle('All 24 targets: weak / strong endpoints and full conditional grid (gray)',fontsize=13)
    fig.savefig(HERE/'all_24_scenarios.pdf');fig.savefig(HERE/'all_24_scenarios.png',dpi=140);plt.close(fig)
    counts=end.groupby(['response_case','hypothesis']).size().to_dict()
    print('Endpoint families:',counts)
    print('Saved ensemble_scenarios.pdf and all_24_scenarios.pdf')


if __name__=='__main__':main()
