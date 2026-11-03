"""Readable pilot figures and source-by-source photometry timing audit."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.time import Time
from paper_protocol import OUT


def fit_figure(name,variant='baseline'):
    results=pd.read_csv(OUT/'fit_results.csv')
    d=results[results.name.eq(name)&results.variant.eq(variant)&results.status.eq('fit')].sort_values('mjd')
    fig,axes=plt.subplots(3,len(d),figsize=(4.4*len(d),8.3),squeeze=False)
    for col,r in enumerate(d.itertuples(index=False)):
        a=np.load(OUT/'fits'/f'{r.epoch_id}__{variant}.npz')
        w,y,e,good,parts=[a[k] for k in ['wave_rest','flux','error','good','components']]
        model=parts.sum(axis=1);scale=r.total5100
        ax=axes[0,col]
        ax.plot(w,np.where(good,y,np.nan)/scale,color='0.45',lw=.6,label='data')
        ax.plot(w,model/scale,'k',lw=1,label='model')
        for v,color,label in [(parts[:,0],'#21866a','AGN'),(parts[:,1:3].sum(axis=1),'#b77725','host'),(parts[:,3],'#8d56ab','Fe II')]:
            ax.plot(w,v/scale,color=color,lw=1,label=label)
        ax.set(title=f'{r.source} {Time(r.mjd,format="mjd").strftime("%Y-%m-%d")}\nEW={r.hb_ew_total_cont:.1f} Å; χ²/ν={r.reduced_chi2:.2f}',xlim=(4200,5500),ylim=(-.2,min(6,np.nanpercentile(y[good]/scale,99.5)*1.25)))
        ax=axes[1,col];continuum=parts[:,:4].sum(axis=1)
        ax.plot(w,np.where(good,y-continuum,np.nan)/scale,color='0.4',lw=.7,label='data − continuum/Fe II')
        ax.plot(w,(parts[:,4:].sum(axis=1))/scale,'k',lw=1)
        ax.plot(w,parts[:,7:9].sum(axis=1)/scale,color='#b83343',lw=1.5,label='broad Hβ')
        ax.plot(w,parts[:,4:7].sum(axis=1)/scale,color='#247ba0',label='narrow Hβ + [O III]')
        ax.set(xlim=(4730,5080),ylim=(-.2,min(6,np.nanpercentile((y-continuum)[good]/scale,99.5)*1.25)))
        ax.axhline(0,color='0.7',lw=.5)
        ax=axes[2,col];ax.plot(w,np.where(good,(y-model)/e,np.nan),color='0.35',lw=.5)
        ax.axhline(0,color='k',lw=.5);ax.axhspan(-3,3,color='0.5',alpha=.12)
        ax.set(xlim=(4200,5500),ylim=(-10,10),xlabel='Rest vacuum wavelength [Å]')
        axes[0,col].set_ylabel('fλ / fitted continuum at 5100 Å')
        axes[1,col].set_ylabel('Continuum-subtracted / C5100')
        axes[2,col].set_ylabel('(data − model) / pixel σ')
    axes[0,0].legend(fontsize=7,ncol=3);axes[1,0].legend(fontsize=7)
    fig.suptitle(f'{name} — {variant}: independent fits; conditional model diagnostics, not discovery significance',fontsize=12)
    fig.tight_layout(rect=(0,0,1,.96));return fig


def envelope_figure(comparison='latest_archive'):
    d=pd.read_csv(OUT/'model_sensitivity_envelope.csv');d=d[d.comparison.eq(comparison)].sort_values(['name','epoch_id'])
    fig,axes=plt.subplots(1,2,figsize=(12,max(5,.36*len(d)+1)),sharey=True)
    for ax,key,title in zip(axes,['hb_ew_total_cont','hb_oiii'],['Broad Hβ equivalent width','Broad Hβ / [O III] λ5008']):
        for i,r in enumerate(d.to_dict('records')):
            lo,hi=r[key+'_ratio_min'],r[key+'_ratio_max']
            supported=r['both_hb_supported'] and (r['both_oiii_supported'] if key=='hb_oiii' else True)
            ax.plot([lo,hi],[i,i],color='#17656b' if supported else '#d18b31',lw=3)
            ax.plot(np.sqrt(lo*hi),i,'o',color='#17656b' if supported else '#d18b31',ms=4)
        ax.axvline(1,color='0.5',ls='--');ax.set_xscale('log');ax.set(title=title,xlabel='Palomar / reference; min–max across three models')
        ax.set_yticks(np.arange(len(d)),d.name+' '+d.epoch_id.str.split('_').str[-1]);ax.grid(axis='x',alpha=.2)
    axes[0].invert_yaxis();fig.suptitle('Pilot model envelope — not confidence intervals; orange = incomplete line-window support')
    fig.tight_layout();return fig


def optical_audit():
    visits=pd.read_csv(OUT/'pilot_epochs.csv');visits=visits[visits.source.eq('NGPS')]
    rows=[]
    for r in visits.itertuples(index=False):
        p=OUT/'optical_alert_refresh'/f'{r.name}_corrected.csv'
        data=pd.read_csv(p) if p.exists() and p.stat().st_size>1 else pd.DataFrame()
        for fid,band in [(1,'g'),(2,'r')]:
            g=data[data.fid.eq(fid)].copy() if len(data) else data
            if len(g):
                g['day']=np.floor(g.mjd)
                near=g[abs(g.mjd-r.mjd)<=30]
                nights=near.groupby(['oid','day']).magpsf_corr.median()
                closest=g.iloc[np.argmin(abs(g.mjd-r.mjd))]
                row=dict(nearest_dt_days=float(closest.mjd-r.mjd),n_points_30d=len(near),
                    n_nights_30d=near.day.nunique(),n_object_ids_30d=near.oid.nunique(),
                    median_mag_30d=float(nights.median()) if len(nights) else np.nan,
                    usable_local_anchor=bool(near.day.nunique()>=3 and near.oid.nunique()==1))
            else:row=dict(nearest_dt_days=np.nan,n_points_30d=0,n_nights_30d=0,n_object_ids_30d=0,median_mag_30d=np.nan,usable_local_anchor=False)
            rows.append(dict(name=r.name,epoch_id=r.epoch_id,band=band,**row))
    result=pd.DataFrame(rows);result.to_csv(OUT/'optical_visit_matches.csv',index=False)
    return result


def optical_figure(names):
    visits=pd.read_csv(OUT/'pilot_epochs.csv')
    fig,axes=plt.subplots(len(names),1,figsize=(12,2.3*len(names)),squeeze=False,sharex=True)
    for name,ax in zip(names,axes[:,0]):
        p=OUT/'optical_alert_refresh'/f'{name}_corrected.csv'
        if p.exists() and p.stat().st_size>1:
            d=pd.read_csv(p)
            if len(d):
                for fid,color,band in [(1,'#1c8a54','g'),(2,'#ba4444','r')]:
                    g=d[d.fid.eq(fid)]
                    ax.errorbar(Time(g.mjd.to_numpy(),format='mjd').datetime,g.magpsf_corr,yerr=g.sigmapsf_corr_ext,fmt='.',ms=3,color=color,alpha=.5,label=band)
        for r in visits[visits.name.eq(name)&visits.source.eq('NGPS')].itertuples(index=False):
            ax.axvline(Time(r.mjd,format='mjd').datetime,color='k',ls='--',alpha=.7)
        ax.set_ylabel(f'{name}\ncorrected alert mag');ax.invert_yaxis();ax.grid(alpha=.2)
        if len(ax.get_legend_handles_labels()[0]):ax.legend(fontsize=8,loc='best')
    axes[-1,0].set_xlim(Time('2025-01-01').datetime,Time('2026-11-15').datetime)
    fig.suptitle('Public ZTF alert detections; dashed lines = Palomar visits. Alert absence is not a flux constraint.')
    fig.tight_layout(rect=(0,0,1,.97));return fig


def main():
    folder=OUT/'figures';folder.mkdir(exist_ok=True)
    sample=pd.read_csv(OUT/'pilot_sample.csv')
    with PdfPages(OUT/'pilot_fit_atlas.pdf') as pdf:
        for name in sample.name:
            fig=fit_figure(name);fig.savefig(folder/f'{name}_decomposition.png',dpi=130);pdf.savefig(fig);plt.close(fig)
    for comparison in ['latest_archive','previous_ngps']:
        fig=envelope_figure(comparison);fig.savefig(folder/f'{comparison}_model_envelope.png',dpi=160);plt.close(fig)
    optical_audit()
    fig=optical_figure(['P21422','P24930','P10381','P11530','P9694']);fig.savefig(folder/'refreshed_optical.png',dpi=140);plt.close(fig)
    print('Wrote 15 source fit figures, PDF atlas, model envelopes and optical audit.')


if __name__=='__main__':main()
