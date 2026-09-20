"""Dated, descriptive evidence for the September review; no CL classifications.

All outputs are local/private. Photometric slopes are archival measurements,
not extrapolations to the observing date. The spectral indices are screening
measurements, not host/Fe-II-subtracted broad-line fits.
"""
from pathlib import Path
import importlib,json,argparse
import numpy as np
import pandas as pd
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'data/reselection_2026-09-20'
WEB=importlib.import_module('16_candidate_webpage')

def trend(points,prefix):
    a=np.asarray(points,float)
    if not len(a):return {}
    a=a[np.isfinite(a).all(axis=1)&(a[:,1]>0)]
    if len(a)<3:return {}
    a=a[np.argsort(a[:,0])]
    mag=-2.5*np.log10(a[:,1]/309540)
    recent=a[:,0]>=a[-1,0]-3*365.25
    # Each WISE visit is one independent epoch here.
    slope=np.polyfit((a[recent,0]-a[-1,0])/365.25,mag[recent],1)[0] if recent.sum()>=4 else np.nan
    return {prefix+'_last_date':Time(a[-1,0],format='mjd').strftime('%Y-%m-%d'),
            prefix+'_slope_mag_year':slope,prefix+'_amplitude':np.percentile(mag,95)-np.percentile(mag,5),
            prefix+'_end_minus_start':np.median(mag[-3:])-np.median(mag[:3])}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--screened',action='store_true');args=ap.parse_args()
    prefix='sep23_screened' if args.screened else 'sep23'
    snapshot=json.loads((OUT/('sep23_screened_snapshot.json' if args.screened else 'sep23_review_snapshot.json')).read_text())
    optics=pd.read_csv(ROOT/'data/optical_change_pool.csv').set_index('name')
    rows=[];objects=[]
    for t in snapshot['targets']:
        records=[]
        for path in [ROOT/'data/spectra_dl'/f"{t['name']}.json",OUT/'sdssv_spectra'/f"{t['name']}.json"]:
            if path.exists():records+=json.loads(path.read_text())
        spec=WEB.spectrum_payload(records)
        t['spec']=spec
        ztf,status=WEB.review_ztf(t['name']);t['ztf']=ztf;t['ztf_status']=status
        t['cut']=None # evidence export does not need embedded images
        row=dict(name=t['name'],ra=t['ra'],dec=t['dec'],z=t['z'],r_historical=t['rmag'],region=t['region'],
                 on_fraction=t['on_fraction'],off_fraction=t['off_fraction'],known_status=t['status'],
                 spectral_dates=t['n_spec'],loaded_dates=len({r['epoch_day'] for r in spec['epochs']}),
                 sdssv_dates=t['sdssv_dates'],ztf_last_date=status['last_date'])
        row.update(trend(t['wise'].get('W1',[]),'unwise'));row.update(trend(t['neo'],'neowise'))
        for band in ['g','r']:
            a=np.asarray(ztf.get(band,[]),float)
            if len(a)<4:continue
            tail=a[a[:,0]>=a[:,0].max()-180]
            recent=a[a[:,0]>=a[:,0].max()-730]
            row[f'ztf_{band}_last180']=float(np.median(tail[:,1]))
            row[f'ztf_{band}_amplitude']=float(np.percentile(a[:,1],95)-np.percentile(a[:,1],5))
            row[f'ztf_{band}_slope_mag_year']=float(np.polyfit((recent[:,0]-recent[-1,0])/365.25,recent[:,1],1)[0]) if len(recent)>=8 and np.ptp(recent[:,0])>180 else np.nan
        if t['name'] in optics.index:
            o=optics.loc[t['name']]
            for key in ['syn_g','syn_r','badspec','mjd_spec']:row[key]=o[key]
            for band in ['g','r']:
                row[f'delta_{band}_latest_minus_spectrum']=row.get(f'ztf_{band}_last180',np.nan)-o[f'syn_{band}']
        rows.append(row);objects.append(t)
    pd.DataFrame(rows).to_csv(OUT/f'{prefix}_photometric_evidence.csv',index=False)
    (OUT/f'{prefix}_evidence_objects.json').write_text(json.dumps(WEB.native(objects),allow_nan=False))
    # Every object is inspected with actual dates and both Balmer regions.
    for page in range((len(objects)+6)//7):
        subset=objects[page*7:(page+1)*7]
        fig,axes=plt.subplots(len(subset),4,figsize=(19,2.7*len(subset)),squeeze=False)
        for axs,t in zip(axes,subset):
            for band,color in [('g','green'),('r','firebrick')]:
                a=np.asarray(t['ztf'].get(band,[]),float)
                if len(a):axs[0].plot(Time(a[:,0],format='mjd').decimalyear,a[:,1],'.-',ms=2,lw=.6,color=color,label=band)
            axs[0].invert_yaxis();axs[0].set_title(f"{t['name']}  z={t['z']:.3f}  r={t['rmag']:.2f}",fontsize=10)
            for points,marker,label in [(t['wise'].get('W1',[]),'o','unWISE'),(t['neo'],'s','NEOWISE')]:
                a=np.asarray(points,float)
                if len(a):
                    good=a[:,1]>0;a=a[good]
                    axs[1].plot(Time(a[:,0],format='mjd').decimalyear,-2.5*np.log10(a[:,1]/309540),marker+'-',ms=2,lw=.7,label=label)
            axs[1].invert_yaxis();axs[1].legend(fontsize=6);axs[1].set_title(t['region'],fontsize=9)
            for e in t['spec']['epochs']:
                w=np.asarray(e['wave'])/(1+t['z']);f=np.asarray(e['flux'],float)
                for ax,(lo,hi) in zip(axs[2:],[(4650,5150),(6300,6850)]):
                    sel=(w>=lo)&(w<=hi)&np.isfinite(f)
                    if sel.sum()<5:continue
                    ax.plot(w[sel],f[sel],lw=.65,label=e['label']+(' !' if e['quality_flag'] else ''))
                for ax in axs[:2]:
                    if e['mjd']:ax.axvline(Time(e['mjd'],format='mjd').decimalyear,color='.7',lw=.4,alpha=.6)
            for ax,(lo,hi) in zip(axs[2:],[(4650,5150),(6300,6850)]):
                ax.set_xlim(lo,hi);ax.legend(fontsize=5,loc='upper right')
                values=np.concatenate([line.get_ydata() for line in ax.lines]) if ax.lines else np.array([])
                if len(values):
                    low,high=np.nanpercentile(values,[2,98]);ax.set_ylim(min(0,low),high*1.1 if high>0 else 1)
            for ax in axs:ax.tick_params(labelsize=7);ax.grid(alpha=.15)
        axes[0,0].set_xlabel('ZTF: calendar year');axes[0,1].set_xlabel('W1: calendar year')
        axes[0,2].set_title('Hβ / [O III] — rest Å');axes[0,3].set_title('Hα — rest Å')
        fig.suptitle('September review — raw archival spectra; ! = metadata flag; no state classifications',fontsize=12)
        fig.tight_layout(rect=(0,0,1,.98));fig.savefig(OUT/f'{prefix}_evidence_{page+1}.png',dpi=115);plt.close(fig)
    print(pd.DataFrame(rows).to_string(index=False),flush=True)

if __name__=='__main__':main()
