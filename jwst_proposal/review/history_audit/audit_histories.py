"""Measured history audit, without assigning nuclear luminosities or dust lags.

Run from any directory. Input caches are read-only. Dates are observed frame;
fluxes include host light. Single-ZTF-object curves avoid mixing field offsets.
"""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
SAMPLE = ROOT/'jwst_proposal/inputs/jwst_sample_cycle6.csv'
SOURCES = {}

def record(path):
    SOURCES[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()

def summarize(d):
    if d.empty:
        return {}
    d=d.sort_values('mjd')
    t=d.mjd.to_numpy();m=d.mag.to_numpy()
    early=np.median(m[t<=t.min()+400]);late=np.median(m[t>=t.max()-400])
    recent=t>=t.max()-3*365.25
    slope=np.polyfit((t[recent]-t[recent].min())/365.25,m[recent],1)[0] if sum(recent)>=4 else np.nan
    return dict(n_visits=len(d),first_year=float(Time(t.min(),format='mjd').decimalyear),
                last_year=float(Time(t.max(),format='mjd').decimalyear),
                baseline_rest_yr=float((t.max()-t.min())/365.25/(1+d.z.iloc[0])),
                early_mag=float(early),late_mag=float(late),
                late_minus_early_mag=float(late-early),
                flux_ratio_late_early=float(10**(-.4*(late-early))),
                p95_p05_mag=float(np.percentile(m,95)-np.percentile(m,5)),
                p95_p05_flux_ratio=float(10**(.4*(np.percentile(m,95)-np.percentile(m,5)))),
                peak_to_trough_mag=float(m.max()-m.min()),
                late_slope_mag_per_observed_yr=float(slope),
                typical_visit_mag_se=float(d.mag_se.median()),
                brightest_year=float(Time(t[np.argmin(m)],format='mjd').decimalyear),
                faintest_year=float(Time(t[np.argmax(m)],format='mjd').decimalyear))

def wise_curves(sample):
    files=[ROOT/'data/neowise_cache/jwstpool/batch_00000.csv',
           ROOT/'data/neowise_cache/zeltyn/batch_00000.csv']
    files += [ROOT/f'data/neowise_cache/poolall/batch_{n:05d}.csv'
              for n in [1900,9000,13900,12600,16300,14600,4100,9800,13700]]
    frames=[]
    wanted=set(sample.internal_id)
    for p in files:
        d=pd.read_csv(p,dtype={'cc_flags':str},low_memory=False)
        d=d[d.name.isin(wanted)].copy()
        if len(d):record(p);d['source']=str(p.relative_to(ROOT));frames.append(d)
    allraw=pd.concat(frames,ignore_index=True)
    rows=[]
    for target in sample.itertuples():
        raw=allraw[allraw.name.eq(target.internal_id)].copy()
        pos=SkyCoord(target.ra*u.deg,target.dec*u.deg)
        raw['sep']=SkyCoord(raw.ra.to_numpy()*u.deg,raw.dec.to_numpy()*u.deg).separation(pos).arcsec
        raw=raw[(raw.qual_frame>0)&(raw.sep<3)].sort_values('mjd').drop_duplicates(['mjd','ra','dec'])
        flags=raw.cc_flags.str.strip().str.zfill(4)
        for j,band in enumerate(['W1','W2']):
            mc=f'w{j+1}mpro';ec=f'w{j+1}sigmpro'
            d=raw[(flags.str[j]=='0')&np.isfinite(raw[mc])&np.isfinite(raw[ec])&(raw[ec]>0)].copy()
            d['season']=(d.mjd.diff().fillna(0)>60).cumsum()
            for _,g in d.groupby('season'):
                if len(g)<3:continue
                rows.append(dict(id=target.id,internal_id=target.internal_id,band=band,z=target.z,
                    mjd=float(g.mjd.median()),mag=float(g[mc].median()),
                    mag_se=float(1.2533*g[mc].std(ddof=1)/np.sqrt(len(g))),
                    exposure_magerr_median=float(g[ec].median()),n_exposures=len(g),
                    median_sep_arcsec=float(g.sep.median()),source=g.source.iloc[0],oid=''))
    return pd.DataFrame(rows)

def optical_curves(sample):
    parent_path=ROOT/'data/master_list_scored.csv';record(parent_path)
    parent=pd.read_csv(parent_path,low_memory=False).dropna(subset=['ra','dec'])
    idx,sep,_=SkyCoord(sample.ra.to_numpy()*u.deg,sample.dec.to_numpy()*u.deg).match_to_catalog_sky(
        SkyCoord(parent.ra.to_numpy()*u.deg,parent.dec.to_numpy()*u.deg))
    rows=[];inventory=[]
    for i,target in enumerate(sample.itertuples()):
        names=list(dict.fromkeys([target.internal_id,target.target]+([parent.iloc[idx[i]]['name']] if sep[i].arcsec<1.5 else [])))
        paths=[OUT/'raw_ztf'/f'{target.id}.csv']
        paths += [ROOT/'data/ztf_cache'/tag/f'{name}.csv' for tag in ['review_dr24_20260920','calib','zeltyn','pool','v2'] for name in names]
        raw=None;path=None
        for p in paths:
            if not p.exists() or p.stat().st_size<10:continue
            try:d=pd.read_csv(p,dtype={'oid':str,'expid':str},low_memory=False)
            except pd.errors.EmptyDataError:continue
            if {'mjd','mag','filtercode','ra','dec'}.issubset(d) and len(d):raw=d;path=p;break
        if raw is None:
            inventory.append(dict(id=target.id,aliases=';'.join(names),status='no local usable optical curve'));continue
        record(path)
        good=np.isfinite(raw.mjd)&np.isfinite(raw.mag)&np.isfinite(raw.magerr)&(raw.magerr>0)&np.isfinite(raw.ra)&np.isfinite(raw.dec)
        raw=raw[good].copy()
        raw=raw[(raw.catflags.fillna(32768).astype(int)&32768)==0]
        raw['sep']=SkyCoord(raw.ra.to_numpy()*u.deg,raw.dec.to_numpy()*u.deg).separation(SkyCoord(target.ra*u.deg,target.dec*u.deg)).arcsec
        raw=raw[raw.sep<1.5].copy()
        raw['season']=np.floor(raw.mjd/90).astype(int)
        for fc,band in [('zg','g'),('zr','r')]:
            b=raw[raw.filtercode.eq(fc)].copy()
            if b.empty:continue
            ranked=[]
            for oid,g in b.groupby('oid'):
                nb=int((g.groupby('season').size()>=3).sum())
                ranked.append((nb,g.mjd.max()-g.mjd.min(),len(g),oid))
            oid=max(ranked)[-1]
            b=b[b.oid.eq(oid)].sort_values('magerr').drop_duplicates('expid').sort_values('mjd')
            inventory.append(dict(id=target.id,aliases=';'.join(names),status='available',band=band,
                source=str(path.relative_to(ROOT)),oid=oid,n_oids=len(ranked),n_exposures=len(b),
                first_mjd=b.mjd.min(),last_mjd=b.mjd.max(),max_sep_arcsec=b.sep.max()))
            for _,g in b.groupby('season'):
                if len(g)<3:continue
                rows.append(dict(id=target.id,internal_id=target.internal_id,band=band,z=target.z,
                    mjd=float(g.mjd.median()),mag=float(g.mag.median()),
                    mag_se=float(1.2533*g.mag.std(ddof=1)/np.sqrt(len(g))),
                    exposure_magerr_median=float(g.magerr.median()),n_exposures=len(g),
                    median_sep_arcsec=float(g.sep.median()),source=str(path.relative_to(ROOT)),oid=oid))
    pd.DataFrame(inventory).to_csv(OUT/'optical_inventory.csv',index=False)
    return pd.DataFrame(rows)

def plot_atlas(sample,curves):
    colors={'W1':'#782ba4','W2':'#c67513','g':'#008562','r':'#c34054'}
    with PdfPages(OUT/'history_atlas.pdf') as pdf:
        for page in range(4):
            fig,axes=plt.subplots(3,2,figsize=(11.7,9.2),sharex=True)
            for ax,target in zip(axes.flat,list(sample.itertuples())[page*6:page*6+6]):
                d=curves[curves.id.eq(target.id)]
                for band in ['g','r','W1','W2']:
                    g=d[d.band.eq(band)].sort_values('mjd')
                    if g.empty:continue
                    # Relative flux to last 400-day median in each band; not host-subtracted.
                    ref=g[g.mjd>=g.mjd.max()-400].mag.median()
                    ax.plot(Time(g.mjd.to_numpy(),format='mjd').decimalyear,
                            10**(-.4*(g.mag-ref)),'.-',color=colors[band],lw=1,ms=3,label=band)
                ax.axhline(1,color='.65',ls=':',lw=.8)
                ax.set_title(f'{target.id}  {target.internal_id}  | selected {target.family}',fontsize=10)
                ax.set_yscale('log');ax.grid(alpha=.15);ax.legend(fontsize=7,ncol=4,loc='best')
                ax.set_ylabel('Total flux / late median');ax.set_xlim(2013.6,2026.2)
            for ax in axes[-1]:ax.set_xlabel('Observed calendar year')
            fig.suptitle('Measured optical and infrared histories — host light included',fontsize=14)
            fig.text(.5,.018,'Bands normalized separately. Lines join observations; gaps are not measurements. No dust-response model is fitted.',ha='center',fontsize=9)
            fig.tight_layout(rect=[0,.035,1,.96]);pdf.savefig(fig);fig.savefig(OUT/f'histories_{page+1}.png',dpi=140);plt.close(fig)

def main():
    OUT.mkdir(parents=True,exist_ok=True);record(SAMPLE)
    sample=pd.read_csv(SAMPLE)
    wise=wise_curves(sample);optical=optical_curves(sample)
    curves=pd.concat([wise,optical],ignore_index=True).sort_values(['id','band','mjd'])
    curves.to_csv(OUT/'binned_measurements.csv',index=False)
    rows=[]
    for target in sample.itertuples():
        row=dict(id=target.id,internal_id=target.internal_id,family=target.family,z=target.z)
        for band in ['W1','W2','g','r']:
            stats=summarize(curves[(curves.id==target.id)&(curves.band==band)])
            row.update({f'{band}_{k}':v for k,v in stats.items()})
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT/'history_summary.csv',index=False)
    plot_atlas(sample,curves)
    (OUT/'provenance.json').write_text(json.dumps(dict(inputs=SOURCES,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        note='Measured total-source photometry; not nuclear luminosity, event dating, dust-lag measurement, or MIRI forecast.'),indent=2)+'\n')
    print(pd.DataFrame(rows)[['id','family','W1_p95_p05_flux_ratio','W1_late_minus_early_mag','W1_late_slope_mag_per_observed_yr','g_late_minus_early_mag','r_late_minus_early_mag']].round(3).to_string(index=False))

if __name__=='__main__':main()
