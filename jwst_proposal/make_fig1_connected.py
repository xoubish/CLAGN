"""Figure 1: P9694 connects its luminosity history and spectra.

The cutout, time series and spectra are saved measurements. The panels retain all
supplied spectral epochs without target-specific exclusions. Light-curve inputs
under inputs/figure_*.csv are downloaded by fetch_fig1_lightcurves.py.
"""
import json
import hashlib
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, FixedLocator, FuncFormatter, NullFormatter
from matplotlib.transforms import Bbox
from astropy.io import fits
from astropy.time import Time
from scipy.ndimage import gaussian_filter1d
from spherex_data import load_spherex, retained_measurements

HERE=Path(__file__).resolve().parent
INPUTS=HERE/'inputs'
year=lambda m:Time(np.asarray(m,float),format='mjd').decimalyear
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman'],
    'font.size':13,'axes.labelsize':13,'axes.titlesize':13,'xtick.labelsize':13,
    'ytick.labelsize':13,'legend.fontsize':12.5,'pdf.fonttype':42,
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.5})
fig=plt.figure(figsize=(7.2,5.4),facecolor='white')
blue,orange,green,pink,purple='#286296','#e36a30','#298b71','#b74680','#7654a2'
W2_COLOR='#a85a16'
ZTF_YEARS=(2017.0,2027.3)
WISE_YEARS=(2010,2027.3)
LIGHTCURVE_YEARS=(min(ZTF_YEARS[0],WISE_YEARS[0]),max(ZTF_YEARS[1],WISE_YEARS[1]))
OPTICAL_FLUX_LIMITS=(.4,12.0)
WISE_ZERO_JY={'W1':309.540,'W2':171.787}
ALERT_MATCH_DAYS=.02
WISE_CALIBRATION='https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/sec4_4h.html'
inventory=[]

# Vacuum wavelengths in microns, verified against the primary line tables.
LINE_SOURCES={
    'optical':'https://classic.sdss.org/dr6/algorithms/linestable.php',
    'infrared':'https://www.stsci.edu/instruments/nicmos/documents/handbooks/instrument/v5/Appendix_26.html',
}
REST={'Mg II':.2799117,'[Ne V]':.3426850,'[O II]':.3727092,'[Ne III]':.3869860,
      'Hδ':.410289,'Hγ':.434168,'He II':.4687020,
      'Hβ':.486268,'[O III] 4959':.4960295,'[O III] 5007':.5008240,
      'He I 5876':.5877290,'[O I]':.6302046,'Hα':.656461,'[N II]':.6585270,
      '[S II] 6716':.671829,'[S II] 6731':.673267,
      '[S III] 9069':.9071100,'[S III] 9531':.9533200,
      'Paδ':1.0052,'He I':1.0833,'Paγ':1.0941,
      'Paβ':1.2822,'Paα':1.8756,'Brγ':2.1661,'Brβ':2.6259}


def mark_lines(ax,z,groups,fontsize=12,rotation=90):
    """Label expected positions, with offsets only for label readability.

    Each group is (label, member lines, label_x or None[, label_y]). Rotated labels
    hang from the top of the axes; horizontal labels may alternate two rows.
    """
    result=[]
    for group in groups:
        label,members,label_x=group[:3]
        label_y=group[3] if len(group)>3 else .975
        waves=[REST[name]*(1+z) for name in members]
        if not all(ax.get_xlim()[0] < wave < ax.get_xlim()[1] for wave in waves):
            continue
        for wave in waves:
            ax.axvline(wave,ymax=.92,color='#a68143',lw=.5,ls=(0,(3,3)),alpha=.65,zorder=1)
        xpos=np.mean(waves) if label_x is None else label_x
        text=ax.text(xpos,label_y,label,rotation=rotation,va='top',ha='center',fontsize=fontsize,
                transform=ax.get_xaxis_transform(),color='#735625',
                bbox=dict(facecolor='white',edgecolor='none',pad=.15,alpha=.88),zorder=5)
        if label_x is not None:
            # Connect displaced labels to their true wavelength markers.
            bottom=ax.transAxes.inverted().transform(
                text.get_window_extent(fig.canvas.get_renderer()).get_points()[0])[1]
            ax.plot([xpos,np.mean(waves)],[bottom-.02,bottom-.065],
                    transform=ax.get_xaxis_transform(),color='#a68143',lw=.5,zorder=1)
        result.append(dict(label=label,members=members,observed_um=waves,
            label_displaced=label_x is not None))
    return result


def load_target(prefix):
    t=json.loads((INPUTS/f'{prefix}_figure_target.json').read_text())
    epochs=sorted(t['spec']['epochs'],key=lambda e:e['mjd'])
    with fits.open(INPUTS/f'{prefix}_ngps_flux_standards.fits') as hdus:
        d,h=hdus[1].data,hdus[1].header
        assert h['FLUXSTD']=='P330E' and h['TARGET']==t['name']
        for e in epochs:
            if e.get('instrument')=='NGPS':
                e['wave']=d['WAVE_VAC_HELIO_A'].copy()
                e['flux']=np.where(d['MASK'],d['FLUX'],np.nan)
                e['arm_breaks']=np.flatnonzero(d['CHANNEL'][1:]!=d['CHANNEL'][:-1])+1
    return t,epochs


def flux_mjy(e):
    w=np.asarray(e['wave'],float);f=np.asarray(e['flux'],float).copy()
    if e.get('instrument')=='NGPS':
        good=np.flatnonzero(np.isfinite(f));breaks=e.get('arm_breaks',[])
        cuts=np.flatnonzero((np.diff(good)>1)|np.isin(good[1:],breaks))+1
        for run in np.split(good,cuts):
            if len(run)>2:
                f[run]=gaussian_filter1d(f[run],6/(2.35482*np.median(np.diff(w[run]))),mode='nearest')
        for b in breaks:f[b]=np.nan
    return w/1e4,f*1e-17*w*w/2.99792458e18/1e-26


def load_wise(name):
    """Visit medians of AllWISE (2010-2011) and NEOWISE (2014-2024) exposures.

    Both bands use one method: quality frames (AllWISE qi_fact>0, NEOWISE
    qual_frame>0), W1/W2 contamination flags '00', moon flag '00', finite magnitudes
    and positive errors in both bands, separation<3 arcsec, visit=round(MJD/180) and
    at least three exposures. Median magnitude and MJD per visit; approximate median
    error 1.2533*sample_std/sqrt(n). Zero magnitudes from the AllWISE explanatory
    supplement; no colour correction or absolute calibration error added.
    """
    allwise=pd.read_csv(INPUTS/'figure_allwise_mep_exposures.csv',
        dtype={'cc_flags':str,'moon_masked':str})
    allwise=allwise.rename(columns={c:c.replace('_ep','') for c in allwise.columns})
    allwise['quality']=allwise['qi_fact']>0
    allwise['dist_x']=0.0;allwise['survey']='AllWISE'
    neo=pd.read_csv(INPUTS/'figure_neowise_exposures.csv',dtype={'cc_flags':str})
    neo=neo.loc[neo['name'].eq(name)].copy()
    neo['quality']=neo['qual_frame']>0;neo['moon_masked']='0000';neo['survey']='NEOWISE'
    columns=['survey','mjd','w1mpro','w1sigmpro','w2mpro','w2sigmpro','quality','cc_flags',
             'moon_masked','dist_x']
    raw=pd.concat([allwise[columns],neo[columns]],ignore_index=True)
    flags=raw['cc_flags'].str.strip().str.zfill(4);moon=raw['moon_masked'].str.strip().str.zfill(4)
    good=raw['quality']&flags.str[:2].eq('00')&moon.str[:2].eq('00')
    good&=np.isfinite(raw[['mjd','w1mpro','w1sigmpro','w2mpro','w2sigmpro']]).all(axis=1)
    good&=(raw['w1sigmpro']>0)&(raw['w2sigmpro']>0)&(raw['dist_x']<3)
    d=raw.loc[good].copy();d['visit']=np.round(d['mjd']/180.0)
    visits=[]
    for visit,g in d.groupby('visit',sort=True):
        if len(g)<3:continue
        entry=dict(visit=int(visit),survey=g['survey'].iloc[0],mjd=float(g['mjd'].median()),n=len(g))
        for band in ['W1','W2']:
            col=band.lower()+'mpro'
            mag=float(g[col].median());magerr=float(1.2533*g[col].std(ddof=1)/np.sqrt(len(g)))
            flux=WISE_ZERO_JY[band]*1e3*10**(-.4*mag)
            entry[band]=dict(mag=mag,mag_err=magerr,flux_mjy=flux,
                flux_err_mjy=flux*np.log(10)/2.5*magerr)
        visits.append(entry)
    values={band:np.array([[v['mjd'],v[band]['flux_mjy'],v[band]['flux_err_mjy']] for v in visits])
            for band in ['W1','W2']}
    provenance=dict(inputs=['figure_allwise_mep_exposures.csv','figure_neowise_exposures.csv'],
        raw_exposures={k:int(v) for k,v in raw['survey'].value_counts().items()},
        accepted_exposures={k:int(v) for k,v in d['survey'].value_counts().items()},
        visits=visits,zero_magnitude_jy=WISE_ZERO_JY,calibration_source=WISE_CALIBRATION,
        method=load_wise.__doc__.strip())
    return values,provenance


def load_ztf():
    """Nightly medians of the ZTF public release, extended by ALeRCE alert photometry.

    Public release: IRSA light-curve service rows with catflags==0, one median per
    night and band. Alerts: ALeRCE detections with corrected magnitudes (reference
    flux added), not flagged dubious, deduplicated by candid. Alert magnitudes are
    shifted by the median (alert - release) difference of same-exposure matches,
    then only nights after the public release are kept.
    """
    dr=pd.read_csv(INPUTS/'figure_ztf_irsa_lightcurve.csv')
    dr=dr.loc[dr['catflags'].eq(0)]
    alerts=pd.read_csv(INPUTS/'figure_ztf_alerce_detections.csv').drop_duplicates('candid')
    alerts=alerts.loc[alerts['corrected'].astype(bool)&~alerts['dubious'].astype(bool)
                      &np.isfinite(alerts['magpsf_corr'])]
    values,provenance={},{}
    for band,code,fid in [('g','zg',1),('r','zr',2)]:
        d=dr.loc[dr['filtercode'].eq(code)]
        nights=d.groupby(np.floor(d['mjd'])).agg(mjd=('mjd','median'),mag=('mag','median'),
            magerr=('magerr','median'),n=('mag','size'))
        a=alerts.loc[alerts['fid'].eq(fid)].copy()
        nearest=np.abs(d['mjd'].values[None,:]-a['mjd'].values[:,None]).argmin(axis=1)
        matched=np.abs(d['mjd'].values[nearest]-a['mjd'].values)<ALERT_MATCH_DAYS
        offset=float(np.median(a['magpsf_corr'].values[matched]-d['mag'].values[nearest][matched]))
        post=a.loc[a['mjd']>d['mjd'].max()].copy();post['mag']=post['magpsf_corr']-offset
        post=post.groupby(np.floor(post['mjd'])).agg(mjd=('mjd','median'),mag=('mag','median'),
            magerr=('sigmapsf_corr_ext','median'),n=('mag','size'))
        values[band]=dict(release=nights[['mjd','mag','magerr']].to_numpy(),
                          alerts=post[['mjd','mag','magerr']].to_numpy())
        provenance[band]=dict(release_exposures=len(d),release_nights=len(nights),
            release_last_mjd=float(d['mjd'].max()),alert_objects=sorted(a['oid'].unique()),
            alert_matched_exposures=int(matched.sum()),alert_minus_release_offset_mag=offset,
            alert_nights_after_release=len(post),
            alert_last_mjd=float(post['mjd'].max()) if len(post) else None)
    provenance['method']=load_ztf.__doc__.strip()
    return values,provenance


def optical_times(epochs):
    """Use one archival color per year in spectra, date keys and time markers."""
    archival_years=sorted({Time(e['mjd'],format='mjd').datetime.year
                          for e in epochs if e.get('instrument')!='NGPS'})
    palette=['#2878a5','#7654a2']
    colors={y:palette[i % len(palette)] for i,y in enumerate(archival_years)}
    times=[]
    for e in epochs:
        date=Time(e['mjd'],format='mjd').datetime
        recent=e.get('instrument')=='NGPS'
        times.append(dict(mjd=float(e['mjd']),date=date.strftime('%Y-%m-%d'),
            label=str(date.year),
            color=orange if recent else colors[date.year],
            in_lightcurve_window=bool(LIGHTCURVE_YEARS[0]<=year(e['mjd'])<=LIGHTCURVE_YEARS[1])))
    return times


def date_key(ax,entries):
    """Place a compact, color-matched date key just above a spectrum."""
    handles=[Line2D([],[],color=e['color'],lw=2,label=e.get('display_label',e['label'])) for e in entries]
    ax.legend(handles=handles,ncol=len(handles),loc='lower left',
        bbox_to_anchor=(0,1.025,1,0),mode='expand',frameon=False,
        handlelength=.9,handletextpad=.3,columnspacing=.7,borderaxespad=0)


def card():
    prefix='P9694'
    x,w=.012,.976
    t,epochs=load_target(prefix)
    optical_dates=optical_times(epochs)
    ax_im=fig.add_axes([x+.010,.79,.095,.1144])
    ax_im.imshow(plt.imread(INPUTS/f'{prefix}_sdss.jpg'),extent=[-20,20,-20,20])
    ax_im.set(xlim=(-7.5,7.5),ylim=(-7.5,7.5));ax_im.set_axis_off()
    ax_im.plot([1,6],[-5.5,-5.5],color='white',lw=1)
    ax_im.text(3.5,-4.8,'5″',ha='center',color='white',fontsize=12)
    ax_g=fig.add_axes([x+.18,.745,.335,.155])
    ax_w=fig.add_axes([x+.615,.745,w-.633,.155])
    ax_g.set_title('ZTF',loc='left',pad=3)
    ax_w.set_title('WISE',loc='left',pad=3)
    ztf,ztf_provenance=load_ztf()
    for band,c in [('g',green),('r',pink)]:
        a=ztf[band]['release']
        ax_g.errorbar(year(a[:,0]),a[:,1],yerr=a[:,2],fmt='.',ms=.8,lw=.25,color=c,alpha=.65)
        a=ztf[band]['alerts']
        if len(a):ax_g.errorbar(year(a[:,0]),a[:,1],yerr=a[:,2],fmt='o',ms=2.2,lw=.3,
            color=c,mfc='white',mew=.5,alpha=.9)
    ax_g.invert_yaxis();ax_g.set_ylabel('mag',labelpad=1)
    ax_g.text(.02,.50,'g',color=green,va='center',transform=ax_g.transAxes)
    ax_g.text(.085,.50,'r',color=pink,va='center',transform=ax_g.transAxes)
    ax_g.set(xlim=ZTF_YEARS,xticks=[2019,2022,2025])
    wise,wise_provenance=load_wise(prefix)
    for band,c,marker in [('W1',purple,'o'),('W2',W2_COLOR,'s')]:
        a=wise[band]
        ax_w.errorbar(year(a[:,0]),a[:,1],yerr=a[:,2],fmt=marker,ms=1.8,lw=.3,color=c,
            mfc='white',mew=.5)
        # Band labels sit in the AllWISE/NEOWISE gap at each band's 2010-2011 level.
        early=a[year(a[:,0])<2012]
        ax_w.text(2012.6,float(np.mean(early[:,1])),band,color=c,ha='center',va='center',
            bbox=dict(facecolor='white',edgecolor='none',pad=.1,alpha=.85))
    ax_w.set_ylabel('mJy',labelpad=1)
    ax_w.set(xlim=WISE_YEARS,xticks=[2010,2015,2020,2025])
    d,groups,_=load_spherex(prefix);keep=retained_measurements(prefix,d)
    spherex_dates=[dict(label=g['label'],display_label=f'S{i+1}: {g["label"]}',
        visit_id=f'S{i+1}',color=g['color'],
        median_mjd=float(np.median(d['mjds'][g['indices']])),
        start_mjd=g['start_mjd'],end_mjd=g['end_mjd']) for i,g in enumerate(groups)]
    for ax in [ax_g,ax_w]:
        ax.yaxis.set_major_locator(MaxNLocator(2,prune='both'))
        ax.tick_params(pad=1,length=2)
        ax.set_xlabel('Year',labelpad=1)
    # Give each timeline one explicit set of epoch links, avoiding six competing
    # markers in both small panels. Exact dates and visit widths remain measured.
    for obs in optical_dates:
        if obs['in_lightcurve_window']:
            date=float(year(obs['mjd']))
            ax_g.axvline(date,color=obs['color'],lw=1,ls='--',alpha=.8,zorder=.8)
            # Left-align a label near the left edge so it clears the panel title.
            near_left=(date-ZTF_YEARS[0])/(ZTF_YEARS[1]-ZTF_YEARS[0])<.2
            ax_g.text(date+(.15 if near_left else 0),1.07,obs['label'],color=obs['color'],
                ha='left' if near_left else 'center',
                fontsize=11.5,transform=ax_g.get_xaxis_transform())
    # Visit ranges are 14-18 days wide: at this scale a shaded span would be
    # narrower than the dotted line itself, so only the median date is drawn.
    for obs,callout_year in zip(spherex_dates,[2024.6,2026.,2027.4]):
        ax_w.axvline(year(obs['median_mjd']),color=obs['color'],lw=1,ls=':',zorder=.9)
        ax_w.annotate(obs['visit_id'],xy=(year(obs['median_mjd']),.98),
            xycoords=ax_w.get_xaxis_transform(),xytext=(callout_year,1.15),
            textcoords=ax_w.get_xaxis_transform(),ha='center',fontsize=11.5,
            color=obs['color'],annotation_clip=False,
            arrowprops=dict(arrowstyle='-',color=obs['color'],lw=.8,shrinkA=1,shrinkB=1))
    ax_o=fig.add_axes([x+.095,.435,w-.118,.195])
    for e,obs in zip(epochs,optical_dates):
        wave,flux=flux_mjy(e);ngps=e.get('instrument')=='NGPS'
        ax_o.plot(wave,flux,color=obs['color'],lw=.95 if ngps else .8)
    ax_o.set(xlim=(.34,.94),yscale='log',ylim=OPTICAL_FLUX_LIMITS,xticks=[.4,.6,.8])
    ax_o.set_ylabel('Optical\n(mJy)',labelpad=1)
    # Show the identical pooled continuum only within measured SPHEREx coverage.
    mean_fit=json.loads((HERE/'review/p9694_hot_dust/mean_fit.json').read_text())
    mean_path=HERE/mean_fit['provenance']['curve_path']
    assert hashlib.sha256(mean_path.read_bytes()).hexdigest()==mean_fit['provenance']['curve_sha256']
    mean_curve=np.genfromtxt(mean_path,delimiter=',',names=True)
    optical_mean_limits=[float(d['wavelength_um'][keep].min()),ax_o.get_xlim()[1]]
    optical_mean_wave=np.linspace(*optical_mean_limits,200)
    ax_o.plot(optical_mean_wave,np.interp(optical_mean_wave,
        mean_curve['wavelength_observed_um'],mean_curve['total_mjy']),
        color='#242424',lw=1.2,ls=(0,(4,2.5)),zorder=4)
    continuum_overlay=dict(type='pooled mean continuum curve',
        curve=mean_fit['provenance']['curve_path'],observed_wavelength_limits_um=optical_mean_limits,
        flux_scale_factor=1.0,measured_points_shown=False,
        interpretation='Identical total model to lower panel, clipped to optical/SPHEREx overlap; no rescaling or extrapolation below SPHEREx coverage')

    optical_groups=[('Mg II',['Mg II'],.353),('[Ne V]',['[Ne V]'],None),
                     ('[O II]',['[O II]'],.454),('[Ne III]',['[Ne III]'],.487),
                     ('Hδ',['Hδ'],None),('Hγ',['Hγ'],None),('He II',['He II'],.571),
                     ('Hβ',['Hβ'],.594),
                     ('[O III]',['[O III] 4959','[O III] 5007'],.632),
                     ('He I',['He I 5876'],None),('[O I]',['[O I]'],None),
                     ('Hα',['Hα'],None),('[S II]',['[S II] 6716','[S II] 6731'],None)]
    optical_markers=mark_lines(ax_o,t['z'],optical_groups,fontsize=10.5)
    optical_key=list({obs['label']:obs for obs in optical_dates}.values())
    date_key(ax_o,optical_key)
    ax_s=fig.add_axes([x+.095,.105,w-.118,.235])
    for group in groups:
        idx=group['indices'];idx=idx[keep[idx]]
        ax_s.errorbar(d['wavelength_um'][idx],d['flux_mjy'][idx],yerr=d['flux_err_mjy'][idx],
            xerr=d['wavelength_half_width_um'][idx],fmt='o',ms=2.8,lw=.4,
            color=group['color'],mec='white',mew=.25,alpha=.95,zorder=3)
    mean_line,=ax_s.plot(mean_curve['wavelength_observed_um'],mean_curve['total_mjy'],
        color='#242424',lw=1.0,ls=(0,(4,2.5)),zorder=2.5,
        label='Mean continuum across visits')
    ax_s.set(xlim=(.7,5.05),yscale='log',ylim=(.8,16),xticks=[1,2,3,4,5])
    ax_s.set_ylabel('SPHEREx\n(mJy)',labelpad=1)
    ax_s.set_xlabel('Observed wavelength (µm)',labelpad=1)
    infrared_groups=[('Hα',['Hα'],None,.975),('[S III]',['[S III] 9069','[S III] 9531'],None,.80),
                     ('He I, Paγ',['He I','Paγ'],None,.975),('Paβ',['Paβ'],None,.80),
                     ('Paα',['Paα'],None,.80),('Brγ',['Brγ'],None,.975)]
    infrared_markers=mark_lines(ax_s,t['z'],infrared_groups,fontsize=11,rotation=0)
    date_key(ax_s,spherex_dates)
    ax_s.add_artist(ax_s.get_legend())
    ax_s.legend(handles=[mean_line],loc='lower right',fontsize=12,
        frameon=True,facecolor='white',edgecolor='none',framealpha=.88,
        handlelength=2.1,handletextpad=.4,borderpad=.2)
    for ax in [ax_o,ax_s]:
        ax.yaxis.set_major_locator(MaxNLocator(2,prune='upper'));ax.tick_params(pad=1,length=2)
        ax.grid(axis='y',lw=.3,color='#dce0e4');ax.set_axisbelow(True)
    ax_o.yaxis.set_major_locator(FixedLocator([.5,1,2,5,10]))
    ax_o.yaxis.set_major_formatter(FuncFormatter(lambda value,pos:f'{value:g}'))
    ax_o.yaxis.set_minor_formatter(NullFormatter())
    ax_o.tick_params(axis='y',which='minor',length=0)
    ax_s.yaxis.set_major_locator(FixedLocator([1,2,5,10]))
    ax_s.yaxis.set_major_formatter(FuncFormatter(lambda value,pos:f'{value:g}'))
    ax_s.yaxis.set_minor_formatter(NullFormatter())
    ax_s.tick_params(axis='y',which='minor',length=0)
    inventory.append(dict(prefix=prefix,name=t['name'],z=t['z'],optical_epochs=len(epochs),
        optical_dates=optical_dates,optical_date_key=[e['label'] for e in optical_key],
        spherex_dates=spherex_dates,spherex_optical_continuum=continuum_overlay,
        spherex_mean_fit=dict(result='review/p9694_hot_dust/mean_fit.json',
            curve=mean_fit['provenance']['curve_path'],
            temperature_K=mean_fit['baseline']['temperature_K'],
            label='Mean continuum across visits',
            model='Pooled host + fixed-slope disc + blackbody; same continuum at all visits',
            n_fit=mean_fit['baseline']['n_used'],
            interpretation=mean_fit['provenance']['interpretation']),
        ztf=ztf_provenance,wise=wise_provenance,wise_visits_displayed=len(wise['W1']),
        spherex_rows=len(keep),spherex_plotted=int(keep.sum()),
        spherex_intervals=[{k:g[k] for k in ['label','start_mjd','end_mjd']} for g in groups],
        optical_line_markers=optical_markers,spherex_line_markers=infrared_markers))


# A single example gives each measured history and spectrum more room.
card()

# Trim the space vacated by the heading without changing the panel scale.
fig.savefig(HERE/'fig1_connected.pdf',dpi=240,
    bbox_inches=Bbox.from_extents(0,.10,7.2,5.10))
(INPUTS/'fig1_connected_provenance.json').write_text(json.dumps(dict(
    examples=inventory,
    spectral_flux='Observed Fnu mJy; no inter-epoch renormalisation; NGPS P330E',
    layout='Full-proposal-width P9694 figure with enlarged type; optical years above ZTF, labelled SPHEREx visits above WISE, and full-width spectra with matching date keys; thirteen optical and six infrared line labels',
    lightcurve_panel_height_fraction=.155,
    manifold_shown=False,
    optical_flux_limits_mjy=list(OPTICAL_FLUX_LIMITS),
    optical_flux_scale='log', optical_linewidth_pt=dict(archival=.8,ngps=.95),
    spherex_flux_scale='log',spherex_flux_limits_mjy=[.8,16],
    spherex_marker_size_pt=2.8,spherex_marker='o',spherex_marker_alpha=.95,
    spherex_marker_edge=dict(color='white',width_pt=.25),spherex_mean_fit_linestyle='dashed',
    source_heading_shown=False,history_heading_shown=False,outer_box_shown=False,
    ztf_year_limits=list(ZTF_YEARS),wise_year_limits=list(WISE_YEARS),wise_bands=['W1','W2'],
    ztf_markers='Filled dots: nightly medians of the public release. Open circles: ALeRCE alert nights after the public release, shifted onto the release system.',
    wise_markers='One visit series per band from the same AllWISE+NEOWISE method; open circles W1, open squares W2.',
    spectral_time_markers='ZTF: dashed optical-epoch lines with explicit year callouts. WISE: dotted SPHEREx visit median dates connected to S1/S2/S3 callouts; visit ranges (14-18 days) are not shaded because they would be narrower than the line. Colors and visit IDs match the spectral date keys.',
    lightcurve_fetch_script='fetch_fig1_lightcurves.py',
    line_marker_interpretation='Expected redshifted positions, not fitted detections; labels may be displaced for readability',
    line_wavelength_sources=LINE_SOURCES),indent=2)+'\n')
print('Figure 1: P9694 alone; all supplied P9694 spectral epochs retained')
