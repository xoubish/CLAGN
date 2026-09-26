"""Figure 1: P1823 connects its luminosity history and spectra.

The cutout, time series and spectra are saved measurements. The panels retain all
supplied spectral epochs, with the previously requested P1823 point exclusion.
"""
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Bbox
from astropy.io import fits
from astropy.time import Time
from scipy.ndimage import gaussian_filter1d
from spherex_data import load_spherex, retained_measurements

HERE=Path(__file__).resolve().parent
INPUTS=HERE/'inputs'
year=lambda m:Time(np.asarray(m,float),format='mjd').decimalyear
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman'],
    'font.size':12,'axes.labelsize':12,'axes.titlesize':12,'xtick.labelsize':12,
    'ytick.labelsize':12,'legend.fontsize':12,'pdf.fonttype':42,
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.5})
fig=plt.figure(figsize=(6.5,5.4),facecolor='white')
blue,orange,green,pink,purple='#286296','#e36a30','#298b71','#b74680','#7654a2'
W2_COLOR='#a85a16'
LIGHTCURVE_YEARS=(2015,2028)
OPTICAL_FLUX_LIMITS=(0,.6)
W2_ZERO_JY=171.787
WISE_CALIBRATION='https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/sec4_4h.html'
inventory=[]

# Vacuum wavelengths in microns, verified against the primary line tables.
LINE_SOURCES={
    'optical':'https://classic.sdss.org/dr6/algorithms/linestable.php',
    'infrared':'https://www.stsci.edu/instruments/nicmos/documents/handbooks/instrument/v5/Appendix_26.html',
}
REST={'Mg II':.2799117,'[O II]':.3727092,'Hγ':.434168,
      'Hβ':.486268,'[O III] 4959':.4960295,'[O III] 5007':.5008240,
      'Hα':.656461,'Paδ':1.0052,'He I':1.0833,'Paγ':1.0941,
      'Paβ':1.2822,'Paα':1.8756,'Brγ':2.1661,'Brβ':2.6259}


def mark_lines(ax,z,groups):
    """Label expected positions, with offsets only for label readability."""
    result=[]
    for label,members,label_x in groups:
        waves=[REST[name]*(1+z) for name in members]
        if not all(ax.get_xlim()[0] < wave < ax.get_xlim()[1] for wave in waves):
            continue
        for wave in waves:
            ax.axvline(wave,ymax=.92,color='#a68143',lw=.5,ls=(0,(3,3)),alpha=.65,zorder=1)
        xpos=np.mean(waves) if label_x is None else label_x
        text=ax.text(xpos,.975,label,rotation=90,va='top',ha='center',
                transform=ax.get_xaxis_transform(),color='#735625',
                bbox=dict(facecolor='white',edgecolor='none',pad=.15,alpha=.88),zorder=5)
        if label_x is not None:
            # Connect displaced labels to their true wavelength markers.
            bottom=ax.transAxes.inverted().transform(
                text.get_window_extent(fig.canvas.get_renderer()).get_points()[0])[1]
            ax.plot([xpos,np.mean(waves)],[bottom-.02,bottom-.065],
                    transform=ax.get_xaxis_transform(),color='#a68143',lw=.5,zorder=1)
        result.append(dict(label=label,members=members,observed_um=waves))
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


def load_w2(name):
    """Visit medians from cached NEOWISE exposures, with scatter-based errors.

    Match the existing W1 visit grouping and quality cuts, additionally requiring
    a finite, positive W2 uncertainty so upper limits are not used as detections.
    """
    raw=pd.read_csv(INPUTS/'figure_neowise_exposures.csv',dtype={'cc_flags':str})
    raw=raw.loc[raw['name'].eq(name)].copy()
    flags=raw['cc_flags'].str.strip().str.zfill(4)
    good=(raw['qual_frame']>0)&flags.str[:2].eq('00')
    good&=np.isfinite(raw[['mjd','w1mpro','w1sigmpro','w2mpro','w2sigmpro']]).all(axis=1)
    good&=(raw['w2sigmpro']>0)&(raw['dist_x']<3)
    d=raw.loc[good].copy();d['visit']=np.round(d['mjd']/180.0)
    visits=[]
    for visit,g in d.groupby('visit',sort=True):
        if len(g)<3:continue
        mag=float(g['w2mpro'].median())
        magerr=float(1.2533*g['w2mpro'].std(ddof=1)/np.sqrt(len(g)))
        flux=W2_ZERO_JY*1e3*10**(-.4*mag)
        visits.append(dict(visit=int(visit),mjd=float(g['mjd'].median()),
            mag=mag,mag_err=magerr,flux_mjy=flux,
            flux_err_mjy=flux*np.log(10)/2.5*magerr,n=len(g)))
    values=np.array([[v['mjd'],v['flux_mjy'],v['flux_err_mjy']] for v in visits])
    return values,dict(raw_exposures=len(raw),accepted_exposures=int(good.sum()),visits=visits)


def optical_times(epochs):
    """Use one archival color per year in spectra, date keys and time markers."""
    archival_years=sorted({Time(e['mjd'],format='mjd').datetime.year
                          for e in epochs if e.get('instrument')!='NGPS'})
    colors={y:to_hex(plt.cm.Blues(.35+.55*i/max(1,len(archival_years)-1)))
            for i,y in enumerate(archival_years)}
    times=[]
    for e in epochs:
        date=Time(e['mjd'],format='mjd').datetime
        recent=e.get('instrument')=='NGPS'
        times.append(dict(mjd=float(e['mjd']),date=date.strftime('%Y-%m-%d'),
            label=date.strftime('%d %b %Y') if recent else str(date.year),
            color=orange if recent else colors[date.year],
            in_lightcurve_window=bool(LIGHTCURVE_YEARS[0]<=year(e['mjd'])<=LIGHTCURVE_YEARS[1])))
    return times


def date_key(ax,entries):
    """Place a compact, color-matched date key just above a spectrum."""
    handles=[Line2D([],[],color=e['color'],lw=1.6,label=e['label']) for e in entries]
    ax.legend(handles=handles,ncol=len(handles),loc='lower left',
        bbox_to_anchor=(0,1.025,1,0),mode='expand',frameon=False,
        handlelength=.9,handletextpad=.3,columnspacing=.7,borderaxespad=0)


def card():
    prefix='P1823'
    x,w=.012,.976
    t,epochs=load_target(prefix)
    optical_dates=optical_times(epochs)
    ax_im=fig.add_axes([x+.024,.755,.105,.1264])
    ax_im.imshow(plt.imread(INPUTS/f'{prefix}_sdss.jpg'),extent=[-20,20,-20,20])
    ax_im.set(xlim=(-7.5,7.5),ylim=(-7.5,7.5));ax_im.set_axis_off()
    ax_im.plot([1,6],[-5.5,-5.5],color='white',lw=1)
    ax_im.text(3.5,-4.8,'5″',ha='center',color='white',fontsize=12)
    ax_g=fig.add_axes([x+.22,.75,.30,.135])
    ax_w=fig.add_axes([x+.64,.75,w-.659,.135])
    ax_g.set_title('ZTF',loc='left',pad=3)
    ax_w.set_title('WISE',loc='left',pad=3)
    for band,c in [('g',green),('r',pink)]:
        a=np.asarray(t['ztf'][band],float)
        ax_g.errorbar(year(a[:,0]),a[:,1],yerr=a[:,2],fmt='.',ms=.8,lw=.25,color=c,alpha=.65)
    ax_g.invert_yaxis();ax_g.set_ylabel('mag',labelpad=1)
    ax_g.text(.02,.08,'g',color=green,transform=ax_g.transAxes)
    ax_g.text(.12,.08,'r',color=pink,transform=ax_g.transAxes)
    ax_g.set(xlim=LIGHTCURVE_YEARS,xticks=[2015,2020,2025])
    for a,face in [(t['wise']['W1'],purple),(t.get('neo',[]),'white')]:
        a=np.asarray(a,float)
        a=a[year(a[:,0])>=LIGHTCURVE_YEARS[0]]
        if len(a):ax_w.errorbar(year(a[:,0]),a[:,1],yerr=a[:,2],fmt='o',ms=1.8,lw=.3,color=purple,mfc=face,mew=.5)
    w2,w2_provenance=load_w2(prefix)
    w2=w2[year(w2[:,0])>=LIGHTCURVE_YEARS[0]]
    ax_w.errorbar(year(w2[:,0]),w2[:,1],yerr=w2[:,2],fmt='s',ms=1.8,lw=.3,
        color=W2_COLOR,mfc='white',mew=.5)
    ax_w.set_ylabel('mJy',labelpad=1)
    ax_w.set(xlim=LIGHTCURVE_YEARS,xticks=[2015,2020,2025])
    for band,c,yy in [('W1',purple,.10),('W2',W2_COLOR,.64)]:
        ax_w.text(.98,yy,band,color=c,ha='right',transform=ax_w.transAxes,
            bbox=dict(facecolor='white',edgecolor='none',pad=.1,alpha=.85))
    d,groups,_=load_spherex(prefix);keep=retained_measurements(prefix,d)
    spherex_dates=[dict(label=g['label'],color=g['color'],
        median_mjd=float(np.median(d['mjds'][g['indices']])),
        start_mjd=g['start_mjd'],end_mjd=g['end_mjd']) for g in groups]
    for ax in [ax_g,ax_w]:
        ax.yaxis.set_major_locator(MaxNLocator(2,prune='both'))
        ax.tick_params(pad=1,length=2)
        ax.set_xlabel('Year',labelpad=1)
        for obs in optical_dates:
            if obs['in_lightcurve_window']:
                ax.axvline(year(obs['mjd']),color=obs['color'],lw=.65,
                    ls='--',alpha=.65,zorder=.8)
        for obs in spherex_dates:
            ax.axvspan(year(obs['start_mjd']),year(obs['end_mjd']),
                color=obs['color'],alpha=.14,lw=0,zorder=.7)
            ax.axvline(year(obs['median_mjd']),color=obs['color'],lw=.85,
                ls=':',zorder=.9)
    ax_o=fig.add_axes([x+.095,.435,w-.118,.195])
    for e,obs in zip(epochs,optical_dates):
        wave,flux=flux_mjy(e);ngps=e.get('instrument')=='NGPS'
        ax_o.plot(wave,flux,color=obs['color'],lw=.55 if ngps else .35)
    ax_o.set(xlim=(.34,.94),ylim=OPTICAL_FLUX_LIMITS,xticks=[.4,.6,.8])
    ax_o.set_ylabel('Optical\n(mJy)',labelpad=1)
    optical_groups=[('Mg II',['Mg II'],None),('[O II]',['[O II]'],None),
                     ('Hγ',['Hγ'],None),('Hβ',['Hβ'],.768),
                     ('[O III]',['[O III] 4959','[O III] 5007'],.84)]
    optical_markers=mark_lines(ax_o,t['z'],optical_groups)
    optical_key=list({obs['label']:obs for obs in optical_dates}.values())
    date_key(ax_o,optical_key)
    ax_s=fig.add_axes([x+.095,.105,w-.118,.235])
    for group in groups:
        idx=group['indices'];idx=idx[keep[idx]]
        ax_s.errorbar(d['wavelength_um'][idx],d['flux_mjy'][idx],yerr=d['flux_err_mjy'][idx],
            xerr=d['wavelength_half_width_um'][idx],fmt='.',ms=1.5,lw=.3,color=group['color'],alpha=.8)
    ax_s.set(xlim=(.7,5.05),ylim=(0,2.8),xticks=[1,2,3,4,5])
    ax_s.set_ylabel('SPHEREx\n(mJy)',labelpad=1)
    ax_s.set_xlabel('Observed wavelength (µm)',labelpad=1)
    infrared_groups=[('Hβ/[O III]',['Hβ','[O III] 4959','[O III] 5007'],.80),
                     ('Hα',['Hα'],None),('Paδ',['Paδ'],1.50),
                     ('He I/Paγ',['He I','Paγ'],1.82),
                     ('Paβ',['Paβ'],None),('Paα',['Paα'],None),
                     ('Brγ',['Brγ'],None),('Brβ',['Brβ'],None)]
    infrared_markers=mark_lines(ax_s,t['z'],infrared_groups)
    date_key(ax_s,spherex_dates)
    for ax in [ax_o,ax_s]:
        ax.yaxis.set_major_locator(MaxNLocator(2,prune='upper'));ax.tick_params(pad=1,length=2)
        ax.grid(axis='y',lw=.3,color='#dce0e4');ax.set_axisbelow(True)
    ax_o.set_yticks([0,.2,.4,.6])
    inventory.append(dict(prefix=prefix,name=t['name'],z=t['z'],optical_epochs=len(epochs),
        optical_dates=optical_dates,optical_date_key=[e['label'] for e in optical_key],
        spherex_dates=spherex_dates,
        w2=w2_provenance,w2_visits_displayed=len(w2),
        spherex_rows=len(keep),spherex_plotted=int(keep.sum()),
        spherex_intervals=[{k:g[k] for k in ['label','start_mjd','end_mjd']} for g in groups],
        optical_line_markers=optical_markers,spherex_line_markers=infrared_markers))


# A single example gives each measured history and spectrum more room.
card()

# Trim the space vacated by the heading without changing the panel scale.
for ext in ['pdf','png']:
    fig.savefig(HERE/f'fig1_connected.{ext}',dpi=240,
        bbox_inches=Bbox.from_extents(0,.10,6.5,5.05))
(INPUTS/'fig1_connected_provenance.json').write_text(json.dumps(dict(
    examples=inventory,
    spectral_flux='Observed Fnu mJy; no inter-epoch renormalisation; NGPS P330E',
    layout='Single P1823 example: light curves above full-width optical and SPHEREx spectra',
    manifold_shown=False,
    optical_flux_limits_mjy=list(OPTICAL_FLUX_LIMITS),
    source_heading_shown=False,history_heading_shown=False,outer_box_shown=False,
    lightcurve_year_limits=list(LIGHTCURVE_YEARS),wise_bands=['W1','W2'],
    spectral_time_markers='Both light curves: dashed lines at every optical epoch within the displayed window; dotted lines at each SPHEREx visit median MJD and shaded visit date ranges. Colors match the date keys above the spectra. The labelled 2002 optical spectrum predates the displayed window.',
    wise_w2_input='figure_neowise_exposures.csv',
    wise_w2_method='NEOWISE visit=round(MJD/180); qual_frame>0; first two cc_flags=00; finite W1/W2 magnitudes and errors; positive W2 error; separation<3 arcsec; at least 3 exposures. Median magnitude and MJD; approximate median error=1.2533*sample_std(mag)/sqrt(n). Convert flux and propagate statistical error; no color correction or absolute calibration error added.',
    wise_w2_zero_mag_jy=W2_ZERO_JY,wise_calibration_source=WISE_CALIBRATION,
    line_marker_interpretation='Expected redshifted positions, not fitted detections; labels may be displaced for readability',
    line_wavelength_sources=LINE_SOURCES),indent=2)+'\n')
print('Figure 1: P1823 alone; all supplied P1823 spectral epochs retained')
