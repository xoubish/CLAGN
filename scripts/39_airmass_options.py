"""Extend the prepared, fully cached sample with exposure-aware airmass options.

No new objects are admitted and no archive selection is overwritten. Backups
may be shared between nights. The separate full-parent search remains intact.
"""
from pathlib import Path
from datetime import timezone
from concurrent.futures import ProcessPoolExecutor
import importlib, json, math, hashlib
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz, get_body
from astropy.utils import iers
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Patch
from zoneinfo import ZoneInfo

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/reselection_2026-09-20'
OBS=importlib.import_module('05_observability')
TZ=ZoneInfo('America/Los_Angeles')
TIERS=[('preferred',1.5,20),('extended',1.8,40),('fallback',2.0,50)]
COLORS={'preferred':'#28618a','extended':'#77aaca','fallback':'#d2e3ee'}
iers.conf.auto_download=False; iers.conf.auto_max_age=None

def exposure(target):
    model=importlib.import_module('22_september_etc')
    ref=model.reference(target)
    if ref is None: return target['name'], {}
    mjd,mag,private,flux=ref
    signature=hashlib.sha256(json.dumps([target['z'],ref,'adopted-2x300-slit1.5-2x3-zenith-seeing-sky18.5-v2']).encode()).hexdigest()
    cache=OUT/'airmass_exposure_cache';cache.mkdir(exist_ok=True)
    path=cache/f"{target['name']}.json"
    if path.exists():
        value=json.loads(path.read_text())
        if value['signature']==signature:
            for plan in value['plans'].values():plan.setdefault('status','scenario')
            return target['name'],value['plans']
    for ch,rn,scale in zip(model.CFG.channels,[2.8,7.8,3.7,4.6],[.193,.189,.186,.186]):
        model.CFG.readnoise[ch]=rn*u.count/u.pix
        model.CFG.platescale[ch]=scale*u.arcsec/u.pix
    wave=4862.7*(1+target['z'])/10
    channels=[ch for ch in ['R','I','G','U'] if model.CFG.channelRange[ch][0].to_value(u.nm)<=wave-4 and model.CFG.channelRange[ch][1].to_value(u.nm)>=wave+4]
    if not channels:return target['name'],{}
    ch=channels[0];plans={}
    bin_A=float((model.CFG.dLambda[ch]*3).to_value(u.AA))
    for tier,x,_ in TIERS:
        # Adopted setting (2026-09-21): 2x300 s, 1.5 arcsec slit, 2x3 binning, 16-minute visit.
        # Report the predicted continuum S/N per Angstrom at the tier's airmass ceiling with
        # seeing 1.3 arcsec at zenith scaled by airmass^0.6 and sky V=18.5.
        cmd=[ch,str(wave-4),str(wave+4),'EXPTIME','300',
             '-slit','SET','1.5','-binspect','3','-binspat','2','-seeing','1.3','500',
             '-airmass',str(x),'-skymag','18.5','-mag',str(mag),
             '-magsystem','AB','-magfilter','match','-noslicer']
        args=model.ETC.parser.parse_args(cmd);model.ETC.check_inputs_add_units(args)
        snr=float(model.ETC.main(args,quiet=True)['SNR'].value)*np.sqrt(2)/np.sqrt(bin_A)
        plans[tier]=dict(exposures=2,seconds_each=300,integration_minutes=10.,visit_minutes=16,airmass=x,
                         slit_arcsec=1.5,binspat=2,binspect=3,predicted_snr_per_angstrom=snr,goal_snr=5,
                         reference_mjd=mjd,reference_private=private,continuum_AB=mag,channel=ch,
                         status='scenario' if snr>=5 else 'below the S/N floor at 2x300 s')
    path.write_text(json.dumps(dict(signature=signature,plans=plans),indent=2))
    return target['name'],plans

def runs(good,edges,times,minimum=30):
    accepted=good[:-1]&good[1:]
    changes=np.diff(np.r_[False,accepted,False].astype(int))
    return [dict(start_utc=times[a].isot,end_utc=times[b].isot,
                 start=times[a].to_datetime(timezone=timezone.utc).astimezone(TZ).strftime('%H:%M'),
                 end=times[b].to_datetime(timezone=timezone.utc).astimezone(TZ).strftime('%H:%M'),
                 minutes=float(edges[b]-edges[a]))
            for a,b in zip(np.where(changes==1)[0],np.where(changes==-1)[0]) if edges[b]-edges[a]>=minimum]

def fits(options,start,plans=None,max_tier='extended'):
    """Return the least-airmass tier accommodating an entire visit."""
    for tier,x,minutes in TIERS:
        if x>dict((t,a) for t,a,_ in TIERS)[max_tier]:break
        duration=(plans or {}).get(tier,{}).get('visit_minutes',minutes+10)
        end=start+pd.Timedelta(minutes=duration)
        for r in options[tier]:
            if pd.Timestamp(r['start_utc'],tz='UTC')<=start+pd.Timedelta(milliseconds=1) and pd.Timestamp(r['end_utc'],tz='UTC')>=end-pd.Timedelta(milliseconds=1):
                return tier,duration
    return None

def main():
    targets=pd.read_csv(OUT/'compact_review_objects.csv')
    with ProcessPoolExecutor(max_workers=4) as pool:
        plans={}
        for i,(name,value) in enumerate(pool.map(exposure,targets.to_dict('records')),1):
            plans[name]=value
            if i%25==0:print('Airmass ETC',i,'/',len(targets),flush=True)
    if not all(plans.values()):raise ValueError('Every prepared target needs an accepted continuum reference')
    coords=SkyCoord(targets.ra.to_numpy()*u.deg,targets.dec.to_numpy()*u.deg)
    windows={n:[] for n in targets.name};meta={};coverage=[]
    for night,(date,part) in OBS.NIGHTS.items():
        t0,t1,_,_=OBS.night_window(date,part);duration=(t1-t0).to_value(u.min)
        # Ten-second edges keep the displayed setting windows consistent with
        # the full-visit packet checks near a tight airmass boundary.
        edges=np.r_[np.arange(0,duration,1/6),duration];times=t0+edges*u.min
        frame=AltAz(obstime=times,location=OBS.PALOMAR.location,pressure=0*u.hPa)
        aa=coords[:,None].transform_to(frame);x=aa.secz.value
        sep=aa.separation(get_body('moon',times,OBS.PALOMAR.location).transform_to(frame)).deg
        nightrows=[]
        for i,r in enumerate(targets.itertuples()):
            tier_ranges={tier:runs((x[i]>=1)&(x[i]<=limit)&(sep[i]>=40),edges,times)
                         for tier,limit,_ in TIERS}
            if not tier_ranges['fallback']:continue
            usable=tier_ranges['extended'];fallback=tier_ranges['fallback']
            good=(x[i]>=1)&(x[i]<=1.8)&(sep[i]>=40)
            fallback_good=(x[i]>=1)&(x[i]<=2)&(sep[i]>=40)
            q=good if usable else fallback_good
            w=dict(night=night,tier_ranges=tier_ranges,ranges=usable,
                   longest_minutes=max((v['minutes'] for v in usable),default=0),
                   fallback_longest_minutes=max(v['minutes'] for v in fallback),
                   min_airmass=float(x[i,q].min()),moon_min=float(sep[i,q].min()),moon_max=float(sep[i,q].max()),
                   minutes_airmass_le1p3=float((((x[i,:-1]<=1.3)&(x[i,:-1]>=1)&good[:-1]&good[1:])*np.diff(edges)).sum()),
                   curve=[[float(times[j].mjd),float(x[i,j]) if 0<x[i,j]<4 else None,float(sep[i,j])] for j in range(0,len(times),30)])
            windows[r.name].append(w);nightrows.append((r,w))
        from astroplan import moon_illumination
        meta[night]=dict(label={'sep23':'Sep 23','oct26':'Oct 26','oct27':'Oct 27'}[night],date=date,
            part='first half' if part=='first' else 'full night',
            window=f"{t0.to_datetime(timezone=timezone.utc).astimezone(TZ):%H:%M} → {t1.to_datetime(timezone=timezone.utc).astimezone(TZ):%H:%M}",
            start_mjd=float(t0.mjd),end_mjd=float(t1.mjd),mjd=float((t0+(t1-t0)/2).mjd),
            moon_percent=round(100*float(moon_illumination(t0+(t1-t0)/2))),count=len(nightrows))
        for minute in np.arange(0,duration-30,30):
            start=pd.Timestamp((t0+minute*u.min).to_datetime(timezone=timezone.utc))
            row=dict(night=night,start_pdt=start.tz_convert(TZ).strftime('%Y-%m-%d %H:%M'))
            for max_tier,_,_ in TIERS:
                for role in ['manifold','reserve']:
                    names=[r.name for r,w in nightrows if r.pool_role==role and fits(w['tier_ranges'],start,plans[r.name],max_tier)]
                    row[f'{max_tier}_{role}']=len(names);row[f'{max_tier}_{role}_names']=';'.join(names)
            coverage.append(row)
        plot_night(night,date,t0,t1,nightrows,plans,coverage)
    value=dict(windows=windows,nights=meta,plans=plans,target_names=sorted(targets.name))
    (OUT/'airmass_options.json').write_text(json.dumps(value))
    pd.DataFrame([dict(name=name,tier=tier,**plan) for name,pp in plans.items() for tier,plan in pp.items()]).to_csv(OUT/'airmass_exposure_plans.csv',index=False)
    pd.DataFrame(coverage).to_csv(OUT/'airmass_visit_coverage.csv',index=False)
    config=json.loads((OUT/'review_selection.json').read_text())
    config.update(airmass_options=True,airmass_preferred=1.5,airmass_max=1.8,airmass_fallback=2.,
                  version='three-nights-airmass-tiers-2026-09-20',shared_night_backups=True)
    (OUT/'review_selection.json').write_text(json.dumps(config,indent=2))
    print(pd.DataFrame(coverage).groupby('night')[[f'{t}_manifold' for t,_,_ in TIERS]].min().to_string(),flush=True)

def plot_night(night,date,t0,t1,rows,plans,coverage):
    # All cached candidates eligible on this night; night sharing is explicit.
    rows=sorted(rows,key=lambda v:(v[1]['tier_ranges']['fallback'][0]['start_utc'],v[0].ra))
    fig,(top,ax)=plt.subplots(2,1,figsize=(14,max(9,.18*len(rows)+4)),gridspec_kw={'height_ratios':[2,max(6,.18*len(rows))]},sharex=True)
    for i,(r,w) in enumerate(rows):
        for tier,_,_ in reversed(TIERS):
            for segment in w['tier_ranges'][tier]:
                left,right=[mdates.date2num(pd.Timestamp(segment[k],tz='UTC').to_pydatetime()) for k in ['start_utc','end_utc']]
                color=COLORS[tier] if r.pool_role=='manifold' else {'preferred':'#a36d31','extended':'#cda370','fallback':'#eee0cd'}[tier]
                ax.barh(i,right-left,left=left,height=.7,color=color)
    ax.set_yticks(range(len(rows)),[f"{r.name}{' [reserve]' if r.pool_role=='reserve' else ''}  r={r.r_planning:.1f}" for r,w in rows],fontsize=7)
    ax.invert_yaxis();ax.set_xlabel('Palomar local time (PDT)');ax.grid(axis='x',alpha=.2)
    ax.xaxis.set_major_locator(mdates.HourLocator(tz=TZ));ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M',tz=TZ))
    c=pd.DataFrame([r for r in coverage if r['night']==night]);xx=[mdates.date2num(pd.Timestamp(v,tz=TZ).to_pydatetime()) for v in c.start_pdt]
    for tier,_,_ in TIERS:top.plot(xx,c[f'{tier}_manifold'],color=COLORS[tier],label=f"X ≤ {dict((t,x) for t,x,_ in TIERS)[tier]}",linewidth=2)
    top.axhline(3,color='#a36d31',linestyle='--',linewidth=1);top.set_ylabel('Manifold choices\nfor a full visit');top.legend(loc='upper right',ncol=3,fontsize=8);top.grid(alpha=.2)
    top.set_title(f'{date}: preferred, extended, and fallback airmass windows\nMoon ≥ 40° throughout; historical r < 19; backups shared between nights')
    ax.set_xlim(mdates.date2num(t0.to_datetime(timezone=timezone.utc)),mdates.date2num(t1.to_datetime(timezone=timezone.utc)))
    fig.legend(handles=[Patch(color=COLORS[t],label=f'{t}: X ≤ {x}') for t,x,_ in TIERS]+[Patch(color='#cda370',label='amber shades: reserves')],loc='lower center',ncol=4,fontsize=9,bbox_to_anchor=(.5,.025))
    fig.text(.02,.009,'Bars: visibility, not a schedule. Top: complete 16-minute visits (2×300 s, 1.5″ slit, 2×3 binning). S/N in the plans is per Å at sky V=18.5 with seeing scaled to airmass; no broad-line guarantee.',fontsize=8)
    fig.tight_layout(rect=(0,.045,1,1))
    for ext in ['png','pdf']:fig.savefig(OUT/'three_night_review'/f'{night}_visibility.{ext}',dpi=150)
    plt.close(fig)

if __name__=='__main__':main()
