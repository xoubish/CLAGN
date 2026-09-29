"""Local geometry, conditional continuum sensitivity and an October queue preview.

Does not modify the existing NGPS queue. Uses 2025/2026 cached brightness,
explicit bright-sky scenarios, and current adopted 1.5 arcsec / 2x3 settings.
"""
from pathlib import Path
import sys, os, importlib, json, warnings
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/clagn-matplotlib')
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz, get_body
from astropy.time import Time
from astropy.utils import iers
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import lil_matrix
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]
sys.path.insert(0,str(ROOT/'scripts'))
OBS=importlib.import_module('05_observability')
MODEL=importlib.import_module('22_september_etc')
iers.conf.auto_download=False; iers.conf.auto_max_age=None
warnings.filterwarnings('ignore',message='Tried to get polar motions')
TZ='America/Los_Angeles'

def local(t):return pd.Timestamp(t.to_datetime()).tz_localize('UTC').tz_convert(TZ)

def load():
    d=pd.read_csv(OUT/'shortlist.csv').fillna('')
    old=pd.read_csv(OUT.parent.parent/'history_audit/binned_measurements.csv')
    old['name']=old.internal_id
    new=pd.read_csv(OUT.parent/'replacement_binned_measurements.csv')
    bins=pd.concat([old[['name','band','mjd','mag']],new[['name','band','mjd','mag']]])
    rows=[]
    for r in d.to_dict('records'):
        b=bins[(bins.name==r['name'])&(bins.band=='r')].sort_values('mjd')
        if len(b):
            end=b.mjd.max();r['r_adopted']=float(b[b.mjd>=end-180].mag.median());r['brightness_year']=float(Time(end,format='mjd').decimalyear)
            r['brightness_source']='cached r 90-day bins, last 180 days'
        else:r.update(r_adopted=r['r_mag'],brightness_year=np.nan,brightness_source='archival catalogue r')
        r['hbeta_A']=4862.7*(1+r['z'])
        rows.append(r)
    return pd.DataFrame(rows)

def geometry(targets):
    coord=SkyCoord(targets.ra.to_numpy()*u.deg,targets.dec.to_numpy()*u.deg)
    grids={};rows=[]
    for night in ['oct26','oct27']:
        dusk,dawn,_,_=OBS.night_window(OBS.NIGHTS[night][0],'full')
        # Reserve the same ten-minute standard blocks at the two ends.
        start=local(dusk).ceil('min')+pd.Timedelta(minutes=10)
        end=local(dawn).floor('min')-pd.Timedelta(minutes=10)
        t=Time(start.to_pydatetime())+np.arange(int((end-start).total_seconds()/30)+1)*30*u.s
        frame=AltAz(obstime=t,location=OBS.PALOMAR.location,pressure=0*u.hPa)
        aa=coord[:,None].transform_to(frame);x=aa.secz.value
        moon=get_body('moon',t,OBS.PALOMAR.location).transform_to(frame)
        sep=aa.separation(moon).deg
        grids[night]=(t,x,sep)
        for i,r in targets.reset_index(drop=True).iterrows():
            for limit in [1.5,1.8,2.0]:
                mask=(x[i]>=1)&(x[i]<=limit)&(sep[i]>=40)
                transitions=np.diff(np.r_[False,mask,False].astype(int))
                for a,b in zip(np.where(transitions==1)[0],np.where(transitions==-1)[0]):
                    if b-a<2:continue
                    rows.append(dict(name=r['name'],night=night,airmass_limit=limit,
                        start_pdt=local(t[a]).isoformat(),end_pdt=local(t[b-1]).isoformat(),
                        minutes=(b-a-1)*.5,best_airmass=float(x[i,a:b].min()),moon_min_deg=float(sep[i,a:b].min())))
    return grids,pd.DataFrame(rows)

def sensitivity(targets,grids):
    for ch,rn,scale in zip(MODEL.CFG.channels,[2.8,7.8,3.7,4.6],[.193,.189,.186,.186]):
        MODEL.CFG.readnoise[ch]=rn*u.count/u.pix;MODEL.CFG.platescale[ch]=scale*u.arcsec/u.pix
    rows=[]
    for i,r in targets.iterrows():
        if not r['shortlist']:continue
        lo,hi=r.hbeta_A/10-4,r.hbeta_A/10+4
        channels=[ch for ch in ['R','I','G','U'] if MODEL.CFG.channelRange[ch][0].to_value(u.nm)<=lo and MODEL.CFG.channelRange[ch][1].to_value(u.nm)>=hi]
        def sn(ch,seconds,mag,sky,x):
            cmd=[ch,str(lo),str(hi),'EXPTIME',str(seconds),'-slit','SET','1.5','-binspect','3','-binspat','2','-seeing','1.3','500','-airmass',str(x),'-skymag',str(sky),'-mag',str(mag),'-magsystem','AB','-magfilter','match','-noslicer']
            a=MODEL.ETC.parser.parse_args(cmd);MODEL.ETC.check_inputs_add_units(a)
            result=MODEL.ETC.main(a,quiet=True)
            return float(result['SNR'].value)*np.sqrt(2)/np.sqrt(float((MODEL.CFG.dLambda[ch]*3).to_value(u.AA)))
        if not channels:raise ValueError(r['name']+' Hbeta falls outside channels')
        options={ch:sn(ch,300,r.r_adopted,18.,1.5) for ch in channels};ch=max(options,key=options.get)
        # Evaluate bright-sky cases; these are not lunar-sky predictions.
        for sky in [18.,18.5]:
            for x in [1.5,2.0]:
                for extra in [0.,.5]:
                    seconds=300;value=sn(ch,seconds,r.r_adopted+extra,sky,x)
                    initial=value
                    while value<5 and seconds<900:
                        seconds+=60;value=sn(ch,seconds,r.r_adopted+extra,sky,x)
                    rows.append(dict(name=r['name'],channel=ch,hbeta_A=r.hbeta_A,r_adopted=r.r_adopted,
                        brightness_year=r.brightness_year,sky_V_mag_arcsec2=sky,airmass=x,extra_faint_mag=extra,
                        snr_per_A_2x300=initial,seconds_each_for_snr5=seconds,
                        verified_snr_per_A=value,reaches_snr5=value>=5))
        print(r['name'],'r=',round(r.r_adopted,2),'base S/N/A=',round(options[ch],2),flush=True)
    result=pd.DataFrame(rows);result.to_csv(OUT/'ngps_sensitivity.csv',index=False)
    return result

def place(targets,grids,durations,label):
    # Twelve-minute placement cells. Longer exposures occupy multiple cells.
    cells=[];edges=[]
    for night,(t,x,sep) in grids.items():
        n=(len(t)-1)//24;offset=len(cells)
        cells.extend([dict(night=night,start_pdt=local(t[k*24]).isoformat(),end_pdt=local(t[(k+1)*24]).isoformat()) for k in range(n)])
        for i,r in targets.iterrows():
            blocks=int(np.ceil(durations[r['name']]/12));limit=2 if r['shortlist'] or r['old_jwst'] else 1.8
            for k in range(n-blocks+1):
                sl=slice(k*24,(k+blocks)*24+1)
                if np.min(x[i,sl])>=1 and np.max(x[i,sl])<=limit and np.min(sep[i,sl])>=40:
                    edges.append((i,offset+k,blocks,float(x[i,sl].max()),float(sep[i,sl].min())))
    nt,nc,ne=len(targets),len(cells),len(edges)
    A=lil_matrix((nt+nc,ne),dtype=float);cost=[]
    for j,(i,k,blocks,x,m) in enumerate(edges):
        A[i,j]=1
        for b in range(blocks):A[nt+k+b,j]=1
        r=targets.iloc[i]
        # First retain as many existing objects as possible; prefer old JWST
        # objects among equivalent-count alternatives. Shortlist is mandatory.
        reward=10000+100*bool(r.old_jwst)
        cost.append(-reward+5*x-.005*m)
    lower=np.zeros(nt+nc);upper=np.ones(nt+nc)
    lower[:nt]=targets['shortlist'].astype(float)
    answer=milp(np.array(cost),integrality=np.ones(ne),bounds=Bounds(0,1),constraints=LinearConstraint(A.tocsr(),lower,upper),options={'time_limit':45,'mip_rel_gap':.0001})
    if answer.x is None:
        return dict(status=answer.message,feasible=False)
    rows=[]
    for j in np.where(answer.x>.5)[0]:
        i,k,blocks,x,m=edges[j];r=targets.iloc[i];last=cells[k+blocks-1]
        rows.append(dict(name=r['name'],shortlist=bool(r['shortlist']),old_jwst=bool(r.old_jwst),night=cells[k]['night'],start_pdt=cells[k]['start_pdt'],end_pdt=last['end_pdt'],reserved_minutes=12*blocks,required_minutes=durations[r['name']],airmass_max=x,moon_min_deg=m))
    frame=pd.DataFrame(rows).sort_values(['night','start_pdt']);frame.to_csv(OUT/f'ngps_queue_preview_{label}.csv',index=False)
    assert frame.name.is_unique and set(targets[targets['shortlist']].name)<=set(frame.name)
    return dict(feasible=True,status=answer.message,selected=len(frame),shortlist_count=int(frame.shortlist.sum()),
        omitted_existing=targets[~targets.name.isin(frame.name)].name.tolist(),
        shortlist_reserved_minutes=int(frame[frame.shortlist].reserved_minutes.sum()),
        shortlist_max_airmass=float(frame[frame.shortlist].airmass_max.max()),
        shortest_moon_separation=float(frame.moon_min_deg.min()),
        solver_gap=float(answer.mip_gap))

def main():
    shortlist=load();shortlist['shortlist']=True
    existing=pd.read_csv(ROOT/'observing/oct26_27_2026/all_80_targets.csv').rename(columns=str.lower)
    existing['shortlist']=existing.name.isin(shortlist.name)
    existing['old_jwst']=existing.jwst_id.fillna('').ne('')
    targets=pd.concat([shortlist,existing[~existing['shortlist']]],ignore_index=True)
    targets['old_jwst']=targets.name.isin(existing[existing.old_jwst].name)
    targets['shortlist']=targets.name.isin(shortlist.name)
    targets.to_csv(OUT/'ngps_planning_targets.csv',index=False)
    grids,windows=geometry(targets);windows[windows.name.isin(shortlist.name)].to_csv(OUT/'ngps_shortlist_windows.csv',index=False)
    # Visibility summary includes failures, rather than silently dropping targets.
    rows=[]
    for _,r in shortlist.iterrows():
        w=windows[(windows.name==r['name'])&(windows.airmass_limit==2)]
        rows.append(dict(name=r['name'],r_adopted=r.r_adopted,brightness_year=r.brightness_year,
            max_window_minutes=float(w.minutes.max()) if len(w) else 0,
            feasible_nights=','.join(sorted(w[w.minutes>=12].night.unique()))))
    pd.DataFrame(rows).to_csv(OUT/'ngps_visibility_summary.csv',index=False)
    sens=sensitivity(targets,grids)
    durations={n:12 for n in targets.name};summary={}
    summary['baseline_2x300']=place(targets,grids,durations,'2x300')
    # Conservative planning at X=2, sky V=18, baseline brightness. Round up
    # to full 12-minute cells for schedule margins; endpoint constraints apply.
    adopted=sens[(sens.sky_V_mag_arcsec2==18)&(sens.airmass==2)&(sens.extra_faint_mag==0)]
    for _,r in adopted.iterrows():durations[r['name']]=2*r.seconds_each_for_snr5/60+2
    summary['snr5_bright_sky']=place(targets,grids,durations,'snr5')
    summary['exposure_case_failures']=adopted[~adopted.reaches_snr5].name.tolist()
    summary['normalization']='Flat fnu near observed Hbeta, AB magnitude equal to last cached r; optimistic for some blue continua/host profiles. No actual 2026 brightness assumed for replacements.'
    summary['settings']='2 exposures; slit 1.5 arcsec; spatial x spectral binning 2x3; no slicer; zenith seeing 1.3 arcsec; 2 min empirical visit overhead.'
    summary['geometry']='Topocentric unrefracted AltAz, 30-second sampling, builtin ephemeris, bundled IERS extrapolation; Moon>=40 deg, X<=2 for shortlist. Two standards/night reserved.'
    (OUT/'ngps_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
