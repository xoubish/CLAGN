"""Build the October 26/27 2026 NGPS proposal from cached science evidence.

Recomputes topocentric geometry; solves both nights jointly; includes every
observable member of the authoritative JWST sample. Does not alter September,
the JWST proposal, or the published observer page. Run with /opt/anaconda3/bin/python.
"""
from pathlib import Path
from datetime import timezone
import hashlib
import importlib
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-october-matplotlib')
import numpy as np
import pandas as pd
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import lil_matrix
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz, get_body
from astropy.time import Time
from astropy.utils import iers
from astroplan import moon_illumination
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from october_review import apply_review, write_review, september_performance, REVIEW_DATE
import october_additions

REVIEW_DATE = october_additions.DATE

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / 'data/reselection_2026-09-20'
DEST = ROOT / 'observing/oct26_27_2026'
OBS = importlib.import_module('05_observability')
PACK = importlib.import_module('40_september_observing_packet')
TZ = 'America/Los_Angeles'
SETTINGS = dict(slit_arcsec=1.5, binspat=2, binspect=3, slitangle='PA')
VISIT = 14
iers.conf.auto_download = False
iers.conf.auto_max_age = None


def truth(v):
    return str(v).lower() == 'true'


def local(t):
    return pd.Timestamp(t.to_datetime(timezone=timezone.utc)).tz_convert(TZ)


def number(v, default=0):
    return float(v) if pd.notna(v) else default


def write_csv(frame, path):
    """Keep calculation keys unchanged; all exported column headers are uppercase."""
    frame.to_csv(path,index=False,header=[str(c).upper() for c in frame.columns])


def load_targets():
    t = pd.read_csv(CACHE / 'compact_review_objects.csv').set_index('name')
    science = pd.read_csv(CACHE / 'three_night_review/science_and_sensitivity.csv').set_index('name')
    jwst = pd.read_csv(ROOT / 'jwst_proposal/inputs/jwst_sample_cycle6.csv')
    observed = set(pd.read_csv(ROOT / 'sep23_data/reduction_20260924/products_p330e/target_summary.csv').NAME)
    for col in science.columns:
        t[col] = science[col]
    t['jwst_id'] = ''
    t['jwst_family'] = ''
    t['jwst_target'] = ''
    t['coordinate_source'] = 'prepared catalogue'
    t['r_adopted'] = t.ztf_r_latest180_mag.fillna(t.r_planning)
    t['brightness_source'] = np.where(t.ztf_r_latest180_mag.notna(), 'cached ZTF last-180-day median', 'archival planning r')
    t['brightness_epoch'] = t.ztf_last_date.fillna('archival; epoch not specified')
    for j in jwst.to_dict('records'):
        name = j['internal_id']
        if name not in t.index:
            t.loc[name, ['ra', 'dec', 'z', 'r_planning', 'r_adopted']] = [j['ra'], j['dec'], j['z'], j['r_mag'], j['r_mag']]
            t.loc[name, ['brightness_source', 'brightness_epoch', 'pool_role', 'field_status', 'field_notes']] = [
                'JWST proposal catalogue r', 'archival; epoch not specified', 'JWST', 'pending',
                'JWST addition: optical acquisition field and nuclear centering need inspection']
            t.loc[name, 'science_question'] = f"JWST {j['id']} {j['family']}: establish optical Balmer state for the MIRI dust-response sample; {j['reference']}"
        else:
            delta = SkyCoord(t.loc[name, 'ra']*u.deg, t.loc[name, 'dec']*u.deg).separation(SkyCoord(j['ra']*u.deg,j['dec']*u.deg)).arcsec
            assert delta < 1, (name, delta)
        t.loc[name, ['ra', 'dec', 'coordinate_source']] = [j['ra'], j['dec'], 'authoritative JWST sample']
        t.loc[name, ['jwst_id', 'jwst_family', 'jwst_target']] = [j['id'], j['family'], j['target']]
    t['observed_sep23'] = t.index.isin(observed)
    t['is_jwst'] = t.jwst_id.fillna('').ne('')
    # September targets are repeated only where requested for JWST coverage.
    t['eligible_science'] = t.is_jwst | (~t.observed_sep23 & t.r_adopted.lt(19) & t.field_status.eq('clear') & t.baseline_quality.map(truth))
    scores, reasons = [], []
    for name, r in t.iterrows():
        optical, ir = truth(r.get('post_spectrum_optical_trigger')), truth(r.get('post_spectrum_ir_flag'))
        manifold = r.pool_role == 'manifold'
        known = r.get('known_state_status') == 'catalog-confirmed CLAGN'
        score = 100*optical + 25*ir + 20*manifold + 10*known
        score += min(10, number(r.get('snr300_x1p3_sky18p5')))
        score += 5*min(1, number(r.get('manifold_cl_neighbor_fraction')))
        # A small field-readiness preference; mandatory JWST objects remain protected.
        score += 2*(r.field_status == 'clear')
        reasons.append('JWST required optical state' if r.is_jwst else
                       'Dated optical '+str(r.get('trigger_direction')) if optical else
                       'Dated IR evolution after reference' if ir else
                       'Known CLAGN current state' if known else
                       'Manifold discovery / archival Balmer comparison' if manifold else
                       'Bright comparison reserve filling a geometry gap')
        scores.append(score)
    t['science_score'] = scores
    t['priority_reason'] = reasons
    t = apply_review(t)
    t = october_additions.apply_additions(t)
    t.index.name = 'name'
    # Reject spatial aliases, including JWST catalogue names differing from pool names.
    coords = SkyCoord(t.ra.to_numpy()*u.deg, t.dec.to_numpy()*u.deg)
    sep = coords[:,None].separation(coords[None,:]).arcsec
    for a, b in zip(*np.where(np.triu(sep < 3, 1))):
        if t.iloc[a].is_jwst and not t.iloc[b].is_jwst:
            t.iloc[b, t.columns.get_loc('eligible_science')] = False
        elif t.iloc[b].is_jwst and not t.iloc[a].is_jwst:
            t.iloc[a, t.columns.get_loc('eligible_science')] = False
        else:
            raise ValueError(f'Unresolved duplicate: {t.index[a]} / {t.index[b]}')
    return t, jwst


def geometry(targets):
    coords = SkyCoord(targets.ra.to_numpy()*u.deg, targets.dec.to_numpy()*u.deg)
    slots, edges, meta, windows = [], [], {}, []
    for night in ['oct26', 'oct27']:
        date = OBS.NIGHTS[night][0]
        dusk, dawn, _, _ = OBS.night_window(date, 'full')
        start = local(dusk).ceil('min')
        end = local(dawn).floor('min')
        # Ten minutes per standard at beginning/end. Remaining 14-minute slots
        # may be left empty by the optimizer and serve as delay/deeper-exposure time.
        first = Time((start+pd.Timedelta(minutes=10)).to_pydatetime())
        last = Time((end-pd.Timedelta(minutes=10)).to_pydatetime())
        n = int(np.floor((last-first).to_value(u.min)/VISIT))
        times = first + np.arange(n*VISIT*2+1)*.5*u.min
        frame = AltAz(obstime=times, location=OBS.PALOMAR.location, pressure=0*u.hPa)
        aa = coords[:,None].transform_to(frame)
        x = aa.secz.value
        sep = aa.separation(get_body('moon',times,OBS.PALOMAR.location).transform_to(frame)).deg
        for k in range(n):
            ix = slice(k*VISIT*2, (k+1)*VISIT*2+1)
            xmax, xmin, moonmin = x[:,ix].max(axis=1), x[:,ix].min(axis=1), sep[:,ix].min(axis=1)
            # Small numerical margin relative to final 20-second checks.
            limits = targets.airmass_limit.to_numpy()-.001
            ok = (xmin>=1)&(xmax<=limits)&(moonmin>=40.02)&targets.eligible_science.to_numpy()
            sid = len(slots)
            slots.append(dict(night=night,slot=k,start=local(times[k*VISIT*2]),end=local(times[(k+1)*VISIT*2])))
            for i in np.where(ok)[0]:
                r = targets.iloc[i]
                if r.protected_followup and night != 'oct27': continue
                if r.requested_night and night != r.requested_night: continue
                stop = k*VISIT*2 + int(r.visit_minutes*2)
                if stop >= len(times): continue
                full = slice(k*VISIT*2, stop+1)
                if x[i,full].min()<1 or x[i,full].max()>limits[i] or sep[i,full].min()<40.02: continue
                edges.append((i,sid,float(x[i,full].max()),float(sep[i,full].min()),float(x[i,full].mean())))
        # Full-night visibility audit, independent of the placement grid.
        ts = dusk + np.r_[np.arange(0,(dawn-dusk).to_value(u.min),.5),(dawn-dusk).to_value(u.min)]*u.min
        fr = AltAz(obstime=ts,location=OBS.PALOMAR.location,pressure=0*u.hPa)
        az = coords[:,None].transform_to(fr)
        xx = az.secz.value
        ms = az.separation(get_body('moon',ts,OBS.PALOMAR.location).transform_to(fr)).deg
        for i, name in enumerate(targets.index):
            for ceiling in [1.5,1.8,2.0]:
                good = (xx[i]>=1)&(xx[i]<=ceiling)&(ms[i]>=40)
                transitions = np.diff(np.r_[False,good,False].astype(int))
                for a,b in zip(np.where(transitions==1)[0],np.where(transitions==-1)[0]):
                    if b-a<2: continue
                    wa,wb=local(ts[a]),local(ts[b-1])
                    duration = int(targets.iloc[i].visit_minutes)
                    if (wb-wa).total_seconds()<duration*60: continue
                    windows.append(dict(name=name,night=night,airmass_limit=ceiling,start_pdt=wa.isoformat(),
                        end_pdt=wb.isoformat(),latest_start_pdt=(wb-pd.Timedelta(minutes=duration)).isoformat(),
                        moon_min=float(ms[i,a:b].min())))
        meta[night] = dict(date=date,dusk_pdt=local(dusk).isoformat(),dawn_pdt=local(dawn).isoformat(),
            minutes=float((dawn-dusk).to_value(u.min)),moon_illumination_percent=100*float(moon_illumination(dusk+(dawn-dusk)/2)),
            standard_start_pdt=start.isoformat(),standard_end_start_pdt=(end-pd.Timedelta(minutes=10)).isoformat(),
            science_minutes=560,standard_minutes=20,remaining_minutes=float((dawn-dusk).to_value(u.min))-580)
    return slots,edges,meta,pd.DataFrame(windows)


def solve(targets, slots, edges, windows):
    nt, ns, ne = len(targets),len(slots),len(edges)
    matrix = lil_matrix((nt+ns+2,ne),dtype=float)
    costs=[]
    covered={e[0] for e in edges}
    previous=october_additions.previous_names()
    for k,(i,s,x,m,xmean) in enumerate(edges):
        r=targets.iloc[i]
        matrix[i,k]=1
        for offset in range(int(np.ceil(r.visit_minutes/VISIT))):
            assert slots[s+offset]['night'] == slots[s]['night']
            matrix[nt+s+offset,k]=1
        matrix[nt+ns+(slots[s]['night']=='oct27'),k]=1
        # Science evidence dominates; use airmass and Moon separation to place
        # a chosen object. No calibrated yield probability is claimed.
        costs.append(-float(r.science_score)-10000*(targets.index[i] in previous)+12*(xmean-1)+20*max(0,x-1.5)-.015*min(m,100)+1e-7*k)
    lower=np.zeros(nt+ns+2); upper=np.ones(nt+ns+2)
    for i,r in enumerate(targets.itertuples()):
        if (r.is_jwst or r.protected_followup or r.requested_addition) and i in covered: lower[i]=1
        if r.protected_followup: assert i in covered, 'Protected follow-up has no feasible slot'
        if r.requested_addition: assert i in covered, f'Requested target has no feasible slot: {r.Index}'
    lower[-2:]=40;upper[-2:]=40
    result=milp(np.array(costs),integrality=np.ones(ne),bounds=Bounds(0,1),
        constraints=LinearConstraint(matrix.tocsc(),lower,upper),options={'time_limit':120,'mip_rel_gap':.000001})
    if result.x is None: raise RuntimeError(f'No feasible 40+40 sequence: {result.message}')
    print('Optimizer:',result.message,'gap',getattr(result,'mip_gap',None),flush=True)
    chosen=[k for k,v in enumerate(result.x) if v>.5]
    rows=[]
    for k in chosen:
        i,s,x,m,mean=edges[k]; slot=dict(slots[s]);r=targets.iloc[i].to_dict();name=targets.index[i]
        slot['end'] = slot['start'] + pd.Timedelta(minutes=int(r['visit_minutes']))
        exact=PACK.geometry(float(r['ra']),float(r['dec']),slot['start'],slot['end'])
        assert exact['moon_min']>=40 and exact['airmass_max_actual']<=r['airmass_limit']
        flags=[]
        if r['r_adopted']>=19:flags.append('JWST faint exception to r<19; assess 2x300 s in Quicklook and use slack for depth')
        if name in ['CLAGN_0987','CLAGN_0969'] or r['r_adopted']<15:flags.append('Resolved galaxy: photometric r is not nuclear slit brightness; check first exposure for saturation and host dilution')
        if r['field_status']!='clear':flags.append(str(r['field_notes']))
        if exact['airmass_max_actual']>1.8:flags.append('JWST high-airmass exception X>1.8')
        if r['observed_sep23']:flags.append('Intentional JWST repeat of September spectrum' if r['is_jwst'] else 'Protected September follow-up; 3x300 s, not the default 2x300 s')
        snr=number(r.get('snr300_x1p3_sky18p5'),np.nan)
        if np.isfinite(snr) and snr<5:flags.append('Cached nominal Hbeta continuum S/N/A below 5; not a slot prediction')
        if not np.isfinite(snr):flags.append('No cached NGPS Hbeta ETC calculation')
        if r['requested_addition']:flags.append('Requested addition; see archive review and exposure-scenario limitations')
        matching=windows[(windows.name==name)&(windows.night==slot['night'])&
            (windows.airmass_limit==r['airmass_limit'])]
        matching=matching[(pd.to_datetime(matching.start_pdt,utc=True)<=slot['start'])&
                          (pd.to_datetime(matching.end_pdt,utc=True)>=slot['end'])]
        assert len(matching)==1,(name,slot)
        rows.append(dict(name=name,night=slot['night'],start_pdt=slot['start'].isoformat(),end_pdt=slot['end'].isoformat(),
            start_utc=slot['start'].tz_convert('UTC').isoformat(),latest_start_pdt=matching.iloc[0].latest_start_pdt,
            ra=r['ra'],dec=r['dec'],z=r['z'],jwst_id=r['jwst_id'],jwst_target=r['jwst_target'],jwst_family=r['jwst_family'],
            observed_sep23=r['observed_sep23'],r_adopted=r['r_adopted'],brightness_source=r['brightness_source'],
            brightness_epoch=r['brightness_epoch'],science_score=r['science_score'],priority_reason=r['priority_reason'],
            science_question=r['science_question'],pool_role=r['pool_role'],protected_followup=r['protected_followup'],
            reference_date=r.get('reference_date'),optical_delta_r=r.get('ztf_r_change_after_reference'),
            ir_delta_w1=r.get('neowise_change_after_reference_mag'),nominal_snr_per_A=snr,
            slit_arcsec=1.5,binspat=2,binspect=3,exposures=int(r['exposures']),seconds_each=300,overhead_minutes=int(r['overhead_minutes']),visit_minutes=int(r['visit_minutes']),airmass_limit=r['airmass_limit'],
            flags='; '.join(flags),**exact))
    df=pd.DataFrame(rows).sort_values(['night','start_pdt']).reset_index(drop=True)
    df.insert(0,'order',df.groupby('night').cumcount()+1)
    assert df.groupby('night').size().tolist()==[40,40]
    assert df.name.nunique()==80
    assert 'P8544' in set(df.name)
    return df,result


def telescope_row(r, note):
    comment=f"{r['start_pdt'][11:16]} PDT; latest start {r['latest_start_pdt'][11:16]}; r={r['r_adopted']:.2f} ({r['brightness_source']}); Moon>={r['moon_min']:.1f}deg; X<={r['airmass_max_actual']:.2f}. {r['priority_reason']}. {r['science_question']}. {r['flags']}"
    return PACK.ngps_row(r,dict(seconds_each=300,exposures=r['exposures'],airmass=r['airmass_limit']),SETTINGS,note,comment)


def write_outputs(targets,jwst,df,slots,edges,meta,result,windows):
    write_csv(df,DEST/'all_80_targets.csv')
    ranked=targets.reset_index().sort_values(['eligible_science','is_jwst','science_score'],ascending=False)
    ranked.insert(0,'review_rank',np.arange(1,len(ranked)+1))
    write_csv(ranked,DEST/'candidate_ranking.csv')
    write_csv(windows,DEST/'visibility_windows.csv')
    standards=pd.read_csv(ROOT/'data/standards_spectrophotometric.csv').set_index('name')
    standards_rows=[];gaps=[];backup_rows=[]
    selected=set(df.name)
    for night in ['oct26','oct27']:
        d=df[df.night==night]
        meta[night]['science_minutes'] = int(d.visit_minutes.sum())
        meta[night]['remaining_minutes'] = meta[night]['minutes'] - meta[night]['science_minutes'] - 20
        write_csv(d,DEST/f'{night}_sequence.csv')
        ngps=[telescope_row(r,f"{r['order']:02d} {r['jwst_id'] or 'SCI'} {r['start_pdt'][11:16]}") for r in d.to_dict('records')]
        (DEST/f'{night}_ngps.csv').write_text(PACK.csv_text(ngps),encoding='ascii')
        (DEST/f'{night}_coordinates.csv').write_text(PACK.csv_text(ngps,True),encoding='ascii')
        cal=[]
        for name,seconds,key in [('BD+28 4211',10,'standard_start_pdt'),('Feige 34',20,'standard_end_start_pdt')]:
            r=standards.loc[name].to_dict();r['name']=name
            start=pd.Timestamp(meta[night][key]);end=start+pd.Timedelta(minutes=10)
            geo=PACK.geometry(r['ra'],r['dec'],start,end)
            assert geo['airmass_max_actual']<1.8
            standards_rows.append(dict(night=night,name=name,start_pdt=start.isoformat(),end_pdt=end.isoformat(),
                seconds_each=seconds,exposures=2,exposure_status='starting estimate; inspect counts before repeat',**geo))
            cal.append(PACK.ngps_row(r,dict(seconds_each=seconds,exposures=2,airmass=1.8),SETTINGS,
                f"STD {start:%H:%M}",f"10 min reserved. Starting exposure estimate; inspect all four channels and adjust before repeat. V={r['V']}; X<={geo['airmass_max_actual']:.2f}"))
        (DEST/f'{night}_standards_ngps.csv').write_text(PACK.csv_text(cal),encoding='ascii')
        occupied=[(pd.Timestamp(r.start_pdt),pd.Timestamp(r.end_pdt)) for r in d.itertuples()]
        occupied += [(pd.Timestamp(r['start_pdt']),pd.Timestamp(r['end_pdt'])) for r in standards_rows if r['night']==night]
        occupied.sort();last=pd.Timestamp(meta[night]['dusk_pdt'])
        for a,b in occupied+[(pd.Timestamp(meta[night]['dawn_pdt']),pd.Timestamp(meta[night]['dawn_pdt']))]:
            assert a>=last
            if (a-last).total_seconds()>1:gaps.append(dict(night=night,start_pdt=last.isoformat(),end_pdt=a.isoformat(),minutes=(a-last).total_seconds()/60))
            last=b
        # Alternatives are matched to exact primary slots, with no selected target
        # reused as a backup and no implicit replacement of a mandatory JWST visit.
        for r in d.to_dict('records'):
            sid=next(s for s,v in enumerate(slots) if v['night']==night and v['start'].isoformat()==r['start_pdt'])
            alternatives=[e for e in edges if e[1]==sid and targets.index[e[0]] not in selected and not targets.iloc[e[0]].is_jwst]
            alternatives.sort(key=lambda e:(-targets.iloc[e[0]].science_score,e[2]))
            # The deeper follow-up has a longer slot: validate alternatives over
            # that entire interval, not just the default placement cell.
            checked = []
            for i,s,x,m,mean in alternatives:
                t=targets.iloc[i]
                g=PACK.geometry(float(t.ra),float(t.dec),pd.Timestamp(r['start_pdt']),pd.Timestamp(r['end_pdt']))
                if g['airmass_max_actual']<=1.8 and g['moon_min']>=40:
                    checked.append((i,s,g['airmass_max_actual'],g['moon_min'],mean))
                if len(checked)==2: break
            alternatives=checked
            for choice,(i,s,x,m,_) in enumerate(alternatives[:2],1):
                t=targets.iloc[i]
                backup_rows.append(dict(night=night,replaces_order=r['order'],replaces_name=r['name'],choice=choice,
                    name=targets.index[i],start_pdt=r['start_pdt'],end_pdt=r['end_pdt'],ra=t.ra,dec=t.dec,r_adopted=t.r_adopted,
                    science_score=t.science_score,priority_reason=t.priority_reason,airmass_max_actual=x,moon_min=m,
                    jwst_replacement_warning='Reschedule JWST target if skipped' if r['jwst_id'] else ''))
        # An exportable figure highlights queue timing and Moon/airmass checks.
        fig,(ax,bx)=plt.subplots(1,2,figsize=(13,12),gridspec_kw={'width_ratios':[2.5,1]})
        base=pd.Timestamp(meta[night]['dusk_pdt'])
        for y,r in enumerate(d.to_dict('records')):
            minute=(pd.Timestamp(r['start_pdt'])-base).total_seconds()/60
            ax.barh(y,r['visit_minutes'],left=minute,color='#c67b22' if r['jwst_id'] else '#256d92')
            ax.text(minute+r['visit_minutes']+3,y,r['start_pdt'][11:16],va='center',fontsize=7)
        labels=[f"{r.order:02d}  {r.name}"+(f" [{r.jwst_id}]" if r.jwst_id else '') for r in d.itertuples()]
        ax.set_yticks(range(40),labels,fontsize=8);ax.set_ylim(39.5,-.5);ax.grid(axis='x',alpha=.25)
        ticks=np.arange(0,meta[night]['minutes'],60)
        ax.set_xticks(ticks,[(base+pd.Timedelta(minutes=float(v))).strftime('%H:%M') for v in ticks]);ax.set_xlabel('PDT; amber = JWST')
        bx.scatter(d.airmass_max_actual,range(40),c=d.moon_min,cmap='viridis',vmin=40,vmax=130)
        bx.axvline(1.8,color='gray',linestyle=':');bx.set_xlim(1,2.05);bx.set_ylim(39.5,-.5)
        bx.set_yticks(range(40),[f"{m:.0f}°" for m in d.moon_min],fontsize=8);bx.yaxis.tick_right();bx.set_xlabel('Maximum visit airmass');bx.set_title('Right labels: minimum Moon separation',fontsize=8)
        fig.suptitle(f"{meta[night]['date']}: 40 targets · normally 2×300 s + 4 min · P8544: 3×300 s; HIRAS-F23290-1252: 4×300 s\n{meta[night]['remaining_minutes']:.0f} min spare after 20 min standards; all Moon separations ≥40°",fontsize=12)
        fig.tight_layout(rect=(0,0,1,.955))
        for ext in ['png','pdf']:fig.savefig(DEST/f'{night}_schedule.{ext}',dpi=160)
        plt.close(fig)
    write_csv(pd.DataFrame(standards_rows),DEST/'standards.csv')
    write_csv(pd.DataFrame(gaps),DEST/'free_time.csv')
    backs=pd.DataFrame(backup_rows);write_csv(backs,DEST/'slot_backups.csv')
    for night in ['oct26','oct27']:
        rows=[]
        for r in backs[backs.night==night].to_dict('records'):
            rows.append(PACK.ngps_row(r,dict(seconds_each=300,exposures=2,airmass=1.8),SETTINGS,
                f"B{r['replaces_order']:02d}.{r['choice']} {r['start_pdt'][11:16]}",
                f"Replacement for {r['replaces_name']} at {r['start_pdt'][11:16]} PDT only. r={r['r_adopted']:.2f}; X<={r['airmass_max_actual']:.2f}; Moon>={r['moon_min']:.1f}deg. {r['priority_reason']}. {r['jwst_replacement_warning']}. Skip if already observed."))
        (DEST/f'{night}_backups_ngps.csv').write_text(PACK.csv_text(rows),encoding='ascii')
    audit=[]
    for j in jwst.to_dict('records'):
        name=j['internal_id'];match=df[df.name==name];w=windows[(windows.name==name)&(windows.airmass_limit==2)]
        audit.append(dict(jwst_id=j['id'],name=name,target=j['target'],family=j['family'],
            observed_sep23=bool(targets.loc[name,'observed_sep23']),
            scheduled_night=match.iloc[0].night if len(match) else '',
            start_pdt=match.iloc[0].start_pdt if len(match) else '',
            status='October scheduled' if len(match) else 'September spectrum retained; no October 14-min window at X<=2 and Moon>=40',
            viable_nights=','.join(sorted(w.night.unique())),proposal_r=j['r_mag'],
            october_max_airmass=match.iloc[0].airmass_max_actual if len(match) else np.nan,
            october_min_moon=match.iloc[0].moon_min if len(match) else np.nan))
    audit=pd.DataFrame(audit);write_csv(audit,DEST/'jwst_coverage.csv')
    assert (audit.scheduled_night.ne('')|audit.observed_sep23).all()
    performance=write_review(DEST,targets,df,windows,write_csv)
    october_additions.write_audit(targets,df,write_csv)
    summary=dict(review_date=REVIEW_DATE,september_performance=performance,
        settings={**SETTINGS,'exposures':2,'seconds_each':300,'overhead_minutes':4,'visit_minutes':14,
        'moon_min_deg':40,'airmass_normal_max':1.8,'airmass_jwst_max':2,'preferred_airmass_max':1.5},nights=meta,
        science_targets=len(df),unique_targets=df.name.nunique(),jwst_total=len(jwst),jwst_october=int(df.jwst_id.ne('').sum()),
        jwst_september_only=audit.loc[audit.scheduled_night.eq(''),'name'].tolist(),
        jwst_repeats=df.loc[df.observed_sep23 & df.jwst_id.ne(''),'name'].tolist(),
        september_followups=df.loc[df.protected_followup,'name'].tolist(),
        requested_additions=df.loc[df.name.str.startswith('H'),'name'].tolist(),
        fainter_than_r19=df.loc[df.r_adopted.ge(19),['name','r_adopted']].to_dict('records'),
        minimum_moon_deg=float(df.moon_min.min()),maximum_airmass=float(df.airmass_max_actual.max()),
        optimizer_message=result.message,optimizer_gap=float(result.mip_gap),
        geometry='Astropy builtin ephemeris; topocentric unrefracted AltAz; 30s placement, 20s final visit checks; bundled IERS extrapolation',
        scope='Prepared pool plus authoritative JWST targets and six requested H identities; preserve prior membership where feasible; not a complete parent-catalogue search',
        inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [
            CACHE/'compact_review_objects.csv',CACHE/'three_night_review/science_and_sensitivity.csv',ROOT/'jwst_proposal/inputs/jwst_sample_cycle6.csv',
            ROOT/'jwst_proposal/proposal.tex',ROOT/'jwst_proposal/review/sample_reassessment/current_target_assessment.csv',
            ROOT/'jwst_proposal/review/history_audit/binned_measurements.csv',
            ROOT/'sep23_data/reduction_20260924/audit/raw_frames.csv',ROOT/'sep23_data/reduction_20260924/products_p330e/target_summary.csv',
            ROOT/'sep23_data/reduction_20260924/science_review/README.md',ROOT/'scripts/october_review.py',
            ROOT/'scripts/53_october_observing_packet.py',ROOT/'scripts/october_additions.py',ROOT/'scripts/64_helou_targets.py',
            ROOT/'data/additional_targets/helou_review/feasibility.csv',ROOT/'data/additional_targets/helou_review/targets.json']})
    (DEST/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    return audit,summary


def report(df,audit,summary):
    night_lines=[]
    for n,m in summary['nights'].items():
        d=df[df.night==n]
        night_lines.append(f"| {m['date']} | {m['dusk_pdt'][11:16]}–{m['dawn_pdt'][11:16]} | 40 | {int(d.jwst_id.ne('').sum())} | {m['remaining_minutes']:.0f} min | {m['moon_illumination_percent']:.0f}% |")
    exceptions=df[df.r_adopted.ge(19)|df.airmass_max_actual.gt(1.8)]
    performance=summary['september_performance']
    decisions=pd.read_csv(DEST/'selection_decisions.csv').rename(columns=str.lower)
    changes=decisions[decisions.decision.isin(['added','removed'])]
    text=f'''# October 26–27, 2026: reviewed Palomar / NGPS queue

Updated 2026-10-22 with the [requested H-target additions](H_target_additions.md). The original October 14 review was against the September reductions, raw timestamps, current `proposal.tex`, authoritative 24-target sample, and later history reassessment. [Observer page](../../docs/index.html#oct26). Dates label the PDT evening; UTC observations fall on October 27/28. Telescope upload is still a separate action.

## Decision

Keep **40 science targets per night**, as requested on October 14. Protect all 23 geometrically feasible JWST objects and add **P8544 on October 27 at 3×300 s** to test the September weakening lead at the same integration depth. HIRAS-F23290-1252 uses **4×300 s** in a 26-minute block at X≤1.5; the other 78 October visits use **2×300 s**. All use a 1.5″ slit, 2×3 spatial×spectral binning and parallactic angle. No object is duplicated across the October nights.

The default block is now **14 minutes: 10 minutes integration + 4 minutes overhead**. P8544 reserves 20 minutes (15 + 5), and HIRAS-F23290-1252 reserves 26 minutes (20 + 6). Four minutes is a base planning allowance, slightly below the measured median; the remaining gaps must absorb the difference, acquisition delays and any extra depth. This is an aggressive 40-target programme, not a promise that all targets will finish in poorer conditions.

| Evening | Astronomical dark (PDT) | Science targets | JWST | Spare after standards | Moon illumination |
| --- | --- | ---: | ---: | ---: | ---: |
{chr(10).join(night_lines)}

## What September actually delivered

18 science targets, 39 exposures, 195 minutes on source across {performance['elapsed_science_minutes']:.1f} elapsed minutes from the first science exposure to the last exposure's end. The 17 measurable start-to-next-target intervals give **{performance['median_overhead_minutes']:.2f} minutes median non-integration time**, with a {performance['p90_overhead_minutes']:.2f}-minute 90th percentile. This includes readout, acquisition and operational pauses; it is not a measurement of readout alone. The previous two-minute total allowance was not supported by these timestamps. See [per-target timing](september_performance.csv).

Median U/G/R/I continuum S/N was approximately 18/22/31/38 **per native pixel**, not per Å. Most objects used 2×300 s; P1823 used four exposures and P8544 three. These measurements support the normal starting exposure, but are not a slot-specific October S/N forecast. P8544's I-arm S/N was only about 12.6 even at 3×300 s, so its repeat must not be shortened to the default.

The September triage found no unambiguous wholesale appearance/disappearance of a broad Balmer line. P9694 is the cleanest strengthening lead; P11530 has a valuable long flare/fading history; P8544 is a weakening candidate. P10381 has a published fading episode and a recovery question. P2921 and P9227 require stronger telluric/continuum checks. [All 18 follow-up decisions](september_followup_decisions.csv).

## JWST and nightly priorities

The final proposal tests warm-dust response to diverse recorded histories, including flares and recoveries. The later sustained-change shortlist is a modelling exercise, not an adopted replacement sample; do not remove proposal targets just because their history reverses. Every JWST row now states its individual history and required state check. Selection-family labels are historical labels, not claims that a source is currently rising or fading.

**October 26:** protect UGC 03601 / CLAGN_0969 (R09), which has no usable October 27 window. Complete the assigned evening JWST targets, new discovery targets and morning JWST block.

**October 27:** protect P9506 (F02), P10381 (F04), and P9694 (R01), which cannot be placed on October 26 under the Moon cut. Obtain the deeper P8544 repeat at better geometry than its short October 26 window. P7837 (F07) receives one deliberate JWST repeat on its assigned night.

**P11530 (F06):** no complete visit at Moon ≥40° and X≤2 on either night. Retain its September spectrum; all 24 proposal objects have September and/or planned October optical coverage, but only {summary['jwst_october']} have new October visits. P1823's tiny setting window is too short for a full default visit and overlaps the opening standard; P2921 has no normal full-visit window. Retain their existing September spectra. Do not automatically repeat the remaining September objects.

The eight candidate replacements from the later modelling exercise remain [development reserves](development_reserves.csv), not ready-to-upload alternatives: their acquisition fields and usable spectral baselines have not been validated in this packet. P20778 also has cached r≈19.37. The submitted 24-target sample is preserved.

## What changed in membership

{changes.to_markdown(index=False)}

[All candidate decisions, including night changes](selection_decisions.csv) · [ranked evidence](candidate_ranking.csv). Ordinary new candidates retain the existing dated-optical/IR, manifold and known-CLAGN ranking, accepted baseline, clear catalogue field screen and cached r<19 requirement. P8544 is an explicit protected repeat. Geometry is solved jointly over both nights; selected targets cannot overlap or repeat. This is a review of the prepared pool, not an exhaustive new parent-catalogue search.

## Downloads and execution

| Night | Ordered science CSV | NGPS upload | Coordinates | Standards | Backups | Timeline |
| --- | --- | --- | --- | --- | --- | --- |
| Oct 26 | [sequence](oct26_sequence.csv) | [NGPS](oct26_ngps.csv) | [coordinates](oct26_coordinates.csv) | [standards](oct26_standards_ngps.csv) | [backups](oct26_backups_ngps.csv) | [PDF](oct26_schedule.pdf) |
| Oct 27 | [sequence](oct27_sequence.csv) | [NGPS](oct27_ngps.csv) | [coordinates](oct27_coordinates.csv) | [standards](oct27_standards_ngps.csv) | [backups](oct27_backups_ngps.csv) | [PDF](oct27_schedule.pdf) |

[All 80](all_80_targets.csv) · [JWST coverage](jwst_coverage.csv) · [visibility windows](visibility_windows.csv) · [free intervals](free_time.csv) · [slot backups](slot_backups.csv) · [validation and input hashes](validation.json).

NGPS CSVs encode settings and order, not absolute starts or gaps. Follow the sequence times and `latest_start_pdt`; advancing a rising/Moon-limited object can violate a constraint. Spare time is fragmented and is not all available at dawn. When late, drop the lowest-priority ordinary target that no longer fits before sacrificing a protected JWST visit or P8544 / HIRAS-F23290-1252 depth. Backups are replacements for their listed slot, not additional targets; skip any already used. A backup does not satisfy a displaced protected JWST, P8544 or requested H observation. Coordinate-only files are headerless.

## Calibration, brightness and limits

Reserve ten minutes at each end for BD+28 4211 and Feige 34, with 2×10 s and 2×20 s starting estimates. Inspect counts in all channels before repeating, record seeing/centering and compare the sensitivity curves before combining the night's data. September's P330E and BD+28 calibrations differed by about 41–53% in flux despite repeatability within each star; statistical S/N does not establish absolute spectrophotometric accuracy. Maintain one documented flux-calibration scale and retain alternative-standard comparisons. The [official NGPS guide](https://caltechopticalobservatories.github.io/NGPS/users-manual/quick-start.html) supports same-setting standards at both ends and Quicklook checks for additional exposure needs (checked October 14).

Adopted brightness now uses the proposal history audit's dated survey bins where available, replacing catalogue totals; all remaining values retain their provenance. This is still cached photometry, usually through 2025, not an October 2026 measurement. Saved 2026 alerts inform the history interpretation but are not silently substituted for total-source r photometry. JWST additions still need acquisition-field/nuclear-centering checks. NGC 2617 and UGC 03601 are resolved: survey brightness is not a nuclear slit magnitude. Inspect first exposures for saturation and host dilution.

Faint/high-airmass exceptions:

{exceptions[['name','jwst_id','r_adopted','airmass_max_actual']].to_markdown(index=False,floatfmt='.2f')}

All science visits satisfy Moon ≥40° for their full duration. Ordinary targets have X≤1.8; protected JWST targets may use X≤2.0. Minimum separation is {summary['minimum_moon_deg']:.2f}°, maximum X is {summary['maximum_airmass']:.3f}. Geometry uses unrefracted topocentric AltAz, 30-second placement and 20-second final checks. Bundled IERS tables require extrapolation; rerun near the observing dates and use live sequencer checks. No fresh photometry or new broad-line decomposition was performed in this review.

Rebuild with `/opt/anaconda3/bin/python scripts/53_october_observing_packet.py`, then `scripts/51_observer_page.py`. The September packet and JWST proposal remain historical/scientific inputs. The review's previous queue is preserved in `before_review_20261014.csv`.
'''
    (DEST/'README.md').write_text(text)


if __name__=='__main__':
    DEST.mkdir(parents=True,exist_ok=True)
    october_additions.snapshot()
    snapshot=DEST/'before_review_20261014.csv'
    if not snapshot.exists():
        previous=DEST/'all_80_targets.csv'
        if previous.exists():
            snapshot.write_bytes(previous.read_bytes())
        else:
            snapshot.write_text('NAME,NIGHT\n')
    targets,jwst=load_targets()
    SLOTS,EDGES,META,WINDOWS=geometry(targets)
    print('Feasible assignments:',len(EDGES),flush=True)
    queue,result=solve(targets,SLOTS,EDGES,WINDOWS)
    audit,summary=write_outputs(targets,jwst,queue,SLOTS,EDGES,META,result,WINDOWS)
    report(queue,audit,summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in ['inputs','nights']},indent=2))
