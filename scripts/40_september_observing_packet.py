"""A reviewed September 23 sequence, time-matched alternatives, and dark cards.

This freezes a reviewable packet from the prepared sample. It never uploads to
NGPS, commits, publishes, or silently changes the sequence as archives update.
"""
from pathlib import Path
from datetime import timezone
from concurrent.futures import ThreadPoolExecutor
import argparse,base64,copy,importlib,json,re
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord,AltAz,get_body
from astropy.utils import iers

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/reselection_2026-09-20'
DEST=ROOT/'observing/sep23';DEST.mkdir(parents=True,exist_ok=True)
OBS=importlib.import_module('05_observability');OPTIONS=importlib.import_module('39_airmass_options')
TZ='America/Los_Angeles'
iers.conf.auto_download=False;iers.conf.auto_max_age=None

# Deliberate science choices, not a spectral-count or brightness ranking.
ORIGINAL_PRIMARY=[
 ('P659','2026-09-23 20:16','A single 2013 broad-Hbeta reference and a candidate manifold location; test whether the broad component has weakened or strengthened.','Discovery test; no post-spectrum optical trigger established.'),
 ('P8548','2026-09-23 20:46','Compare the 2014 and 2021 Hbeta/Halpha profiles with the new spectrum; W1 spans about 0.43 mag.','Low-redshift line comparison; the IR amplitude alone does not time a transition.'),
 ('P9227','2026-09-23 21:16','W1 became about 0.35 mag fainter after the 2017 reference; test whether broad Hbeta and Halpha are weaker.','Dated IR trigger; allow for a dust lag and the end of NEOWISE in 2024.'),
 ('P9506','2026-09-23 21:46','W1 became about 0.41 mag fainter after the 2014 reference; compare broad Balmer emission with that baseline.','Dated IR trigger; no automatic prediction of the present spectral state.'),
 ('P10381','2026-09-23 22:16','Use the repeated broad-line spectra to measure the present state and check the unverified literature association.','Whole-spectrum flux offsets occur among epochs; these do not establish a changing-look event.'),
 ('P11113','2026-09-23 22:46','Compare broad Hbeta and Halpha with the 2018 spectrum in a candidate manifold region.','Manifold test; the recent IR brightening slope is not a dated transition prediction.'),
 ('P10596','2026-09-23 23:16','Compare the broad-Hbeta profile across the 2000-2024 history; test persistence or further change relative to the recent baseline.','Raw overlays suggest a weaker broad component than in 2000; aperture and flux matching are still required.'),
 ('P12457','2026-09-23 23:56','Search for a new broad Balmer component against the narrow-line-dominated 2000 and 2021 spectra.','Turn-on hypothesis; host dilution makes continuum S/N an incomplete measure of broad-line sensitivity.'),
]
HEADERS=['name','RA','DECL','slitwidth','exptime','nexp','binspat','binspect','slitangle','airmass_max','Note','Comment']
FIELD_NOTES={
 'P659':'Faint neighbours at 8.2arcsec PA149deg and 11.9arcsec PA111deg; check slit and sky apertures.',
 'P8548':'Resolved host; center on the compact nucleus. Archival aperture light includes the host.',
 'P9227':'Compact central source; no prominent close companion in the 40arcsec image.',
 'P9506':'r=19.35 star at 19.2arcsec PA218deg; verify nuclear centering and sky apertures.',
 'P10381':'r=19.35 extended neighbour at 14.9arcsec PA234deg; inspect slit orientation.',
 'P11113':'Compact central nucleus; no catalogue neighbour flag. Check the wider field before acquisition.',
 'P10596':'Compact central source; no prominent close companion in the 40arcsec image.',
 'P12457':'Resolved host and narrow-line-dominated reference; center on nucleus and inspect host subtraction.',
 'P3642':'Resolved host with a close source south of the nucleus and another source near the west edge of the 40arcsec field; inspect slit PA and sky apertures.',
 'P7281':'Compact central source in the 40arcsec image; check the wider field and acquisition centering.',
 'P7837':'Central source has visible extended light; center on the nucleus and check host contribution.',
 'P11010':'A field source lies north-east of the nucleus; keep it out of the extraction and sky apertures.',
 'P11082':'Compact central source in the 40arcsec image; inspect the wider field before acquisition.',
 'P10961':'Compact central source in the 40arcsec image; inspect the wider field before acquisition.',
 'P10983':'Visible field sources north-west and south-west of the nucleus; verify slit PA and uncontaminated sky apertures.',
}

def stamp(s):return pd.Timestamp(s,tz=TZ)
def utc(t):return t.tz_convert('UTC').strftime('%Y-%m-%dT%H:%M:%SZ')
def local(t):return t.tz_convert(TZ).strftime('%Y-%m-%d %H:%M')
def payload(path):return json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>',path.read_text(),re.S).group(1))

def geometry(ra,dec,start,end):
    times=Time(pd.date_range(start,end,freq='20s').union(pd.DatetimeIndex([end])).to_pydatetime())
    frame=AltAz(obstime=times,location=OBS.PALOMAR.location,pressure=0*u.hPa)
    aa=SkyCoord(ra*u.deg,dec*u.deg).transform_to(frame)
    moon=get_body('moon',times,OBS.PALOMAR.location).transform_to(frame)
    assert (aa.alt.deg>0).all()
    return dict(airmass_start=float(aa.secz.value[0]),airmass_end=float(aa.secz.value[-1]),
                airmass_max_actual=float(aa.secz.value.max()),moon_min=float(aa.separation(moon).deg.min()))

def csv_text(rows):
    # The documented NGPS parser does not interpret quoted CSV fields.
    lines=[','.join(HEADERS)]
    for row in rows:
        values=[str(row[h]) for h in HEADERS]
        assert len(row['Note'])<=24 and len(row['Comment'])<=1024
        assert all(not any(c in v for c in [',','\n','\r','"']) and v.isascii() for v in values)
        lines.append(','.join(values))
    return '\n'.join(lines)+'\n'

def ngps_row(target,plan,note,comment):
    c=SkyCoord(target['ra']*u.deg,target['dec']*u.deg)
    return dict(name=target['name'],RA=c.ra.to_string(unit=u.hourangle,sep=':',precision=3,pad=True),
        DECL=c.dec.to_string(unit=u.deg,sep=':',precision=2,pad=True,alwayssign=True),
        slitwidth='SET 1.0',exptime=f"SET {plan['seconds_each']}",nexp=plan['exposures'],
        binspat=2,binspect=2,slitangle='PA',airmass_max=plan['airmass'],Note=note,Comment=comment)

def window_notes(w,plan,tier):
    result=[];starts=[]
    for r in w['tier_ranges'][tier]:
        begin=pd.Timestamp(r['start_utc'],tz='UTC').tz_convert(TZ).ceil('min')
        end=pd.Timestamp(r['end_utc'],tz='UTC').tz_convert(TZ).floor('min')
        latest=(end-pd.Timedelta(minutes=plan['visit_minutes'])).floor('min')
        if latest<begin:continue
        starts.append((begin,latest))
        result.append(dict(start=local(begin),end=local(end),latest_visit_start=local(latest)))
    return result,starts

def make_visit(target,w,plans,start,max_minutes=30):
    for tier,limit,_ in OPTIONS.TIERS:
        if limit>1.8:break
        plan=plans.get(tier)
        if plan is None:continue
        if plan['visit_minutes']>max_minutes:continue
        end=start+pd.Timedelta(minutes=plan['visit_minutes'])
        if not any(pd.Timestamp(r['start_utc'],tz='UTC')<=start and pd.Timestamp(r['end_utc'],tz='UTC')>=end for r in w['tier_ranges'][tier]):continue
        g=geometry(target['ra'],target['dec'],start,end)
        assert g['airmass_max_actual']<=limit+1e-7 and g['moon_min']>=40-1e-7
        windows,starts=window_notes(w,plan,tier)
        latest=next(b for a,b in starts if a<=start<=b)
        return dict(name=target['name'],tier=tier,plan=plan,start_pdt=local(start),end_pdt=local(end),
                    start_utc=utc(start),end_utc=utc(end),latest_start_pdt=local(latest),windows=windows,**g)
    return None

def build():
    opts=json.loads((OUT/'airmass_options.json').read_text())
    revision=json.loads((DEST/'snr5_plan.json').read_text())
    assert revision['settings']['goal_snr']==5
    opts['plans']=revision['plans']
    targets=pd.read_csv(OUT/'compact_review_objects.csv').set_index('name',drop=False)
    science=pd.read_csv(OUT/'three_night_review/science_and_sensitivity.csv').set_index('name')
    local_data=payload(OUT/'candidate_review_local.html');public_data=payload(ROOT/'docs/index.html')
    local_targets={t['name']:t for t in local_data['targets']};public_names={t['name'] for t in public_data['targets']}
    public_targets={t['name']:t for t in public_data['targets']}
    original={v[0]:v for v in ORIGINAL_PRIMARY}
    sequence_choices=revision['sequence']
    assert all(v['name'] in public_names for v in sequence_choices)
    names={v['name'] for v in sequence_choices};primaries=[];primary_csv=[];backups=[];backup_csv=[];used_backups=set()
    for index,choice in enumerate(sequence_choices,1):
        name=choice['name'];start=choice['start_pdt']
        question=original[name][2] if name in original else str(science.loc[name].science_question)+'. Compare the broad Balmer profile and its strength with the dated archival spectra.'
        caution=original[name][3] if name in original else 'Promoted from the backup list for a complete observing window. No confirmed state change is inferred from the manifold position.'
        t=targets.loc[name].to_dict();w=next(w for w in opts['windows'][name] if w['night']=='sep23')
        visit=make_visit(t,w,{choice['tier']:opts['plans'][name][choice['tier']]},stamp(start),choice['duration']);assert visit,(name,'primary does not fit')
        tag=f'P{index:02}';visit.update(rank=index,role='primary',science_question=question,caution=caution,field_note=FIELD_NOTES.get(name,str(t['field_notes'])),promoted=name not in original)
        s=science.loc[name];assert pd.isna(s.ztf_r_latest180_mag) or s.ztf_r_latest180_mag<19
        latest_date=max((e['date'] for e in public_targets[name]['epochs']),default='unknown')
        ref=Time(visit['plan']['reference_mjd'],format='mjd').strftime('%Y-%m-%d')
        comment=(f"PRIMARY {tag}; visit {visit['start_pdt']} to {visit['end_pdt']} PDT; UTC {visit['start_utc']} to {visit['end_utc']}; "
                 f"latest visit start {visit['latest_start_pdt']} PDT; X<={visit['plan']['airmass']}; Moon>={visit['moon_min']:.1f}deg during planned visit; "
                 f"{visit['plan']['exposures']}x{visit['plan']['seconds_each']}s; combined continuum SNR>=5; sky V=18; 10min overhead plus slot rounding; "
                 f"window "+' / '.join(f"{v['start']} to {v['end']} PDT" for v in visit['windows'])+
                 f"; archival r={t['r_planning']:.2f}; continuum reference {ref}; latest public spectral date {latest_date}; "
                 f"{question} {visit['field_note']} Weak broad-line nondetection needs deeper data.")
        primary_csv.append(ngps_row(t,visit['plan'],f"{tag} by {visit['latest_start_pdt'][-5:]}",comment))
        choices=[]
        for other in targets.to_dict('records'):
            n=other['name']
            if n in names:continue
            ss=science.loc[n]
            if pd.notna(ss.ztf_r_latest180_mag) and ss.ztf_r_latest180_mag>=19:continue
            other_w=next((v for v in opts['windows'][n] if v['night']=='sep23'),None)
            if not other_w:continue
            v=make_visit(other,other_w,opts['plans'].get(n,{}),stamp(start),visit['plan']['visit_minutes'])
            if not v:continue
            score=float(ss.review_order_score)+5*bool(other['balmer_pair_in_range'])+5*float(other['manifold_cl_neighbor_fraction'])
            score-=1000 if other['pool_role']=='reserve' else 0
            score-=3 if v['tier']=='extended' else 0
            choices.append((score,n,v,other,ss))
        choices.sort(key=lambda x:(x[1] in used_backups,-x[0],x[1]));assert len(choices)>=2,(tag,'fewer than two alternatives fitting the full slot')
        chosen=choices[:3]
        visit['backups']=[n for _,n,_,_,_ in chosen]
        for rank,(_,n,v,other,ss) in enumerate(chosen,1):
            used_backups.add(n);v.update(role='backup',replaces=tag,backup_rank=rank,pool_role=other['pool_role'],science_question=ss.science_question)
            backups.append(v)
            comment=(f"BACKUP ONLY for {tag}; choice {rank}; replacement visit {v['start_pdt']} to {v['end_pdt']} PDT; "
                     f"UTC {v['start_utc']} to {v['end_utc']}; latest visit start {v['latest_start_pdt']} PDT; "
                     f"X<={v['plan']['airmass']}; Moon>={v['moon_min']:.1f}deg during replacement visit; "
                     f"window "+' / '.join(f"{a['start']} to {a['end']} PDT" for a in v['windows'])+
                     f"; {v['plan']['exposures']}x{v['plan']['seconds_each']}s; continuum SNR>=5; sky V=18; 10min overhead plus rounding; archival r={other['r_planning']:.2f}; "
                     f"{ss.science_question}. Replace the primary; never run the whole backup list. Names can recur for different slots; skip objects already observed. Inspect field and PA.")
            backup_csv.append(ngps_row(other,v['plan'],f"B{index:02}{rank} for {tag}",comment))
        primaries.append(visit)
    # Two explicit standard visits. Short exposure settings require quicklook
    # saturation checks; the ten-minute blocks include acquisition and repeats.
    stds=pd.read_csv(ROOT/'data/standards_spectrophotometric.csv').set_index('name');standards=[];std_csv=[]
    for name,label,start,seconds in [('P330E','P330E','2026-09-23 20:06',30),('BD+28 4211','BD284211','2026-09-24 00:28',5)]:
        s=stds.loc[name];a=stamp(start);b=a+pd.Timedelta(minutes=10);g=geometry(s.ra,s.dec,a,b)
        assert g['airmass_max_actual']<1.8 and g['moon_min']>=40
        plan=dict(exposures=2,seconds_each=seconds,airmass=1.8)
        note='STD check saturation'
        comment=(f"CALSPEC {name}; visit {local(a)} to {local(b)} PDT; UTC {utc(a)} to {utc(b)}; "
                 f"2x{seconds}s INITIAL setting only; inspect peak counts and adjust to avoid saturation; same 1arcsec slit and 2x2 binning as science; "
                 f"X<1.8; Moon>={g['moon_min']:.1f}deg; calspec identifier {s.calspec}; allow acquisition and repeat checks within this ten-minute block.")
        row=ngps_row(dict(name=label,ra=s.ra,dec=s.dec),plan,note,comment);std_csv.append(row)
        standards.append(dict(name=name,csv_name=label,role='standard',start_pdt=local(a),end_pdt=local(b),start_utc=utc(a),end_utc=utc(b),plan=plan,**g))
    sequence=[standards[0],*primaries,standards[1]]
    for a,b in zip(sequence,sequence[1:]):assert a['end_utc']<=b['start_utc']
    t0,t1,_,_=OBS.night_window('2026-09-23','first')
    assert Time(sequence[0]['start_utc'])>=t0 and Time(sequence[-1]['end_utc'])<=t1
    all_csv=[std_csv[0],*primary_csv,std_csv[1]]
    files={'sep23_primaries_ngps.csv':csv_text(all_csv),'sep23_backups_ngps.csv':csv_text(backup_csv)}
    for filename,text in files.items():(DEST/filename).write_text(text,encoding='ascii')
    pd.DataFrame([dict(name=v['name'],role=v['role'],start_pdt=v['start_pdt'],end_pdt=v['end_pdt'],
        start_utc=v['start_utc'],end_utc=v['end_utc'],exposures=v['plan']['exposures'],seconds_each=v['plan']['seconds_each'],
        airmass_max=v['airmass_max_actual'],moon_min=v['moon_min'],backups=';'.join(v.get('backups',[]))) for v in sequence]).to_csv(DEST/'sep23_sequence.csv',index=False)
    packet=dict(night='2026-09-23',timezone='PDT (UTC-07:00)',primaries=primaries,backups=backups,sequence=sequence,
        conditional=[],settings=revision['settings'],reserved=revision['reserved'],promoted=revision['promoted'],
        files=files,full_page='candidate_review_local.html',status='Continuum SNR 5 screening sequence; inspect fields and quicklook data on the night.')
    (DEST/'packet.json').write_text(json.dumps(packet,indent=2))
    print(pd.read_csv(DEST/'sep23_sequence.csv').to_string(index=False),flush=True)
    print('Backups',len(backups),'unique',len(used_backups),flush=True)
    return packet,local_data,public_data

def images(packet):
    image_module=importlib.import_module('18_fetch_review_images')
    targets=pd.read_csv(OUT/'compact_review_objects.csv');names={v['name'] for v in packet['primaries']}
    with ThreadPoolExecutor(max_workers=4) as executor:
        results=list(executor.map(image_module.one,targets[targets.name.isin(names)].itertuples(index=False)))
    assert all(r['status']=='available' and r['field_arcsec']==40 for r in results),results
    (DEST/'cutout_manifest.json').write_text(json.dumps(results,indent=2))

def render(packet,source,private):
    OLD=importlib.import_module('07_make_webpage')
    wanted={v['name']:v for v in packet['primaries']};value=copy.deepcopy(source)
    public_names={t['name'] for t in source['targets']}
    value['targets']=[t for t in value['targets'] if t['name'] in wanted]
    assert len(value['targets'])==len(wanted)
    value['targets'].sort(key=lambda t:wanted[t['name']]['rank'])
    for t in value['targets']:
        v=copy.deepcopy(wanted[t['name']]);v['plan']={k:v['plan'][k] for k in ['exposures','seconds_each','visit_minutes','airmass','goal_snr'] if k in v['plan']}
        t['visit']=v
        path=ROOT/'data/cutouts'/f"{t['name']}_sdss.jpg"
        metadata=json.loads(path.with_suffix('.json').read_text());assert metadata['field_arcsec']==40
        t['cut40']='data:image/jpeg;base64,'+base64.b64encode(path.read_bytes()).decode()
        if not private:
            t.pop('science',None);t.pop('exposure_plans',None)
            t['visit']['backups']=[n for n in t['visit']['backups'] if n in public_names]
    value['packet']=copy.deepcopy(packet)
    if not private:
        # All primary identities, brightness selections and continuum
        # references are public. The backup packet can contain internal data.
        assert all(not wanted[n]['plan'].get('reference_private',False) for n in wanted)
        value['packet']['files']={'sep23_primaries_ngps.csv':packet['files']['sep23_primaries_ngps.csv']}
        value['packet'].pop('backups',None);value['packet'].pop('conditional',None)
        for v in value['packet']['primaries']+value['packet']['sequence']:
            v['plan']={k:val for k,val in v['plan'].items() if k in ['exposures','seconds_each','visit_minutes','airmass','goal_snr']}
            if 'backups' in v:v['backups']=[n for n in v['backups'] if n in public_names]
    charts=OLD.TEMPLATE[OLD.TEMPLATE.index('function mjdToYear'):OLD.TEMPLATE.index('/* ---------- manifold thumbnail')]
    charts+='\n'+(ROOT/'web/candidate_spectra.js').read_text()
    template=(ROOT/'web/sep23_primaries_template.html').read_text()
    encoded=json.dumps(value,separators=(',',':'),allow_nan=False).replace('<','\\u003c')
    page=template.replace('__PAYLOAD__',encoded).replace('__CHART_FUNCTIONS__',charts)
    if private:(DEST/'sep23_primaries_local.html').write_text(page)
    else:
        (ROOT/'docs/sep23_primaries.html').write_text(page)
        (ROOT/'web/sep23_primaries.html').write_text(page.replace("'index.html'","'clagn_night_sheet.html'"))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fetch-images',action='store_true');args=ap.parse_args()
    packet,local_data,public_data=build()
    if args.fetch_images:images(packet)
    render(packet,local_data,True);render(packet,public_data,False)
    report=f"# September 23 observing packet\n\n{len(packet['primaries'])} science primaries in observing order, each with two target-specific exposures for combined continuum S/N >=5 per 2-pixel spectral bin near observed Hbeta. Includes two ten-minute standard visits.\n\n"
    report+='The original eight primaries remain; promoted backups: '+', '.join(packet['promoted'])+'.\n\n'
    report+='Each visit includes ten minutes for acquisition/readout plus rounding up to five-minute slots. Exposures are rounded up to 30 seconds with a 60-second minimum per exposure. Twenty minutes from Sep24 00:08 to 00:28 PDT are reserved for extra depth or delays, before the final standard. This reserve is not a target row or automatic wait command.\n\n'
    report+='ETC assumptions: archival local Hbeta continuum brightness; sky V=18 mag/arcsec2; zenith seeing 1.3 arcsec at 500 nm; 1arcsec central slice; 2x2 binning; optimal point-source extraction; airmass evaluated conservatively at the selected 1.5 or 1.8 ceiling. The forward calculation includes both reads. A source 0.5 mag fainter can fall below S/N 5; the scenario is recorded per visit in packet.json. No broad-line detection significance is promised.\n\n'
    report+='Use the primary CSV in order. Backups are replacement choices that fit their associated primary slot. A target can appear for several slots with different exposure settings; choose the row for the slot being replaced and skip any target already observed. Never append the entire backup list to an automatic run. Standard exposure settings require saturation checks. Inspect the slit field and Quicklook data; deepen ambiguous potential turn-offs before classifying them.\n\n'
    report+='Public primary exposure references and identities are retained for a consistent shared page. Candidates needing private-only identity or continuum information remain eligible for the complete local backup packet. SDSS-V spectra remain available on the local page. The manifold is a selection prior, not a forecast of the current state.\n\n'
    report+='[Dark primary page](sep23_primaries_local.html) · [NGPS primary sequence](sep23_primaries_ngps.csv) · [NGPS backups](sep23_backups_ngps.csv) · [Detailed timing](sep23_sequence.csv)\n\n'
    report+=pd.read_csv(DEST/'sep23_sequence.csv').to_markdown(index=False)
    report+='\n\nFormat checked against https://caltechopticalobservatories.github.io/NGPS/users-manual/target-lists.html and https://caltechopticalobservatories.github.io/NGPS/users-manual/quick-start.html. CSV imports have not been exercised on the observatory installation.\n'
    (DEST/'README.md').write_text(report)

if __name__=='__main__':main()
