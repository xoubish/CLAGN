"""A reviewed September 23 sequence, time-matched alternatives, and dark cards.

Renders the plan written by 41_september_snr5_plan.py: one instrument setting for every
science target and visits placed by predicted slot S/N. It never uploads to NGPS, commits,
publishes, or silently changes the sequence as archives update.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,base64,copy,importlib,json,re
import numpy as np
import pandas as pd
import sep23_backups as backup_selection
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord,AltAz,get_body
from astropy.utils import iers

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/reselection_2026-09-20'
DEST=ROOT/'observing/sep23';DEST.mkdir(parents=True,exist_ok=True)
OBS=importlib.import_module('05_observability')
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
HOST_CONTAMINATED={'P8548','P12457'}
FIELD_NOTES={
 'P659':'Faint neighbours at 8.2arcsec PA149deg and 11.9arcsec PA111deg; check slit and sky apertures.',
 'P8548':'Resolved host; center on the compact nucleus. ETC S/N uses AGN plus host light and overstates AGN-only S/N.',
 'P9227':'Compact central source; no prominent close companion in the 40arcsec image.',
 'P9506':'r=19.35 star at 19.2arcsec PA218deg; verify nuclear centering and sky apertures.',
 'P10381':'r=19.35 extended neighbour at 14.9arcsec PA234deg; inspect slit orientation.',
 'P11113':'Compact central nucleus; no catalogue neighbour flag. Check the wider field before acquisition.',
 'P10596':'Compact central source; no prominent close companion in the 40arcsec image.',
 'P12457':'Host-dominated narrow-line reference; center on nucleus. ETC S/N uses AGN plus host light and overstates AGN-only S/N.',
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

def csv_text(rows, coordinates_only=False):
    # The documented NGPS parser does not interpret quoted CSV fields.
    columns=['name','RA','DECL'] if coordinates_only else HEADERS
    lines=[] if coordinates_only else [','.join(h.upper() for h in columns)]
    for row in rows:
        values=[str(row[h]) for h in columns]
        if not coordinates_only:
            assert len(row['Note'])<=24 and len(row['Comment'])<=1024,(row['name'],len(row['Comment']))
        assert all(not any(c in v for c in [',','\n','\r','"']) and v.isascii() for v in values)
        lines.append(','.join(values))
    return ''.join(line+'\n' for line in lines)

ASCII={'\u03b2':'beta','\u03b1':'alpha','\u03b3':'gamma','\u2013':'-','\u2014':'-','\u2033':'arcsec','\u2032':'arcmin','\u00c5':'A','\u2248':'~','\u2265':'>=','\u2264':'<=','\u00b0':'deg','\u00b1':'+/-','\u2026':'...','\u2019':"'",'\u201c':'"','\u201d':'"','\u00d7':'x','\u03bc':'u'}
def ascii_field(text,limit):
    """The documented NGPS CSV parser takes plain ASCII without commas or quotes."""
    text=str(text)
    for k,v in ASCII.items():text=text.replace(k,v)
    text=text.encode('ascii','ignore').decode().replace(',',';').replace('"',"'").replace('\n',' ').replace('\r',' ')
    return text[:limit].rstrip()

def ngps_row(target,plan,settings,note,comment):
    note=ascii_field(note,24);comment=ascii_field(comment,1024)
    c=SkyCoord(target['ra']*u.deg,target['dec']*u.deg)
    return dict(name=target['name'],RA=c.ra.to_string(unit=u.hourangle,sep=':',precision=3,pad=True),
        DECL=c.dec.to_string(unit=u.deg,sep=':',precision=2,pad=True,alwayssign=True),
        slitwidth=f"SET {settings['slit_arcsec']}",exptime=f"SET {plan['seconds_each']}",nexp=plan['exposures'],
        binspat=settings['binspat'],binspect=settings['binspect'],slitangle=settings.get('slitangle','PA'),
        airmass_max=plan['airmass'],Note=note,Comment=comment)

def make_visit(target,plans,slots,windows,start,require_goal=True):
    """A visit at a planned start; backups must meet the S/N floor, chosen primaries need only a feasible slot."""
    slot=slots.get(start)
    if not slot or (require_goal and not slot['meets_goal']):return None
    plan=plans.get(slot['tier'])
    if plan is None:return None
    a=stamp(start);b=a+pd.Timedelta(minutes=plan['visit_minutes'])
    g=geometry(target['ra'],target['dec'],a,b)
    assert g['airmass_max_actual']<=plan['airmass']+1e-6 and g['moon_min']>=40-1e-6,(target['name'],start)
    tier_windows=windows.get(slot['tier'],[])
    latest=next((w['latest_visit_start'] for w in tier_windows if w['start']<=start<=w['latest_visit_start']),start)
    return dict(name=target['name'],tier=slot['tier'],plan=plan,slot=slot,start_pdt=local(a),end_pdt=local(b),
                start_utc=utc(a),end_utc=utc(b),latest_start_pdt=latest,windows=tier_windows,
                snr_per_angstrom=slot['snr_per_angstrom'],sky_V=slot['sky_V'],seeing_arcsec=slot['seeing_arcsec'],
                airmass_mean=slot['airmass_mean'],**g)

def standards_comparison(primaries, targets, standards):
    """Same-band archival flux ratios for the actual primary list, using public spectra."""
    from spectral_utils import accepted_reference, archival_records
    from synphot import SourceSpectrum, SpectralElement, Observation, units
    from synphot.models import Empirical1D
    band = SpectralElement.from_filter('johnson_v')
    vega = SourceSpectrum.from_vega()
    rows = []
    for visit in primaries:
        name = visit['name']
        row = dict(name=name, archival_V=None, reference_mjd=None,
                   P330E_flux_ratio=None, BD284211_flux_ratio=None,
                   reference_access='public', status='No accepted public reference')
        result = accepted_reference(targets.loc[name].to_dict(),
                                    [r for r in archival_records(name) if not r.get('proprietary')])
        if result:
            record, _, _ = result
            w, f = np.asarray(record['wave'], float), np.asarray(record['flux'], float)
            good = np.isfinite(w) & np.isfinite(f)
            support = band.waveset.to_value(u.AA)[band(band.waveset).value > .001]
            coverage = good & (w >= support.min()) & (w <= support.max())
            if (coverage.sum() > 20 and w[good].min() <= support.min()
                    and w[good].max() >= support.max() and np.max(np.diff(w[coverage])) < 13):
                spectrum = SourceSpectrum(Empirical1D, points=w[good]*u.AA,
                                          lookup_table=f[good]*1e-17*units.FLAM)
                value = float(Observation(spectrum, band).effstim(units.VEGAMAG, vegaspec=vega).value)
                row.update(archival_V=value, reference_mjd=record['mjd'], status='Synthetic Johnson V; archival aperture flux',
                           P330E_flux_ratio=10**(.4*(value-float(standards.loc['P330E'].V))),
                           BD284211_flux_ratio=10**(.4*(value-float(standards.loc['BD+28 4211'].V))))
            else:
                row['status'] = 'Insufficient V-band coverage'
        rows.append(row)
    pd.DataFrame(rows).to_csv(DEST/'standards_comparison.csv', index=False)

def build():
    revision=json.loads((DEST/'snr5_plan.json').read_text());S=revision['settings']
    assert S['exposures']==2 and S['seconds_each']==300 and S['goal_snr']==5
    plans=revision['plans'];slots=revision['slots'];windows=revision['windows']
    targets=pd.read_csv(OUT/'compact_review_objects.csv').set_index('name',drop=False)
    science=pd.read_csv(OUT/'three_night_review/science_and_sensitivity.csv').set_index('name')
    local_data=payload(OUT/'candidate_review_local.html');public_data=payload(ROOT/'docs/index.html')
    public_names={t['name'] for t in public_data['targets']};public_targets={t['name']:t for t in public_data['targets']}
    original={v[0]:v for v in ORIGINAL_PRIMARY};protected=set(revision['protected']);added=set(revision['promoted']);user_pick=bool(revision.get('user_selection'))
    sequence_choices=revision['sequence']
    assert all(v['name'] in public_names for v in sequence_choices)
    names={v['name'] for v in sequence_choices};primaries=[];primary_csv=[];backups=[];backup_csv=[];used_backups=set()
    backup_candidates=backup_selection.candidates(targets, science, public_names, names)
    planner=importlib.import_module('41_september_snr5_plan')
    backup_geometry=planner.night_geometry(targets)
    backup_indices={name:i for i,name in enumerate(targets.index)}
    for index,choice in enumerate(sequence_choices,1):
        name=choice['name'];start=choice['start_pdt']
        if name in original:question,caution=original[name][2],original[name][3]
        else:
            question=str(science.loc[name].science_question)+'. Compare the broad Balmer profile and its strength with the dated archival spectra.'
            caution=((f"PI-selected primary; selection record updated {revision['user_selection'].get('decided', 'date unavailable')}." if user_pick else 'Promoted from the backup list for a complete observing window.') if name in protected else
                     'Filled automatically from the reviewed September pool by review score and slot S/N.')+' No confirmed state change is inferred from the manifold position.'
        t=targets.loc[name].to_dict()
        visit=make_visit(t,plans[name],slots[name],windows[name],start,require_goal=False);assert visit,(name,'primary does not fit')
        if name in HOST_CONTAMINATED:caution+=' ETC continuum includes host light; AGN-only S/N is lower and has not been estimated.'
        if visit['plan'].get('below_floor') or not visit['slot'].get('meets_goal',True):caution+=f" Predicted continuum S/N {visit['snr_per_angstrom']:.1f} per Angstrom is below the floor of 5 even with {visit['plan']['exposures']} exposures; a weak broad line will not be constrained."
        if t['pool_role']=='reserve':caution+=' Comparison reserve: a bright quasar outside the selected manifold regions, observed as a control.'
        tag=f'P{index:02}'
        visit.update(rank=index,role='primary',science_question=question,caution=caution,field_note=FIELD_NOTES.get(name,str(t['field_notes'])),
                     promoted=name not in original,added=name in added,protected=name in protected,host_contaminated=name in HOST_CONTAMINATED,
                     science_value=choice['science_value'],slot_quality=choice['slot_quality'])
        s=science.loc[name];assert pd.isna(s.ztf_r_latest180_mag) or s.ztf_r_latest180_mag<19
        latest_date=max((e['date'] for e in public_targets[name]['epochs']),default='unknown')
        ref='collaboration data' if visit['plan'].get('reference_private') else Time(visit['plan']['reference_mjd'],format='mjd').strftime('%Y-%m-%d')
        csv_question=original[name][2] if name in original else str(s.science_question).rstrip('.')+'.'
        comment=(f"Visit {visit['start_pdt']} to {visit['end_pdt'][-5:]} PDT; Moon min {visit['moon_min']:.1f}deg; "
                 f"model continuum S/N {visit['snr_per_angstrom']:.1f}/A near Hbeta; "
                 f"archival r={t['r_planning']:.2f}; continuum ref {ref}; latest public spectrum {latest_date}; "
                 f"{csv_question} {visit['field_note'].rstrip('.')}."
                 +(" Below S/N floor: weak broad lines unconstrained." if visit['snr_per_angstrom']<S['goal_snr'] else '')
                 +(" Setting target: skip after deadline in Note." if visit['airmass_end']>visit['airmass_start'] else ''))
        primary_csv.append(ngps_row(t,visit['plan'],S,f"{tag} by {visit['latest_start_pdt'][-5:]}",comment))
        chosen=backup_selection.select(backup_candidates, visit, backup_geometry, backup_indices,
                                       make_visit, OUT/'sep23_backup_etc')
        visit['backups']=[v['name'] for v,_ in chosen]
        for rank,(v,other) in enumerate(chosen,1):
            n=v['name']
            used_backups.add(n);v.update(role='backup',replaces=tag,backup_rank=rank,pool_role=other['pool_role'])
            backups.append(v)
            r_source='ZTF 180-day median' if v['eligibility']['r_source'].startswith('ZTF') else 'archival'
            field=FIELD_NOTES.get(n,str(other['field_notes'])).rstrip('.')
            comment=(f"Replacement visit {v['start_pdt']} to {v['end_pdt'][-5:]} PDT; "
                     f"latest start {v['latest_start_pdt']} PDT; Moon min {v['moon_min']:.1f}deg; "
                     f"model continuum S/N {v['snr_per_angstrom']:.1f}/A near Hbeta; "
                     f"r={v['eligibility']['r_mag']:.2f} ({r_source}); public ref MJD {v['eligibility']['reference_mjd']:.0f}; "
                     f"{v['eligibility']['region']}. {field}. Replace only the slot in Note; skip if already observed.")
            backup_csv.append(ngps_row(other,v['plan'],S,f"B{index:02}{rank} for {tag}",comment))
        primaries.append(visit)
    # Two explicit standard visits. Short exposure settings require quicklook
    # saturation checks; the ten-minute blocks include acquisition and repeats.
    stds=pd.read_csv(ROOT/'data/standards_spectrophotometric.csv').set_index('name');standards=[];std_csv=[]
    # Exposure times sized for the 1.5arcsec slit and 2x3 binning: roughly 5,000-8,000 peak counts per binned pixel,
    # inside the manual's 1,000 to 40,000 range and near its preferred 10,000. Still initial values: inspect the first frame.
    for name,label,start,seconds in [('P330E','P330E','2026-09-23 20:06',60),('BD+28 4211','BD284211','2026-09-24 00:28',10)]:
        s=stds.loc[name];a=stamp(start);b=a+pd.Timedelta(minutes=10);g=geometry(s.ra,s.dec,a,b)
        assert g['airmass_max_actual']<1.8 and g['moon_min']>=40
        plan=dict(exposures=2,seconds_each=seconds,airmass=1.8)
        note='STD check saturation'
        comment=(f"CALSPEC {s.calspec}; visit {local(a)} to {local(b)[-5:]} PDT incl. acquisition/checks. "
                 "Initial exposure: inspect all arms before repeating; adjust to avoid saturation. "
                 "Quicklook spatial-profile counts: >=1000; aim ~10000; keep <40000.")
        row=ngps_row(dict(name=label,ra=s.ra,dec=s.dec),plan,S,note,comment);std_csv.append(row)
        standards.append(dict(name=name,csv_name=label,role='standard',start_pdt=local(a),end_pdt=local(b),start_utc=utc(a),end_utc=utc(b),plan=plan,**g))
    sequence=[standards[0],*primaries,standards[1]]
    standards_comparison(primaries, targets, stds)
    for a,b in zip(sequence,sequence[1:]):assert a['end_utc']<=b['start_utc'],(a['name'],b['name'])
    t0,t1,_,_=OBS.night_window('2026-09-23','first')
    assert Time(sequence[0]['start_utc'])>=t0 and Time(sequence[-1]['end_utc'])<=t1
    all_csv=[std_csv[0],*primary_csv,std_csv[1]]
    files={'sep23_primaries_ngps.csv':csv_text(all_csv),'sep23_backups_ngps.csv':csv_text(backup_csv)}
    public_csv=copy.deepcopy(all_csv)
    private_names={v['name'] for v in primaries if v['plan'].get('reference_private')}
    for row in public_csv:
        if row['name'] in private_names:
            row['Comment']=re.sub(r'model continuum S/N .*?; ', 'Continuum estimate available on local page; ', row['Comment'])
            row['Comment']=row['Comment'].replace(' Below S/N floor: weak broad lines unconstrained.', '')
    public_files={'sep23_primaries_ngps.csv':csv_text(public_csv), 'sep23_backups_ngps.csv':csv_text(backup_csv)}
    for exports, primary_rows in [(files, all_csv), (public_files, public_csv)]:
        science_rows=[row for row in primary_rows if row['name'] in {v['name'] for v in primaries}]
        exports['sep23_primaries_ngps_coordinates.csv']=csv_text(science_rows, coordinates_only=True)
        exports['sep23_backups_ngps_coordinates.csv']=csv_text(backup_csv, coordinates_only=True)
    for filename,text in files.items():(DEST/filename).write_text(text,encoding='ascii')
    pd.DataFrame([dict(name=v['name'],role=v['role'],start_pdt=v['start_pdt'],end_pdt=v['end_pdt'],
        start_utc=v['start_utc'],end_utc=v['end_utc'],exposures=v['plan']['exposures'],seconds_each=v['plan']['seconds_each'],
        airmass_max=round(v['airmass_max_actual'],3),moon_min=round(v['moon_min'],1),
        snr_per_angstrom=round(v['snr_per_angstrom'],1) if 'snr_per_angstrom' in v else np.nan,
        sky_V=round(v['sky_V'],2) if 'sky_V' in v else np.nan,backups=';'.join(v.get('backups',[]))) for v in sequence]).to_csv(DEST/'sep23_sequence.csv',index=False)
    # Run-level information for the observer's page: nights, calibrations, standards, procedure, science aim.
    from astroplan import moon_illumination
    nights=[]
    for key,(date,part) in OBS.NIGHTS.items():
        a,b,dusk,dawn=OBS.night_window(date,part)
        mid=a+(b-a)/2
        nights.append(dict(key=key,date=date,part='first half' if part=='first' else 'full night',
            window_pdt=f"{pd.Timestamp(a.utc.iso,tz='UTC').tz_convert(TZ):%H:%M} to {pd.Timestamp(b.utc.iso,tz='UTC').tz_convert(TZ):%H:%M}",
            twilight_pdt=f"{pd.Timestamp(dusk.utc.iso,tz='UTC').tz_convert(TZ):%H:%M} to {pd.Timestamp(dawn.utc.iso,tz='UTC').tz_convert(TZ):%H:%M} (18 deg)",
            moon_percent=int(round(100*float(moon_illumination(mid)))),
            status=('Sequence set: this page.' if key=='sep23' else 'Sequence not chosen yet; use the candidate explorer with this night selected.')))
    run=dict(nights=nights,
        setting=dict(text=f"{S['slit_arcsec']} arcsec slit, {S['binspat']}x{S['binspect']} binning (spatial x spectral), {S['seconds_each']} s sub-exposures, {S['exposures']} per target unless noted, slit angle at the parallactic angle, single slit (no slicer).",
                     why=[f"Under the 93 to 99 percent Moon and the 1.5 to 1.8 arcsec seeing expected at the slit, the ETC slit scan gives 1.0 arcsec 77 to 82 percent of the best S/N, 1.5 arcsec 89 to 94 percent, 2.0 arcsec 96 to 99 percent; 1.5 arcsec keeps R about 1650 and half the wavelength zero-point sensitivity of 2.0 arcsec.",
                          "2x3 is the documented binning for a 1.5 arcsec slit; it gains 23 to 25 percent S/N per Angstrom in the G channel (read noise 7.8 e) and keeps 2.3 bins per resolution element. Quicklook cosmic-ray defaults are tuned for 2x3.",
                          f"{S['exposures']}x{S['seconds_each']} s reaches continuum S/N 5 per Angstrom at Hbeta to about AB 18.6 to 19.1 depending on Moon distance; faint favourites get a third or fourth exposure instead of a different setting."]),
        calibrations=['Afternoon, at 2x3 binning: 3 ThAr and 3 FeAr arcs per channel (slit 1.5 arcsec is within the 2 arcsec limit for arcs), 7 biases per channel.',
                      'Dome flats at 1.5 arcsec slit + 2x3 binning: at least 5 per channel, 7 to 10 for U and G.',
                      'One configuration is used for science and standards, so no other slit or binning needs calibrating.',
                      'Internal arcs can be repeated between targets without moving the telescope if the wavelength solution drifts.'],
        standards=[dict(name=v['name'],calspec=stds.loc[v['name']].calspec,exposures=v['plan']['exposures'],seconds_each=v['plan']['seconds_each'],start_pdt=v['start_pdt'],
                        note='HST CALSPEC; Quicklook builds the sensitivity function only from CALSPEC stars and uses the standard closest in time.') for v in sequence if v['role']=='standard'],
        procedure=['Load the primary CSV through the GUI target-list import. The sequencer runs rows in stored order; rows with nexp above 1 expand into that many exposures; slitangle PA is recomputed per exposure.',
                   'The Note column gives the latest start for each primary. The sequencer waits for a target to drop below airmass_max, so a SETTING target that is late never runs: skip it and continue.',
                   'Acquisition is automatic on the ACAM (about 90 s). Check the slit view, then guiding, before the first exposure; inspect the wide field for the neighbours listed on each card.',
                   'Standards: take the first exposure, check counts in the Quicklook spatial-profile display in all four channels (at least about 1,000, preferably about 10,000, below about 40,000), then repeat or adjust.',
                   'Behind schedule: drop the last automatic fill first (P9506, then P9227), never a PI choice; the packed sequence has an 8-minute buffer before the closing standard at 00:28.',
                   'Ahead of schedule or a target fails acquisition: take the backup listed for that slot from the backup CSV (same setting, same start time); skip names already observed.',
                   'Quicklook writes spec1d and spec2d under the reduced-data directory within tens of seconds; look at Hbeta (and Halpha where in range) against the archival spectrum on the card before moving on. An ambiguous faint broad line needs deeper data; only add an exposure after checking the remaining schedule and acquisition/readout time.'],
        data=['Import Quicklook CSVs with scripts/12_ngps_ingest.py and supply --mjd for the exposure epoch. It screens EW changes; it cannot classify a broad-line transition. Run scripts/16_candidate_webpage.py then scripts/51_observer_page.py to show new spectra on the local page.',
              'Record start time, seeing, sky and any deviation from the sequence in the observing log; the Comment column is copied into the NGPS log automatically.'])
    # The science text of the observer page lives in web/observer_page_template.html (fillAbout), where it can
    # quote the live pool and sequence counts; the packet carries only the run logistics.
    packet=dict(backup_policy=backup_selection.POLICY, night='2026-09-23',timezone='PDT (UTC-07:00)',run=run,primaries=primaries,backups=backups,sequence=sequence,
        conditional=[],settings=S,reserved=revision['reserved'],promoted=revision['promoted'],protected=revision['protected'],
        user_selection=revision.get('user_selection'),demoted=revision.get('demoted',[]),waived_rules=revision.get('waived_rules',{}),
        original_primaries=revision['original_primaries'],files=files,public_files=public_files,full_page='observer_page_local.html',
        status=f"Fixed {S['exposures']}x{S['seconds_each']}s science sequence with a {S['slit_arcsec']}arcsec slit and {S['binspat']}x{S['binspect']} binning; existing order preserved and revalidated with corrected slot S/N; inspect fields and quicklook data on the night.")
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

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fetch-images',action='store_true');args=ap.parse_args()
    packet,local_data,public_data=build()
    if args.fetch_images:images(packet)
    # Pages retired on 2026-09-21: 16 builds the payloads and 51 renders the single
    # observer page. Run both after this script.
    S=packet['settings'];gaps=packet['reserved']['gaps'];revision_demoted=packet.get('demoted',[])
    gap_text=', '.join(f"{g['start_pdt'][11:]}-{g['end_pdt'][11:]} PDT ({g['minutes']} min)" for g in gaps) or 'none'
    added=', '.join(packet['promoted']) or 'none'
    deeper=[v for v in packet['primaries'] if v['plan']['exposures']!=S['exposures']]
    deeper_text=('; '+', '.join(f"{v['name']} {v['plan']['exposures']}x{S['seconds_each']} s ({v['plan']['visit_minutes']} min)" for v in deeper)+' for depth') if deeper else ''
    report=(f"# September 23 observing packet\n\n{len(packet['primaries'])} science primaries in observing order, all with one instrument setting: "
            f"{S['slit_arcsec']}arcsec slit, {S['binspat']}x{S['binspect']} binning (spatial x spectral), {S['exposures']}x{S['seconds_each']}-second exposures{deeper_text}. "
            f"A standard visit is {S['visit_minutes']} minutes: {S['exposures']*S['seconds_each']//60} minutes of integration plus {S['overhead_minutes']} minutes including slew, acquisition and normal two-exposure readouts. Extra exposures add 0.6 minute readout each; visits round up. "
            f"Two ten-minute standard visits bookend the sequence. Unscheduled time inside the science block: {gap_text}; total {packet['reserved']['minutes']} minutes. It absorbs delays or takes a backup and is not a target row.\n\n")
    chosen_text=(f"Current PI choices (selection updated {packet['user_selection'].get('decided', 'date unavailable')}): {', '.join(packet['protected'])}. Filled automatically: {added}. Demoted to backups from the previous packet: {', '.join(revision_demoted) or 'none'}. " if packet.get('user_selection') else f"Protected from the previous packet: {', '.join(packet['protected'])}. Added on 2026-09-21 to use the shorter visits: {added}. ")
    report+=(chosen_text+
             "The current primaries and their start times are locked under user_selection.json preserve_sequence; approved replacements are recorded in its revisions and geometry is revalidated before writing. The original order came from an integer program. "
             "Every eligible target's S/N in every candidate slot is in snr5_slot_table.csv.\n\n")
    report+=(f"ETC assumptions: official NGPS ETC, single slit, optimal point-source extraction, May 2026 read noise and plate scales; archival local Hbeta continuum brightness; "
             f"zenith seeing {S['seeing_zenith_500nm']} arcsec at 500 nm scaled by airmass^{S['seeing_airmass_power']} to the target; sky brightness from the Krisciunas and Schaefer (1991) moonlight model at the mid-visit Moon geometry (93 percent illumination). "
             f"The floor is a combined continuum S/N >= {S['goal_snr']} per Angstrom, close to the earlier 2-pixel-bin criterion; both reads are included. "
             "P8548 and P12457 include host light: these are total-continuum estimates and overstate AGN-only S/N; no numerical nuclear S/N is available. No broad-line detection significance is promised.\n\n")
    report+=('Use the primary CSV in order. Up to two backups per primary: public spectroscopic quasars with r<19 (latest available ZTF median, otherwise archival), Moon separation >40 degrees and airmass <1.5 throughout the full visit at Palomar. Each inherits the exact exposure count and duration of its primary, including overhead. Literature/Zeltyn regions rank first, followed by labelled on/off neighbour fraction, Balmer coverage and brightness. A public spectral reference and predicted S/N >=5 are required. Qualifying backups and their CSV appear on both pages. If fewer than two qualify, show the shortfall without relaxing the constraints. A target can appear for several slots; choose the row for the slot being replaced and skip any target already observed. '
             'Never append the entire backup list to an automatic run. Standard exposure settings require saturation checks. Inspect the slit field and Quicklook data; only add exposures if the remaining schedule permits; ambiguous broad-line states remain unclassified.\n\n')
    report+='Public primary identities and telescope settings are retained, including private-reference primaries. Private-reference S/N and science derivatives are omitted from the public payload and download. Local CSVs retain the complete information. SDSS-V spectra remain available on the local page. The manifold is a selection prior, not a forecast of the current state.\n\n'
    report+='[Observer page](../../data/reselection_2026-09-20/observer_page_local.html) · [NGPS primary sequence](sep23_primaries_ngps.csv) · [NGPS backups](sep23_backups_ngps.csv) · [Detailed timing](sep23_sequence.csv) · [Per-slot S/N table](snr5_slot_table.csv)\n\n'
    report+=pd.read_csv(DEST/'sep23_sequence.csv').to_markdown(index=False)
    report+='\n\nFormat checked against https://caltechopticalobservatories.github.io/NGPS/users-manual/target-lists.html and https://caltechopticalobservatories.github.io/NGPS/users-manual/quick-start.html. CSV imports have not been exercised on the observatory installation.\n'
    (DEST/'README.md').write_text(report)

if __name__=='__main__':main()
