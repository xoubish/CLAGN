"""Attach the October queue and uppercase CSV downloads to the observer page."""
import csv
import io
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PACKET=ROOT/'observing/oct26_27_2026'


def uppercase_header(text):
    rows=list(csv.reader(io.StringIO(text)))
    if not rows:return text
    rows[0]=[c.upper() for c in rows[0]]
    out=io.StringIO();writer=csv.writer(out,lineterminator='\n');writer.writerows(rows)
    return out.getvalue()


def read_table(name):
    with (PACKET/name).open() as f:
        return [{k.lower():v for k,v in row.items()} for row in csv.DictReader(f)]


def csv_table(rows):
    out=io.StringIO();keys=list(rows[0])
    writer=csv.writer(out,lineterminator='\n')
    writer.writerow([k.upper() for k in keys])
    writer.writerows([[r[k] for k in keys] for r in rows])
    return out.getvalue()


def attach_october(payload,private,public_reference_names):
    if not (PACKET/'validation.json').exists():return
    rows=read_table('all_80_targets.csv')
    safe_columns=['order','name','night','start_pdt','end_pdt','start_utc','latest_start_pdt',
        'ra','dec','z','observed_sep23','r_adopted',
        'brightness_source','brightness_epoch','pool_role','slit_arcsec','binspat','binspect',
        'exposures','seconds_each','overhead_minutes','visit_minutes','protected_followup','airmass_start',
        'airmass_end','airmass_max_actual','moon_min','airmass_limit']
    if not private:
        # The requested operational identities/settings are shown for all 80.
        # Reference-dependent measurements retain the existing local-only boundary.
        sanitized=[]
        for r in rows:
            s={k:r[k] for k in safe_columns}
            s['priority_reason']=r['priority_reason'] if r['name'] in public_reference_names or r['jwst_id'] else 'Prepared candidate; optical-state comparison'
            s['science_question']=r['science_question'] if r['name'] in public_reference_names or r['jwst_id'] else ''
            s['flags']='; '.join(f for f in r['flags'].split('; ') if 'S/N' not in f and 'ETC' not in f)
            s['target_label']=r['jwst_target']
            if r['jwst_id']:
                s['priority_reason']='Optical-state measurement'
                s['science_question']='Compare broad Balmer lines and continuum with archival spectra to establish the current optical state.'
            s['pool_role']='additional' if s['pool_role']=='JWST' else s['pool_role']
            s['brightness_source']=s['brightness_source'].replace('JWST history audit: ', '').replace('JWST proposal catalogue r','archival catalogue r')
            s['flags']=s['flags'].replace('JWST addition:', 'Additional target:').replace('JWST ', '')
            sanitized.append(s)
        rows=sanitized
    files=payload.setdefault('files',{})
    files['all_80_targets.csv']=csv_table(rows)
    files.pop('jwst_coverage.csv',None)
    for name in (['jwst_coverage.csv'] if private else [])+['standards.csv','free_time.csv']:
        files[name]=uppercase_header((PACKET/name).read_text())
    backups=read_table('slot_backups.csv')
    if not private:
        for b in backups:
            b.pop('science_score',None)
            b.pop('jwst_replacement_warning',None)
            if b['name'] not in public_reference_names:b['priority_reason']='Prepared alternative; optical-state comparison'
    files['slot_backups.csv']=csv_table(backups)
    for night in ['oct26','oct27']:
        selected=[r for r in rows if r['night']==night]
        files[f'{night}_sequence.csv']=csv_table(selected)
        for suffix in ['ngps','coordinates','standards_ngps','backups_ngps']:
            name=f'{night}_{suffix}.csv';text=(PACKET/name).read_text()
            if not private and suffix in ['ngps','backups_ngps']:
                parsed=list(csv.DictReader(io.StringIO(text)))
                by_name={r['name']:r for r in selected}
                for record in parsed:
                    if suffix=='ngps':
                        r=by_name[record['NAME']]
                        record['NOTE']=f"{int(r['order']):02d} SCI {r['start_pdt'][11:16]}"
                        comment=f"{r['start_pdt'][11:16]} PDT; latest start {r['latest_start_pdt'][11:16]}; r={float(r['r_adopted']):.2f} archival; Moon>={float(r['moon_min']):.1f}deg; X<={float(r['airmass_max_actual']):.2f}. {r['priority_reason']}. {r['flags']}"
                    else:
                        # Preserve the slot and geometry from the operational file;
                        # remove the reference-dependent prioritization sentence.
                        comment=record['COMMENT'].split('deg.')[0]+'deg. Prepared slot alternative; skip if already observed; review priorities before replacing a target.'
                    record['COMMENT']=comment.encode('ascii','ignore').decode().replace(',',';').replace('"',"'")[:1024].rstrip()
                text=csv_table(parsed)
            files[name]=text
        coords={r['NAME']:r for r in csv.DictReader(io.StringIO(files[f'{night}_ngps.csv']))}
        for r in selected:
            r['ra_hms']=coords[r['name']]['RA'];r['dec_dms']=coords[r['name']]['DECL']
    selected_names={r['name'] for r in rows}
    windows=read_table('visibility_windows.csv')
    files['october_visibility_windows.csv']=csv_table([w for w in windows if w['name'] in selected_names])
    summary=json.loads((PACKET/'validation.json').read_text())
    settings=dict(summary['settings'])
    if not private:
        settings['airmass_exception_max']=settings.pop('airmass_jwst_max')
    payload['october']=dict(rows=rows,backups=backups,
        standards=read_table('standards.csv'),nights=summary['nights'],settings=settings,
        review_date=summary.get('review_date'),september_performance=summary.get('september_performance',{}),
        membership_changes=[{k:r[k] for k in ['name','decision','previous_night','scheduled_night']}
            for r in read_table('selection_decisions.csv') if r['decision'] in ['added','removed']])
    if private:payload['october']['jwst']=read_table('jwst_coverage.csv')
    for name,text in list(files.items()):
        if not name.endswith('_coordinates.csv'):
            files[name]=uppercase_header(text)
    if not private:
        dest=ROOT/'docs/october';dest.mkdir(exist_ok=True)
        # Remove the former public audit as well as omitting it from new builds.
        (dest/'jwst_coverage.csv').unlink(missing_ok=True)
        names=['all_80_targets.csv','standards.csv','free_time.csv','slot_backups.csv','october_visibility_windows.csv']
        names += [f'{n}_{suffix}.csv' for n in ['oct26','oct27'] for suffix in ['sequence','ngps','coordinates','standards_ngps','backups_ngps']]
        for name in names:(dest/name).write_text(files[name])
