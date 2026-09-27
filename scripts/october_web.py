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
        'ra','dec','z','jwst_id','jwst_target','jwst_family','observed_sep23','r_adopted',
        'brightness_source','brightness_epoch','pool_role','slit_arcsec','binspat','binspect',
        'exposures','seconds_each','overhead_minutes','visit_minutes','airmass_start',
        'airmass_end','airmass_max_actual','moon_min']
    if not private:
        # The requested operational identities/settings are shown for all 80.
        # Reference-dependent measurements retain the existing local-only boundary.
        sanitized=[]
        for r in rows:
            s={k:r[k] for k in safe_columns}
            s['priority_reason']=r['priority_reason'] if r['name'] in public_reference_names or r['jwst_id'] else 'Prepared candidate; optical-state comparison'
            s['science_question']=r['science_question'] if r['name'] in public_reference_names or r['jwst_id'] else ''
            s['flags']='; '.join(f for f in r['flags'].split('; ') if 'S/N' not in f and 'ETC' not in f)
            sanitized.append(s)
        rows=sanitized
    files=payload.setdefault('files',{})
    files['all_80_targets.csv']=csv_table(rows)
    for name in ['jwst_coverage.csv','standards.csv','free_time.csv']:
        files[name]=uppercase_header((PACKET/name).read_text())
    backups=read_table('slot_backups.csv')
    if not private:
        for b in backups:
            b.pop('science_score',None)
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
                        comment=f"{r['start_pdt'][11:16]} PDT; latest start {r['latest_start_pdt'][11:16]}; r={float(r['r_adopted']):.2f} archival; Moon>={float(r['moon_min']):.1f}deg; X<={float(r['airmass_max_actual']):.2f}. {r['priority_reason']}. {r['flags']}"
                    else:
                        # Preserve the slot and geometry from the operational file;
                        # remove the reference-dependent prioritization sentence.
                        comment=record['COMMENT'].split('deg.')[0]+'deg. Prepared slot alternative; skip if already observed; reschedule any displaced JWST target.'
                    record['COMMENT']=comment.encode('ascii','ignore').decode().replace(',',';').replace('"',"'")[:1024]
                text=csv_table(parsed)
            files[name]=text
        coords={r['NAME']:r for r in csv.DictReader(io.StringIO(files[f'{night}_ngps.csv']))}
        for r in selected:
            r['ra_hms']=coords[r['name']]['RA'];r['dec_dms']=coords[r['name']]['DECL']
    selected_names={r['name'] for r in rows}
    windows=read_table('visibility_windows.csv')
    files['october_visibility_windows.csv']=csv_table([w for w in windows if w['name'] in selected_names])
    summary=json.loads((PACKET/'validation.json').read_text())
    payload['october']=dict(rows=rows,backups=backups,jwst=read_table('jwst_coverage.csv'),
        standards=read_table('standards.csv'),nights=summary['nights'],settings=summary['settings'])
    for name,text in list(files.items()):
        if not name.endswith('_coordinates.csv'):
            files[name]=uppercase_header(text)
    if not private:
        dest=ROOT/'docs/october';dest.mkdir(exist_ok=True)
        names=['all_80_targets.csv','jwst_coverage.csv','standards.csv','free_time.csv','slot_backups.csv','october_visibility_windows.csv']
        names += [f'{n}_{suffix}.csv' for n in ['oct26','oct27'] for suffix in ['sequence','ngps','coordinates','standards_ngps','backups_ngps']]
        for name in names:(dest/name).write_text(files[name])
