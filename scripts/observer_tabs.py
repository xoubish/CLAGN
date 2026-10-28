"""Catalogue download under About plus one target list per observing tab."""
import csv
import io
from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord

ROOT=Path(__file__).resolve().parents[1]


def csv_text(rows):
    out=io.StringIO()
    writer=csv.DictWriter(out,fieldnames=list(rows[0]),lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
    return out.getvalue()


def attach_tabs(payload):
    targets=payload['targets'];october=payload.get('october',{}).get('rows',[])
    inventory={t['name']:dict(name=t['name'],ra=t['ra'],dec=t['dec'],z=t.get('z'),
        r=t.get('science',{}).get('ztf_r_latest180_mag') or t.get('rmag'),
        target_label=t.get('jname',''),scheduled_night='',observed=bool(t.get('observed')),
        sep23=t['name'] in payload.get('sep23_sequence',{}),has_detail=True) for t in targets}
    for r in october:
        name=r['name']
        if name not in inventory:
            inventory[name]=dict(name=name,ra=float(r['ra']),dec=float(r['dec']),z=float(r['z']),
                r=float(r['r_adopted']),observed=False,sep23=False,has_detail=False)
        inventory[name].update(target_label=r.get('target_label','') or inventory[name].get('target_label',''),scheduled_night=r['night'])
    field_cards={r['name']:r for r in payload.get('field_cards',[])}
    for name,r in field_cards.items():
        inventory[name]=dict(name=name,ra=r['ra'],dec=r['dec'],z=None,r=None,
            target_label=r['label'],scheduled_night=r['night'],observed=False,sep23=False,has_detail=True)
    inventory=sorted(inventory.values(),key=lambda t:t['ra'])
    pool=[]
    for r in inventory:
        c=SkyCoord(r['ra']*u.deg,r['dec']*u.deg)
        r['ra_hms']=c.ra.to_string(unit=u.hourangle,sep=':',precision=3,pad=True)
        r['dec_dms']=c.dec.to_string(unit=u.deg,sep=':',precision=2,pad=True,alwayssign=True)
        pool.append(dict(NAME=r['name'],RA=r['ra_hms'],DECL=r['dec_dms'],SLITWIDTH='SET 1.5',
            EXPTIME='SET 300',NEXP=2,BINSPAT=2,BINSPECT=3,SLITANGLE='PA',AIRMASS_MAX=1.8,
            NOTE='POOL',COMMENT='Candidate catalogue; not an observing sequence. Use the selected night tab for exposure settings and order.'))
    for r in pool:
        if r['NAME'] in field_cards:
            f=field_cards[r['NAME']]
            r.update(EXPTIME=f"SET {f['seconds_each']}",NEXP=f['exposures'],SLITANGLE=f['pa_deg'])
    files=payload['files']
    files['all_targets.csv']=csv_text(pool)
    primaries=list(csv.DictReader(io.StringIO(files['sep23_primaries_ngps.csv'])))
    files['sep23_targets.csv']=csv_text([r for r in primaries if r['NAME'] in payload['sep23_sequence']])
    by_name={t['name']:t for t in inventory}
    observed=[]
    for o in payload['observed']:
        t=by_name[o['name']]
        observed.append(dict(ORDER=o['order'],NAME=o['name'],RA=t['ra_hms'],DECL=t['dec_dms'],
            START_UTC=o['start_utc'],END_UTC=o['end_utc'],START_PDT=o['start_pdt'],END_PDT=o['end_pdt'],
            NEXP=o['exposures'],TOTAL_EXPTIME_S=o['exposure_seconds'],FLUX_STANDARD=o['flux_standard'],
            **{f'SNR_{ch}':o['snr'][ch] for ch in 'UGRI'}))
    files['observed_targets.csv']=csv_text(observed)
    for night in ['oct26','oct27']:
        rows=[]
        for o in payload.get(f'observed_{night}',[]):
            t=by_name[o['name']];t['observed']=True
            rows.append(dict(ORDER=o['order'],NAME=o['name'],RA=t['ra_hms'],DECL=t['dec_dms'],
                START_UTC=o['start_utc'],END_UTC=o['end_utc'],START_PDT=o['start_pdt'],END_PDT=o['end_pdt'],
                NEXP=o['exposures'],TOTAL_EXPTIME_S=o['exposure_seconds'],FLUX_STANDARD=o['flux_standard'],
                **{f'SNR_{ch}':o.get('snr',{}).get(ch,'') for ch in 'UGRI'}))
        files[f'observed_{night}_targets.csv']=csv_text(rows)
    payload['all_targets']=inventory
    payload['tab_downloads']=dict(all='all_targets.csv',sep23='sep23_targets.csv',observed='observed_targets.csv',
        oct26='oct26_ngps.csv',observed_oct26='observed_oct26_targets.csv',oct27='oct27_ngps.csv',
        observed_oct27='observed_oct27_targets.csv')
    dest=ROOT/'docs/targets';dest.mkdir(exist_ok=True)
    for name in payload['tab_downloads'].values():(dest/name).write_text(files[name])
