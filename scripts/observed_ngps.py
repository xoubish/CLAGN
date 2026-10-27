"""Package the reduced NGPS spectra of each observed night for the observer page.

Only the delivered science products are exported (September 23: the 18 P330E-calibrated
spectra; October 26: the 43 Feige 34-calibrated spectra plus the 16 Cyg relative spectrum),
never raw frames or unrelated private archival spectra. Each night's manifest supports
rebuilding without the reduction directories.
"""
import csv
import json
import math
import os
import shutil
from datetime import datetime, timedelta, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

ROOT=Path(__file__).resolve().parents[1]
NIGHTS={
    'sep23':dict(reduction=ROOT/'sep23_data/reduction_20260924',products='products_p330e',
        assets=ROOT/'docs/observed/sep23_p330e',standard='P330E',standard_label='P330E',
        standard_note='P330E · 60 s standards',telluric='BD284211',targets=18,exposures=39,rename={},
        quality_note='P330E flux calibration; central slice; vacuum heliocentric wavelengths. Masked pixels omitted. Absolute slit losses are not included in statistical errors.'),
    'oct26':dict(reduction=ROOT/'oct26_data/reduction_20261027',products='products_feige34',
        assets=ROOT/'docs/observed/oct26_feige34',standard='FEIGE34',standard_label='Feige 34',
        standard_note='Feige 34 · 20 s standard (the opening BD+28 4211 lost light)',telluric='BD284211',
        targets=43,exposures=85,
        # The frames were labelled P8828 at the telescope; their coordinates are P8823's.
        rename={'P8828':'P8823'},
        quality_note='Feige 34 flux calibration; central slice; vacuum heliocentric wavelengths. Masked pixels omitted. Absolute slit losses are not included in statistical errors.'),
}
STELLAR=dict(name='Z16CygABC',products='products_16cyg',stem='16CYG_STAR',frames=['0051','0052','0053'],
    note='A and B fell along the central slice (one at the slice end, saturated in R/I); C cannot be separated from their PSF wings. Spectrum of the unsaturated central-slice star from the 1 s frames, fluxed with Feige 34, not telluric-corrected, about 2.3 mag of slit loss: use for shape and lines, not absolute flux.')


def utc(mjd):
    return datetime(1858,11,17,tzinfo=timezone.utc)+timedelta(days=mjd)


def visit_times(start,end):
    first=utc(start);last=utc(end);la=ZoneInfo('America/Los_Angeles')
    return dict(start_mjd=start,start_utc=first.isoformat(),end_utc=last.isoformat(),
        start_pdt=first.astimezone(la).strftime('%H:%M'),end_pdt=last.astimezone(la).strftime('%H:%M'),
        local_date=first.astimezone(la).strftime('%b %-d'),utc_date=first.strftime('%B %-d'))


def page_names():
    # Header names use the display form of the JWST-sample names (JCLAGN_0808, JP7837).
    with (ROOT/'jwst_proposal/inputs/jwst_sample_cycle6.csv').open() as f:
        return {(r['internal_id'] if r['internal_id'].startswith('J') else 'J'+r['internal_id']):r['internal_id']
                for r in csv.DictReader(f)}


def build_observed(night='sep23'):
    cfg=NIGHTS[night];reduction=cfg['reduction'];products=reduction/cfg['products'];assets=cfg['assets']
    manifest=assets/'manifest.json'
    if not (products/'target_summary.csv').exists():
        return json.loads(manifest.read_text()) if manifest.exists() else []
    with (products/'target_summary.csv').open() as f:
        targets={r['NAME']:r for r in csv.DictReader(f)}
    with (reduction/'audit/raw_frames.csv').open() as f:
        frames=list(csv.DictReader(f))
    summary=reduction/'audit/audit_summary.json'
    excluded=set(json.loads(summary.read_text()).get('excluded_files',[])) if summary.exists() else set()
    aliases=page_names()
    visits=[]
    assets.mkdir(parents=True,exist_ok=True)
    for product,row in targets.items():
        assert row['FLUX_STANDARD']==cfg['standard']
        name=cfg['rename'].get(product,product)
        name=aliases.get(name,name)
        exposures=[f for f in frames if f['name']==product and f['type']=='SCI' and f['file'] not in excluded]
        assert len(exposures)==int(row['NEXP'])
        start=min(float(f['mjd']) for f in exposures)
        end=max(float(f['mjd'])+float(f['shutter_seconds'])/86400 for f in exposures)
        times=visit_times(start,end)
        with (products/f'spectra/{product}.csv').open() as f:
            data=list(csv.DictReader(f))
        wave=[float(r['WAVE_VAC_HELIO_A']) for r in data]
        flux=[round(float(r['FLUX']),6) if r['MASK'].lower()=='true' and math.isfinite(float(r['FLUX'])) else None for r in data]
        assert all(b>a for a,b in zip(wave,wave[1:]))
        record=dict(name=name,night=night,**times,
            exposures=int(row['NEXP']),exposure_seconds=round(float(row['EXPTIME_S']),3),
            snr={ch:round(float(row[f'SNR_{ch}']),1) for ch in 'UGRI'},flux_standard=cfg['standard'],
            standard_label=cfg['standard_label'],standard_note=cfg['standard_note'],
            telluric_standard=cfg['telluric'],plot=f'{name}.png',csv=f'{name}.csv',fits=f'{name}.fits',
            epoch=dict(label=f"{utc(start).date()} · NGPS · {cfg['standard_label']}",date=str(utc(start).date()),mjd=start,
                epoch_day=math.floor(start),instrument='NGPS',coadd=True,quality_note=cfg['quality_note'],
                wave=wave,flux=flux,
                arm_breaks=[i for i in range(1,len(data)) if data[i]['CHANNEL']!=data[i-1]['CHANNEL']]))
        if product!=name and product in cfg['rename']:
            record['header_name']=product
        visits.append(record)
        for folder,extension in [('plots','png'),('spectra','csv'),('spectra','fits')]:
            shutil.copy2(products/f'{folder}/{product}.{extension}',assets/f'{name}.{extension}')
    assert len(visits)==cfg['targets'] and sum(v['exposures'] for v in visits)==cfg['exposures']
    if night=='oct26':
        visits.append(stellar_visit(reduction,frames,assets))
    visits.sort(key=lambda v:v['start_mjd'])
    for rank,v in enumerate(visits,1):v['order']=rank
    manifest.write_text(json.dumps(visits,separators=(',',':'),allow_nan=False)+'\n')
    return visits


def stellar_visit(reduction,frames,assets):
    products=reduction/STELLAR['products']
    exposures=[f for f in frames if f['file'][-9:-5] in STELLAR['frames']]
    start=min(float(f['mjd']) for f in exposures)
    end=max(float(f['mjd'])+float(f['shutter_seconds'])/86400 for f in exposures)
    name=STELLAR['name']
    for source,target in [(f"{STELLAR['stem']}.png",f'{name}.png'),(f"{STELLAR['stem']}.csv",f'{name}.csv'),
                          (f"{STELLAR['stem']}.fits",f'{name}.fits'),('spatial_profiles.png',f'{name}_profiles.png')]:
        shutil.copy2(products/source,assets/target)
    return dict(name=name,night='oct26',stellar=True,**visit_times(start,end),exposures=len(exposures),
        exposure_seconds=round(sum(float(f['shutter_seconds']) for f in exposures),3),
        flux_standard='FEIGE34',standard_label='Feige 34',standard_note='Feige 34 · relative spectrum only',
        plot=f'{name}.png',csv=f'{name}.csv',fits=f'{name}.fits',profiles=f'{name}_profiles.png',note=STELLAR['note'])


def add_observed(payload,visits,dest,night='sep23'):
    """Attach one night's visits; September keeps the original `observed` keys."""
    assets=NIGHTS[night]['assets']
    key='observed' if night=='sep23' else f'observed_{night}'
    by_name={t['name']:t for t in payload['targets']}
    photo=ROOT/'docs/observed/palomar-sep23.jpg'
    if photo.exists():
        payload['observed_photo']=Path(os.path.relpath(photo,dest.parent)).as_posix()
    payload[key]=[]
    for visit in visits:
        record={k:v for k,v in visit.items() if k!='epoch'}
        for field in ['plot','csv','fits','profiles']:
            if field in record:
                record[field]=Path(os.path.relpath(assets/visit[field],dest.parent)).as_posix()
        t=by_name.get(visit['name'])
        record['detail']=t is not None and not visit.get('stellar')
        if night=='sep23':
            assert t is not None,visit['name']
        if record['detail']:
            t[key]=record
            if not t.get('spec'):t['spec']={}
            epochs=t['spec'].setdefault('epochs',[])
            epochs.append(visit['epoch'])
            epochs.sort(key=lambda e:e.get('mjd') or 0)
            t.setdefault('epochs',[]).append(dict(mjd=visit['start_mjd'],date=visit['epoch']['date'],src='NGPS',file_available=True))
            t['n_spec']=len({math.floor(e['mjd']) for e in t['epochs'] if e.get('mjd') is not None})
        payload[key].append(record)
    return payload
