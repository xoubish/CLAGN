"""Package the explicitly requested September NGPS spectra for the observer page.

Only the 18 P330E-calibrated science products are exported, never raw frames or
unrelated private archival spectra. The manifest supports rebuilding without raw data.
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
REDUCTION=ROOT/'sep23_data/reduction_20260924'
ASSETS=ROOT/'docs/observed/sep23_p330e'


def build_observed():
    products=REDUCTION/'products_p330e'
    manifest=ASSETS/'manifest.json'
    if not (products/'target_summary.csv').exists():
        return json.loads(manifest.read_text()) if manifest.exists() else []
    with (products/'target_summary.csv').open() as f:
        targets={r['NAME']:r for r in csv.DictReader(f)}
    with (REDUCTION/'audit/raw_frames.csv').open() as f:
        frames=list(csv.DictReader(f))
    visits=[]
    for name,row in targets.items():
        assert row['FLUX_STANDARD']=='P330E'
        exposures=[f for f in frames if f['name']==name and f['type']=='SCI']
        assert len(exposures)==int(row['NEXP'])
        start=min(float(f['mjd']) for f in exposures)
        end=max(float(f['mjd'])+float(f['shutter_seconds'])/86400 for f in exposures)
        def time_at(mjd):return datetime(1858,11,17,tzinfo=timezone.utc)+timedelta(days=mjd)
        first=time_at(start);last=time_at(end)
        with (products/f'spectra/{name}.csv').open() as f:
            data=list(csv.DictReader(f))
        wave=[float(r['WAVE_VAC_HELIO_A']) for r in data]
        flux=[round(float(r['FLUX']),6) if r['MASK'].lower()=='true' and math.isfinite(float(r['FLUX'])) else None for r in data]
        assert all(b>a for a,b in zip(wave,wave[1:]))
        record=dict(name=name,start_mjd=start,start_utc=first.isoformat(),end_utc=last.isoformat(),
            start_pdt=first.astimezone(ZoneInfo('America/Los_Angeles')).strftime('%H:%M'),
            end_pdt=last.astimezone(ZoneInfo('America/Los_Angeles')).strftime('%H:%M'),
            exposures=int(row['NEXP']),exposure_seconds=round(float(row['EXPTIME_S']),3),
            snr={ch:round(float(row[f'SNR_{ch}']),1) for ch in 'UGRI'},flux_standard='P330E',
            telluric_standard='BD284211',plot=f'{name}.png',csv=f'{name}.csv',fits=f'{name}.fits',
            epoch=dict(label=f'{first.date()} · NGPS · P330E',date=str(first.date()),mjd=start,
                epoch_day=math.floor(start),instrument='NGPS',coadd=True,
                quality_note='P330E flux calibration; central slice; vacuum heliocentric wavelengths. Masked pixels omitted. Absolute slit losses are not included in statistical errors.',
                wave=wave,flux=flux))
        visits.append(record)
        ASSETS.mkdir(parents=True,exist_ok=True)
        for folder,extension in [('plots','png'),('spectra','csv'),('spectra','fits')]:
            shutil.copy2(products/f'{folder}/{name}.{extension}',ASSETS/f'{name}.{extension}')
    visits.sort(key=lambda v:v['start_mjd'])
    assert len(visits)==18 and sum(v['exposures'] for v in visits)==39
    for rank,v in enumerate(visits,1):v['order']=rank
    manifest.write_text(json.dumps(visits,separators=(',',':'),allow_nan=False)+'\n')
    return visits


def add_observed(payload,visits,dest):
    by_name={t['name']:t for t in payload['targets']}
    payload['observed']=[]
    for visit in visits:
        assert visit['name'] in by_name,visit['name']
        record={k:v for k,v in visit.items() if k!='epoch'}
        for field in ['plot','csv','fits']:
            record[field]=Path(os.path.relpath(ASSETS/visit[field],dest.parent)).as_posix()
        t=by_name[visit['name']]
        t['observed']=record
        epochs=t.setdefault('spec',{}).setdefault('epochs',[])
        epochs.append(visit['epoch'])
        epochs.sort(key=lambda e:e.get('mjd') or 0)
        t.setdefault('epochs',[]).append(dict(mjd=visit['start_mjd'],date=visit['epoch']['date'],src='NGPS',file_available=True))
        t['n_spec']=len({math.floor(e['mjd']) for e in t['epochs'] if e.get('mjd') is not None})
        payload['observed'].append(record)
    return payload
