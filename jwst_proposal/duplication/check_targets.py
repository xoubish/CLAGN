"""Coordinate-based JWST duplication screening, including planned CAOM rows.

Follows STScI's MAST duplication notebook. No public-data, date, instrument,
or calibration-level restriction is applied. Matches require subsequent review.
"""
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import warnings

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from astroquery.mast import ObservationsClass

HERE = Path(__file__).resolve().parent
SAMPLE = HERE.parent/'inputs/jwst_sample_cycle6.csv'
CACHE = HERE/'cache'
RADIUS_ARCSEC = 120.


def clean(value):
    if np.ma.is_masked(value):
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def query(row):
    coord = SkyCoord(row['ra_hms'], row['dec_dms'], unit=(u.hourangle, u.deg))
    path = CACHE/f"{row['id']}.ecsv"
    if path.exists():
        obs = Table.read(path, format='ascii.ecsv')
        queried = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()
    else:
        client = ObservationsClass()
        client.TIMEOUT = 60
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', AstropyWarning)
            obs = client.query_criteria(obs_collection='JWST', coordinates=coord,
                                        radius=RADIUS_ARCSEC*u.arcsec)
        obs.write(path, format='ascii.ecsv', overwrite=True)
        queried = datetime.now(timezone.utc).isoformat()
    matches = []
    for r in obs:
        m = {k: clean(r[k]) for k in obs.colnames}
        m['sample_id'], m['sample_target'] = row['id'], row['target']
        if m.get('s_ra') is not None and m.get('s_dec') is not None:
            center = SkyCoord(m['s_ra'], m['s_dec'], unit=u.deg)
            m['center_separation_arcsec'] = float(coord.separation(center).arcsec)
        matches.append(m)
    summary = dict(id=row['id'], target=row['target'], family=row['family'],
                   ra_deg=float(coord.ra.deg), dec_deg=float(coord.dec.deg),
                   radius_arcsec=RADIUS_ARCSEC, queried_utc=queried,
                   mast_rows=len(matches),
                   planned_rows=sum(m.get('calib_level') == -1 for m in matches),
                   miri_rows=sum('MIRI' in str(m.get('instrument_name', '')).upper() for m in matches),
                   programs=';'.join(sorted({str(m['proposal_id']) for m in matches})),
                   instruments=';'.join(sorted({str(m['instrument_name']) for m in matches})),
                   review='Pending inspection of matches' if matches else 'No JWST matches in 120-arcsec search',
                   cache_file=str(path.relative_to(HERE)))
    (CACHE/f"{row['id']}_query.json").write_text(json.dumps(summary, indent=2)+'\n')
    return summary, matches


def main():
    CACHE.mkdir(exist_ok=True)
    with SAMPLE.open() as handle:
        sample = list(csv.DictReader(handle))
    summaries, matches, errors = [], [], []
    with ThreadPoolExecutor(max_workers=4) as pool:
        jobs = {pool.submit(query, row): row for row in sample}
        for job in as_completed(jobs):
            row = jobs[job]
            try:
                summary, found = job.result()
                summaries.append(summary)
                matches.extend(found)
                print(row['id'], row['target'], 'rows=', summary['mast_rows'],
                      'MIRI=', summary['miri_rows'], 'planned=', summary['planned_rows'],
                      'programs=', summary['programs'], flush=True)
            except Exception as exc:
                errors.append(dict(id=row['id'], error=repr(exc)))
                print('ERROR', row['id'], repr(exc), flush=True)
    summaries.sort(key=lambda s:s['id'])
    if summaries:
        with (HERE/'target_screening.csv').open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
            writer.writeheader(); writer.writerows(summaries)
    (HERE/'matches.json').write_text(json.dumps(matches, indent=2)+'\n')
    report = dict(queried_utc=datetime.now(timezone.utc).isoformat(),
                  method='astroquery.mast.Observations.query_criteria, obs_collection=JWST, 120 arcsec.',
                  source='https://spacetelescope.github.io/mast_notebooks/notebooks/JWST/duplication_checking/duplication_checking.html',
                  coordinate_columns=['ra_hms', 'dec_dms'],
                  input_sha256=hashlib.sha256(SAMPLE.read_bytes()).hexdigest(),
                  includes_planned_calib_level_minus1=True, public_only_filter=False,
                  targets_requested=len(sample), targets_queried=len(summaries),
                  matched_targets=sum(s['mast_rows'] > 0 for s in summaries),
                  miri_matched_targets=sum(s['miri_rows'] > 0 for s in summaries),
                  errors=errors, reviewed=False)
    (HERE/'screening_summary.json').write_text(json.dumps(report, indent=2)+'\n')
    if errors:
        raise RuntimeError('Incomplete query results; inspect screening_summary.json')


if __name__ == '__main__':
    main()
