"""Read-only NGPS raw-data inventory and calibration QA for a nightly reduction."""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits


def section(text):
    match = re.fullmatch(r'\[(\d+):(\d+),(\d+):(\d+)\]', text)
    if match is None:
        return None
    x1, x2, y1, y2 = map(int, match.groups())
    if min(x1, y1) < 1 or x2 < x1 or y2 < y1:
        return None
    return np.s_[y1-1:y2, x1-1:x2]


def audit(raw, out):
    out.mkdir(parents=True, exist_ok=True)
    manifest, frames, detectors = [], [], []
    for path in sorted(raw.rglob('*')):
        if not path.is_file() or path.name == '.DS_Store':
            continue
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1024*1024), b''):
                digest.update(block)
        manifest.append(dict(path=str(path.relative_to(raw)), bytes=path.stat().st_size,
                             sha256=digest.hexdigest()))
        if path.suffix.lower() != '.fits':
            continue
        with fits.open(path, memmap=False, checksum=True) as hdus:
            h = hdus[0].header
            name = str(h.get('NAME', '')).strip()
            setup = path.name.startswith('focus') or h.get('IMGTYPE') == 'TEST'
            record = dict(file=path.name, name_raw=name, name=name.upper(),
                          type=h.get('IMGTYPE'), excluded_setup=setup,
                          shutter_seconds=h.get('SHUTTIME'), exptime_ms=h.get('EXPTIME'),
                          mjd=h.get('MJD'), slit_arcsec=h.get('SLITW'),
                          ra=h.get('RA'), dec=h.get('DECL'), airmass=h.get('AIRMASS'),
                          channels=','.join(x.name for x in hdus[1:]),
                          checksum_ok=all(x.verify_checksum() != 0 and x.verify_datasum() != 0 for x in hdus),
                          checksum_missing_hdus=','.join(x.name for x in hdus if x.verify_checksum() == 2))
            frames.append(record)
            for image in hdus[1:]:
                if image.name not in 'UGRI' or image.data is None:
                    continue
                ih = image.header
                data = np.asarray(image.data, dtype=float)
                data_section, bias_section = section(ih['DATASEC']), section(ih['BIASSEC'])
                science = data if data_section is None else data[data_section]
                overscan = None if bias_section is None else data[bias_section]
                detectors.append(dict(file=path.name, name=name.upper(), type=h.get('IMGTYPE'),
                    excluded_setup=setup, channel=image.name, binspat=ih.get('BINSPAT'),
                    binspec=ih.get('BINSPEC'), rows=data.shape[0], cols=data.shape[1],
                    datasec=ih['DATASEC'], biassec=ih['BIASSEC'],
                    valid_header_sections=data_section is not None and bias_section is not None,
                    overscan_median=float(np.median(overscan)) if overscan is not None else np.nan,
                    median=float(np.median(science)), p99=float(np.percentile(science, 99)),
                    maximum=float(np.max(science)),
                    pixels_ge_40000=int(np.sum(science >= 40000)),
                    pixels_ge_60000=int(np.sum(science >= 60000)),
                    pixels_eq_65535=int(np.sum(science == 65535))))
    pd.DataFrame(manifest).to_csv(out/'raw_sha256_manifest.csv', index=False)
    f = pd.DataFrame(frames)
    d = pd.DataFrame(detectors)
    f.to_csv(out/'raw_frames.csv', index=False)
    d.to_csv(out/'raw_detector_qa.csv', index=False)
    d[~d.excluded_setup].groupby(['channel', 'binspat', 'binspec', 'type']).size().rename('count').to_csv(out/'calibration_counts.csv')
    science = f[(f.type == 'SCI') & ~f.name.isin(['P330E', 'BD284211'])]
    summary = science.groupby('name').agg(exposures=('file','count'),
        integration_seconds=('shutter_seconds','sum'), airmass_min=('airmass','min'),
        airmass_max=('airmass','max'))
    summary.to_csv(out/'observed_targets.csv')
    (out/'audit_summary.json').write_text(json.dumps(dict(raw_root=str(raw.resolve()),
        files_hashed=len(manifest), fits_files=len(f), checksum_failures=f.loc[~f.checksum_ok,'file'].tolist(),
        science_targets=len(summary), science_exposures=len(science),
        standards=f[f.name.isin(['P330E','BD284211'])].file.tolist(),
        note='Raw pixels only: high-count totals include arcs, sky lines and cosmic rays; inspect extracted traces before excluding data.'), indent=2))
    print(summary.to_string())
    print('Checksum failures:', f.loc[~f.checksum_ok,'file'].tolist())
    print('Standard detector statistics:')
    print(d[d.name.isin(['P330E','BD284211'])][['file','channel','p99','maximum','pixels_ge_40000','pixels_ge_60000']].to_string(index=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw', type=Path, default=Path('sep23_data/shemmati'))
    p.add_argument('--out', type=Path, default=Path('sep23_data/reduction_20260924/audit'))
    args = p.parse_args()
    audit(args.raw, args.out)
