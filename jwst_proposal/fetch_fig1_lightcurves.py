"""Download the P9694 light-curve inputs used by make_fig1_connected.py.

Writes three CSV files under inputs/ (all public archives, no credentials):
  figure_allwise_mep_exposures.csv   AllWISE multi-epoch photometry, 2010-2011
  figure_ztf_irsa_lightcurve.csv     ZTF public data release, IRSA light-curve service
  figure_ztf_alerce_detections.csv   ZTF alert-stream detections via ALeRCE (extends the
                                     public release to the most recent alerts)
The cached NEOWISE single-exposure file (figure_neowise_exposures.csv) is unchanged;
a fresh IRSA download on 2026-09-28 was identical row for row.
"""
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

HERE = Path(__file__).resolve().parent
INPUTS = HERE / 'inputs'
RA, DEC = 11.5333, 9.1225          # P9694 = J004607.98+090720.9
RADIUS_DEG = 3 / 3600
SOURCES = dict(
    allwise_mep='https://irsa.ipac.caltech.edu/TAP/sync',
    ztf_irsa='https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves',
    alerce='https://api.alerce.online/ztf/v1/objects',
)


def get(url, params=None):
    if params:
        url = url + '?' + urlencode(params)
    with urlopen(Request(url, headers={'User-Agent': 'clagn-jwst-fig1/1.0'}), timeout=300) as r:
        return r.read().decode()


def allwise_mep():
    query = ('SELECT mjd,w1mpro_ep,w1sigmpro_ep,w2mpro_ep,w2sigmpro_ep,qi_fact,cc_flags,'
             'saa_sep,moon_masked,ra,dec FROM allwise_p3as_mep WHERE '
             f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA},{DEC},{RADIUS_DEG:.6f}))=1")
    text = get(SOURCES['allwise_mep'], dict(QUERY=query, RESPONSEFORMAT='csv'))
    if text.startswith('<?xml'):
        raise RuntimeError('IRSA TAP returned an error:\n' + text)
    path = INPUTS / 'figure_allwise_mep_exposures.csv'
    path.write_text(text)
    return path, len(pd.read_csv(path))


def ztf_irsa():
    text = get(SOURCES['ztf_irsa'], dict(POS=f'CIRCLE {RA} {DEC} {RADIUS_DEG:.5f}',
                                         BANDNAME='g,r', FORMAT='csv'))
    path = INPUTS / 'figure_ztf_irsa_lightcurve.csv'
    path.write_text(text)
    return path, len(pd.read_csv(path))


def alerce():
    objects = json.loads(get(SOURCES['alerce'] + '/', dict(ra=RA, dec=DEC, radius=5, page_size=20)))
    frames = []
    for item in objects['items']:
        rows = json.loads(get(f"{SOURCES['alerce']}/{item['oid']}/detections"))
        frame = pd.DataFrame(rows)
        frame.insert(0, 'oid', item['oid'])
        frames.append(frame)
    table = pd.concat(frames, ignore_index=True)
    keep = ['oid', 'candid', 'mjd', 'fid', 'magpsf', 'sigmapsf', 'magpsf_corr', 'sigmapsf_corr',
            'sigmapsf_corr_ext', 'corrected', 'dubious', 'isdiffpos', 'rb', 'drb', 'distnr',
            'ra', 'dec', 'step_id_corr']
    table = table[[c for c in keep if c in table.columns]].sort_values(['oid', 'mjd'])
    path = INPUTS / 'figure_ztf_alerce_detections.csv'
    table.to_csv(path, index=False)
    return path, len(table)


if __name__ == '__main__':
    for fn in (allwise_mep, ztf_irsa, alerce):
        path, n = fn()
        print(f'{path.name}: {n} rows')
