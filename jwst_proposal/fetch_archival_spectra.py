"""Download public SDSS spectra (all epochs with a public sas_url) for the five NGPS-observed
targets, matched by position in data/spectra_epochs_*.csv, into inputs/archival_spectra/<id>/.
Saves one CSV per epoch (wave_A, flux_1e-17, ivar) plus an index CSV."""
import glob, os, io
from pathlib import Path
import numpy as np, pandas as pd
from astropy.io import fits
from urllib.request import urlopen, Request

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
OUT = HERE/'inputs/archival_spectra'; OUT.mkdir(exist_ok=True)
IDS = ['F02', 'F04', 'F06', 'F07', 'R01']

def main():
    s = pd.read_csv(HERE/'inputs/jwst_sample_cycle6.csv'); five = s[s.id.isin(IDS)]
    ep = pd.concat([pd.read_csv(f, low_memory=False) for f in glob.glob(str(ROOT/'data/spectra_epochs*.csv'))], ignore_index=True)
    ep = ep[ep['sas_url'].astype(str).str.startswith('https://data.sdss.org')]
    index = []
    for _, r in five.iterrows():
        d = np.hypot((ep.ra-r.ra)*np.cos(np.radians(r.dec)), ep.dec-r.dec)*3600
        m = ep[d < 2].drop_duplicates('sas_url').sort_values('mjd')
        tdir = OUT/r.id; tdir.mkdir(exist_ok=True)
        for _, e in m.iterrows():
            name = f"{r.id}_mjd{int(e.mjd)}_{str(e.get('survey','')).lower()}.csv"
            path = tdir/name
            if not path.exists():
                try:
                    data = urlopen(Request(e.sas_url, headers={'User-Agent': 'clagn/1.0'}), timeout=120).read()
                    with fits.open(io.BytesIO(data)) as h:
                        t = h[1].data
                        pd.DataFrame(dict(wave_A=10**t['loglam'], flux=t['flux'], ivar=t['ivar'])).to_csv(path, index=False)
                except Exception as exc:
                    print('FAILED', r.id, e.mjd, e.sas_url, str(exc)[:80]); continue
            index.append(dict(id=r.id, internal_id=r.internal_id, mjd=float(e.mjd), survey=e.get('survey', ''),
                              instrument=e.get('instrument', ''), file=str(path.relative_to(HERE)), url=e.sas_url))
            print('ok', r.id, int(e.mjd), e.get('survey', ''))
    pd.DataFrame(index).to_csv(OUT/'index.csv', index=False); print('index rows', len(index))

if __name__ == '__main__':
    main()
