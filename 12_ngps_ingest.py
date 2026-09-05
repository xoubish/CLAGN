"""
12_ngps_ingest.py  --  bring NGPS Quicklook spectra onto the night sheet and call the state change on the spot.

The Quicklook DRP (C. Fremling, manual 2025-07-02) writes /media/data_archive/<UTDATE>_reduced/ on the Quicklook machine:
  spec1d/   1-d spectra saved with the GUI's 'Extract/Save' button, "wavelength, flux CSV format", one file per channel
            (filenames carry the object NAME and the channel, e.g. *_R*, *_I*), flux-calibrated with the CALSPEC sensitivity function
  spec2d/   multi-extension FITS (data, sky model, wavelength solution, bad-pixel map, illumination flat)
Copy the directory off the mountain (scp -r ... <UTDATE>_reduced) and point this script at it.

For every CSV whose name contains a target name from data/targets_*.csv (or the NGPS lists), the channels are merged onto the
6 A grid used for the archival spectra, written to data/ngps_spectra/<name>.csv (wave_A, flux in 1e-17 erg/s/cm2/A, err if present),
and the same model-free rest-frame EW indices as 03d_fetch_spectra.py are measured (Hb 4800-4930 vs 4700-4790/5090-5150;
Ha 6480-6650 vs 6350-6450/6700-6800).  The result is compared with the latest archival epoch in data/spectra_lines.csv and a
verdict is printed per target; data/ngps_lines.csv keeps the numbers.  Then run 07_make_webpage.py (or ./rescore.sh) so the
spectra appear in each card's NGPS slot, overlaid on the archival epochs.

Usage: /opt/anaconda3/bin/python 12_ngps_ingest.py <UTDATE_reduced dir | spec1d dir | file.csv ...> [--name P5554 file.csv]
       --name pins one file (or several channel files) to a target when the filename does not carry the target name.
"""
import os, re, sys, glob
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, 'data')
OUT = os.path.join(DATA, 'ngps_spectra')
GRID = np.arange(3000.0, 10500.0, 6.0)           # observed-frame 6 A bins, same as 03d (which starts at 3600; NGPS reaches 3050)
EW_WINDOWS = {'Hb': (4800, 4930, 4700, 4790, 5090, 5150), 'Ha': (6480, 6650, 6350, 6450, 6700, 6800)}


def read_spec1d(path):
    """Two or three numeric columns (wavelength, flux[, error]); header line optional; comma, whitespace or tab separated."""
    raw = open(path, errors='ignore').read().splitlines()
    rows = []
    for line in raw:
        parts = re.split(r'[,\s;]+', line.strip())
        try:
            vals = [float(p) for p in parts if p != '']
        except ValueError:
            continue                                  # header / comment line
        if len(vals) >= 2:
            rows.append(vals[:3])
    if not rows:
        raise ValueError('no numeric rows')
    a = np.array([r + [np.nan] * (3 - len(r)) for r in rows], float)
    w, f, e = a[:, 0], a[:, 1], a[:, 2]
    if np.nanmedian(w) < 100:                          # microns
        w = w * 1e4
    return w, f, e


def to_grid(w, f, e):
    idx = np.searchsorted(GRID, w) - 1
    fb = np.full(len(GRID), np.nan); eb = np.full(len(GRID), np.nan)
    good = np.isfinite(f) & (idx >= 0) & (idx < len(GRID))
    for i in np.unique(idx[good]):
        m = good & (idx == i)
        fb[i] = np.median(f[m])
        if np.isfinite(e[m]).any():
            eb[i] = np.sqrt(np.nanmean(e[m] ** 2) / max(1, m.sum()))
    return fb, eb


def ew_indices(wave, flux, z):
    out = {}
    if z is None or not np.isfinite(z):
        return out
    rest = wave / (1 + z)
    for nm, (l0, l1, c0, c1, c2, c3) in EW_WINDOWS.items():
        ml = (rest > l0) & (rest < l1) & np.isfinite(flux)
        mc = (((rest > c0) & (rest < c1)) | ((rest > c2) & (rest < c3))) & np.isfinite(flux)
        if ml.sum() > 5 and mc.sum() > 5:
            cont = np.polyfit(rest[mc], flux[mc], 1); cfit = np.polyval(cont, rest[ml])
            out[nm] = float(np.nansum((flux[ml] / cfit - 1) * np.gradient(rest[ml]))) if np.all(cfit > 0) else np.nan
    return out


def verdict(name, ew_now, last):
    """Compare with the latest archival epoch: factor-of-two changes in the broad-line EW are called."""
    parts = []; call = 'no change'
    for ln, col in (('Hb', 'EW_Hb_rest'), ('Ha', 'EW_Ha_rest')):
        now = ew_now.get(ln); ref = last.get(col) if last is not None else None
        if now is None or not np.isfinite(now):
            continue
        if ref is None or not np.isfinite(ref):
            parts.append(f'{ln} EW {now:.0f} A (no archival index)'); continue
        ratio = (now + 1e-3) / (ref + 1e-3) if ref > 1 else np.nan
        parts.append(f'{ln} EW {ref:.0f} -> {now:.0f} A')
        if ref > 5 and now < 0.5 * ref:
            call = 'DIMMED: broad line lost >50% -> turn-off candidate'
        elif ref > 0 and now > 2 * max(ref, 5):
            call = 'BRIGHTENED: broad line gained >2x -> turn-on candidate'
        elif np.isfinite(ratio) and call == 'no change' and (now > 1.4 * ref or now < ref / 1.4):
            call = 'changed ~40%: re-observe or check flux calibration'
    return f'{name}: ' + '; '.join(parts) + f'  =>  {call}'


def main():
    args = sys.argv[1:]
    pins = {}
    while '--name' in args:
        i = args.index('--name'); pins.setdefault(args[i + 1], []).append(args[i + 2]); del args[i:i + 3]
    files = []
    for a in args:
        if os.path.isdir(a):
            files += glob.glob(os.path.join(a, 'spec1d', '*')) + glob.glob(os.path.join(a, '*.csv')) + glob.glob(os.path.join(a, '*.txt'))
        else:
            files += glob.glob(a)
    files = sorted(set(f for f in files if os.path.isfile(f)))
    tg = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(DATA, 'targets_*.csv'))]).drop_duplicates('name').set_index('name')
    names = sorted(tg.index, key=len, reverse=True)
    groups = {n: list(v) for n, v in pins.items()}
    for f in files:
        if any(f in v for v in groups.values()):
            continue
        base = os.path.basename(f)
        hit = next((n for n in names if re.search(r'(?<![A-Za-z0-9])' + re.escape(n) + r'(?![0-9])', base)), None)
        if hit:
            groups.setdefault(hit, []).append(f)
        else:
            print(f'   unmatched file (use --name <target> <file>): {base}')
    if not groups:
        sys.exit('no spec1d files matched a target name')
    lines = pd.read_csv(os.path.join(DATA, 'spectra_lines.csv')) if os.path.exists(os.path.join(DATA, 'spectra_lines.csv')) else None
    os.makedirs(OUT, exist_ok=True); rows = []
    for name, fl in sorted(groups.items()):
        W, F, E = [], [], []
        for f in fl:
            try:
                w, fx, e = read_spec1d(f)
            except Exception as ex:
                print(f'   {name}: cannot read {os.path.basename(f)}: {ex}'); continue
            W.append(w); F.append(fx); E.append(e)
        if not W:
            continue
        w = np.concatenate(W); fx = np.concatenate(F); e = np.concatenate(E)
        scale = 1e17 if np.nanmedian(np.abs(fx)) < 1e-8 else 1.0          # Quicklook writes cgs; the page uses 1e-17 cgs
        fb, eb = to_grid(w, fx * scale, e * scale)
        pd.DataFrame(dict(wave_A=GRID, flux=np.round(fb, 3), err=np.round(eb, 3))).dropna(subset=['flux']).to_csv(os.path.join(OUT, f'{name}.csv'), index=False)
        z = float(tg.at[name, 'z']) if name in tg.index and pd.notna(tg.at[name, 'z']) else np.nan
        ew = ew_indices(GRID, fb, z)
        last = None
        if lines is not None and (lines.name == name).any():
            L = lines[(lines.name == name) & lines.mjd.notna()].sort_values('mjd')
            last = L.iloc[-1].to_dict() if len(L) else None
        v = verdict(name, ew, last); print(v)
        rows.append(dict(name=name, files=len(fl), z=z, lam_min=float(np.nanmin(w)), lam_max=float(np.nanmax(w)), EW_Hb_rest=ew.get('Hb', np.nan),
                         EW_Ha_rest=ew.get('Ha', np.nan), EW_Hb_last=(last or {}).get('EW_Hb_rest', np.nan), EW_Ha_last=(last or {}).get('EW_Ha_rest', np.nan),
                         last_mjd=(last or {}).get('mjd', np.nan), verdict=v.split('=>')[-1].strip()))
    pd.DataFrame(rows).to_csv(os.path.join(DATA, 'ngps_lines.csv'), index=False)
    print(f'\n{len(rows)} targets -> data/ngps_spectra/*.csv and data/ngps_lines.csv; now run 07_make_webpage.py to update the cards')


if __name__ == '__main__':
    main()
