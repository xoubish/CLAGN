"""Ingest NGPS quicklook CSVs and screen model-free EW changes (not a CL classification).

Usage: python scripts/12_ngps_ingest.py <directory | files...>
       [--name TARGET file.csv] [--mjd EXPOSURE_MJD]
Channels for a target are combined into six-Angstrom bins. Supply spectra from
one observing epoch per invocation. The optional MJD is the exposure epoch,
never the ingestion time. Rebuild with scripts/16_candidate_webpage.py followed
by scripts/51_observer_page.py; NGPS traces remain on the local page only.
"""
import os, re, sys, glob, json
from pathlib import Path
from spectral_utils import bin_indices, accepted_reference
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(HERE, 'data')
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
    idx = bin_indices(GRID, w)
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
        blue = (rest > c0) & (rest < c1) & np.isfinite(flux)
        red = (rest > c2) & (rest < c3) & np.isfinite(flux)
        step = np.median(np.diff(rest))
        covered = ml.sum() > 5 and np.ptp(rest[ml]) >= (l1-l0)-2*step and np.max(np.diff(rest[ml])) < 1.6*step
        if covered and blue.sum() >= 3 and red.sum() >= 3:
            cont = np.polyfit(rest[mc], flux[mc], 1); cfit = np.polyval(cont, rest[ml])
            out[nm] = float(np.nansum((flux[ml] / cfit - 1) * np.gradient(rest[ml]))) if np.all(cfit > 0) else np.nan
    return out


def verdict(name, ew_now, last):
    """Report each finite EW comparison separately; no broad-line claims from these indices."""
    parts = []
    changes = []
    for ln, col in (('Hb', 'EW_Hb_rest'), ('Ha', 'EW_Ha_rest')):
        now = ew_now.get(ln)
        ref = last.get(col) if last is not None else None
        if now is None or not np.isfinite(now) or ref is None or not np.isfinite(ref):
            continue
        parts.append(f'{ln} EW {ref:.1f} -> {now:.1f} A')
        if ref > 5 and now < 0.5*ref:
            changes.append(f'{ln} decreased >50%')
        elif now > 2*max(ref, 5):
            changes.append(f'{ln} increased >2x')
        elif ref > 1 and (now > 1.4*ref or now < ref/1.4):
            changes.append(f'{ln} changed >40%')
    call = ('unclassified: no measurable archival comparison' if not parts else
            '; '.join(changes) + '; inspect calibration and decompose lines' if changes else
            'no large EW change in measured indices; broad-line state unclassified')
    return f'{name}: ' + '; '.join(parts) + f'  =>  {call}'


def target_catalog():
    paths = [Path(DATA)/'reselection_2026-09-20/compact_review_objects.csv']
    paths += sorted(Path(DATA).glob('targets_*.csv'))
    tables = [pd.read_csv(p) for p in paths if p.exists()]
    if not tables:
        raise ValueError('No target catalog available')
    return pd.concat(tables, ignore_index=True).drop_duplicates('name', keep='first').set_index('name')


def main():
    args = sys.argv[1:]
    observed_mjd = None
    if '--mjd' in args:
        i = args.index('--mjd')
        observed_mjd = float(args[i+1]); del args[i:i+2]
        if not np.isfinite(observed_mjd) or observed_mjd <= 40000:
            sys.exit('Invalid exposure MJD')
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
    tg = target_catalog()
    unknown = set(pins) - set(tg.index)
    if unknown:
        sys.exit(f'Unknown target names: {sorted(unknown)}')
    if any(pd.isna(tg.at[n, 'z']) for n in pins):
        sys.exit('Pinned target lacks a valid catalog redshift')
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
        Path(OUT, f'{name}.json').write_text(json.dumps(dict(source='NGPS', mjd=observed_mjd, grid_version=2), indent=2))
        last = None
        reference = accepted_reference(dict(name=name, z=z))
        if reference:
            record, _, _ = reference
            indices = ew_indices(np.asarray(record['wave'], float), np.asarray(record['flux'], float), z)
            last = dict(mjd=record['mjd'], EW_Hb_rest=indices.get('Hb'), EW_Ha_rest=indices.get('Ha'))
        v = verdict(name, ew, last); print(v)
        rows.append(dict(name=name, files=len(fl), z=z, lam_min=float(np.nanmin(w)), lam_max=float(np.nanmax(w)), EW_Hb_rest=ew.get('Hb', np.nan),
                         EW_Ha_rest=ew.get('Ha', np.nan), EW_Hb_last=(last or {}).get('EW_Hb_rest', np.nan), EW_Ha_last=(last or {}).get('EW_Ha_rest', np.nan),
                         last_mjd=(last or {}).get('mjd', np.nan), verdict=v.split('=>')[-1].strip()))
    result = pd.DataFrame(rows)
    previous = Path(DATA)/'ngps_lines.csv'
    if previous.exists() and rows:
        old = pd.read_csv(previous)
        result = pd.concat([old[~old.name.isin(result.name)], result], ignore_index=True)
    if rows:
        result.to_csv(previous, index=False)
    print(f'\n{len(rows)} targets -> data/ngps_spectra/*.csv and data/ngps_lines.csv; now run scripts/16_candidate_webpage.py then scripts/51_observer_page.py to update local cards')


if __name__ == '__main__':
    main()
