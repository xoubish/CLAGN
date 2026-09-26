"""Load the supplied P1823 SPHEREx dataset used in Figure 1.

Preserves the exported errors, spectral widths, epoch grouping and the single
previously requested P1823 exclusion. Input provenance is in inputs/provenance.json.
"""
from pathlib import Path
import csv
from functools import lru_cache
import numpy as np
from astropy.time import Time

HERE = Path(__file__).resolve().parent / 'inputs'
CSV_PATHS = {'P1823': HERE / 'spaxel_scryer_Shooby_AGN_ra245p1631_dec43p1373_cleaned.csv'}
PASS_COLORS = ['#2878a5', '#a14c7d', '#27836c']


def retained_measurements(prefix, data):
    """Exclude only the explicitly identified P1823 point; keep CSV unchanged."""
    keep = np.ones(len(data['flux_mjy']), dtype=bool)
    if prefix == 'P1823':
        rejected = (np.isclose(data['wavelength_um'], 4.06996, rtol=0, atol=1e-6)
                    & np.isclose(data['mjds'], 61075.067882, rtol=0, atol=1e-6))
        if rejected.sum() != 1:
            raise ValueError('Expected exactly one explicitly excluded P1823 point.')
        keep &= ~rejected
    return keep


@lru_cache(maxsize=1)
def load_spherex(prefix='P1823'):
    """Preserve every supplied row; group observations only by gaps >45 days."""
    with CSV_PATHS[prefix].open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    name = {'P1823': 'Shooby_AGN'}.get(prefix, prefix)
    if not rows or any(r['object_name'] != name or r['spectrum_type'] != 'cleaned'
                       or r.get('spectrum_role', 'source') != 'source' for r in rows):
        raise ValueError(f'Expected the cleaned source spectrum of {name}.')
    data = {key: np.array([float(r[key]) for r in rows]) for key in
            ['wavelength_um', 'wavelength_half_width_um', 'flux_mjy', 'flux_err_mjy', 'mjds']}
    if not all(np.isfinite(a).all() for a in data.values()):
        raise ValueError('Non-finite SPHEREx data require explicit review.')
    if np.any(data['flux_err_mjy'] <= 0) or np.any(data['wavelength_half_width_um'] < 0):
        raise ValueError('Invalid supplied uncertainty or spectral width.')
    if 'is_binned' in rows[0]:
        if any(r['is_binned'] != 'false' or r['n_images'] != '1' for r in rows):
            raise ValueError('Grouping assumes individual, unbinned measurements.')
    elif any(r['measurement_method'] != 'native_raw_mef_batch_exact_fractional'
             or not r['image_names'].endswith('.fits')
             or ';' in r['image_names'] or ',' in r['image_names'] for r in rows):
        raise ValueError('Expected one native image per exported measurement.')
    order = np.argsort(data['mjds'])
    cuts = np.flatnonzero(np.diff(data['mjds'][order]) > 45) + 1
    groups = []
    for i, idx in enumerate(np.split(order, cuts)):
        idx = idx[np.argsort(data['wavelength_um'][idx])]
        start, end = data['mjds'][idx].min(), data['mjds'][idx].max()
        first, last = Time(start, format='mjd').to_datetime(), Time(end, format='mjd').to_datetime()
        if (first.year, first.month) == (last.year, last.month):
            label = first.strftime('%b %Y')
        elif first.year == last.year:
            label = first.strftime('%b') + '–' + last.strftime('%b %Y')
        else:
            label = first.strftime('%b %Y') + '–' + last.strftime('%b %Y')
        groups.append(dict(indices=idx, start_mjd=float(start), end_mjd=float(end),
                           label=label, color=PASS_COLORS[i % len(PASS_COLORS)]))
    return data, groups, rows
