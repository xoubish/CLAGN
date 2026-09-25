"""SPHEREx spectra with expected redshifted line positions.

All three targets use their supplied cleaned numerical CSVs.
Rest wavelengths are vacuum values in microns. References and qualifications
are recorded in spherex_line_identifications.md.
"""
from pathlib import Path
import csv
import json
import hashlib
from functools import lru_cache
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time

HERE = Path(__file__).resolve().parent
A3_CSV = HERE / 'spaxel_scryer_P9694_ra11p5333_dec9p1225_cleaned.csv'
CSV_PATHS = {'A3': A3_CSV,
             'P2190': HERE / 'spaxel_scryer_P2190_ra240p2575_dec36p9388_cleaned.csv',
             'P1823': HERE / 'spaxel_scryer_Shooby_AGN_ra245p1631_dec43p1373_cleaned.csv'}
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
def load_spherex(prefix='A3'):
    """Preserve every supplied row; group observations only by gaps >45 days."""
    with CSV_PATHS[prefix].open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    name = {'A3': 'P9694', 'P1823': 'Shooby_AGN'}.get(prefix, prefix)
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


def load_a3_spherex():
    return load_spherex('A3')


def write_data_summary(prefix='A3'):
    data, groups, rows = load_spherex(prefix)
    keep = retained_measurements(prefix, data)
    source = CSV_PATHS[prefix]
    flags = ['incomplete_aperture', 'incomplete_annulus', 'asymmetric_annulus',
             'saturated_or_nonlinear']
    summary = dict(
        source_file=source.name, source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        target=rows[0]['object_name'], rows_supplied=len(rows), rows_plotted=int(keep.sum()),
        excluded_measurements=[dict(csv_row_index=int(i), wavelength_um=float(data['wavelength_um'][i]),
            mjd=float(data['mjds'][i]), flux_mjy=float(data['flux_mjy'][i]),
            reason='User-requested exclusion of isolated high point; input CSV preserved')
            for i in np.flatnonzero(~keep)],
        flux_column='flux_mjy', error_column='flux_err_mjy',
        horizontal_bars='wavelength_half_width_um: exported spectral half-widths, not wavelength errors',
        grouping='Chronological groups separated by gaps >45 days; no smoothing, rebinning or epoch averaging',
        observation_intervals=[dict(label=g['label'], count=len(g['indices']),
            start_utc=Time(g['start_mjd'], format='mjd').isot,
            end_utc=Time(g['end_mjd'], format='mjd').isot,
            wavelength_min_um=float(data['wavelength_um'][g['indices']].min()),
            wavelength_max_um=float(data['wavelength_um'][g['indices']].max())) for g in groups],
        flags_true={flag: sum(r[flag].lower() == 'true' for r in rows)
                    for flag in flags if flag in rows[0]},
        raw_flag_names=sorted({r.get('flag_names', '') for r in rows}),
        extraction_methods=sorted({r['measurement_method'] for r in rows}),
        uncertainty_methods=sorted({r.get('uncertainty_method', 'exported variance_total_mjy2') for r in rows}),
        local_background_values=sorted({r.get('local_background', r.get('background_mode', '')) for r in rows}),
        aperture_radius_arcsec=[min(float(r['aperture_radius_arcsec']) for r in rows),
                                max(float(r['aperture_radius_arcsec']) for r in rows)],
        interpretation='Expected line markers are not fitted detections. Epochs are plotted separately; this figure does not establish significant variability.')
    (HERE / f'inputs/{prefix}_spherex_csv_provenance.json').write_text(json.dumps(summary, indent=2) + '\n')


def write_a3_data_summary():
    write_data_summary('A3')


def line_groups(z):
    groups = [('Hα', ['Hα']), ('Paδ', ['Paδ']),
              ('He I + Paγ', ['He I', 'Paγ']), ('Paβ', ['Paβ']),
              ('Paα', ['Paα']), ('Brγ', ['Brγ']), ('Brβ', ['Brβ'])]
    if z > .6:
        groups.insert(0, ('Hβ / [O III]', ['Hβ', '[O III] 4959', '[O III] 5007']))
    return groups


def draw_csv(ax, prefix, z, fontsize):
    data, groups, _ = load_spherex(prefix)
    keep = retained_measurements(prefix, data)
    compact = fontsize < 8
    for group in groups:
        idx = group['indices'][keep[group['indices']]]
        ax.errorbar(data['wavelength_um'][idx], data['flux_mjy'][idx],
                    xerr=data['wavelength_half_width_um'][idx], yerr=data['flux_err_mjy'][idx],
                    fmt='o', ms=1.8 if compact else 3.2, color=group['color'],
                    ecolor=group['color'], elinewidth=.35 if compact else .65,
                    markeredgewidth=0, capsize=0, alpha=.85, label=group['label'], zorder=3)
    ymax = (18 if compact else 16) if prefix == 'A3' else (1.9 if compact else 1.8)
    if prefix == 'P1823':
        ymax = 2.2
    ax.set(xlim=(.70, 5.06), ylim=(0, ymax),
           xlabel='Observed wavelength (µm)', ylabel=r'$F_\nu$ (mJy)')
    ax.set_xticks(np.arange(1, 5.1, .5))
    ax.set_yticks([0, 5, 10, 15] if prefix == 'A3' else [0, .5, 1, 1.5])
    if prefix == 'P1823':
        ax.set_yticks([0, .5, 1, 1.5, 2])
    ax.grid(axis='y', color='#e3e7eb', lw=.5)
    ax.set_axisbelow(True)
    ax.tick_params(length=2.5, pad=2, labelsize=7 if compact else 9)
    ax.set_xlabel('Observed wavelength (µm)', fontsize=8 if compact else 10, labelpad=2)
    ax.set_ylabel(r'$F_\nu$ (mJy)', fontsize=8 if compact else 10, labelpad=3)
    legend_options = (dict(loc='lower right', bbox_to_anchor=(1, 1.035), ncol=3)
                      if prefix == 'P1823' and compact else
                      dict(loc='upper left', bbox_to_anchor=(.18, .67))
                      if prefix == 'P1823' else
                      dict(loc='upper right', bbox_to_anchor=(.88, 1) if prefix == 'P2190' else (1, 1)))
    ax.legend(frameon=False, fontsize=6.1 if compact else 8.5,
              handletextpad=.4, borderaxespad=.25, labelspacing=.3,
              **legend_options)
    for label, members in line_groups(z):
        waves = [REST[name] * (1 + z) for name in members]
        if not all(.75 <= wave <= 5 for wave in waves):
            continue
        for wave in waves:
            ax.axvline(wave, color='#a68143', lw=.55, ls=(0, (3, 3)), alpha=.7, zorder=1)
        ax.text(np.mean(waves) + .025, .965, label, rotation=90, va='top', ha='left',
                transform=ax.get_xaxis_transform(), fontsize=fontsize, color='#735625',
                bbox=dict(facecolor='white', edgecolor='none', pad=.3, alpha=.9))
    if not compact:
        top = ax.secondary_xaxis('top', functions=(lambda w: w / (1 + z), lambda w: w * (1 + z)))
        top.set_xlabel('Rest wavelength (µm)', fontsize=10, labelpad=5)
        top.set_xticks(np.arange(1, 4.1, .5))
        top.tick_params(labelsize=8, length=3)
LINES = [
    ('Hβ', .486268, 'SDSS'),
    ('[O III] 4959', .4960295, 'SDSS'),
    ('[O III] 5007', .5008240, 'SDSS'),
    ('Hα', .656461, 'SDSS'),
    ('Paδ', 1.0052, 'STScI NICMOS'),
    ('He I', 1.0833, 'STScI NICMOS'),
    ('Paγ', 1.0941, 'STScI NICMOS'),
    ('Paβ', 1.2822, 'STScI NICMOS'),
    ('Paα', 1.8756, 'STScI NICMOS'),
    ('Brδ', 1.9451, 'STScI NICMOS'),
    ('Brγ', 2.1661, 'STScI NICMOS'),
    ('Brβ', 2.6259, 'STScI NICMOS'),
    # Derived from the tabulated vacuum wavenumber; the old handbook's
    # wavelength column has a transposition typo (4.5225 instead of 4.05225).
    ('Brα', 10000 / 2467.765, 'STScI NICMOS vacuum wavenumber'),
]
REST = {name: wave for name, wave, _ in LINES}

# Pixel centres of plot borders, verified against the 1.0/5.0-micron gridlines.
# Both original screenshots use linear wavelength limits of 0.7 to 5.1 microns.
REGISTRATION = {
    'A3': dict(left=175.5, right=2303.5, top=2487.5, bottom=3335.5,
               view_top=.699, width=2344, height=3480),
    'P2190': dict(left=175.5, right=2195.5, top=2379.5, bottom=3227.5,
                  view_top=.69, width=2236, height=3372),
}


def draw_spherex_panel(ax, prefix, z, fontsize=8):
    if prefix in CSV_PATHS:
        draw_csv(ax, prefix, z, fontsize)
        return
    im = plt.imread(HERE / 'inputs' / f'{prefix}_spherex_inspector.png')
    reg = REGISTRATION[prefix]
    height, width = im.shape[:2]
    if (width, height) != (reg['width'], reg['height']):
        raise ValueError('Screenshot dimensions changed; recalibrate wavelength registration.')
    ax.imshow(im, origin='upper', interpolation='none', aspect='auto')
    ax.set_xlim(-.5, width - .5)
    ax.set_ylim(height - .5, reg['view_top'] * height)
    ax.axis('off')

    def pixel(wave):
        return reg['left'] + (wave - .7) / 4.4 * (reg['right'] - reg['left'])

    for label, members in line_groups(z):
        waves = [REST[name] * (1 + z) for name in members]
        if not all(.75 <= wave <= 5 for wave in waves):
            continue
        for wave in waves:
            ax.plot([pixel(wave)] * 2, [reg['top'] + 10, reg['bottom']],
                    color='#e6bc6a', lw=.6, ls=(0, (3, 3)), alpha=.8)
        label_x = pixel(sum(waves) / len(waves))
        # Narrow labels can be vertical without covering the underlying peaks.
        ax.text(label_x + 4, reg['top'] + 14, label, rotation=90,
                ha='left', va='top', color='#ffe1a6', fontsize=fontsize,
                bbox=dict(facecolor='#080808', edgecolor='none', alpha=.75, pad=.4))


def write_line_table(prefix, z):
    with (HERE / f'spherex_{prefix}_line_positions.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['line', 'rest_vacuum_um', 'redshift', 'observed_um',
                         'within_nominal_0p75_5um', 'interpretation', 'rest_wavelength_source'])
        for name, rest, source in LINES:
            observed = rest * (1 + z)
            writer.writerow([name, f'{rest:.7f}', z, f'{observed:.6f}',
                             .75 <= observed <= 5, 'expected position; not a measured detection', source])
