"""Single-target Figure 1 layouts: A3 / P9694, P2190 and P1823.

Archival optical/WISE data are the bundled proposal snapshot. Both NGPS epochs use
the saved P330E calibration in the bundled FITS file. The SDSS JPEG
has its retrieval provenance alongside it. SPHEREx data use the supplied
cleaned CSVs, with per-measurement dates and uncertainties. A3 is the active proposal Figure 1; the legacy preview filename is retained.
"""
from pathlib import Path
import json
import os
import argparse

os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FixedLocator, FuncFormatter, NullFormatter
import numpy as np
from astropy.time import Time
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d
from spherex_line_labels import (draw_spherex_panel, write_line_table,
                                load_spherex, write_data_summary, line_groups, REST,
                                retained_measurements)

HERE = Path(__file__).resolve().parent
INPUTS = HERE / 'inputs'
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--target', choices=['A3', 'P2190', 'P1823'], default='A3')
args = parser.parse_args()
prefix = args.target
target = (next(t for t in json.loads((INPUTS / 'figure_targets.json').read_text())['targets']
               if t['name'] == 'P9694') if prefix == 'A3'
          else json.loads((INPUTS / f'{prefix}_figure_target.json').read_text()))
epochs = sorted(target['spec']['epochs'], key=lambda e: e['mjd'])
with fits.open(INPUTS / f'{prefix}_ngps_flux_standards.fits') as hdus:
    data, header = hdus[1].data, hdus[1].header
    assert header['TARGET'] == target['name'] and header['FLUXSTD'] == 'P330E'
    for epoch in epochs:
        if epoch.get('instrument') == 'NGPS':
            assert abs(epoch['mjd'] - header['MJD']) < 1e-6
            epoch['wave'] = data['WAVE_VAC_HELIO_A'].copy()
            epoch['flux'] = np.where(data['MASK'], data['FLUX'], np.nan)
            epoch['arm_breaks'] = (np.flatnonzero(data['CHANNEL'][1:] !=
                                                   data['CHANNEL'][:-1]) + 1).tolist()
            epoch['flux_standard'] = 'P330E'
z = target['z']
COLORS = ['#75a6c9', '#244c78', '#d25a2f']
if prefix == 'P1823':
    archival = [e for e in epochs if e.get('instrument') != 'NGPS']
    shades = iter(plt.cm.Blues(np.linspace(.32, .88, len(archival))))
    COLORS = ['#d25a2f' if e.get('instrument') == 'NGPS' else next(shades) for e in epochs]
    seen_years = set()
    for e in archival:
        yr = e['date'][:4]
        e['figure_label'] = yr if yr not in seen_years else '_nolegend_'
        seen_years.add(yr)
GREEN, PINK, PURPLE = '#178169', '#b44675', '#725397'
INK, MUTED, GRID = '#172733', '#596571', '#e3e7eb'

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['Times New Roman', 'DejaVu Serif'],
    'font.size': 8, 'axes.titlesize': 8.5, 'axes.labelsize': 8,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': '#9ba5ad', 'axes.linewidth': .6,
    'text.color': INK, 'axes.labelcolor': INK,
    'xtick.color': MUTED, 'ytick.color': MUTED,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})


def year(mjd):
    return Time(np.asarray(mjd, float), format='mjd').decimalyear


def spectrum(epoch):
    """Same NGPS display smoothing as the existing Figure 1; no rescaling."""
    wave = np.asarray(epoch['wave'], float)
    flux = np.asarray(epoch['flux'], float).copy()
    if epoch.get('instrument') == 'NGPS':
        good = np.flatnonzero(np.isfinite(flux))
        breaks = epoch.get('arm_breaks', [])
        cuts = np.flatnonzero((np.diff(good) > 1) | np.isin(good[1:], breaks)) + 1
        for run in np.split(good, cuts):
            if len(run) > 2:
                sigma = 6 / (2.35482 * np.median(np.diff(wave[run])))
                flux[run] = gaussian_filter1d(flux[run], sigma, mode='nearest')
        for b in breaks:
            if 0 <= b < len(flux):
                flux[b] = np.nan
    return wave, flux


def style(ax):
    ax.grid(axis='y', color=GRID, lw=.5)
    ax.set_axisbelow(True)
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.tick_params(length=2.5, pad=2)


fig = plt.figure(figsize=(6.5, 5.65), facecolor='white')
fig.text(.97, .971, f'z = {z:.3f}', ha='right', fontsize=9)

# Small identification image alongside two aligned histories.
ax_image = fig.add_axes([.065, .643, .242, .278])
ax_image.imshow(plt.imread(INPUTS / f'{prefix}_sdss.jpg'), extent=[-20, 20, -20, 20])
ax_image.axis('off')
ax_image.set_title('(a) SDSS gri', loc='left', pad=4)
half_field = 7.5
ax_image.set_xlim(-half_field, half_field)
ax_image.set_ylim(-half_field, half_field)
ax_image.plot([1, 6], [-6, -6],
              color='white', lw=1.5)
ax_image.text(3.5, -5.4, '5″',
              color='white', ha='center', fontsize=7)

ax_opt = fig.add_axes([.407, .788, .565, .139])
ax_ir = fig.add_axes([.407, .642, .565, .117], sharex=ax_opt)
for band, color in [('g', GREEN), ('r', PINK)]:
    values = np.asarray(target['ztf'][band], float)
    ax_opt.errorbar(year(values[:, 0]), values[:, 1], yerr=values[:, 2],
                    fmt='.', ms=1.9, elinewidth=.35, color=color, alpha=.8,
                    label=f'ZTF {band}')
ax_opt.invert_yaxis()
ax_opt.set_title('(b) Light curves', loc='left', pad=5)
ax_opt.set_ylabel('AB mag')
ax_opt.legend(loc='upper left', frameon=False, ncol=2, handletextpad=.2,
              columnspacing=.7, borderaxespad=.15)
plt.setp(ax_opt.get_xticklabels(), visible=False)

for key, markerface, label in [('wise', PURPLE, 'unWISE'), ('neo', 'white', 'NEOWISE')]:
    values = np.asarray(target['wise']['W1'] if key == 'wise' else target['neo'], float)
    ax_ir.errorbar(year(values[:, 0]), values[:, 1], yerr=values[:, 2],
                   fmt='o', ms=2.6, elinewidth=.5, color=PURPLE,
                   mfc=markerface, mew=.65, label=label)
ax_ir.set_ylabel('W1 (mJy)')
ax_ir.set_xlabel('Year', labelpad=1)
ax_ir.legend(loc='upper left', frameon=False, ncol=2, handletextpad=.3,
             columnspacing=.7, borderaxespad=.15, fontsize=6.5)
ax_ir.set_xlim(2009.6, 2027.0)
ax_ir.set_xticks([2010, 2014, 2018, 2022, 2026])
for ax in (ax_opt, ax_ir):
    style(ax)
    for epoch, color in zip(epochs, COLORS):
        label = epoch.get('label', '')
        if '–' in label and len(label.split('–')) == 2:
            start, end = label.split('–')
            ax.axvspan(Time(start).decimalyear, Time(end).decimalyear,
                       color=color, alpha=.12, lw=0)
        if 2009.6 <= year(epoch['mjd']) <= 2027:
            ax.axvline(year(epoch['mjd']), color=color, lw=.8, ls='--', alpha=.85, zorder=0)
    for group in load_spherex(prefix)[1]:
        ax.axvspan(year(group['start_mjd']), year(group['end_mjd']),
                   facecolor=group['color'], alpha=.20, lw=.6, edgecolor=group['color'])

OPTICAL_LINES = ([(3728.48, '[O II]'), (3869.86, '[Ne III]'),
                  (4102.89, 'Hδ'), (4341.68, 'Hγ'), (4862.68, 'Hβ'),
                  (5008.240, '[O III]'), (6302.046, '[O I]'),
                  (6564.61, 'Hα + [N II]'), (6725.48, '[S II]')]
                 if prefix == 'A3' else
                 [(2799.12, 'Mg II'), (3728.48, '[O II]'), (3869.86, '[Ne III]'),
                  (4102.89, 'Hδ'), (4341.68, 'Hγ'), (4862.68, 'Hβ'),
                  (5008.240, '[O III]')])


def optical_spectrum_mjy(epoch):
    wave, flux = spectrum(epoch)
    # F_lambda is in 1e-17 erg/s/cm^2/Angstrom; c is in Angstrom/s.
    # 1 mJy = 1e-26 erg/s/cm^2/Hz. Both axes are observed-frame.
    flux_mjy = flux * 1e-17 * wave**2 / 2.99792458e18 / 1e-26
    label = (epoch['date'] + ' NGPS' if epoch.get('instrument') == 'NGPS'
             else epoch.get('figure_label', epoch.get('label', epoch['date'])))
    if '–' in label:
        label = '2021–2022 coadd'
    return wave / 1e4, flux_mjy, label


def draw_combined_spectrum(ax, compact=True):
    """Observed wavelengths and F_nu in mJy, without inter-dataset scaling."""
    optical_handles = []
    for epoch, color in zip(epochs, COLORS):
        wave, flux_mjy, label = optical_spectrum_mjy(epoch)
        handle, = ax.plot(wave, flux_mjy, color=color,
                          lw=.65 if compact else .85, label=label, zorder=3)
        if label != '_nolegend_':
            optical_handles.append(handle)

    data, groups, _ = load_spherex(prefix)
    retained = retained_measurements(prefix, data)
    infrared_handles = []
    for group in groups:
        idx = group['indices'][retained[group['indices']]]
        handle = ax.errorbar(data['wavelength_um'][idx], data['flux_mjy'][idx],
                            xerr=data['wavelength_half_width_um'][idx],
                            yerr=data['flux_err_mjy'][idx], fmt='o',
                            ms=1.8 if compact else 2.8, color=group['color'],
                            elinewidth=.35 if compact else .55,
                            markeredgewidth=0, alpha=.85, label=group['label'], zorder=4)
        infrared_handles.append(handle)
    ax.set_xscale('log')
    ax.set(xlim=(.34, 5.08), ylim=(0, 16 if prefix == 'A3' else 2.1),
           xlabel='Observed wavelength (µm)', ylabel=r'$F_\nu$ (mJy)')
    ax.xaxis.set_major_locator(FixedLocator([.4, .5, .7, 1, 1.5, 2, 3, 4, 5]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks([0, 4, 8, 12, 16] if prefix == 'A3' else [0, .5, 1, 1.5, 2])
    if prefix == 'P1823':
        ax.set_ylim(0, 1.25 * max(np.max((data['flux_mjy'] + data['flux_err_mjy'])[retained]),
            max(np.nanmax(optical_spectrum_mjy(e)[1]) for e in epochs)))
        ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.grid(axis='y', color=GRID, lw=.5)
    ax.set_axisbelow(True)
    ax.tick_params(length=2.5, pad=2, labelsize=7 if compact else 9)
    ax.set_xlabel('Observed wavelength (µm)', fontsize=8 if compact else 10, labelpad=4)
    ax.set_ylabel(r'$F_\nu$ (mJy)', fontsize=8 if compact else 10, labelpad=3)
    optical_legend = ax.legend(handles=optical_handles, loc='upper left',
                               title='Optical spectra', title_fontsize=7 if compact else 9,
                               fontsize=6.5 if compact else 8, frameon=False,
                               borderaxespad=.4, labelspacing=.25)
    ax.add_artist(optical_legend)
    ax.legend(handles=infrared_handles, loc='upper right', title='SPHEREx',
              title_fontsize=7 if compact else 9, fontsize=6.5 if compact else 8,
              frameon=False, borderaxespad=.4, labelspacing=.25)

    for rest, label in OPTICAL_LINES:
        wave = rest * (1 + z) / 1e4
        ax.axvline(wave, ymax=.72, color='#b4bdc4', lw=.5, zorder=0)
        ax.text(wave * 1.009, .72, label, transform=ax.get_xaxis_transform(),
                rotation=90, va='top', fontsize=6 if compact else 8, color=MUTED)
    for label, members in line_groups(z):
        # Optical markers already identify these same observed wavelengths.
        if label == 'Hβ / [O III]' or (prefix == 'A3' and label == 'Hα'):
            continue
        waves = [REST[name] * (1 + z) for name in members]
        if not all(.75 <= w <= 5 for w in waves):
            continue
        for wave in waves:
            ax.axvline(wave, ymax=.72, color='#a68143', lw=.55,
                       ls=(0, (3, 3)), alpha=.7, zorder=1)
        ax.text(np.mean(waves) * 1.009, .72, label, transform=ax.get_xaxis_transform(),
                rotation=90, va='top', fontsize=6 if compact else 8, color='#735625')


# Separate spectral panels preserve the optical detail and infrared continuum.
ax_spec = fig.add_axes([.09, .362, .882, .174])
optical_flux = []
optical_xlim = (.34, .94) if prefix == 'P1823' else (.34, 1.04)
for epoch, color in zip(epochs, COLORS):
    wave, flux, label = optical_spectrum_mjy(epoch)
    ax_spec.plot(wave, flux, color=color, lw=.65, label=label)
    optical_flux.extend(flux[(wave >= optical_xlim[0]) & (wave <= optical_xlim[1]) & np.isfinite(flux)])
ax_spec.set(xlim=optical_xlim, ylim=(0, max(optical_flux) * 1.12),
            xlabel='Observed wavelength (µm)', ylabel=r'$F_\nu$ (mJy)')
ax_spec.set_xticks(np.arange(.4, optical_xlim[1], .1))
ax_spec.set_title('(c) Optical spectra', loc='left', pad=5)
ax_spec.set_xlabel('Observed wavelength (µm)', labelpad=2)
ax_spec.set_ylabel(r'$F_\nu$ (mJy)', labelpad=3)
legend_options = (dict(loc='lower right', bbox_to_anchor=(1, 1.04), ncol=3)
                  if prefix == 'P1823' else dict(loc='upper right'))
ax_spec.legend(frameon=False, fontsize=6.4, handlelength=1.4,
               borderaxespad=.2, labelspacing=.2, **legend_options)
for rest, label in OPTICAL_LINES:
    wave = rest * (1 + z) / 1e4
    ax_spec.axvline(wave, color='#b4bdc4', lw=.5, zorder=0)
    ax_spec.text(wave + .003, .96, label, transform=ax_spec.get_xaxis_transform(),
                 fontsize=6.4, rotation=90, va='top')
style(ax_spec)

fig.text(.065, .267, '(d) SPHEREx', fontsize=8.5)
ax_spherex = fig.add_axes([.09, .092, .882, .162])
draw_spherex_panel(ax_spherex, prefix, z, fontsize=6.0)

for suffix in ('pdf', 'png'):
    stem = 'fig1_single_example_preview' if prefix == 'A3' else f'fig1_{prefix}_example_preview'
    path = HERE / f'{stem}.{suffix}'
    fig.savefig(path, dpi=240)
    print(path)
plt.close(fig)

# Standalone, wide version of the unified observed-frame spectrum.
fig = plt.figure(figsize=(10, 4.5), facecolor='white')
ax = fig.add_axes([.075, .145, .9, .76])
draw_combined_spectrum(ax, compact=False)
ax.set_title(f'Optical–infrared spectra  |  z = {z:.3f}', loc='left', fontsize=12, pad=10)
for suffix in ('pdf', 'png'):
    path = HERE / f'optical_spherex_{prefix}_observed.{suffix}'
    fig.savefig(path, dpi=240)
    print(path)
plt.close(fig)

# Larger companion plot makes the line labels and supplied error bars readable.
fig = plt.figure(figsize=(8, 4.8), facecolor='white')
heading = 'SPHEREx spectrum'
fig.text(.055, .952, f'{heading}  |  z = {z:.6f}', fontsize=12, weight='bold')
ax = fig.add_axes([.085, .155, .88, .615])
draw_spherex_panel(ax, prefix, z, fontsize=9)
for suffix in ('pdf', 'png'):
    path = HERE / f'spherex_{prefix}_lines.{suffix}'
    fig.savefig(path, dpi=240)
    print(path)
plt.close(fig)
write_line_table(prefix, z)
write_data_summary(prefix)
