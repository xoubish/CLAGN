"""Redraw the published NGC 7469 spectrum from calibrated native figure pixels.

Digitized display traces are not the original spectral samples or new fits.
"""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
from pypdf import PdfReader
from scipy.signal import savgol_filter

BASE = Path(__file__).resolve().parents[1]
OUT = BASE/'review/fig2_three_panel/panel_a'
OUT.mkdir(parents=True, exist_ok=True)
BUILD = BASE/'review/fig2_three_panel/build'
BUILD.mkdir(parents=True, exist_ok=True)
SOURCE = BASE/'review/miri_literature/donnan_2402.17479.pdf'
DATA = BASE/'inputs/fig2_ngc7469_digitized.csv'
protected = [BASE/n for n in ['proposal.tex', 'fig2_miri.pdf']]
before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
page = PdfReader(SOURCE).pages[12]
original = next(v.image for v in page.images if v.image.size == (4026, 3355))
pixels = np.asarray(original.convert('RGBA'))

# Calibrate the logarithmic axes from all visible major tick centers.
wave_ticks = np.array([1.5, 3., 5., 10., 25.])
x_ticks = np.array([3307.5, 3462.5, 3576.5, 3731.5, 3936.5])
flux_ticks = np.array([1., .1, .01])
y_ticks = np.array([248.5, 424.5, 599.5])
xslope, xintercept = np.polyfit(np.log(wave_ticks), x_ticks, 1)
yslope, yintercept = np.polyfit(np.log10(flux_ticks), y_ticks, 1)
xres = float(np.max(abs(xslope*np.log(wave_ticks)+xintercept-x_ticks)))
yres = float(np.max(abs(yslope*np.log10(flux_ticks)+yintercept-y_ticks)))
assert xres < .6 and yres < .6
def wave(x): return np.exp((np.asarray(x)-xintercept)/xslope)
def flux(y): return 10**((np.asarray(y)-yintercept)/yslope)

# Native plot interior, excluding the original legend rectangle.
x0, x1, y0, y1 = 3303, 3962, 140, 626
rgb = pixels[y0:y1, x0:x1, :3].astype(float)
yy, xx = np.indices(rgb.shape[:2]); xx, yy = xx+x0, yy+y0
interior = (pixels[y0:y1, x0:x1, 3] > 200) & ~((xx < 3580) & (yy < 307))
colors = {'continuum': (218, 165, 32), 'warm': (31, 119, 180), 'hot': (148, 103, 189)}
traces = {}
for name, color in colors.items():
    mask = interior & (np.max(abs(rgb-np.asarray(color)), axis=2) < 30)
    columns = np.where(mask.any(axis=0))[0]
    rows = np.array([np.median(np.where(mask[:, x])[0])+y0 for x in columns])
    traces[name] = (columns+x0, rows)

# Column-median positions of the black observed-spectrum symbols. A generous
# envelope around the continuum rejects axes without removing line peaks.
cx, cy = traces['continuum']
cy_grid = np.interp(np.arange(x0, x1), cx, cy)
black = (rgb.max(axis=2) < 65) & interior
black &= (yy > cy_grid[None, :]-85) & (yy < cy_grid[None, :]+24)
black &= (xx >= cx.min()) & (xx <= cx.max())
columns = np.where(black.any(axis=0))[0]
rows = np.array([np.median(np.where(black[:, x])[0])+y0 for x in columns])
traces['observed'] = (columns+x0, rows)
counts = {name: len(x) for name, (x, y) in traces.items()}
assert counts['continuum'] > 600 and counts['observed'] > 500
assert counts['warm'] > 250 and counts['hot'] > 450

# Preserve actual extracted columns; NaNs retain missing pixels/dash gaps.
grid = np.arange(x0, x1)
table = pd.DataFrame({'source_x_pixel': grid, 'rest_um': wave(grid)})
for name, (x, y) in traces.items():
    values = np.full(len(grid), np.nan); values[x-x0] = flux(y)
    table[name+'_jy'] = values
table.to_csv(DATA, index=False)

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'pdf.fonttype': 42, 'axes.linewidth': .8})
fig = plt.figure(figsize=(6.5, 2.65))
ax = fig.add_axes([.085, .18, .90, .59])
ax.set_xscale('log'); ax.set_yscale('log')
ax.set(xlim=(1.35, 29.), ylim=(.007, 5.0),
       xlabel='Rest wavelength (µm)', ylabel=r'$F_\nu$ (Jy)')
ax.xaxis.set_major_locator(FixedLocator([1.5, 3, 5, 7, 10, 15, 20, 28]))
ax.xaxis.set_major_formatter(FixedFormatter(['1.5', '3', '5', '7', '10', '15', '20', '28']))
ax.xaxis.set_minor_formatter(NullFormatter())
ax.tick_params(which='both', direction='in', right=True, labelsize=9)
for name, color, label, style, width in [
        ('observed', '#252525', 'NIRSpec + MRS', 'None', 0),
        ('continuum', '#ba8618', 'Dust continuum', '-', 1.7),
        ('hot', '#9467bd', 'Hot dust', '--', 1.35),
        ('warm', '#1f77b4', 'Warm dust', '--', 1.35)]:
    x, y = traces[name]
    if name == 'observed':
        # Unconnected dots avoid artificial teeth from joining digitized
        # symbol-column medians. Thin only the display, preserving the CSV.
        ax.plot(wave(x[::2]), flux(y[::2]), color=color, marker='o', ms=1.9,
                markeredgewidth=0, ls='None', label=label, zorder=5)
        continue
    if name in ['hot', 'warm']:
        # Suppress pixel/dash-edge stair steps in the published smooth model
        # components. Preserve the raw extracted columns in the CSV.
        display_x = np.arange(x.min(), x.max()+1)
        y = savgol_filter(np.interp(display_x, x, y), 17, 2)
        x = display_x
    ax.plot(wave(x), flux(y), color=color, lw=width, ls=style, label=label,
            zorder=4 if name=='continuum' else 3)
ax.axvspan(8, 13, facecolor='#e8b446', alpha=.12, zorder=0)
for w in [9.7, 18.]:
    ax.axvline(w, color='#996516', lw=.8, ls=':', zorder=1)
    ax.text(w, .011, f'{w:g} µm', ha='center', fontsize=9, color='#79540c',
            bbox=dict(facecolor='white', alpha=.9, edgecolor='none', pad=1))
ax.text(6.0, .0155, 'Silicates:', fontsize=8.5, color='#79540c', ha='right')
for w in [6.2, 7.7, 11.3, 12.7, 17.]:
    ax.plot([w, w], [.88, .93], transform=ax.get_xaxis_transform(), color='#82549a', lw=.85)
ax.text(5.7, .905, 'PAH', transform=ax.get_xaxis_transform(), ha='right',
        va='center', color='#82549a', fontsize=8.5)

# Convert the observed-frame instrument limits to NGC 7469's rest frame.
# SPHEREx is comparison coverage, not the source of the plotted near-IR data.
z_ngc7469 = .0163
channels = [(4.9, 7.65), (7.51, 11.7), (11.55, 17.98), (17.7, 27.9)]
bar_colors = ['#247ba0', '#5264a5', '#ac6c32', '#a84352']
for (lo, hi), color in zip(channels, bar_colors):
    lo, hi = lo/(1+z_ngc7469), hi/(1+z_ngc7469)
    ax.plot([lo, hi], [1.08, 1.08], transform=ax.get_xaxis_transform(), color=color,
            lw=5, solid_capstyle='butt', clip_on=False)
lo, hi = 1.35, 5/(1+z_ngc7469)
ax.plot([lo, hi], [1.08, 1.08], transform=ax.get_xaxis_transform(), color='#34836b',
        lw=5, solid_capstyle='butt', clip_on=False)
ax.text(np.sqrt(lo*hi), 1.17, 'SPHEREx', transform=ax.get_xaxis_transform(),
        ha='center', fontsize=10, color='#28715d', weight='bold')
mrs_lo, mrs_hi = channels[0][0]/(1+z_ngc7469), channels[-1][1]/(1+z_ngc7469)
ax.text(np.sqrt(mrs_lo*mrs_hi), 1.17, 'JWST MIRI / MRS',
        transform=ax.get_xaxis_transform(), ha='center', fontsize=10,
        color='#283542', weight='bold')
fig.text(.535, .971, 'Separating hot and warm dust', ha='center', va='top',
         fontsize=11, weight='bold', color='#283542')
ax.legend(loc='upper left', fontsize=8, frameon=False, handlelength=2.0,
          markerscale=1.8,
          labelspacing=.22, borderaxespad=.4)
ax.text(.025, .49, 'NGC 7469, z = 0.0163', transform=ax.transAxes,
        ha='left', va='center', fontsize=9, color='#283542', weight='bold')
ax.text(.025, .40, 'Donnan et al. (2024)', transform=ax.transAxes,
        ha='left', va='center', fontsize=8, color='#283542')
fig.savefig(OUT/'literature_with_our_leverage.pdf')
fig.savefig(OUT/'literature_with_our_leverage.png', dpi=220)
plt.close(fig)

# Trace overlay on native pixels for visual verification of the digitization.
fig, ax = plt.subplots(figsize=(10, 7)); ax.imshow(original)
for name, (x, y) in traces.items():
    ax.scatter(x[::8], y[::8], s=7, marker='+', label=name)
ax.set(xlim=(3260, 3998), ylim=(640, 125)); ax.legend(loc='lower left')
fig.savefig(BUILD/'digitization_audit.png', dpi=160, bbox_inches='tight'); plt.close(fig)
assert before == {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
report = dict(source='https://arxiv.org/abs/2402.17479', figure='Figure 8, NGC 7469 nucleus',
    source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    method='Native pixels; logarithmic-axis calibration; colour traces and column-median black spectrum trace.',
    data=str(DATA.relative_to(BASE)), data_sha256=hashlib.sha256(DATA.read_bytes()).hexdigest(),
    x_tick_pixels=x_ticks.tolist(), wavelength_ticks_um=wave_ticks.tolist(),
    y_tick_pixels=y_ticks.tolist(), flux_ticks_jy=flux_ticks.tolist(),
    x_tick_max_residual_pixels=xres, y_tick_max_residual_pixels=yres,
    traced_columns=counts, original_data_samples=False,
    observed_display='Dark unconnected dots, every second extracted column; raw CSV unchanged.',
    channel_bars_observed_um=channels, channel_bars_source_redshift=z_ngc7469,
    channel_bars_rest_um=[[lo/(1+z_ngc7469),hi/(1+z_ngc7469)] for lo,hi in channels],
    redshift_source='https://www.aanda.org/articles/aa/pdf/2021/05/aa39291-20.pdf',
    observed_spectrum_instruments=['JWST NIRSpec IFU', 'JWST MIRI MRS'],
    spherex_bar='Comparison access at the same source redshift; not an input to the published spectrum.',
    channel_source='https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-observing-modes/miri-medium-resolution-spectroscopy',
    interpretation='Digitized published example, not target predictions. Hot/warm dash gaps interpolated and pixel steps smoothed with a 17-pixel quadratic Savitzky-Golay filter; no extrapolation. PAH ticks mark positions.',
    protected_sha256=before, active_proposal_unchanged=True)
(OUT/'checks.json').write_text(json.dumps(report, indent=2)+'\n')
(BASE/'inputs/fig2_ngc7469_digitization.json').write_text(json.dumps(report, indent=2)+'\n')
print(f'Digitized NGC 7469: {counts}; axis residuals {xres:.3f}/{yres:.3f} pixels.')
