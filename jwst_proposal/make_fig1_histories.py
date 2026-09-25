"""Figure 1 (v0.3): dated optical spectra, ZTF and W1 histories for the three anchors,
with the Cycle 6 window. Data from the bundled observer-page snapshot; no rescaling.
NGPS gets display-only 6-A FWHM smoothing per run/arm (as in v0.2)."""
from pathlib import Path
import json
import numpy as np
from astropy.time import Time
from scipy.ndimage import gaussian_filter1d
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator

HERE = Path(__file__).resolve().parent
SRC = HERE / 'inputs/figure_targets.json'
data = json.loads(SRC.read_text(encoding='utf-8'))
T = {t['name']: t for t in data['targets']}
RAMP = LinearSegmentedColormap.from_list('blue', ['#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])
ORANGE, AQUA, MAGENTA, VIOLET, INK, INK2, GRID, AXIS, MUTED = '#eb6834', '#1baf7a', '#d55181', '#4a3aa7', '#0b0b0b', '#52514e', '#e1e0d9', '#c3c2b7', '#898781'
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Times New Roman', 'Times'], 'font.size': 8.5,
                     'axes.labelsize': 8.5, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'pdf.fonttype': 42,
                     'axes.edgecolor': AXIS, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .6})
yr = lambda m: Time(np.asarray(m, float), format='mjd').decimalyear
def arrays(e):
    w = np.asarray(e['wave'], float); f = np.asarray(e['flux'], float)
    if e.get('instrument') == 'NGPS':
        good = np.flatnonzero(np.isfinite(f)); breaks = set(e.get('arm_breaks', []))
        cuts = np.flatnonzero((np.diff(good) > 1) | np.isin(good[1:], list(breaks))) + 1
        for run in np.split(good, cuts):
            if len(run) > 2: f[run] = gaussian_filter1d(f[run], 6/(2.35482*np.median(np.diff(w[run]))), mode='nearest')
        for b in breaks:
            if 0 <= b < len(f): f[b] = np.nan
    return w, f

NAMES = [('P11530', 'A1  faded 2019–21'), ('P10381', 'A2  faded 2014–19, recovering'), ('P9694', 'A3  brightening')]
fig = plt.figure(figsize=(6.5, 3.4), facecolor="white")
gs = fig.add_gridspec(3, 3, height_ratios=[1.35, 1, 1], hspace=.62, wspace=.42, left=.085, right=.985, top=.89, bottom=.115)
for c, (name, role) in enumerate(NAMES):
    t = T[name]; z = t['z']
    eps = sorted(t['spec']['epochs'], key=lambda e: e['mjd'])
    arch = [e for e in eps if e.get('instrument') != 'NGPS' and not e.get('quality_flag')]
    ngps = [e for e in eps if e.get('instrument') == 'NGPS']
    mj = np.array([e['mjd'] for e in arch]); span = max(1., mj.max()-mj.min())
    ax = fig.add_subplot(gs[0, c]); ymax = 0
    for e in arch:
        w, f = arrays(e); m = (w > 3800) & (w < 10000); ymax = max(ymax, np.nanpercentile(f[m], 99.5))
        ax.plot(w, f, color=RAMP((e['mjd']-mj.min())/span), lw=.5, alpha=.9)
    for e in ngps:
        w, f = arrays(e); m = (w > 3800) & (w < 10000); ymax = max(ymax, np.nanpercentile(f[m], 99.5))
        ax.plot(w, f, color=ORANGE, lw=.8, zorder=5)
    ax.set_xlim(3500, 10300); ax.set_ylim(0, 1.1*ymax); ax.set_xticks([4000, 7000, 10000])
    ax.yaxis.set_major_locator(MaxNLocator(3)); ax.grid(axis='y', color=GRID, lw=.5); ax.set_axisbelow(True)
    for ln in t.get('lines', []):
        if ln['name'] in ('Hβ', 'Hα'):
            ax.axvline(ln['angstrom'], color=AXIS, lw=.5, zorder=0)
            ax.text(ln['angstrom'], .97, ln['name'], transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=7.5, color=INK2,
                    bbox=dict(facecolor='white', edgecolor='none', pad=.3, alpha=.9))
    ax.set_title(f"{name.replace('P','')}: {t['jname'][:10]}  z={z:.3f}\n{role}", fontsize=8.5, color=INK, loc='left', pad=3)
    if c == 0: ax.set_ylabel(r'$F_\lambda$ (10$^{-17}$ cgs/Å)')
    ax.set_xlabel('Observed wavelength (Å)', labelpad=1)
    ax.text(.98, .78, f'{len(arch)} archival', color=INK2, transform=ax.transAxes, ha='right', fontsize=7.5)
    ax.text(.98, .62, '2026 NGPS', color=ORANGE, transform=ax.transAxes, ha='right', fontsize=7.5)

    ax = fig.add_subplot(gs[1, c])
    for band, col in [('g', AQUA), ('r', MAGENTA)]:
        v = np.asarray(t['ztf'].get(band, []), float)
        if len(v): ax.plot(yr(v[:, 0]), v[:, 1], '.', ms=2, color=col, alpha=.8)
    ax.invert_yaxis(); ax.set_xlim(2009.5, 2028.7); ax.set_xticks([2010, 2016, 2022, 2028]); ax.yaxis.set_major_locator(MaxNLocator(3))
    for e in arch: ax.plot([yr(e['mjd'])]*2, [.9, 1.0], transform=ax.get_xaxis_transform(), lw=.7, color='#2a78d6')
    ax.axvline(yr(ngps[0]['mjd']), color=ORANGE, lw=.9, ls='--')
    ax.axvspan(2027.5, 2028.5, color=GRID, alpha=.6, lw=0)
    if c == 0: ax.set_ylabel('ZTF mag')
    ax.text(.03, .1, 'g', color=AQUA, transform=ax.transAxes, weight='bold', fontsize=8); ax.text(.1, .1, 'r', color=MAGENTA, transform=ax.transAxes, weight='bold', fontsize=8)

    ax = fig.add_subplot(gs[2, c])
    wise = np.asarray(t['wise'].get('W1', []), float); neo = np.asarray(t.get('neo', []), float)
    if len(wise): ax.errorbar(yr(wise[:, 0]), wise[:, 1], yerr=wise[:, 2], fmt='o', ms=2.6, lw=.5, color=VIOLET, mfc=VIOLET)
    if len(neo): ax.errorbar(yr(neo[:, 0]), neo[:, 1], yerr=neo[:, 2], fmt='o', ms=2.6, lw=.5, color=VIOLET, mfc='white', mew=.7)
    ax.axvline(yr(ngps[0]['mjd']), color=ORANGE, lw=.9, ls='--'); ax.axvspan(2027.5, 2028.5, color=GRID, alpha=.6, lw=0)
    ax.set_xlim(2009.5, 2028.7); ax.set_xticks([2010, 2016, 2022, 2028]); ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.set_ylim(0, None); ax.grid(axis='y', color=GRID, lw=.5); ax.set_axisbelow(True)
    if c == 0: ax.set_ylabel('W1 (mJy)')
    ax.set_xlabel('Year', labelpad=1)
    if c == 2: ax.text(2028.0, .06, 'Cycle 6', rotation=90, transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=7, color=INK2)
    ax.text(2024.7, .06, 'NEOWISE ends', transform=ax.get_xaxis_transform(), ha='right', va='bottom', fontsize=6.5, color=MUTED)
fig.savefig(HERE/'fig1_histories.pdf'); fig.savefig(HERE/'fig1_histories.png', dpi=220)
print('ok')
