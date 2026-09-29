"""Figure: what one MIRI/MRS epoch measures for each target's own history.

Left: one example target. Points are the measured total W1 history relative
to the 2010 AllWISE epoch; curves are the model 12-um nuclear flux under the
instant (no memory), delayed (fixed dust) and evolving hypotheses, all anchored
to the AllWISE 2010 W3 measurement. The MIRI 2028 error bar is the adopted
cross-epoch precision. Right: the same three predictions for all 24 targets.
Inputs come from forecast_warm_response.py (inputs/warm_response_forecast*.{csv,json}).
"""
from pathlib import Path
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
import numpy as np
import pandas as pd

import forecast_warm_response as F

HERE = Path(__file__).resolve().parent
INPUTS = HERE/'inputs'
BLUE, RED, INK, GOLD = '#24658a', '#bc5039', '#283542', '#bb8c35'
GREY, PURPLE, GREEN = '#7d878f', '#72518b', '#4c826e'
EXAMPLES = ('F03',)
SPECTRUM_WAVES = np.geomspace(4.5, 24., 70)


def style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_color('#9ca4a9')
    ax.tick_params(length=3, color='#9ca4a9', labelsize=10)


def spectrum_at_miri(tid):
    """Model rest-frame spectra at the MIRI epoch for one target, three hypotheses."""
    sample = pd.read_csv(INPUTS/'jwst_sample_cycle6.csv').set_index('id')
    binned = pd.read_csv(HERE/'review/history_audit/binned_measurements.csv')
    allwise = pd.read_csv(INPUTS/'sample_allwise_psd.csv')
    row = sample.loc[tid].copy(); row['id'] = tid
    aw = allwise.loc[allwise['id'].eq(tid)].iloc[0]
    host = F.host_fraction(aw)
    t, l, t0 = F.history(row, binned, allwise, host)
    tau = row['tau_dust_yr'] if np.isfinite(row['tau_dust_yr']) else F.koshida_lag_years(row['r_mag'], row['z'])
    grains = {g['name']: g for g in F.shell.GRAINS}
    rin = grains['Sil_21.gz']['inner_radius']/grains['Gra_21.gz']['inner_radius']*tau
    out = {}
    for h in ['instant', 'delayed', 'evolving']:
        out[h] = np.array([F.model_flux(np.array([F.MIRI_EPOCH]), t, l, w, rin, h, row['z'])[0] for w in SPECTRUM_WAVES])
    return out


def main(examples=EXAMPLES):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5,
                         'axes.labelsize': 10.5, 'axes.labelcolor': INK,
                         'text.color': INK, 'pdf.fonttype': 42})
    table = pd.read_csv(INPUTS/'warm_response_forecast.csv')
    metadata = json.loads((INPUTS/'warm_response_forecast_summary.json').read_text())
    F.DRIVER_GAMMA = metadata['driver_gamma']  # same heating conversion as the saved forecast
    series = json.loads((INPUTS/'warm_response_forecast_series.json').read_text())['series']
    miri = 2028.0
    fig = plt.figure(figsize=(7.6, 3.9), facecolor='white')
    axes_left = [fig.add_axes([.085, .665, .39, .27])]
    ax_sp = fig.add_axes([.085, .215, .39, .31])
    ax_r = fig.add_axes([.575, .125, .41, .72])

    for ax, tid in zip(axes_left, examples):
        s = series[tid]; row = table.set_index('id').loc[tid]
        years = np.array(s['years'])
        ax.axvspan(2010.8, miri-.15, color='#f2f3f5', zorder=0)
        # history_l is the inferred heating proxy. Plot actual total W1 photometry.
        aw = pd.read_csv(INPUTS/'sample_allwise_psd.csv').set_index('id').loc[tid]
        measured = pd.read_csv(HERE/'review/history_audit/binned_measurements.csv')
        measured = measured[measured.id.eq(tid) & measured.band.eq('W1')].sort_values('mjd')
        measured_years = np.r_[F.decimal_year(aw['w1mjdmean']), F.decimal_year(measured.mjd.to_numpy())]
        measured_ratio = np.r_[1., 10**(-.4*(measured.mag.to_numpy()-aw['w1mpro']))]
        ax.plot(measured_years, measured_ratio, 'o', ms=2.6, color=PURPLE, mec='white', mew=.3, zorder=4)
        ax.plot(years, s['instant'], color=GREY, lw=1.3, ls=':', zorder=2)
        ax.plot(years, s['delayed'], color=INK, lw=1.9, zorder=3)
        ax.plot(years, s['evolving'], color=RED, lw=1.5, ls='--', zorder=3)
        ax.plot([s['history_t'][0]], [1.], 's', ms=5, color=GOLD, mec='white', mew=.4, zorder=5)
        y = s['delayed'][-1]; sig = float(row['sigma_dex'])
        ax.errorbar([miri], [y], yerr=[[y*(1-10**-sig)], [y*(10**sig-1)]], fmt='D', ms=4,
                    color=INK, mfc='white', mew=1, capsize=2.5, lw=1, zorder=6)
        ax.set(xlim=(2009.6, 2030.4), yscale='log')
        lo, hi = min(min(measured_ratio), min(s['delayed']), min(s['evolving']), min(s['instant'])), \
                 max(max(measured_ratio), max(s['delayed']), max(s['evolving']), max(s['instant']))
        ax.set_ylim(lo/1.12, hi*1.28)
        # Bracket the spread of the three predictions at the MIRI epoch.
        ends = [s['instant'][-1], s['delayed'][-1], s['evolving'][-1]]
        spread = np.log10(max(ends)/min(ends))
        ax.plot([miri+.9, miri+.9], [min(ends), max(ends)], color='#5d6870', lw=.9)
        ax.plot([miri+.7, miri+.9], [min(ends)]*2, color='#5d6870', lw=.9)
        ax.plot([miri+.7, miri+.9], [max(ends)]*2, color='#5d6870', lw=.9)
        ax.text(miri+1.15, np.sqrt(min(ends)*max(ends)), f'{spread:.2f}\ndex', fontsize=7.5, color='#5d6870',
                ha='left', va='center', linespacing=1.0)
        ticks = [t for t in [.5, .7, 1, 1.5, 2, 3, 4] if lo/1.12 < t < hi*1.28]
        ax.yaxis.set_major_locator(FixedLocator(ticks))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis='y', which='minor', length=0)
        style(ax)
        label = f"{tid}   z = {row['z']:.2f}   {'fading' if row['family'] == 'fade' else 'rising'} selection"
        ax.text(.02, .95, label, transform=ax.transAxes, fontsize=8.5, va='top', fontweight='bold')
    axes_left[0].set_xlabel('Year', labelpad=1)
    axes_left[0].set_ylabel('Flux / 2010', labelpad=2)
    axes_left[0].text(2019.4, axes_left[0].get_ylim()[0]*1.08, '2010 reference → MIRI epoch',
                      ha='center', va='bottom', fontsize=8, color='#5d6870')
    handles = [Line2D([], [], marker='o', color=PURPLE, ls='', ms=3.5, label='Total W1 (measured)'),
               Line2D([], [], marker='s', color=GOLD, ls='', ms=5, label='2010 W3 reference'),
               Line2D([], [], color=GREY, ls=':', lw=1.3, label='12 µm, no memory'),
               Line2D([], [], color=INK, lw=1.9, label='12 µm, delayed'),
               Line2D([], [], color=RED, ls='--', lw=1.5, label='12 µm, evolving dust'),
               Line2D([], [], marker='D', color=INK, mfc='white', ls='', ms=4, label='Forecast 2028 ±1σ')]
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(.02, .005), fontsize=7.4, frameon=False,
               ncol=3, handlelength=1.6, handletextpad=.5, columnspacing=1.0, borderaxespad=.2)

    # Spectrum panel: the same target's MRS spectrum at the MIRI epoch, relative to no memory.
    tid = examples[0]
    spec = spectrum_at_miri(tid)
    z = float(table.set_index('id').loc[tid, 'z'])
    ax_sp.axvspan(8, 13, color='#eef0f2', zorder=0)
    ax_sp.axhline(0, color=GREY, ls=':', lw=1.2)
    # Scale every hypothesis to the same 5-um flux: this is what the spectrum alone constrains.
    i5 = np.argmin(abs(SPECTRUM_WAVES-5.))
    norm = {h: v/v[i5] for h, v in spec.items()}
    ax_sp.plot(SPECTRUM_WAVES, 100*(norm['delayed']/norm['instant']-1), color=INK, lw=1.9)
    ax_sp.plot(SPECTRUM_WAVES, 100*(norm['evolving']/norm['instant']-1), color=RED, lw=1.5, ls='--')
    mrs = (4.9/(1+z), 27.9/(1+z))
    ax_sp.set(xscale='log', xlim=(max(4.5, mrs[0]), min(24., mrs[1])), xlabel='Rest wavelength (µm)',
              ylabel='% vs no memory\n(scaled at 5 µm)')
    ax_sp.xaxis.set_major_locator(FixedLocator([5, 7, 10, 15, 20]))
    ax_sp.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
    ax_sp.xaxis.set_minor_formatter(NullFormatter())
    peak = max(100*(norm['delayed']/norm['instant']-1).max(), 100*(norm['evolving']/norm['instant']-1).max())
    top = peak*1.9+8
    ax_sp.set_ylim(-10, top)
    # Diagnostics the same spectrum delivers, marked in the free band above the curves.
    row1, row2 = peak*1.25+4, peak*1.6+6
    for w in [9.7, 18.]:
        ax_sp.axvline(w, color=GOLD, lw=.8, ls=':', ymax=(row1+2-ax_sp.get_ylim()[0])/(top-ax_sp.get_ylim()[0]))
    ax_sp.text(13.2, row1, 'silicates: geometry', color=GOLD, fontsize=7.2, ha='center', va='center')
    for w in [6.2, 7.7, 11.3]:
        ax_sp.plot([w, w], [row1-3, row1+3], color=GREEN, lw=1.2)
    ax_sp.text(6.9, row1-5, 'PAH: host', color=GREEN, fontsize=7.2, ha='center', va='top')
    for w in [7.65, 14.32, 15.56]:
        ax_sp.plot([w, w], [row2-3, row2+3], color=PURPLE, lw=1.2)
    ax_sp.text(15.0, row2+5, '[Ne V], [Ne III]: ionising power', color=PURPLE, fontsize=7.2, ha='center', va='bottom')
    ax_sp.text(7.65, row2+5, '[Ne VI]', color=PURPLE, fontsize=7.2, ha='center', va='bottom')
    ax_sp.text(10.5, -8, '8–13 µm', fontsize=8, color='#5d6870', ha='center', va='bottom')
    ax_sp.set_title(f'{tid}: MRS spectrum at the MIRI epoch', loc='left', fontsize=8.5, fontweight='bold', pad=3)
    style(ax_sp)

    # Right: the three predictions per target in one space, sorted by their spread.
    summary = json.loads((INPUTS/'warm_response_forecast_summary.json').read_text())['summary']
    t = table.copy()
    cols = ['log_ratio_instant', 'log_ratio_delayed', 'log_ratio_evolving']
    t['spread'] = t[cols].max(axis=1)-t[cols].min(axis=1)
    t['sigma_meas'] = np.sqrt(t['sigma_dex']**2-t['sigma_host_dex']**2)
    t = t.sort_values(t['spread']/t['sigma_meas'], ascending=True).reset_index(drop=True) \
        if False else t.assign(signif=t['spread']/t['sigma_meas']).sort_values('signif', ascending=True).reset_index(drop=True)
    y = np.arange(len(t))
    ax_r.axvline(0, color='#c9cfd4', lw=.8)
    # Two-tone error bars: light = with the current host-fraction uncertainty, dark = measurement terms only.
    ax_r.errorbar(t['log_ratio_delayed'], y, xerr=t['sigma_dex'], fmt='none', ecolor='#e3e6e9', elinewidth=2.6,
                  capsize=0, zorder=1)
    ax_r.errorbar(t['log_ratio_delayed'], y, xerr=t['sigma_meas'], fmt='none', ecolor='#b9c0c6', elinewidth=2.6,
                  capsize=0, zorder=1.5)
    n3m = int((t['spread'] > 3*t['sigma_meas']).sum()); n3t = int((t['spread'] > 3*t['sigma_dex']).sum())
    for i, r in t.iterrows():
        ax_r.plot([t.loc[i, cols].min(), t.loc[i, cols].max()], [i, i], color='#8f979e', lw=.7, zorder=2)
    ax_r.plot(t['log_ratio_instant'], y, 'o', ms=4.2, color='white', mec=GREY, mew=1.1, zorder=3)
    ax_r.plot(t['log_ratio_delayed'], y, 'o', ms=4.2, color=INK, zorder=4)
    ax_r.plot(t['log_ratio_evolving'], y, 'x', ms=5, color=RED, mew=1.2, zorder=5)
    base = summary['baseline']; err = base['mean_error_dex']
    ax_r.set_yticks(y)
    ax_r.set_yticklabels([f"{r['id']}" for _, r in t.iterrows()], fontsize=7)
    for lab, (_, r) in zip(ax_r.get_yticklabels(), t.iterrows()):
        lab.set_color(RED if r['family'] == 'fade' else BLUE)
    ax_r.set(xlabel='log$_{10}$ (12 µm flux / 2010)', ylim=(-.8, len(t)-.2))
    style(ax_r)
    ax_r.tick_params(axis='y', length=0)
    ax_r.set_title('24 histories: W3-proxy forecasts', fontsize=9.5,
                   loc='left', pad=7, fontweight='bold')
    handles = [Line2D([], [], marker='o', color='white', mec=GREY, mew=1.1, ls='', ms=4.2, label='no memory'),
               Line2D([], [], marker='o', color=INK, ls='', ms=4.2, label='delayed'),
               Line2D([], [], marker='x', color=RED, ls='', ms=5, mew=1.2, label='evolving dust'),
               Line2D([], [], color='#b9c0c6', lw=2.6, label='±1σ measurement'),
               Line2D([], [], color='#e3e6e9', lw=2.6, label='+ host-fraction term')]
    ax_r.set_xlim(-.42, .74)
    ax_r.legend(handles=handles, loc='lower right', fontsize=7.4, frameon=False, handletextpad=.4, borderaxespad=.1)
    n3 = base['n_memory_over_3sigma']; n3e = base['n_evolution_over_3sigma']
    fig.savefig(HERE/'fig2_forecast.pdf', dpi=220)
    plt.close(fig)
    print(f'fig2_forecast.pdf: examples {examples}; >3 sigma memory {n3}/24, evolving {n3e}/24')


if __name__ == '__main__':
    main()
