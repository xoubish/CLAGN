"""Plots shared by the executed exploration notebooks and source atlas."""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time
from analysis import ROOT, OUT, LINES, line_index

COLORS = {'sep23': '#d97706', 'oct26': '#d63d4f', 'oct27': '#8d41c7',
          'SDSS': '#2573a6', 'DESI': '#269280'}


def load_tables():
    names = ['objects', 'visits', 'epochs', 'archive_audit', 'line_indices', 'comparisons',
             'photometry', 'photometry_summary', 'photometry_matches', 'jwst_coverage']
    return {name: pd.read_csv(OUT/f'{name}.csv', low_memory=False) for name in names}


def year(mjd):
    return Time(np.asarray(mjd, float), format='mjd').decimalyear


def configure():
    plt.rcParams.update({'figure.dpi': 110, 'font.size': 10, 'axes.spines.top': False,
        'axes.spines.right': False, 'axes.grid': True, 'grid.alpha': .16,
        'savefig.bbox': 'tight', 'axes.titlesize': 11})


def dossier(name, spectra, tables):
    obj = tables['objects'].set_index('name').loc[name]
    ss = sorted([s for s in spectra if s['name'] == name], key=lambda s: s['mjd'])
    fig, axs = plt.subplots(2, 3, figsize=(17, 9), layout='constrained')
    axfull, axhb, axother, axopt, axir, axhist = axs.flat
    other = 'Ha' if 6800*(1+obj.z) < 10200 else 'MgII'
    for i, s in enumerate(ss):
        ngps = s['source'] == 'NGPS'
        old = [a for a in ss if a['source'] != 'NGPS']
        highlighted = ngps or (len(old) and (s is old[0] or s is old[-1]))
        color = COLORS.get(s['night'] if ngps else s['source'], '#444444') if highlighted else '#b5bbc3'
        label = f"{Time(s['mjd'], format='mjd').isot[:10]} {s['source']}"
        if s['max_mjd']-s['min_mjd'] > 1:
            label += f" ({s['max_mjd']-s['min_mjd']:.0f}d coadd)"
        rest = s['wave']/(1+s['z'])
        for ax, line in [(axfull, 'full'), (axhb, 'Hb'), (axother, other)]:
            metric = line_index(s, 'Hb' if line == 'full' else line)
            norm = metric.get('continuum', np.nan)
            if not np.isfinite(norm) or norm <= 0:
                valid = np.isfinite(s['flux']) & (rest > 4500) & (rest < 5500)
                norm = np.nanmedian(s['flux'][valid]) if valid.any() else np.nan
            if not np.isfinite(norm) or norm <= 0:
                continue
            y = s['flux']/norm
            bounds = (max(2200, np.nanmin(rest)), min(7500, np.nanmax(rest))) if line == 'full' else LINES[line]['panel']
            sel = (rest > bounds[0]) & (rest < bounds[1])
            ax.plot(rest[sel], y[sel], color=color, lw=.65 if highlighted else .5,
                    alpha=.9 if highlighted else .35, label=label if highlighted and line == 'full' else None)
            if ngps and line != 'full':
                ax.fill_between(rest[sel], y[sel]-s['error'][sel]/norm, y[sel]+s['error'][sel]/norm,
                                color=color, alpha=.08)
            ax.set_xlabel('Rest wavelength (Å)')
            ax.set_ylabel('Flux / local pseudo-continuum')
            if line != 'full':
                ax.set_xlim(*bounds)
                for a, b in LINES[line]['wings']:
                    ax.axvspan(a, b, color='#169b8e', alpha=.035)
                for a, b in LINES[line]['sides']:
                    ax.axvspan(a, b, color='#777777', alpha=.025)
    for ax in [axfull, axhb, axother]:
        vals = np.concatenate([np.asarray(l.get_ydata()) for l in ax.lines]) if ax.lines else np.array([])
        vals = vals[np.isfinite(vals)]
        if len(vals):
            ax.set_ylim(min(-.1, np.percentile(vals, 1)), max(1.5, min(15, np.percentile(vals, 99.3)*1.15)))
    axfull.legend(fontsize=6, loc='upper right')
    axfull.set_title('All retained epochs; intermediate archives in grey')
    axhb.set_title('Hβ + [O III]; shaded wings and continuum windows')
    axother.set_title(other+' profile; no resolution matching')
    photo = tables['photometry'].query('name == @name')
    for (band, stream), g in photo.groupby(['band', 'stream']):
        if band in ['g', 'r']:
            color = '#239267' if band == 'g' else '#d84949'
            alert = 'alerts' in stream
            # Display medians to retain long baselines without overplotting raw cadence.
            b = g.groupby(np.floor(g.mjd/30)).agg(mjd=('mjd', 'median'), value=('value', 'median'))
            axopt.plot(year(b.mjd), b.value, 'o' if not alert else 'v', ms=3.3,
                color=color, mfc='none' if alert else color, label=band+' '+stream, alpha=.8)
        if band in ['W1', 'W2']:
            axir.errorbar(year(g.mjd), g.value, yerr=g.error, fmt='o' if 'NEOWISE' in stream else '.',
                ms=3, lw=.5, alpha=.7, label=band+' '+stream)
    for s in ss:
        if s['source'] == 'NGPS':
            for ax in [axopt, axir]:
                ax.axvline(year([s['mjd']])[0], color=COLORS[s['night']], ls='--', alpha=.6, lw=.8)
    axopt.set_title('Optical history · 30-day display bins; alerts open')
    axir.set_title('Infrared history · separate reduction streams')
    for ax in [axopt, axir]:
        ax.invert_yaxis(); ax.set_xlabel('Calendar year')
        if ax.lines:
            ax.legend(fontsize=6, loc='best')
    axopt.set_ylabel('AB magnitude (total source)')
    axir.set_ylabel('Vega magnitude (total source)')
    metrics = tables['line_indices']
    h = metrics[(metrics.name == name) & (metrics.line == 'Hb') & (metrics.variant == 'nominal') & metrics.valid]
    for source, g in h.groupby('source'):
        color = '#b84c41' if source == 'NGPS' else COLORS.get(source, 'grey')
        axhist.scatter(year(g.mjd), g.wing_ew, label=source, color=color, s=22)
        er = g.wing_ew_err_diagonal
        mask = np.isfinite(er)
        axhist.errorbar(year(g.loc[mask, 'mjd']), g.loc[mask, 'wing_ew'], yerr=er[mask], fmt='none', color=color, lw=.7)
    axhist.set_title('Hβ wing pseudo-EW · statistical errors where available')
    axhist.set_ylabel('Rest-frame signed wing index (Å)'); axhist.set_xlabel('Calendar year')
    if len(h):
        axhist.legend(fontsize=7)
    else:
        axhist.text(.5, .5, 'No accepted Hβ wing measurement', transform=axhist.transAxes, ha='center')
    fig.suptitle(f"{name}  |  z={obj.z:.4f}  |  {obj.n_archive} archive + {obj.n_NGPS} NGPS epochs\n"
                 'Exploratory: host, Fe II, aperture, telluric and flux-calibration effects remain', fontsize=14)
    return fig


def export_atlas(spectra, tables):
    from matplotlib.backends.backend_pdf import PdfPages
    import html
    directory = OUT/'dossiers'
    directory.mkdir(parents=True, exist_ok=True)
    rows = []
    with PdfPages(OUT/'all_sources_atlas.pdf') as pdf:
        for i, name in enumerate(tables['objects'].name, 1):
            fig = dossier(name, spectra, tables)
            filename = f'{name}.png'
            fig.savefig(directory/filename, dpi=110)
            pdf.savefig(fig, dpi=100)
            plt.close(fig)
            rows.append(dict(name=name, figure=f'dossiers/{filename}'))
            if i % 20 == 0:
                print(f'Atlas: {i}/{len(tables["objects"])}', flush=True)
    pd.DataFrame(rows).to_csv(OUT/'dossier_index.csv', index=False)
    # This is a local report index, with no scripts, server, or publication step.
    links = '\n'.join(f'<li><a href="{html.escape(r["figure"], quote=True)}">{html.escape(r["name"])}</a></li>' for r in rows)
    (OUT/'source_index.html').write_text('<!doctype html><meta charset="utf-8"><title>Palomar source atlas</title>'
        '<h1>Palomar 2026: all 99 sources</h1><p>Local exploratory products; may contain collaboration spectra.</p>'
        '<p><a href="all_sources_atlas.pdf">Full PDF atlas</a></p><ul>'+links+'</ul>')
    return pd.DataFrame(rows)
