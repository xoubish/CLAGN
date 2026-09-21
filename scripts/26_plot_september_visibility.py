"""Export September visibility plots from the validated public review snapshot.

Uses the stored topocentric Moon/airmass calculations and accepted intervals.
These are visibility windows, not scheduled exposure blocks.
"""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import json
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'plots/sep23_2026'
TZ = ZoneInfo('America/Los_Angeles')
INK = '#23384a'


def date_number(mjd):
    return mdates.date2num(Time(mjd, format='mjd').to_datetime())


def interval_number(value):
    return mdates.date2num(datetime.fromisoformat(value).replace(tzinfo=ZoneInfo('UTC')))


def format_time(value):
    return mdates.num2date(value, tz=TZ).strftime('%H:%M')


def setup_time_axis(ax, start, end):
    ax.set_xlim(start, end)
    ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0, 30], tz=TZ))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M', tz=TZ))
    ax.grid(axis='x', color='#dce2e8', linewidth=.7, zorder=0)
    midnight = mdates.date2num(datetime(2026, 9, 24, tzinfo=TZ))
    ax.axvline(midnight, color='#907aad', linewidth=1.2, alpha=.7, zorder=1)
    ax.tick_params(length=0, labelsize=10, colors=INK, pad=7)
    for spine in ax.spines.values():
        spine.set_color('#dce2e8')


def save(fig, stem):
    for suffix in ['png', 'pdf', 'svg']:
        fig.savefig(DEST / f'{stem}.{suffix}', dpi=180, facecolor='white')
    plt.close(fig)


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    html = (ROOT / 'docs/index.html').read_text()
    data = json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>', html, re.S).group(1))
    night = data['nights']['sep23']
    start, end = date_number(night['start_mjd']), date_number(night['end_mjd'])
    targets = []
    for target in data['targets']:
        window = next((w for w in target['nights'] if w['night'] == 'sep23'), None)
        if window and window['longest_minutes'] >= 120 and target['field_status'] == 'clear':
            targets.append(dict(target, window=window))
    targets.sort(key=lambda t: (t['pool_role'] == 'reserve', t['window']['ranges'][0]['start_utc'], t['ra']))
    main_targets = [t for t in targets if t['pool_role'] == 'manifold']
    reserve_targets = [t for t in targets if t['pool_role'] == 'reserve']
    assert main_targets and reserve_targets
    assert all(t['window']['moon_min'] >= 40 for t in targets)
    nmain, nreserve = len(main_targets), len(reserve_targets)
    rows = list(range(nmain)) + list(range(nmain + 1, nmain + 1 + nreserve))
    color_map = plt.get_cmap('Blues_r')
    norm = Normalize(1, 1.5)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'text.color': INK, 'axes.labelcolor': INK,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})

    fig = plt.figure(figsize=(14, 11.5))
    ax = fig.add_axes([.19, .20, .66, .65])
    export = []
    for target, row in zip(targets, rows):
        window = target['window']
        curve = np.asarray(window['curve'], float)
        x = date_number(curve[:, 0])
        x[0], x[-1] = start, end
        midpoints = (x[:-1] + x[1:]) / 2
        airmass = (curve[:-1, 1] + curve[1:, 1]) / 2
        moon = (curve[:-1, 2] + curve[1:, 2]) / 2
        accepted = np.zeros(len(midpoints), dtype=bool)
        for interval in window['ranges']:
            a, b = interval_number(interval['start_utc']), interval_number(interval['end_utc'])
            accepted |= (midpoints >= a) & (midpoints <= b)
            export.append(dict(name=target['name'], pool=target['pool_role'], historical_r=target['rmag'],
                               start_pdt=format_time(a), end_pdt=format_time(b),
                               start_utc=interval['start_utc'], end_utc=interval['end_utc'],
                               minutes=interval['minutes'], minimum_moon_deg=window['moon_min']))
        ax.barh(row, end-start, left=start, height=.76, color='#eef1f4', zorder=1)
        ax.pcolormesh(x, [row-.38, row+.38], np.ma.masked_where(~accepted, airmass)[None, :],
                      cmap=color_map, norm=norm, shading='flat', zorder=2, rasterized=False)
        # Explicitly show where the Moon, rather than elevation, rules out a slot.
        excluded = (~accepted) & (moon < 40)
        for left, right in zip(x[:-1][excluded], x[1:][excluded]):
            ax.barh(row, right-left, left=left, height=.76, color='#f7d9c8', linewidth=0, zorder=2)
        label = ' / '.join(f"{format_time(interval_number(r['start_utc']))}–{format_time(interval_number(r['end_utc']))}" for r in window['ranges'])
        ax.text(1.015, row, label, transform=ax.get_yaxis_transform(), va='center', fontsize=9, color='#526575')
    ax.set_ylim(rows[-1] + .7, -1.05)
    ax.set_yticks(rows)
    ax.set_yticklabels([f"{t['name']}   r {t['rmag']:.2f}" for t in targets], fontsize=10)
    setup_time_axis(ax, start, end)
    ax.text(-.01, -.82, f'{nmain} MAIN CANDIDATES · ≥2-HOUR WINDOWS', transform=ax.get_yaxis_transform(), fontsize=10, fontweight='bold')
    ax.text(-.01, nmain, f'{nreserve} BRIGHT QUASAR RESERVES', transform=ax.get_yaxis_transform(), fontsize=10, fontweight='bold', va='center')
    ax.text(1.015, -.82, 'Usable window · PDT', transform=ax.get_yaxis_transform(), fontsize=9, color='#526575')
    ax.set_xlabel('Palomar local time · PDT (UTC−7) · midnight begins September 24', labelpad=13)
    fig.text(.055, .955, 'September 23 · when each target is observable', fontsize=22, fontweight='bold')
    fig.text(.055, .915, f"Palomar / NGPS · first half-night, {format_time(start)}–{format_time(end)} PDT · Moon ≈{night['moon_percent']}% illuminated", fontsize=12)
    fig.text(.055, .882, 'Colored bars: airmass ≤1.5 and Moon separation ≥40° simultaneously. Darker blue = lower airmass.', fontsize=11)
    color_ax = fig.add_axes([.19, .105, .29, .018])
    colorbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=color_map), cax=color_ax, orientation='horizontal')
    colorbar.set_ticks([1, 1.1, 1.2, 1.3, 1.4, 1.5])
    colorbar.set_label('Airmass within the usable window', fontsize=10)
    colorbar.ax.tick_params(labelsize=9)
    fig.legend(handles=[Patch(facecolor='#eef1f4', label='Outside accepted window'),
                        Patch(facecolor='#f7d9c8', label='Moon closer than 40°')],
               loc='lower left', bbox_to_anchor=(.53, .085), frameon=False, fontsize=10)
    fig.text(.055, .028, 'Candidate visibility; exposure sequence pending. Historical r shown. Screened fields still need slit-PA review.\n'
             'Half-night ends at the midpoint of astronomical dusk and dawn; final handover time may differ.', fontsize=9, color='#627381')
    save(fig, 'visibility')

    fig, axes = plt.subplots(2, 1, figsize=(13.5, 9), sharex=True, gridspec_kw={'height_ratios': [1.3, 1]})
    fig.subplots_adjust(left=.09, right=.76, top=.86, bottom=.11, hspace=.24)
    for ax, group, title in zip(axes, [main_targets, reserve_targets], ['Main candidates', 'Bright quasar reserves']):
        ax.axhspan(1, 1.3, color='#e5f2ec', zorder=0)
        ax.axhspan(1.3, 1.5, color='#f0f6f2', zorder=0)
        for i, target in enumerate(group):
            curve = np.asarray(target['window']['curve'], float)
            x = date_number(curve[:, 0])
            color = plt.get_cmap('tab20')(i / 20)
            ax.plot(x, curve[:, 1], color=color, linewidth=1.1, alpha=.4)
            accepted = np.zeros(len(x), dtype=bool)
            for interval in target['window']['ranges']:
                a, b = interval_number(interval['start_utc']), interval_number(interval['end_utc'])
                accepted |= (x >= a-1/86400) & (x <= b+1/86400)
            ax.plot(x, np.where(accepted, curve[:, 1], np.nan), color=color, linewidth=2.4, label=target['name'])
        ax.axhline(1.5, linestyle='--', color='#86969e', linewidth=1)
        ax.set_ylim(2.7, .97)
        ax.set_yticks([1, 1.3, 1.5, 2, 2.5])
        ax.set_ylabel('Airmass · lower is better')
        ax.set_title(title, loc='left', fontsize=13, fontweight='bold', pad=9)
        ax.legend(loc='upper left', bbox_to_anchor=(1.025, 1.02), frameon=False, fontsize=9, ncol=2 if len(group)>8 else 1)
        setup_time_axis(ax, start, end)
    axes[-1].set_xlabel('Palomar local time · PDT (UTC−7) · midnight begins September 24', labelpad=13)
    fig.text(.09, .955, 'September 23 · visibility through the half-night', fontsize=21, fontweight='bold')
    fig.text(.09, .91, 'Bold curves satisfy the accepted Moon + airmass windows. Faint extensions show the rest of each track.', fontsize=11)
    fig.text(.09, .025, 'Airmass 1.0 is overhead. Green shading marks airmass ≤1.5; deeper green marks ≤1.3. Exposure sequence pending.', fontsize=10, color='#627381')
    save(fig, 'airmass')
    pd.DataFrame(export).to_csv(DEST / 'visibility_windows.csv', index=False)
    (DEST / 'README.md').write_text(
        '# September 23, 2026 visibility\n\n'
        f'Snapshot: {data["generated"]}. {nmain} screened main candidates and {nreserve} reserves, each with at least a two-hour accepted window.\n\n'
        'All times on the plots are Palomar local time, PDT (UTC−7); the night crosses into September 24. '
        'The first half-night runs from astronomical dusk to the midpoint of dusk and dawn. '
        'The plot uses the public review snapshot’s topocentric Moon separations and geometric sec(z) airmasses. '
        'Accepted windows are conservative five-minute intervals; the time labels are truncated to minutes. '
        'These are availability windows, not assigned exposures or an optimized observing sequence.\n\n'
        'PNG, PDF and SVG versions are provided for both the target-window plot and the airmass curves. '
        'The CSV preserves exact UTC boundaries. Regenerate with `26_plot_september_visibility.py`.\n')
    print(f'Wrote visibility and airmass PNG/PDF/SVG plots for {nmain} main candidates and {nreserve} reserves to {DEST}')


if __name__ == '__main__':
    main()
