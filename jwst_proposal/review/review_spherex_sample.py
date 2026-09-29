"""Build an exploratory, reproducible review of all 24 SPHEREx targets.

Reads a snapshot of individually saved results, including incomplete targets.
Does not replace proposal inputs or interpret continuum colors as temperatures.
"""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
import collections
import hashlib
import html
import json
import sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.time import Time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RAW = ROOT / 'local_data/spherex/sample'
OUT = HERE / 'spherex_sample_review'
COLORS = ['#247ba0', '#d36936', '#8b5a9b', '#278568', '#8b7735']
COLS = ['wavelength_um', 'flux_mjy', 'flux_err_mjy', 'mjds']
WINDOWS = {'short': (2.1, 2.4), 'long': (3.6, 3.9)}


def year(mjd):
    return Time(mjd, format='mjd').decimalyear


def read_target(tid):
    original = (ROOT / 'inputs/spaxel_scryer_P9694_ra11p5333_dec9p1225_cleaned.csv'
                if tid == 'R01' else ROOT / f'inputs/spherex_extracted/{tid}.csv')
    counts = collections.Counter()
    digest = hashlib.sha256()
    if tid in ('F04', 'F06', 'R01'):
        data = original.read_bytes(); digest.update(data)
        df = pd.read_csv(original)
        if 'spectrum_type' in df:
            df = df[df.spectrum_type.eq('cleaned') & df.spectrum_role.eq('source')]
        complete, requested, source = True, len(df), str(original.relative_to(ROOT))
        counts['accepted'] = len(df)
    else:
        records = []
        for path in sorted((RAW / tid).glob('*.json')):
            data = path.read_bytes(); digest.update(data)
            try:
                rec = json.loads(data)
                counts[rec['status']] += 1
                if rec['status'] == 'accepted': records.append(rec['measurement'])
            except (ValueError, KeyError):
                counts['unreadable'] += 1
        df = pd.DataFrame(records) if records else pd.DataFrame(columns=COLS)
        index = pd.read_csv(RAW / f'{tid}_index.csv')
        if 'dataproduct_subtype' in index:
            index = index[index.dataproduct_subtype.eq('science')]
        requested = index.access_url.nunique()
        complete = counts['accepted'] + counts['quality_rejected'] == requested
        source = f'local_data/spherex/sample/{tid}/*.json'
    valid = (np.isfinite(df[COLS].astype(float)).all(axis=1) & (df.flux_err_mjy > 0))
    counts['invalid_accepted'] = int((~valid).sum())
    df = df[valid].sort_values('mjds').copy()
    if len(df):
        df['visit'] = np.r_[0, np.cumsum(np.diff(df.mjds) > 45)]
    else:
        df['visit'] = pd.Series(dtype=int)
    return df, dict(id=tid, source=source, snapshot_sha256=digest.hexdigest(),
                    download_complete=bool(complete), requested=int(requested), **counts)


def window_mean(group, lo, hi):
    """Equal-weight three fixed wavelength cells; no interpolation over gaps.

    Inverse-variance mean within each cell; propagated supplied errors only.
    This is a descriptive continuum proxy, not synthetic calibrated photometry.
    """
    means, variances = [], []
    for a, b in zip(np.linspace(lo, hi, 4)[:-1], np.linspace(lo, hi, 4)[1:]):
        g = group[(group.rest_um >= a) & (group.rest_um < b)]
        if not len(g): return np.nan, np.nan
        w = 1 / g.flux_err_mjy.to_numpy() ** 2
        means.append(float(np.sum(w * g.flux_mjy) / w.sum()))
        variances.append(1 / w.sum())
    return float(np.mean(means)), float(np.sqrt(np.sum(variances)) / 3)


def savefig(fig, name):
    fig.savefig(OUT / f'{name}.png', dpi=150, bbox_inches='tight')
    fig.savefig(OUT / f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    sample = pd.read_csv(ROOT / 'inputs/jwst_sample_cycle6.csv')
    history = pd.read_csv(HERE / 'data_anchored_scenarios/sample_measurements.csv').set_index('id')
    datasets, inventory, visits = {}, [], []
    for target in sample.itertuples():
        df, info = read_target(target.id)
        df['rest_um'] = df.wavelength_um / (1 + target.z)
        datasets[target.id] = df
        info.update(z=target.z, valid_points=len(df), visits=int(df.visit.nunique()))
        inventory.append(info)
        for visit, group in df.groupby('visit'):
            rec = dict(id=target.id, visit=int(visit), n=len(group),
                       first_mjd=float(group.mjds.min()), last_mjd=float(group.mjds.max()),
                       mean_year=float(year(group.mjds.median())),
                       rest_min_um=group.rest_um.min(), rest_max_um=group.rest_um.max(),
                       download_complete=info['download_complete'])
            for name, bounds in WINDOWS.items():
                rec[f'{name}_mjy'], rec[f'{name}_err_mjy'] = window_mean(group, *bounds)
            a, b, ea, eb = (rec[k] for k in ['short_mjy', 'long_mjy', 'short_err_mjy', 'long_err_mjy'])
            rec['color_eligible'] = bool(np.isfinite([a,b,ea,eb]).all() and a > 3*ea and b > 3*eb)
            rec['color'] = b/a if rec['color_eligible'] else np.nan
            rec['color_err'] = rec['color'] * np.hypot(ea/a, eb/b) if rec['color_eligible'] else np.nan
            visits.append(rec)
        df.to_csv(OUT / f'{target.id}_snapshot.csv', index=False)
    inv = pd.DataFrame(inventory).fillna(0)
    epochs = pd.DataFrame(visits)
    inv.to_csv(OUT / 'inventory.csv', index=False)
    epochs.to_csv(OUT / 'visit_measurements.csv', index=False)
    # Atlas: all finite fluxes, including negative values, on a linear scale.
    for page in range(4):
        fig, axes = plt.subplots(3, 2, figsize=(12, 10), constrained_layout=True)
        for ax, target in zip(axes.flat, list(sample.itertuples())[page*6:(page+1)*6]):
            df = datasets[target.id]; row = inv[inv.id.eq(target.id)].iloc[0]
            for visit, g in df.groupby('visit'):
                label = f'{year(g.mjds.min()):.2f}–{year(g.mjds.max()):.2f} (n={len(g)})'
                ax.errorbar(g.rest_um, g.flux_mjy, yerr=g.flux_err_mjy, fmt='.', ms=3,
                            elinewidth=.5, alpha=.65, color=COLORS[int(visit)%len(COLORS)], label=label)
            for lo, hi in WINDOWS.values(): ax.axvspan(lo, hi, color='grey', alpha=.09)
            ax.axhline(0, color='grey', lw=.5)
            ax.set(xlabel='Rest wavelength (µm)', ylabel='Flux density (mJy)', xlim=(.5, 5.1),
                   title=f'{target.id} · z={target.z:.3f} · {len(df)} points · '+
                   ('existing input' if target.id in ('F04','F06','R01') else ('processed' if row.download_complete else 'PARTIAL')))
            if len(df): ax.legend(fontsize=7, loc='best')
            else: ax.text(.5,.5,'No usable extracted measurements',ha='center',transform=ax.transAxes)
        fig.suptitle(f'All-target SPHEREx review · {page+1}/4\nUnrescaled visits; shaded continuum windows; supplied errors only', fontsize=14)
        savefig(fig, f'spectra_{page+1}')
    # Measured history amplitude/direction versus latest available continuum color.
    eligible = epochs[epochs.color_eligible].sort_values('mean_year')
    latest = eligible.groupby('id').tail(1)
    fig, axes = plt.subplots(1, 2, figsize=(12,5), constrained_layout=True)
    for ax, xcol, xlabel in zip(axes, ['W1_last3_over_first3','W1_p95_over_p05'],
            ['W1 late / early flux (median last / first 3 visits)', 'W1 variability amplitude (95th / 5th percentile)']):
        for r in latest.itertuples():
            x = history.loc[r.id,xcol]
            ax.errorbar(x, r.color, yerr=r.color_err, fmt='o', color=COLORS[0],
                        mfc=COLORS[0] if r.download_complete else 'white', ms=5, lw=.7)
            ax.annotate(r.id,(x,r.color),xytext=(3,3),textcoords='offset points',fontsize=7)
        ax.set(xlabel=xlabel, ylabel='Latest eligible SPHEREx Fν(3.6–3.9) / Fν(2.1–2.4 µm)')
        if xcol == 'W1_last3_over_first3': ax.axvline(1,color='grey',ls=':',lw=.8)
    fig.suptitle(f'History versus present continuum shape · {len(latest)}/24 eligible targets\nHost + AGN light; no temperature inference · open points: incomplete downloads', fontsize=12)
    savefig(fig,'history_color')
    # Paired differences: flux in one window versus color, errors are correlated.
    changes=[]
    for tid, group in eligible.groupby('id'):
        if len(group)<2: continue
        a,b=group.iloc[0],group.iloc[-1]
        ratio=b.long_mjy/a.long_mjy
        err=ratio*np.hypot(a.long_err_mjy/a.long_mjy,b.long_err_mjy/b.long_mjy)
        color=b.color/a.color
        colorerr=color*np.hypot(a.color_err/a.color,b.color_err/b.color)
        changes.append(dict(id=tid,flux_ratio=ratio,flux_ratio_err=err,color_ratio=color,color_ratio_err=colorerr,
                            first_visit=int(a.visit),last_visit=int(b.visit),download_complete=bool(b.download_complete)))
    pd.DataFrame(changes).to_csv(OUT/'visit_changes.csv',index=False)
    fig, ax=plt.subplots(figsize=(8,6),constrained_layout=True)
    for r in changes:
        ax.errorbar(r['flux_ratio'],r['color_ratio'],xerr=r['flux_ratio_err'],yerr=r['color_ratio_err'],
                    fmt='o',color=COLORS[1],mfc=COLORS[1] if r['download_complete'] else 'white',lw=.7)
        ax.annotate(r['id'],(r['flux_ratio'],r['color_ratio']),xytext=(4,4),textcoords='offset points',fontsize=8)
    ax.axhline(1,color='grey',lw=.7,ls=':');ax.axvline(1,color='grey',lw=.7,ls=':')
    ax.set(xlabel='Long-window brightness: last / first eligible visit',
           ylabel='Continuum color: last / first eligible visit',
           title=f'Visit-to-visit changes · {len(changes)}/24 eligible targets\nProvisional errors; shared-flux covariance and calibration require review')
    savefig(fig,'visit_changes')
    # Coverage: overview and zoom; catalogued coverage is distinct from quality.
    fig,axes=plt.subplots(1,2,figsize=(13,9),sharey=True,constrained_layout=True)
    for i,target in enumerate(sample.itertuples()):
        h=history.loc[target.id]
        axes[0].plot([h.W1_first_year,h.W1_last_year],[i,i],color='#aaa',lw=3)
        for r in epochs[epochs.id.eq(target.id)].itertuples():
            lo,hi=year(r.first_mjd),year(r.last_mjd)
            for ax in axes: ax.plot([lo,hi],[i,i],color=COLORS[r.visit%len(COLORS)],lw=4,marker='|')
        if not len(datasets[target.id]): axes[1].text(2025.1,i,'No accepted extraction',va='center',fontsize=8,color='#b34b35')
    axes[0].set_yticks(range(24),sample.id);axes[0].invert_yaxis()
    axes[0].set(xlim=(2013.5,2027),xlabel='Observed year',title='Grey: W1 span · colors: usable SPHEREx visits')
    axes[1].set(xlim=(2025,2027),xlabel='Observed year',title='SPHEREx zoom (a visit need not cover all wavelengths)')
    savefig(fig,'coverage')
    ngps=ROOT.parent/'sep23_data/reduction_20260924/products_p330e/spectra'
    archive=pd.read_csv(ROOT/'inputs/archival_spectra/index.csv')
    inventory_summary=dict(created_utc=datetime.now(timezone.utc).isoformat(),
        complete_additional_targets=int(inv[~inv.id.isin(['F04','F06','R01'])].download_complete.sum()),
        targets_with_valid_points=int((inv.valid_points>0).sum()),
        targets_with_eligible_color=len(latest),targets_with_paired_colors=len(changes),
        targets_with_local_archival_spectra=int(sample.id.isin(archive.id.unique()).sum()),
        targets_with_local_ngps_spectra=sum((ngps/f'{r.internal_id}.csv').exists() for r in sample.itertuples()),
        visits_by_target={r.id:int(r.visits) for r in inv.itertuples()},
        limitations=['Snapshot includes incomplete downloads; rerun after completion.',
                    'QR2 extractions and supplied R01 data have not had their calibration errors homogenized.',
                    'Color error bars use an independent-error approximation; calibration covariance is not propagated.',
                    'Colors require all three fixed cells in both rest-frame windows and >3 sigma window means.',
                    'Ineligible targets/visits remain in atlas and CSVs; color eligibility is not sample membership.',
                    'Color contains stellar, disc and dust emission; not a dust temperature.',
                    'Changes use earliest/latest eligible visits, which may differ between targets.',
                    'Existing WISE history statistics are descriptive and include host light.',
                    'No sample-wide physical-model recovery or correlation significance has been established.',
                    'Some spectra show sharp steps or isolated discrepant points; assess detector/aperture effects before physical fitting.',
                    'Accepted downloader status does not guarantee finite photometry: invalid values are explicitly counted and excluded here.'])
    (OUT/'summary.json').write_text(json.dumps(inventory_summary,indent=2)+'\n')
    sections=[('Spectra: every target',[(f'spectra_{i}',f'Spectral atlas {i}/4') for i in range(1,5)]),
              ('Candidate sample-wide plots',[('history_color','History and continuum shape'),('visit_changes','Changes between visits'),('coverage','Time coverage')])]
    parts=['<!doctype html><meta charset="utf-8"><title>SPHEREx full-sample review</title>',
           '<style>body{font:16px system-ui;max-width:1200px;margin:40px auto;padding:0 20px;color:#24313a}img{width:100%}table{border-collapse:collapse}td,th{padding:6px;border-bottom:1px solid #ddd}a{color:#247ba0}</style>',
           '<h1>SPHEREx: all 24 targets</h1>',f'<p>Snapshot: {inventory_summary["created_utc"]}. Exploratory review, including incomplete downloads.</p>',
           '<p>Fluxes retain the original scale. Colors describe total extracted light. Formal error bars are provisional until calibration and extraction quality are checked.</p>',
           '<ul>'+''.join(f'<li>{html.escape(x)}</li>' for x in inventory_summary['limitations'])+'</ul>']
    for title, figures in sections:
        parts.append(f'<h2>{title}</h2>')
        for name,label in figures: parts.append(f'<h3>{label}</h3><a href="{name}.pdf">Figure PDF</a><img src="{name}.png" alt="{label}">')
    parts += ['<h2>Inventory</h2>',inv.drop(columns=['snapshot_sha256','source']).to_html(index=False),
              '<p><a href="visit_measurements.csv">Visit measurements CSV</a> · <a href="inventory.csv">Full inventory CSV</a></p>']
    (OUT/'index.html').write_text('\n'.join(parts))
    print(json.dumps(inventory_summary,indent=2))


if __name__=='__main__': main()
