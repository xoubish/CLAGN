"""Two science tests on a shared wavelength axis, using the saved shell physics.

Illustrative fourfold step responses, not target fits or detection forecasts.
The comparison in B fits a single hot-dust amplitude, not the entire hot SED.
"""
from pathlib import Path
import csv
import hashlib
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import make_fig2_diagnostics as shell

HERE = Path(__file__).resolve().parent
BLUE, RED, INK, GOLD = '#24658a', '#bc5039', '#283542', '#a37a25'
AGE = 10.0  # initial silicate inner-radius light-crossing times


def main():
    sample_path = HERE/'inputs/jwst_sample_cycle6.csv'
    sample = list(csv.DictReader(sample_path.open()))
    z = np.array([float(row['z']) for row in sample])
    common = [4.9/(1+z.min()), 27.9/(1+z.max())]
    hot_end = 5.0/(1+z.max())
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.labelsize': 11.5, 'text.color': INK,
                         'axes.labelcolor': INK, 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(7.6, 6.45), facecolor='white')
    fig.text(.115, .98, 'What MIRI adds: two tests of the dust response',
             fontsize=13, weight='bold', va='top')
    fig.text(.115, .934, '24 histories  →  one warm spectrum per target  →  two physical tests',
             fontsize=11)
    coverage = fig.add_axes([.115, .843, .86, .056])
    coverage.set(xlim=(2, common[1]), ylim=(0, 1))
    coverage.axis('off')
    coverage.axvspan(2, hot_end, color='#e8e5ef')
    coverage.text((2+hot_end)/2, .5, 'Hot\ndust', ha='center', va='center', fontsize=9.5)
    coverage.axvspan(*common, color='#dce9ec')
    coverage.text(np.mean(common), .5, 'MIRI/MRS: warm continuum + silicate pair',
                  ha='center', va='center', fontsize=11, weight='bold', color=BLUE)
    a = fig.add_axes([.115, .525, .86, .263])
    b = fig.add_axes([.115, .168, .86, .263])
    a.set_title('A   How strongly does warm dust remember the history?',
                loc='left', fontsize=12, weight='bold', pad=12)
    b.set_title('B   Does the dust boundary move, or only its temperature?',
                loc='left', fontsize=12, weight='bold', pad=12)
    for ax in (a, b):
        ax.set_xlim(2, common[1])
        ax.axvspan(*common, color='#f1f6f7', zorder=0)
        ax.axvspan(8, 13, color='#dce7eb', zorder=0)
        for w in (9.7, 18.):
            ax.axvline(w, color=GOLD, lw=.85, ls=':', zorder=1)
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['left', 'bottom']].set_color('#a7b0b6')
        ax.tick_params(color='#a7b0b6', labelsize=10.5)
        ax.set_xticks([2, 5, 8, 10, 13, 15, 18, 21])
    a.set(ylim=(.64, 1.70), yticks=[.75, 1, 1.25, 1.5],
          ylabel='Flux / current-power\nequilibrium')
    a.tick_params(labelbottom=False)
    a.axhline(1, color='#6e7981', lw=1, ls=(0, (3, 2)))
    a.text(12.6, 1.018, 'Current-power equilibrium', fontsize=10, va='bottom', color='#5e6971')
    a.text(14.5, 1.58, 'Fading: warm excess', color=RED, fontsize=11.5, weight='bold')
    a.text(14.5, .68, 'Rising: warm deficit', color=BLUE, fontsize=11.5, weight='bold')
    a.text(10.5, 1.635, '8–13 µm', ha='center', fontsize=10, color='#52636d')
    b.set(ylim=(-44, 53), yticks=[-40, -20, 0, 20, 40],
          ylabel='Change from delayed\nfixed dust (%)', xlabel='Rest wavelength (µm)')
    b.axhline(0, color='#6e7981', lw=1, ls=(0, (3, 2)))
    b.text(5.1, 2.4, 'Delayed fixed dust', fontsize=10, color='#5e6971')
    b.text(14.5, 44, 'Rising: evolving boundary', color=BLUE, fontsize=10.5, weight='bold')
    b.text(14.5, -40, 'Fading: evolving boundary', color=RED, fontsize=10.5, weight='bold')
    for w in (9.7, 18.):
        b.text(w, 51 if w == 9.7 else 16, f'{w:g} µm silicate', ha='center', va='top', fontsize=9.5, color=GOLD)
    b.text(2.3, 47, 'Equal 2–4 µm luminosity', fontsize=9, color='#5e6971')
    rows, calculations = [], []
    for ratio, color in [(4., BLUE), (.25, RED)]:
        fixed, evolving, equilibrium_error = shell.spectra(ratio, age=AGE)
        equilibrium = shell.spectra(ratio, age=1e5)[0]
        _, matched, scale = shell.match_hot_luminosity(fixed, evolving)
        memory = fixed/equilibrium
        residual = 100*(matched/fixed-1)
        a.plot(shell.WAVE, memory, color=color, lw=2.4)
        b.plot(shell.WAVE, residual, color=color, lw=2.4)
        assert np.all(np.isfinite(memory)) and np.all(np.isfinite(residual))
        rows.extend(zip([ratio]*len(shell.WAVE), shell.WAVE, memory, residual))
        calculations.append(dict(luminosity_ratio=ratio, hot_luminosity_amplitude_scale=scale,
                                 equilibrium_relative_error=equilibrium_error,
                                 memory_ratio_at_12um=float(np.interp(12, shell.WAVE, memory)),
                                 evolution_percent_at_12um=float(np.interp(12, shell.WAVE, residual))))
    fig.text(.115, .059, 'Continuum + silicate profiles', fontsize=11, weight='bold', color=BLUE)
    fig.text(.115, .030, 'Measure memory; test structural change', fontsize=10)
    fig.text(.575, .059, 'Spatial profiles + PAH bands', fontsize=11, weight='bold', color=BLUE)
    fig.text(.575, .030, 'Separate nuclear and host emission', fontsize=10)
    fig.savefig(HERE/'fig2_miri.pdf')
    fig.savefig('/tmp/clagn-fig2-miri.png', dpi=160)
    plt.close(fig)
    with (HERE/'inputs/fig2_miri_curves.csv').open('w') as handle:
        writer = csv.writer(handle)
        writer.writerow(['luminosity_ratio', 'rest_um', 'delayed_over_equilibrium',
                         'hot_matched_evolving_vs_fixed_percent'])
        writer.writerows(rows)
    provenance = dict(
        figure='fig2_miri.pdf', generator='make_fig2_miri.py',
        status='Illustrative shell calculation; no observed target fits or detection significance.',
        physics='Unchanged make_fig2_diagnostics.py; optically thin silicate/graphite shells.',
        age_initial_silicate_Rin_over_c=AGE,
        panel_A='Fixed distribution with delayed heating divided by its current-power equilibrium.',
        panel_B='100*(evolving/fixed-1); evolving amplitude matched to fixed integrated rest 2–4 micron luminosity. Boundary adjusts immediately as changed illumination arrives. Signs depend on this prescription and normalization.',
        coverage=dict(mrs_common_rest_um=common, spherex_common_rest_upper_um=hot_end,
                      note='Intersection across actual sample redshifts, nominal wavelength limits only.'),
        normalization_note='Matching one hot-band luminosity does not establish identical near-IR spectra.',
        host_note='Spatial/PAH payoff is labelled; host emission is not simulated in these curves.',
        sample_sha256=hashlib.sha256(sample_path.read_bytes()).hexdigest(),
        model_sha256=hashlib.sha256((HERE/'make_fig2_diagnostics.py').read_bytes()).hexdigest(),
        calculations=calculations)
    (HERE/'inputs/fig2_miri_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    main()
