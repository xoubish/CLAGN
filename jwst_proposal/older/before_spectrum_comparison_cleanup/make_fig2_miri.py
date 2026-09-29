"""MIRI spectrum and two compact science tests, from unchanged shell physics.

The spectrum is a physical illustration; PAH/ionic ticks locate diagnostics
and do not synthesize line emission or imply detections.
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
BLUE, RED, INK = '#24658a', '#bc5039', '#283542'
PURPLE, GREEN, GOLD = '#805298', '#42806c', '#a37a25'
AGE = 10.0


def main():
    sample_path = HERE/'inputs/jwst_sample_cycle6.csv'
    sample = list(csv.DictReader(sample_path.open()))
    z = np.array([float(row['z']) for row in sample])
    common = [4.9/(1+z.min()), 27.9/(1+z.max())]
    hot_end = 5.0/(1+z.max())
    models, rows, calculations = {}, [], []
    for ratio in (4., .25):
        fixed, evolving, error = shell.spectra(ratio, age=AGE)
        equilibrium = shell.spectra(ratio, age=1e5)[0]
        _, matched, scale = shell.match_hot_luminosity(fixed, evolving)
        memory = fixed/equilibrium
        residual = 100*(matched/fixed-1)
        models[ratio] = (fixed, matched, memory, residual)
        assert np.all(np.isfinite(memory)) and np.all(np.isfinite(residual))
        rows.extend(zip([ratio]*len(shell.WAVE), shell.WAVE, memory, residual))
        calculations.append(dict(luminosity_ratio=ratio, hot_luminosity_amplitude_scale=scale,
                                 equilibrium_relative_error=error,
                                 memory_ratio_at_12um=float(np.interp(12, shell.WAVE, memory)),
                                 evolution_percent_at_12um=float(np.interp(12, shell.WAVE, residual))))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5,
                         'axes.labelsize': 10.5, 'text.color': INK,
                         'axes.labelcolor': INK, 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(7.6, 6.45), facecolor='white')
    fig.text(.11, .983, 'Recorded histories → MIRI spectra → dust physics',
             fontsize=13, weight='bold', va='top')
    coverage = fig.add_axes([.11, .883, .86, .045])
    coverage.set(xlim=(2, common[1]), ylim=(0, 1)); coverage.axis('off')
    coverage.axvspan(2, hot_end, color='#e8e5ef')
    coverage.text(2.05, .55, 'SPHEREx', ha='left', va='center', fontsize=8.5)
    coverage.axvspan(*common, color='#dce9ec')
    coverage.text(np.mean(common), .5, 'MIRI/MRS: the missing warm-dust spectrum',
                  ha='center', va='center', fontsize=11, weight='bold', color=BLUE)
    a = fig.add_axes([.11, .535, .86, .305])
    b = fig.add_axes([.11, .135, .365, .245])
    c = fig.add_axes([.605, .135, .365, .245])
    for ax in (a, b, c):
        ax.set_xlim(2, common[1])
        ax.axvspan(*common, color='#f1f6f7', zorder=0)
        ax.axvspan(8, 13, color='#dce7eb', zorder=0)
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['left', 'bottom']].set_color('#a7b0b6')
        ax.tick_params(color='#a7b0b6', labelsize=10)
        ax.set_xticks([2, 5, 10, 15, 20])
    a.set_title('A   Resolve the continuum and both silicate features',
                loc='left', fontsize=11.5, weight='bold', pad=10)
    fixed, matched, memory, residual = models[4.]
    norm = np.interp(5., shell.WAVE, fixed)
    a.plot(shell.WAVE, fixed/norm, color=INK, lw=2, label='Delayed fixed dust')
    a.plot(shell.WAVE, matched/norm, color=PURPLE, lw=2, ls='--', label='Evolving boundary')
    a.set(ylim=(0, 4.15), yticks=[0, 1, 2, 3], ylabel='Model flux density\n(relative to fixed dust at 5 µm)',
          xlabel='Rest wavelength (µm)')
    a.legend(loc='upper right', bbox_to_anchor=(1, .94), frameon=False, fontsize=9.5)
    # Diagnostic ticks are markers only: these features are not added to the model.
    for w in (6.2, 7.7, 11.3):
        a.plot([w, w], [3.25, 3.44], color=GREEN, lw=1.3)
        a.text(w, 3.47, 'PAH', ha='center', va='bottom', fontsize=8.5, color=GREEN)
    for w, label in [(7.65, '[Ne VI]'), (14.32, '[Ne V]')]:
        a.plot([w, w], [3.8, 3.92], color=PURPLE, lw=1.3)
        a.text(w, 3.94, label, ha='center', va='bottom', fontsize=8.5, color=PURPLE)
    for w, xytext, label in [(9.7, (12.3, 2.55), '9.7 µm silicate'),
                             (18., (17.1, 1.5), '18 µm silicate')]:
        a.annotate(label, xy=(w, np.interp(w, shell.WAVE, matched/norm)), xytext=xytext,
                   ha='center', color=GOLD, fontsize=10,
                   arrowprops=dict(arrowstyle='-', color=GOLD, lw=.9))
    a.text(2.35, 2.35, 'Rising example\nEqual 2–4 µm luminosity', fontsize=9, color='#5e6971')
    a.text(10.5, .13, '8–13 µm', ha='center', fontsize=8.5, color='#52636d')
    fig.text(.11, .451, 'PAH bands + spatial profiles: isolate the nucleus', fontsize=9.4, color=GREEN)
    fig.text(.65, .451, 'Ionic lines: test excitation', fontsize=9.4, color=PURPLE)
    b.set_title('B   Measure spectral memory', loc='left', fontsize=11, weight='bold', pad=12)
    c.set_title('C   Test boundary evolution', loc='left', fontsize=11, weight='bold', pad=12)
    b.set(ylim=(.63, 1.70), yticks=[.75, 1, 1.25, 1.5],
          ylabel='Flux / current-power\nequilibrium', xlabel='Rest wavelength (µm)')
    b.axhline(1, color='#6e7981', lw=1, ls=(0, (3, 2)))
    b.text(9, 1.028, 'Current-power equilibrium', fontsize=8, color='#5e6971')
    for ratio, color in [(4., BLUE), (.25, RED)]:
        b.plot(shell.WAVE, models[ratio][2], color=color, lw=2.2)
    b.text(10, 1.60, 'Fading: excess', color=RED, fontsize=10, weight='bold')
    b.text(10, .68, 'Rising: deficit', color=BLUE, fontsize=10, weight='bold')
    c.set(ylim=(-4, 53), yticks=[0, 20, 40],
          ylabel='Change from fixed dust (%)', xlabel='Rest wavelength (µm)')
    c.axhline(0, color='#6e7981', lw=1, ls=(0, (3, 2)))
    use = shell.WAVE >= common[0]
    c.plot(shell.WAVE[use], residual[use], color=PURPLE, lw=2.2, ls='--')
    c.text(3, 2, 'Delayed fixed dust', fontsize=8.5, color='#5e6971')
    c.text(3, 46, 'Same rising example as A', fontsize=9, color=PURPLE)
    for w in (9.7, 18.):
        c.axvline(w, color=GOLD, lw=.8, ls=':', zorder=1)
    fig.text(.11, .025, '24 recorded histories + warm spectra: measure memory and test dust evolution.',
             fontsize=10, weight='bold')
    fig.savefig(HERE/'fig2_miri.pdf')
    fig.savefig('/tmp/clagn-fig2-miri.png', dpi=160)
    plt.close(fig)
    with (HERE/'inputs/fig2_miri_curves.csv').open('w') as handle:
        writer = csv.writer(handle)
        writer.writerow(['luminosity_ratio', 'rest_um', 'delayed_over_equilibrium',
                         'hot_matched_evolving_vs_fixed_percent'])
        writer.writerows(rows)
    np.savetxt(HERE/'inputs/fig2_miri_spectrum.csv',
               np.column_stack([shell.WAVE, fixed/norm, matched/norm]), delimiter=',',
               header='rest_um,fixed_relative_Fnu,evolving_relative_Fnu', comments='')
    provenance = dict(
        figure='fig2_miri.pdf', generator='make_fig2_miri.py',
        status='Illustrative shell calculation; no observed target fits or detection significance.',
        physics='Unchanged make_fig2_diagnostics.py; optically thin silicate/graphite shells.',
        age_initial_silicate_Rin_over_c=AGE,
        panel_A='Fourfold rising event: delayed fixed and evolving spectra, normalized by fixed Fnu at rest 5 micron, after matching integrated rest 2–4 micron luminosity.',
        panel_B='Fixed distribution with delayed heating divided by its current-power equilibrium, for fourfold rises and fades.',
        panel_C='100*(evolving/fixed-1) for the same rising example as A; displayed over common MRS coverage. Boundary adjusts immediately as changed illumination arrives. Sign depends on this prescription and normalization.',
        diagnostic_markers=dict(PAH_rest_um=[6.2, 7.7, 11.3], NeVI_rest_um=7.65, NeV_rest_um=14.32,
                                interpretation='Positions only; no PAH or ionic emission synthesized, no detection claim.'),
        coverage=dict(mrs_common_rest_um=common, spherex_common_rest_upper_um=hot_end,
                      note='Intersection across actual sample redshifts, nominal wavelength limits only. SPHEREx extends to shorter wavelengths than displayed.'),
        normalization_note='Matching one hot-band luminosity does not establish identical near-IR spectra.',
        host_note='Spatial/PAH payoff is labelled; host emission is not simulated in these curves.',
        sample_sha256=hashlib.sha256(sample_path.read_bytes()).hexdigest(),
        model_sha256=hashlib.sha256((HERE/'make_fig2_diagnostics.py').read_bytes()).hexdigest(),
        calculations=calculations)
    (HERE/'inputs/fig2_miri_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    main()
