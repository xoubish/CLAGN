"""Figure 2: MRS extraction -> observables -> history comparisons -> science.

The spatial profile is an explicitly labelled schematic. The spectrum uses
unchanged shell physics to locate observables; it is not a target fit, a
forecast, or a simulated observation. No model-discrimination claim is plotted.
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
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import make_fig2_diagnostics as shell

HERE = Path(__file__).resolve().parent
BLUE, INK, GREEN = '#24658a', '#283542', '#42806c'
PURPLE, GOLD, GREY = '#805298', '#a37a25', '#727e86'


def box(fig, x, y, w, h, color, edge):
    patch = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.006,rounding_size=0.008',
                           transform=fig.transFigure, facecolor=color, edgecolor=edge,
                           linewidth=.7, zorder=-1)
    fig.add_artist(patch)


def arrow(fig, start, end, color=GREY):
    fig.add_artist(FancyArrowPatch(start, end, transform=fig.transFigure,
                                  arrowstyle='-|>', mutation_scale=11,
                                  linewidth=1, color=color))


def text(fig, x, y, message, size=10, color=INK, weight='normal', **kwargs):
    return fig.text(x, y, message, fontsize=size, color=color, weight=weight, **kwargs)


def main():
    sample_path = HERE/'inputs/jwst_sample_cycle6.csv'
    sample = list(csv.DictReader(sample_path.open()))
    redshift = np.array([float(row['z']) for row in sample])
    common = [4.9/(1+redshift.min()), 27.9/(1+redshift.max())]
    fixed = shell.spectra(4., age=10.)[0]
    normalized = fixed/np.interp(5., shell.WAVE, fixed)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.labelsize': 9.5, 'text.color': INK,
                         'axes.labelcolor': INK, 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(7.6, 6.9), facecolor='white')
    text(fig, .055, .983, 'What we will measure — and what those measurements decide',
         size=12.5, weight='bold', va='top')
    text(fig, .055, .929, '1  Separate nucleus and host', size=11, weight='bold')
    text(fig, .428, .929, '2  Extract the spectral diagnostics', size=11, weight='bold')

    spatial = fig.add_axes([.09, .747, .23, .152])
    x = np.linspace(-3, 3, 300)
    nuclear = np.exp(-.5*(x/.43)**2)
    host = .42*np.exp(-.5*(x/1.55)**2)
    spatial.plot(x, nuclear+host, color=INK, lw=1.6)
    spatial.plot(x, nuclear, color=BLUE, lw=1.8)
    spatial.plot(x, host, color=GREEN, lw=1.6, ls='--')
    spatial.set(xlim=(-3, 3), ylim=(0, 1.68), xticks=[], yticks=[],
                xlabel='Position across nucleus', ylabel='Surface brightness')
    spatial.text(-2.8, 1.49, 'MRS slice: schematic', fontsize=8.3, color=GREY)
    spatial.text(1.0, .94, 'Nucleus', fontsize=8.4, color=BLUE)
    spatial.text(1.4, .37, 'Host', fontsize=8.4, color=GREEN)
    spatial.spines[['top', 'right']].set_visible(False)
    arrow(fig, (.344, .821), (.396, .821), BLUE)

    spectrum = fig.add_axes([.428, .747, .535, .152])
    spectrum.plot(shell.WAVE, normalized, color=BLUE, lw=1.9)
    spectrum.axvspan(8, 13, color='#e2edf0', zorder=0)
    spectrum.set(xlim=common, ylim=(0, 4), yticks=[], xticks=[5, 10, 15, 20],
                 xlabel='Rest wavelength (µm)')
    spectrum.tick_params(labelsize=9, length=3)
    spectrum.spines[['top', 'right', 'left']].set_visible(False)
    for w in (6.2, 7.7, 11.3):
        spectrum.plot([w, w], [3.0, 3.23], color=GREEN, lw=1.1)
    spectrum.text(6.9, 3.34, 'PAHs', fontsize=8.3, color=GREEN, ha='center')
    for w in (7.65, 14.32):
        spectrum.plot([w, w], [3.57, 3.77], color=PURPLE, lw=1.1)
    spectrum.text(14.32, 3.81, '[Ne V]', fontsize=8, color=PURPLE, ha='center')
    spectrum.text(7.65, 3.81, '[Ne VI]', fontsize=8, color=PURPLE, ha='center')
    for w, xytext, label in [(9.7, (12.8, 2.35), '9.7 µm'), (18., (18., 1.5), '18 µm')]:
        spectrum.annotate(label, xy=(w, np.interp(w, shell.WAVE, normalized)), xytext=xytext,
                          ha='center', fontsize=8.7, color=GOLD,
                          arrowprops=dict(arrowstyle='-', color=GOLD, lw=.8))
    spectrum.text(10.5, .12, '8–13 µm', ha='center', fontsize=8, color=BLUE)
    spectrum.text(20.9, 3.25, 'Illustrative\nspectrum', ha='right', fontsize=8, color=GREY)

    # Inventory: measurements, not independently inferred physical parameters.
    cards = [(.055, BLUE, '#eff5f7', 'Warm continuum',
              ['8–13 µm nuclear luminosity', 'Continuum levels + slopes', 'Near 5, 12 and 18 µm']),
             (.365, GOLD, '#faf6ed', 'Paired silicate features',
              ['9.7 / 18 µm strengths', 'Peak, width and full profile', 'Jointly fitted continuum']),
             (.675, GREEN, '#f0f6f2', 'Host and gas emission',
              ['PAH fluxes + equivalent widths', '[Ne V/VI], H₂ fluxes or limits', 'Spatial extent + nuclear fraction'])]
    for left, color, fill, title, lines in cards:
        box(fig, left, .527, .285, .154, fill, '#c8d2d6')
        text(fig, left+.012, .653, title, size=10.2, color=color, weight='bold')
        for y, line in zip([.618, .586, .554], lines):
            text(fig, left+.012, y, line, size=9.1)
    # Every observable feeds the joint historical analysis.
    for center in (.1975, .5075, .8175):
        arrow(fig, (center, .519), (center, .505))
    box(fig, .055, .426, .905, .073, '#f4f4f5', '#c8d2d6')
    text(fig, .507, .472, 'Connect every spectrum to its own recorded history', size=10.5,
         weight='bold', ha='center')
    text(fig, .507, .442, '2010 W3/W4  •  optical + WISE light curves  •  SPHEREx + optical spectra',
         size=9.4, ha='center')
    for center in (.275, .742):
        arrow(fig, (center, .423), (center, .400), BLUE)

    box(fig, .055, .179, .438, .216, '#eff5f7', '#b4c9d2')
    box(fig, .522, .179, .438, .216, '#f5f0f7', '#cbbbd3')
    text(fig, .070, .369, 'Q1  How much spectral memory?',
         size=10, weight='bold', color=BLUE)
    text(fig, .070, .332, 'Fit the warm change relative to 2010.', size=9.5)
    text(fig, .070, .299, 'Measure excess / deficit versus current power.', size=9.5)
    text(fig, .070, .266, 'Compare direction, amplitude and event age.', size=9.5)
    text(fig, .070, .211, '→ Memory strength and persistence',
         size=9.6, weight='bold', color=BLUE)
    text(fig, .537, .369, 'Q2  Is dust evolution required?',
         size=10, weight='bold', color=PURPLE)
    text(fig, .537, .332, 'Fit continuum + silicates to the full history.', size=9.5)
    text(fig, .537, .299, 'Test delayed, fixed dust across geometries.', size=9.5)
    text(fig, .537, .266, 'Assess coherent spectral residuals.', size=9.5)
    text(fig, .537, .211, '→ Fixed dust / structural evolution',
         size=9.6, weight='bold', color=PURPLE)

    text(fig, .055, .128, 'Complementary tests from the same spectra', size=10.5, weight='bold')
    text(fig, .055, .089, 'Silicates + optical state + ionic lines → assess obscuration and intrinsic change', size=9.7)
    text(fig, .055, .052, 'PAHs + H₂ + spatial structure → characterize host star formation and gas excitation', size=9.7)
    fig.savefig(HERE/'fig2_miri.pdf')
    fig.savefig('/tmp/clagn-fig2-miri.png', dpi=165)
    plt.close(fig)

    np.savetxt(HERE/'inputs/fig2_miri_workflow_spectrum.csv',
               np.column_stack([shell.WAVE, normalized]), delimiter=',',
               header='rest_um,illustrative_nuclear_relative_Fnu', comments='')
    provenance = dict(
        figure='fig2_miri.pdf', generator='make_fig2_miri.py',
        status='Measurement-to-inference workflow, not a forecast or simulated MRS observation.',
        spatial_schematic='Dimensionless Gaussian nucleus (sigma=.43) plus extended host (amplitude=.42, sigma=1.55); arbitrary illustrative shape, not a calibrated MRS PSF or fitted target.',
        spectrum='Unchanged optically thin shell calculation; delayed fixed dust after a fourfold rise at age 10 initial silicate Rin/c; normalized at rest 5 micron. Used only to locate continuum and feature observables.',
        diagnostic_markers=dict(PAH_rest_um=[6.2, 7.7, 11.3], NeVI_rest_um=7.65, NeV_rest_um=14.32,
                                interpretation='Positions only; no PAH or ionic line emission synthesized.'),
        mrs_common_rest_um=common,
        science_tests=dict(memory='Warm change and excess/deficit relative to current power, versus direction/amplitude/event age.',
                           evolution='Joint continuum/silicate residuals relative to delayed fixed-dust models across allowed geometries.'),
        measurements=['nuclear 8–13 micron continuum luminosity', 'continuum levels and slopes',
                      'paired silicate strengths and profiles', 'PAH fluxes and equivalent widths',
                      'Ne V/VI and H2 fluxes or limits', 'spatial extent and nuclear fraction'],
        analysis_plan='review/mrs_measurement_plan.md',
        sample_sha256=hashlib.sha256(sample_path.read_bytes()).hexdigest(),
        model_sha256=hashlib.sha256((HERE/'make_fig2_diagnostics.py').read_bytes()).hexdigest())
    (HERE/'inputs/fig2_miri_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    main()
