"""Figure 2: the spectral imprint of past activity, from the shared shell model.

Panels B and C use the optically thin, delayed-heating shell calculation in
make_fig2_diagnostics.py (fixed dust distribution only) so that Figures 2 and 3
rest on one model family. Panel A is a cartoon; panel D uses the published MRS
point-spread-function width against hosts of fixed physical size at the sample
redshifts. Nothing here is fitted to a target: amplitudes are illustrative.
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
from matplotlib.patches import Circle, Ellipse
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
from matplotlib.transforms import Affine2D
import numpy as np
from astropy.cosmology import Planck18
from scipy.integrate import trapezoid

import make_fig2_diagnostics as shell

HERE = Path(__file__).resolve().parent
BLUE, RED, GOLD = '#24658a', '#bc5039', '#bb8c35'
PURPLE, GREEN, INK = '#72518b', '#4c826e', '#283542'
BAND = (8., 13.)                      # rest-frame diagnostic band, microns
AGES_B = (3., 10.)                    # panel B snapshots, initial silicate R_in/c
AGE_GRID = np.geomspace(.5, 60., 36)  # panel C event ages, same units
RATIOS_C = (4., 2., .5, .25)
EQUILIBRIUM_AGE = 1e5
SCATTER_DEX, GROUP_SIZE = .20, 12     # planning numbers quoted in the text
PSF_WAVELENGTH_UM = 10.
HOST_HALF_LIGHT_KPC = 3.
CH2_FIELD_ARCSEC = (4.0, 4.8)         # MIRI MRS channel 2 field of view (JDox)


def mrs_psf_fwhm(wave_um):
    """Empirical MRS PSF FWHM, arcsec (Law et al. 2023, AJ 166, 45, eq. for FWHM(lambda))."""
    return .033*wave_um + .106


def panel(ax, title):
    ax.set_title(title, loc='left', fontsize=11, fontweight='bold', pad=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_color('#9ca4a9')
    ax.tick_params(length=3, color='#9ca4a9', labelsize=10.5)


def agn_cartoon(ax, center, bright_outer, vertical=0):
    """Identical engines and geometry; only the outer dust emissivity differs.

    Visual shorthand for delayed response, not an echo surface or a resolved
    prediction. Both inner hot components have responded to the current power.
    """
    initial_patches = len(ax.patches)
    outer = '#d58a32' if bright_outer else '#9baab6'
    for scale, alpha in [(1.16, .045), (1.08, .07), (1., .15)]:
        ax.add_patch(Ellipse((center, 0), 1.85*scale, 1.02*scale,
                             facecolor=outer, edgecolor='none', alpha=alpha))
    ax.add_patch(Ellipse((center, 0), 1.06, .57, facecolor='white', edgecolor='none'))
    for radius, count, phase in [(.74, 22, .05), (.92, 28, .17)]:
        theta = np.linspace(0, 2*np.pi, count, endpoint=False)+phase
        for angle in theta:
            ax.add_patch(Circle((center+radius*np.cos(angle), .54*radius*np.sin(angle)),
                                .041, facecolor=outer, edgecolor='white', lw=.25,
                                alpha=.9 if bright_outer else .65, zorder=2))
    theta = np.linspace(0, 2*np.pi, 22, endpoint=False)
    for angle in theta:
        ax.add_patch(Circle((center+.52*np.cos(angle), .27*np.sin(angle)),
                            .029, facecolor='#e8b051', edgecolor='none', zorder=3))
    for width, height, color in [(.49, .23, '#fff0c5'), (.38, .16, '#e9b044'), (.26, .10, '#dc813b')]:
        ax.add_patch(Ellipse((center, 0), width, height, facecolor=color, edgecolor='none', zorder=5))
    ax.add_patch(Circle((center, 0), .068, facecolor=INK, edgecolor='white', lw=.45, zorder=6))
    for patch in list(ax.patches)[initial_patches:]:
        patch.set_transform(Affine2D().translate(0, vertical) + ax.transData)


def band_luminosity(spectrum):
    use = (shell.WAVE >= BAND[0]) & (shell.WAVE <= BAND[1])
    nu = shell.C/(shell.WAVE[use]*1e-6)
    return -trapezoid(spectrum[use], nu)


def memory_curves():
    """Fixed-distribution spectra divided by equilibrium at the current power."""
    equilibrium = {r: shell.spectra(r, age=EQUILIBRIUM_AGE)[0] for r in set(RATIOS_C)}
    spectra_b = {}
    for ratio in (4., .25):
        for age in AGES_B:
            fixed = shell.spectra(ratio, age=age)[0]
            spectra_b[(ratio, age)] = fixed/equilibrium[ratio]
    offsets = {}
    for ratio in RATIOS_C:
        eq_band = band_luminosity(equilibrium[ratio])
        offsets[ratio] = np.array([np.log10(band_luminosity(shell.spectra(ratio, age=a)[0])/eq_band)
                                   for a in AGE_GRID])
    return spectra_b, offsets


def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.labelsize': 11, 'axes.labelcolor': INK,
                         'text.color': INK, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    sample_path = HERE/'inputs/jwst_sample_cycle6.csv'
    with sample_path.open() as handle:
        rows = list(csv.DictReader(handle))
    zkey = next(key for key in rows[0] if key.lower() in ['z', 'redshift'])
    redshift = np.array([float(row[zkey]) for row in rows])
    common = (4.9/(1+redshift.min()), 27.9/(1+redshift.max()))
    assert common[0] < 9.7 < 18 < common[1]
    spectra_b, offsets = memory_curves()
    mean_sigma = SCATTER_DEX/np.sqrt(GROUP_SIZE)

    fig = plt.figure(figsize=(7.6, 2.55), facecolor='white')
    a = fig.add_axes([.012, .235, .165, .61])
    b = fig.add_axes([.268, .235, .175, .61])
    c = fig.add_axes([.548, .235, .160, .61])
    d = fig.add_axes([.805, .235, .165, .61])

    # A: same current power, different histories.
    a.set_title('A  AGN history', loc='left', fontsize=11, fontweight='bold', pad=10)
    a.set(xlim=(-1.18, 1.18), ylim=(-1.65, 1.65))
    a.axis('off')
    for vertical, bright, colour, label in [(.76, False, BLUE, 'Rising'), (-.88, True, RED, 'Fading')]:
        agn_cartoon(a, 0, bright, vertical)
        a.text(-1.05, vertical+.66, label, color=colour, fontsize=10.5, fontweight='bold')
        now = vertical+.62
        past = now+(.30 if bright else -.30)
        a.plot([.10, .46, .46, .96], [past, past, now, now], color=colour, lw=1.8,
               ls='--' if bright else '-', solid_capstyle='round')
        a.scatter([.96], [now], s=16, color=INK, zorder=7)
    a.text(.5, -.15, 'Same current power', ha='center', fontsize=10, transform=a.transAxes)

    # B: warm continuum relative to equilibrium at the current power.
    panel(b, 'B  Continuum memory')
    b.axvspan(common[0], common[1], color='#eef0f2', zorder=0)
    b.axvspan(*BAND, color='#dfe3e7', zorder=0)
    b.axhline(1, color='#69747b', ls=':', lw=1.2)
    for centre in [9.7, 18.]:
        b.axvline(centre, color=GOLD, ls=':', lw=.8, alpha=.8)
    styles = {AGES_B[0]: dict(ls='-', lw=2.1), AGES_B[1]: dict(ls='--', lw=1.7)}
    for (ratio, age), curve in spectra_b.items():
        b.plot(shell.WAVE, curve, color=RED if ratio < 1 else BLUE, **styles[age])
    b.set(xscale='log', yscale='log', xlim=(3, common[1]), ylim=(.33, 3.0),
          xlabel='Rest λ (µm)', ylabel='Continuum / equilibrium')
    b.xaxis.set_major_locator(FixedLocator([5, 10, 20]))
    b.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
    b.xaxis.set_minor_formatter(NullFormatter())
    b.yaxis.set_major_locator(FixedLocator([.5, 1, 2]))
    b.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
    b.yaxis.set_minor_formatter(NullFormatter())
    b.text(3.15, 2.40, 'Fading: excess', color=RED, fontsize=9.5, va='center')
    b.text(3.15, .375, 'Rising: deficit', color=BLUE, fontsize=9.5, va='center')
    b.text(20.8, 1.03, 'equilibrium', ha='right', va='bottom', color='#465159', fontsize=8.5)
    at = lambda ratio, age, w: spectra_b[(ratio, age)][np.argmin(abs(shell.WAVE-w))]
    b.text(20.8, at(4., AGES_B[1], 20.8)*1.07, f't = {AGES_B[1]:g}', color=BLUE,
           fontsize=8.5, ha='right', va='bottom')
    b.text(20.8, at(4., AGES_B[0], 20.8)/1.07, f't = {AGES_B[0]:g}', color=BLUE,
           fontsize=8.5, ha='right', va='top')
    b.text(10.3, 2.85, '8–13 µm', color='#5d6870', fontsize=8.5, ha='center', va='center')

    # C: persistence of the 8-13 um offset against event age, with the planned precision.
    panel(c, 'C  Persistence')
    c.axhspan(-mean_sigma, mean_sigma, color='#e4e7ea', zorder=0)
    c.axhline(0, color='#69747b', ls=':', lw=1.2)
    for age in AGES_B:
        c.axvline(age, color='#b9c0c6', lw=.7, ls=(0, (2, 2)), zorder=0)
    for ratio, curve in offsets.items():
        colour = RED if ratio < 1 else BLUE
        big = ratio in (4., .25)
        c.plot(AGE_GRID, curve, color=colour, lw=2.1 if big else 1.2, alpha=1 if big else .75)
        label = {4.: 'L×4', 2.: 'L×2', .5: 'L/2', .25: 'L/4'}[ratio]
        c.text(AGE_GRID[0]*1.08, curve[0]+(.03 if ratio < 1 else -.03), label, color=colour,
               fontsize=8.5, ha='left', va='bottom' if ratio < 1 else 'top')
    c.set(xscale='log', xlim=(.4, 60), ylim=(-.52, .62),
          xlabel='Event age ($R_{\\rm in}/c$)', ylabel='8–13 µm offset (dex)')
    c.xaxis.set_major_locator(FixedLocator([1, 3, 10, 30]))
    c.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
    c.xaxis.set_minor_formatter(NullFormatter())
    c.set_yticks([-.4, -.2, 0, .2, .4])
    c.text(58, -mean_sigma-.02, '12-target mean ±1σ', ha='right', va='top',
           color='#5d6870', fontsize=8.5)
    c.text(AGES_B[0], .58, 'B', ha='center', va='top', color='#8a9299', fontsize=8.5)
    c.text(AGES_B[1], .58, 'B', ha='center', va='top', color='#8a9299', fontsize=8.5)

    # D: unresolved nucleus against hosts of fixed physical size at the sample redshifts.
    panel(d, 'D  Nucleus, host')
    x = np.linspace(-3.2, 3.2, 800)
    fwhm = mrs_psf_fwhm(PSF_WAVELENGTH_UM)
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    nucleus = np.exp(-.5*(x/sigma)**2)
    d.fill_between(x, 0, nucleus, color=PURPLE, alpha=.22)
    d.plot(x, nucleus, color=PURPLE, lw=1.7)
    hosts = []
    for z, style in [(float(np.median(redshift)), dict(ls='-', lw=1.6)),
                     (float(redshift.max()), dict(ls='--', lw=1.4))]:
        kpc_per_arcsec = Planck18.kpc_proper_per_arcmin(z).value/60
        r_e = HOST_HALF_LIGHT_KPC/kpc_per_arcsec
        host = np.exp(-1.678*np.abs(x)/r_e)
        d.plot(x, host, color=GREEN, **style)
        hosts.append(dict(z=z, kpc_per_arcsec=kpc_per_arcsec, half_light_arcsec=r_e))
    half_field = CH2_FIELD_ARCSEC[0]/2
    for sign in (-1, 1):
        d.axvline(sign*half_field, color='#8b9197', lw=.8, ls=(0, (3, 2)), ymax=.36)
    d.set(xlim=(-3.2, 3.2), ylim=(0, 1.24), xlabel='Offset (arcsec)', ylabel='Relative brightness')
    d.set_xticks([-2, 0, 2]); d.set_yticks([1], ['1'])
    d.text(-3.05, 1.0, f'Nucleus\n(PSF, {PSF_WAVELENGTH_UM:g} µm)', color=PURPLE, fontsize=8.5,
           ha='left', va='top', linespacing=1.0)
    d.text(1.95, .74, f'{HOST_HALF_LIGHT_KPC:g} kpc host\nz = {hosts[0]["z"]:.2f}', color=GREEN,
           fontsize=8.5, ha='center', va='center', linespacing=1.0)
    d.annotate(f'z = {hosts[1]["z"]:.2f}', xy=(.45, float(np.exp(-1.678*.45/hosts[1]['half_light_arcsec']))),
               xytext=(1.25, .47), color=GREEN, fontsize=8.5, ha='left', va='center',
               arrowprops=dict(arrowstyle='-', color=GREEN, lw=.6, shrinkB=2))
    d.text(-half_field-.1, .44, 'Ch. 2\nfield', color='#5d6870', fontsize=8, ha='right',
           va='center', linespacing=1.0)
    fig.savefig(HERE/'fig2_diagnostics.pdf', dpi=220)
    plt.close(fig)

    # Save the plotted model curves and provenance.
    with (HERE/'inputs/fig2_memory_response.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['panel', 'luminosity_ratio', 'age_Rin_over_c', 'rest_um',
                         'fixed_over_equilibrium', 'log10_band_offset'])
        for (ratio, age), curve in spectra_b.items():
            writer.writerows(('B', ratio, age, w, v, '') for w, v in zip(shell.WAVE, curve))
        for ratio, curve in offsets.items():
            writer.writerows(('C', ratio, age, '', '', v) for age, v in zip(AGE_GRID, curve))
    offset_table = {str(r): {f'{a:g}': float(np.interp(a, AGE_GRID, curve))
                             for a in [1, 3, 10, 30]} for r, curve in offsets.items()}
    provenance = dict(
        figure='fig2_diagnostics.pdf', generator='make_fig2_memory.py',
        scope='Illustrative delayed-heating calculation for a fixed dust distribution; not target spectra, measured lags or a population-power forecast. Panels B and C share the shell model, grain properties and step history of Figure 3 (make_fig2_diagnostics.py).',
        presentation='A: cartoons of equal current power with different histories. B: fixed-distribution spectrum divided by equilibrium at the current power at two event ages. C: log10 ratio of 8-13 um luminosity to equilibrium against event age for four step amplitudes, with the planned 12-target mean precision. D: MRS PSF against exponential hosts of fixed physical size.',
        actual_targets_fitted=[], sample_z_range=[float(redshift.min()), float(redshift.max())],
        common_mrs_rest_um=list(common),
        panel_a=dict(cartoons='Same central-engine brightness, hot inner rim, outer dust geometry and clump count; outer dust brightness alone differs to depict delayed illumination. Step histories are qualitative; the dot marks the observation.'),
        panel_b=dict(model='make_fig2_diagnostics.spectra(ratio, age)[0] (fixed distribution) divided by the same at age 1e5 (equilibrium at the new luminosity)',
                     luminosity_ratios=[4., .25], ages_initial_silicate_Rin_over_c=list(AGES_B),
                     ordinate='Fnu divided by the equilibrium Fnu at the current power; unity is equilibrium.',
                     shading='Light: common nominal rest-frame MRS coverage for the sample. Dark: the 8-13 um diagnostic band. Gold dotted: 9.7 and 18 um.',
                     assumptions='Optically thin isotropic shells, dM/dlnr proportional to r to 100 initial silicate inner radii, 0.1 um silicate/graphite grains, uniform echo delays on [0, 2r/c]. Amplitudes are illustrative.'),
        panel_c=dict(quantity='log10 of the 8-13 um luminosity of the delayed fixed-distribution spectrum over its equilibrium value at the current power',
                     luminosity_ratios=list(RATIOS_C), age_grid=[float(a) for a in AGE_GRID],
                     offsets_dex=offset_table,
                     precision_band=dict(half_width_dex=float(mean_sigma), formula='scatter/sqrt(n)',
                                         scatter_dex=SCATTER_DEX, n=GROUP_SIZE,
                                         meaning='±1 sigma of one 12-target group mean for the planning scatter quoted in the text; the rising-fading difference has 0.408*s = 0.082 dex.'),
                     age_units='Initial silicate inner-radius light-crossing times; the conversion to years depends on each target luminosity and is not drawn.'),
        panel_d=dict(psf_fwhm_arcsec=float(fwhm), psf_relation='FWHM = 0.033*lambda_um + 0.106 arcsec',
                     psf_reference='Law et al. 2023, AJ 166, 45, https://doi.org/10.3847/1538-3881/acdddc',
                     wavelength_um=PSF_WAVELENGTH_UM, nucleus='Gaussian with that FWHM, peak-normalised',
                     hosts=hosts, host_profile='Exponential, I proportional to exp(-1.678 r/r_e), peak-normalised; half-light radius 3 kpc proper',
                     cosmology='astropy Planck18', field='MIRI MRS channel 2 field of view 4.0 x 4.8 arcsec (JDox); dashed lines at ±2.0 arcsec',
                     meaning='Shows the angular contrast between an unresolved nucleus and a typical host, not a measured MRS PSF or resolved torus imaging.'),
        limitations=['No unique accretion-history recovery.',
                     'Fixed-distribution model only; Figure 3 addresses dust evolution.',
                     'Amplitudes depend on the assumed radial dust distribution and grain mix; they are not target predictions.',
                     'Ages are in model units; the light-year scale of R_in varies across the sample.',
                     'Long-wavelength response may average over epochs before the observed light curves.'],
        references=['https://doi.org/10.1051/0004-6361/201117750', 'https://doi.org/10.3847/1538-4357/ab6aa1',
                    'https://doi.org/10.1086/172149', 'https://doi.org/10.3847/1538-3881/acdddc',
                    'https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-observing-modes/miri-medium-resolution-spectroscopy'],
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sample_sha256=hashlib.sha256(sample_path.read_bytes()).hexdigest(),
        previous_figure='review/local_backups/fig2_schematic_20260928/ (all-schematic version replaced 2026-09-28)')
    (HERE/'inputs/fig2_diagnostics_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print('Figure 2 from the shell model; 8-13 um offsets (dex):')
    for ratio, table in offset_table.items():
        print(f'  L ratio {ratio}: ' + ', '.join(f'age {k}: {v:+.2f}' for k, v in table.items()))


if __name__ == '__main__':
    main()
