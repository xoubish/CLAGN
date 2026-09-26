"""Illustrate the proposed dust-response test, without inventing target fits.

Reproducible optically thin shell calculation with published grain efficiencies.
This is a controlled physical illustration, not a posterior predictive spectrum.
See inputs/fig2_diagnostics_provenance.json for assumptions and numerical checks.
"""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Patch, Wedge
import numpy as np
from scipy.integrate import trapezoid

HERE = Path(__file__).resolve().parent
INPUTS = HERE / 'inputs'
H, C, K = 6.62607015e-34, 299792458., 1.380649e-23
BLUE, ORANGE = '#24658a', '#be513b'
WAVE = np.geomspace(1., 30., 1000)


def efficiency(name, wave):
    lines = gzip.decompress((INPUTS/'dust_opacity'/name).read_bytes()).decode().splitlines()
    start = next(i for i, line in enumerate(lines)
                 if '= radius(micron)' in line and np.isclose(float(line.split()[0]), .1))
    data = np.array([[float(v) for v in line.split()] for line in lines[start+2:start+243]])
    data = data[np.argsort(data[:, 0])]
    assert data.shape == (241, 4) and np.all(data[:, 1] > 0)
    return np.exp(np.interp(np.log(wave), np.log(data[:, 0]), np.log(data[:, 1])))


def planck_nu(wave, temperature):
    nu = C/(np.asarray(wave)*1e-6)
    x = H*nu/(K*np.asarray(temperature)[..., None])
    return 2*H*nu**3/C**2/np.expm1(np.minimum(x, 700))


def setup_grains():
    # Equilibrium cooling integral in frequency, using the full table's IR range.
    wave = np.geomspace(.01, 1000, 6000)
    nu = C/(wave*1e-6)
    temps = np.geomspace(10, 6000, 1400)
    grains = []
    for name, rho, mass_fraction, edge_temperature in [
            ('Sil_21.gz', 3.5, .6, 1000), ('Gra_21.gz', 2.24, .4, 1500)]:
        q = efficiency(name, wave)
        cooling = -trapezoid(q*planck_nu(wave, temps), nu, axis=-1)
        q_heat = float(efficiency(name, np.array([.1]))[0])
        grains.append(dict(name=name, rho=rho, mass_fraction=mass_fraction,
                           q=efficiency(name, WAVE), q_heat=q_heat,
                           log_t=np.log(temps), log_cool=np.log(cooling),
                           edge_temperature=edge_temperature))
    # Incident flux is monochromatic at 0.1 um, fixed spectral shape, varying amplitude.
    sil = grains[0]
    heating = np.exp(np.interp(np.log(sil['edge_temperature']), sil['log_t'],
                               sil['log_cool']))/sil['q_heat']
    for g in grains:
        edge_cooling = np.exp(np.interp(np.log(g['edge_temperature']),
                                       g['log_t'], g['log_cool']))
        g['inner_radius'] = np.sqrt(heating*g['q_heat']/edge_cooling)
    return grains, heating


GRAINS, HEATING = setup_grains()


def spectra(luminosity_ratio, age=20., n_shells=2000):
    """Age and radius in initial silicate inner-radius light-crossing units.

    For isotropic thin shells, observer-frame echo delays are uniform on [0,2r/c].
    The evolved boundary changes immediately when the new illumination arrives;
    this limiting case supplies no measured dust-formation/destruction timescale.
    Shell dust mass per log radius is proportional to r, from the inner edge to 100.
    """
    r = np.geomspace(.1, 100., n_shells)
    dlogr = np.log(r[1]/r[0])
    fraction = np.clip(age/(2*r), 0, 1)
    fixed, evolved = np.zeros_like(WAVE), np.zeros_like(WAVE)
    errors = []
    for g in GRAINS:
        base_power = HEATING*g['q_heat']/r**2
        power_new = base_power*luminosity_ratio
        t_old = np.exp(np.interp(np.log(base_power), g['log_cool'], g['log_t']))
        t_new = np.exp(np.interp(np.log(power_new), g['log_cool'], g['log_t']))
        b_old = planck_nu(WAVE, t_old)
        b_new = planck_nu(WAVE, t_new)
        # Integrate boundary cells fractionally to reduce grid-edge errors.
        old = np.clip((np.log(r)+dlogr/2-np.log(g['inner_radius']))/dlogr, 0, 1)
        new = np.clip((np.log(r)+dlogr/2-np.log(g['inner_radius']*np.sqrt(luminosity_ratio)))/dlogr, 0, 1)
        weight = r*dlogr*g['mass_fraction']/g['rho']
        f = old[:, None]*((1-fraction[:, None])*b_old + fraction[:, None]*b_new)
        e = (old*(1-fraction))[:, None]*b_old + (new*fraction)[:, None]*b_new
        fixed += g['q']*np.sum(weight[:, None]*f, axis=0)
        evolved += g['q']*np.sum(weight[:, None]*e, axis=0)
        use = (old > 0) & (fraction > 0)
        reconstructed = np.exp(np.interp(np.log(t_new[use]), g['log_t'], g['log_cool']))
        errors.append(float(np.max(np.abs(reconstructed/power_new[use]-1))) if np.any(use) else 0.)
    return fixed, evolved, max(errors)


def match_hot_luminosity(fixed, evolved):
    """Fit one nuisance amplitude: equal integrated rest 2--4 um luminosities.

    This controls for the total hot-dust flux; it is not a fit to observed spectra.
    All grain temperatures, opacities, and relative radial weights stay unchanged.
    """
    use = (WAVE >= 2.) & (WAVE <= 4.)
    nu = C/(WAVE[use]*1e-6)
    scale = trapezoid(fixed[use], nu)/trapezoid(evolved[use], nu)
    return fixed, evolved*scale, float(scale)


def draw_black_hole(ax, x, bright):
    """Schematic black hole and tilted accretion disk, deliberately not to scale."""
    disk = '#d98b32' if bright else '#a1917d'
    rim = '#f0b958' if bright else '#b8aa95'
    ax.add_patch(Ellipse((x, 0), .17, .077, angle=-25, facecolor=disk,
                         edgecolor=rim, lw=.6, zorder=5))
    ax.add_patch(Ellipse((x, 0), .12, .040, angle=-25, facecolor='#312922',
                         edgecolor='none', zorder=6))
    ax.add_patch(Circle((x, 0), .030, facecolor='#101216',
                        edgecolor='#101216', lw=.4, zorder=7))


def draw_cartoon(ax, ratio):
    """Vector cross-sections illustrating boundaries, not a simulated dust image."""
    ax.set(xlim=(-1.38, 1.38), ylim=(-.60, .80), aspect='equal')
    ax.axis('off')
    original, outer = .2, .52
    for x, moving, color in [(-.73, False, BLUE), (.73, True, ORANGE)]:
        inner = original*np.sqrt(ratio) if moving else original
        ax.add_patch(Wedge((x, 0), outer, 0, 360, width=outer-inner,
                           facecolor=color, alpha=.09, edgecolor='none'))
        # A deterministic texture indicates dust; dot counts carry no mass information.
        rng = np.random.default_rng(8)
        radii = outer*np.sqrt(rng.uniform(0, 1, 110))
        theta = rng.uniform(0, 2*np.pi, 110)
        keep = radii > inner
        ax.scatter(x+radii[keep]*np.cos(theta[keep]), radii[keep]*np.sin(theta[keep]),
                   s=3, color=color, alpha=.4, linewidths=0)
        ax.add_patch(Circle((x, 0), outer, fill=False, edgecolor=color, lw=.65, alpha=.35))
        ax.add_patch(Circle((x, 0), original, fill=False, edgecolor=BLUE, lw=1.25,
                            alpha=.65 if moving else 1))
        if moving:
            ax.add_patch(Circle((x, 0), inner, fill=False, edgecolor=ORANGE, lw=1.7,
                                linestyle=(0, (3, 2))))
            # Outward after a rise; inward after a fade, as prescribed in the model.
            direction = np.array([np.cos(np.pi/4), np.sin(np.pi/4)])
            start = np.array([x, 0])+original*direction
            end = np.array([x, 0])+inner*direction
            ax.annotate('', xy=end, xytext=start,
                        arrowprops=dict(arrowstyle='-|>', color=ORANGE, lw=1.25,
                                        mutation_scale=10, shrinkA=0, shrinkB=0))
        # The same central luminosity is used in both hypotheses of each column.
        draw_black_hole(ax, x, bright=ratio > 1)
        ax.text(x, .64, 'Evolving' if moving else 'Fixed', ha='center',
                color=color, fontsize=10.5)


def main():
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Times New Roman'],
                         'font.size': 11, 'axes.labelsize': 11, 'axes.titlesize': 12,
                         'xtick.labelsize': 10, 'ytick.labelsize': 10,
                         'pdf.fonttype': 42, 'axes.spines.top': False,
                         'axes.spines.right': False})
    fig = plt.figure(figsize=(6.5, 4.65))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.65], hspace=.18, wspace=.12,
                           left=.105, right=.99, bottom=.22, top=.90)
    calculations, checks, rows = [], [], []
    raw_models = [spectra(4.), spectra(.25)]
    models = [match_hot_luminosity(f, e) for f, e, _ in raw_models]
    for i, (ratio, title) in enumerate([(4., r'Rising: $L\,\times\,4$'), (.25, r'Fading: $L\,/\,4$')]):
        cartoon = fig.add_subplot(grid[0, i])
        draw_cartoon(cartoon, ratio)
        cartoon.set_title(title, fontsize=12, fontweight='bold', pad=6)
        residual = fig.add_subplot(grid[1, i])
        f, e, mass_scale = models[i]
        raw_fixed, raw_evolved, equilibrium_error = raw_models[i]
        norm = np.interp(12., WAVE, f)
        # Warm-continuum diagnostic and broad silicate bands, all at rest wavelength.
        residual.axvspan(2., 4., color='#dbe8ef', alpha=.5, lw=0, zorder=0)
        residual.axvspan(8., 13., color='#ddd4b6', alpha=.23, lw=0, zorder=0)
        for center in [9.7, 18.]:
            residual.axvline(center, color='#9b782d', lw=.8, ls=':', ymax=.87)
            residual.text(center, .98, f'{center:g} µm', transform=residual.get_xaxis_transform(),
                          ha='center', va='top', color='#87631c', fontsize=10)
        # PAH ticks identify contaminants to fit, not synthetic PAH detections.
        for w in [6.2, 7.7, 11.3]:
            residual.plot([w, w], [.02, .055], transform=residual.get_xaxis_transform(),
                    color='#54806a', lw=1.1)
        log_ratio = np.log10(e/f)
        percent_difference = 100*(e/f-1)
        residual.axhline(0, color=BLUE, lw=.8)
        residual.plot(WAVE, percent_difference, color=ORANGE, lw=1.8, ls='--')
        residual.fill_between(WAVE, 0, percent_difference, color=ORANGE, alpha=.17, linewidth=0)
        residual.set(xlim=(2, 24.5), ylim=(-48, 67), yticks=[-40, -20, 0, 20, 40], xticks=[3, 5, 10, 15, 20])
        contrast18 = float(np.interp(18., WAVE, percent_difference))
        residual.plot(18., contrast18, 'o', color=ORANGE, ms=3)
        residual.annotate(f'{contrast18:+.0f}%', xy=(18., contrast18),
                          xytext=(14., 48 if i == 0 else -42), fontsize=12, fontweight='bold',
                          color=ORANGE, arrowprops=dict(arrowstyle='-', lw=.7, color=ORANGE))
        if i == 1:
            residual.annotate('Overlap', xy=(9.7, float(np.interp(9.7, WAVE, percent_difference))),
                              xytext=(6., 24), fontsize=9, color='#555555',
                              arrowprops=dict(arrowstyle='-', lw=.7, color='#555555'))
        if i == 0:
            residual.set_ylabel(r'Difference / fixed (%)')
        else:
            residual.tick_params(labelleft=False)
        f_hi, e_hi, _ = spectra(ratio, n_shells=4000)
        convergence = float(max(np.max(np.abs(f_hi/raw_fixed-1)), np.max(np.abs(e_hi/raw_evolved-1))))
        checks.append(dict(luminosity_ratio=ratio, maximum_equilibrium_relative_error=equilibrium_error,
                           maximum_relative_change_doubling_radial_grid=convergence))
        calculations.append(dict(luminosity_ratio=ratio, label=title,
                                 evolved_dust_mass_amplitude_scale=mass_scale,
                                 percent_difference_at_rest_18um=contrast18,
                                 normalization='same scale for both curves: fixed Fnu at 12um = 1'))
        rows.extend(zip([ratio]*len(WAVE), WAVE, f/norm, e/norm, log_ratio))
    fig.text(.55, .125, r'Rest wavelength ($\mu$m)', ha='center', fontsize=11)
    fig.legend(handles=[
        Patch(facecolor='#dbe8ef', alpha=.5, label='Hot-dust anchor'),
        Patch(facecolor='#ddd4b6', alpha=.23, label='Warm-dust region'),
        Patch(facecolor=ORANGE, alpha=.17, label='Model difference'),
    ], loc='lower center', bbox_to_anchor=(.55, .025), ncol=3,
        frameon=False, fontsize=10, handlelength=1.35, handleheight=1.,
        handletextpad=.45, columnspacing=1.2)
    fig.savefig(HERE/'fig2_diagnostics.pdf')
    fig.savefig(HERE/'fig2_diagnostics.png', dpi=220)
    plt.close(fig)
    with (INPUTS/'fig2_diagnostics_spectra.csv').open('w') as handle:
        writer = csv.writer(handle)
        writer.writerow(['luminosity_ratio', 'rest_um', 'fixed_Fnu_relative',
                         'evolving_Fnu_relative', 'log10_evolving_over_fixed'])
        writer.writerows(rows)
    nochange = spectra(1.)
    before = spectra(4., age=0.)
    assert np.allclose(nochange[0], nochange[1], rtol=1e-12, atol=0)
    assert np.allclose(before[0], before[1], rtol=1e-12, atol=0)
    assert all(c['maximum_relative_change_doubling_radial_grid'] < .01 for c in checks)
    provenance = dict(
        status='Illustrative calculation; not fitted to any observed target or SPHEREx spectrum.',
        model='Optically thin, isotropic shells; no self-absorption, host emission, or stochastic heating.',
        heating='Monochromatic 0.1 um radiation; same spectral shape and step history in each hypothesis.',
        grains='0.1 um spheres; 60/40 silicate/graphite mass, densities 3.5/2.24 g cm^-3.',
        grain_temperatures='Separate equilibrium integrals for each species; inner-edge initial T=1000/1500 K.',
        distribution='dM/dlnr proportional to r; fixed outer radius 100 initial silicate inner radii.',
        echo='Uniform shell delay from 0 to 2r/c, evaluated at t=20 initial silicate inner-radius/c.',
        alternative='Each species inner boundary scales with sqrt(L/L_initial), immediately after illumination arrives.',
        amplitude_fit='Each evolving spectrum is rescaled to match the fixed model integrated 2--4 um luminosity. This represents a free dust mass amplitude; it does not fit colors or any actual data.',
        limitations=['Not a clumpy or disk-wind radiative transfer calculation.',
                     'No fit to actual light curves or near-IR spectra; cannot establish indistinguishable existing-data fits.',
                     'No posterior predictive uncertainty or population power calculation.',
                     'Residual signs are specific to this boundary prescription, not universal rising/fading signatures.',
                     'PAH positions marked only; no invented PAH fluxes.',
                     'Cartoons show boundary motion only; dust dot counts, outer radii and brightness are schematic.',
                     'No ETC error bars are plotted on this figure.'],
        opacities={name:dict(url='https://www.astro.princeton.edu/~draine/dust/diel/'+name,
                           sha256=hashlib.sha256((INPUTS/'dust_opacity'/name).read_bytes()).hexdigest())
                   for name in ['Sil_21.gz', 'Gra_21.gz']},
        opacity_reference='Laor & Draine 1993, ApJ 402, 441; Draine & Lee 1984, ApJ 285, 89.',
        wavelength_axis='Rest-frame model wavelengths; no target-specific coverage claimed.',
        actual_targets_fitted=[], calculations=calculations,
        numerical_checks=checks, no_change_and_pre_echo_equal=True,
        presentation='Compact cartoons and fractional differences, with panel labels only; explanation in proposal caption.')
    (INPUTS/'fig2_diagnostics_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print(json.dumps(dict(checks=checks), indent=2))


if __name__ == '__main__':
    main()
