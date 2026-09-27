"""Exploratory, reproducible P9694 SPHEREx continuum decomposition.

Run with python jwst_proposal/fit_spherex_hot_dust.py. Fits use supplied errors,
nonnegative component amplitudes, a constant host and a separate disc/dust
amplitude and dust temperature per visit. See output report for limitations.
"""
from pathlib import Path
import csv
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares, nnls
from scipy.stats import chi2

from spherex_data import CSV_PATHS, load_spherex

HERE = Path(__file__).resolve().parent
OUT = HERE / 'review' / 'p9694_hot_dust'
Z = 0.2377224
TBOUNDS = (600., 2400.)
# Fixed, flux-independent rest-frame windows; conservative broad-line masks.
LINES = np.array([.656461, 1.0052, 1.0833, 1.0941, 1.2822,
                  1.8756, 2.1661, 2.6259, 4.0523])


def blackbody(wave, temperature, beta=0):
    """F_nu shape normalized to unity at rest 2 microns."""
    return ((2 / wave)**(3 + beta) * np.expm1(14387.76877 / (2 * temperature))
            / np.expm1(14387.76877 / (wave * temperature)))


class Fit:
    def __init__(self, template='Ell5', alpha=1/3, beta=0, quadrature=32,
                 mask_velocity=.02):
        self.data, self.groups, _ = load_spherex('P9694')
        self.alpha, self.beta, self.template = alpha, beta, template
        d = self.data
        lo = (d['wavelength_um'] - d['wavelength_half_width_um']) / (1 + Z)
        hi = (d['wavelength_um'] + d['wavelength_half_width_um']) / (1 + Z)
        self.keep = np.ones(len(lo), dtype=bool)
        self.reasons = [[] for _ in lo]
        windows = [(l*(1-mask_velocity), l*(1+mask_velocity), f'line_{l:g}')
                   for l in LINES] + [(3.20, 3.35, 'PAH_3.3')]
        for a, b, label in windows:
            overlap = (lo <= b) & (hi >= a)
            self.keep &= ~overlap
            for i in np.flatnonzero(overlap):
                self.reasons[i].append(label)
        self.epoch = np.empty(len(lo), dtype=int)
        for e, g in enumerate(self.groups):
            self.epoch[g['indices']] = e
        nodes, weights = np.polynomial.legendre.leggauss(quadrature)
        self.waves = (lo[:, None]+hi[:, None])/2 + (hi-lo)[:, None]*nodes/2
        self.weights = weights/2
        path = HERE / 'inputs' / 'host_templates' / f'{template}_template_norm.sed'
        tab = np.loadtxt(path)
        self.hostwave = tab[:, 0] / 1e4
        # SWIRE supplies F_lambda, NOT F_nu. Overall normalization is free.
        self.hostflux = tab[:, 1] * self.hostwave**2
        self.hostflux /= np.interp(1.6, self.hostwave, self.hostflux)
        self.host = self.host_at(self.waves) @ self.weights
        self.y = d['flux_mjy'][self.keep]
        self.err = d['flux_err_mjy'][self.keep]

    def host_at(self, wave):
        return np.interp(wave, self.hostwave, self.hostflux)

    def design(self, nonlinear):
        temps = nonlinear[:3]
        alpha = nonlinear[3] if self.alpha is None else self.alpha
        disk = self.waves**(-alpha) @ self.weights
        matrix = np.zeros((len(self.keep), 7))
        matrix[:, 0] = self.host
        for e in range(3):
            take = self.epoch == e
            matrix[take, 1+2*e] = disk[take]
            matrix[take, 2+2*e] = (blackbody(self.waves[take], temps[e], self.beta)
                                       @ self.weights)
        return matrix

    def solve_linear(self, nonlinear, y=None):
        a = self.design(nonlinear)[self.keep]
        target = self.y if y is None else y
        amps, _ = nnls(a/self.err[:, None], target/self.err)
        return (a @ amps-target)/self.err, amps

    def optimize(self, fixed=None, start=None, y=None):
        bounds = [TBOUNDS]*3 + ([(-1.5, 1.)] if self.alpha is None else [])
        fixed = fixed or {}
        free = [j for j in range(len(bounds)) if j not in fixed]
        def expand(x):
            p = np.empty(len(bounds))
            for j, v in fixed.items():
                p[j] = v
            p[free] = x
            return p
        def residual(x):
            return self.solve_linear(expand(x), y=y)[0]
        starts = [start] if start is not None else [
            [t]*3 + ([a] if self.alpha is None else [])
            for t in (850., 1400., 2000.)
            for a in ((-.8, 1/3) if self.alpha is None else (1/3,))]
        fits = []
        for initial in starts:
            r = least_squares(residual, np.array(initial)[free],
                              bounds=np.array(bounds)[free].T, x_scale='jac',
                              ftol=1e-9, xtol=1e-9, gtol=1e-9, max_nfev=400)
            fits.append(r)
        r = min(fits, key=lambda r: np.sum(r.fun**2))
        p = expand(r.x)
        residuals, amps = self.solve_linear(p, y=y)
        npar = len(bounds)+7
        statistic = float(residuals @ residuals)
        return dict(parameters=p.tolist(), amplitudes_mjy=amps.tolist(),
                    chi2=statistic, dof=len(self.y)-npar,
                    reduced_chi2=statistic/(len(self.y)-npar),
                    p_value=float(chi2.sf(statistic, len(self.y)-npar)),
                    converged=bool(r.success), template=self.template,
                    alpha=float(p[3] if self.alpha is None else self.alpha),
                    beta=self.beta, n_used=len(self.y),
                    n_per_visit=[int(np.sum(self.keep & (self.epoch == e))) for e in range(3)])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    variants = [('baseline', dict()), ('free_disc_slope', dict(alpha=None)),
                ('young_host', dict(template='Ell2')),
                ('S0_host', dict(template='S0')),
                ('modified_blackbody_beta1', dict(beta=1)),
                ('wider_line_masks', dict(mask_velocity=.03))]
    models, results = {}, {}
    for name, kwargs in variants:
        model = Fit(**kwargs)
        result = model.optimize()
        models[name], results[name] = model, result
        print(name, json.dumps(result), flush=True)
    model, best = models['baseline'], results['baseline']
    # Numerical checks: band integration convergence and noiseless recovery.
    fine = Fit(quadrature=64)
    pred = model.design(best['parameters']) @ best['amplitudes_mjy']
    pred_fine = fine.design(best['parameters']) @ best['amplitudes_mjy']
    fractional_quadrature_change = float(np.max(abs(pred_fine-pred)/pred_fine))
    synthetic_p = [1100., 1450., 1700.]
    synthetic_a = [1.5, .8, 2., 1.2, 3., 1., 2.5]
    synthetic_y = model.design(synthetic_p)[model.keep] @ synthetic_a
    recovered = model.optimize(y=synthetic_y)
    recovery_error = float(np.max(abs(np.array(recovered['parameters'])-synthetic_p)))
    assert recovery_error < 1., recovered
    assert fractional_quadrature_change < .005, fractional_quadrature_change
    # Profile intervals conditional on baseline, independent Gaussian errors.
    profiles = []
    for e in range(3):
        grid = np.unique(np.r_[np.linspace(*TBOUNDS, 121),
                               np.arange(best['parameters'][e]-60, best['parameters'][e]+61, 1)])
        grid = grid[(grid >= TBOUNDS[0]) & (grid <= TBOUNDS[1])]
        values = [model.optimize(fixed={e: float(t)}, start=best['parameters'])['chi2']
                  for t in grid]
        delta = np.array(values)-best['chi2']
        def interval(threshold):
            # Interpolate threshold crossings surrounding the global minimum.
            imin = np.argmin(delta)
            left, right = imin, imin
            while left > 0 and delta[left-1] <= threshold:
                left -= 1
            while right < len(grid)-1 and delta[right+1] <= threshold:
                right += 1
            low = grid[0] if left == 0 else np.interp(
                threshold, [delta[left], delta[left-1]], [grid[left], grid[left-1]])
            high = grid[-1] if right == len(grid)-1 else np.interp(
                threshold, [delta[right], delta[right+1]], [grid[right], grid[right+1]])
            return [float(low), float(high)]
        profiles.append(dict(epoch=model.groups[e]['label'], temperature_K=best['parameters'][e],
                             conditional_68pct_K=interval(1), conditional_95pct_K=interval(3.841459),
                             grid_K=grid.tolist(), delta_chi2=delta.tolist()))
    results['profiles'] = profiles
    results['constant_temperature_comparison'] = {}
    for name in ('baseline', 'free_disc_slope'):
        m, b = models[name], results[name]
        def expand(p):
            return [p[0]]*3 + (list(p[1:]) if m.alpha is None else [])
        trial = least_squares(lambda p: m.solve_linear(expand(p))[0],
                              [np.mean(b['parameters'][:3])] + ([b['alpha']] if m.alpha is None else []),
                              bounds=np.array([TBOUNDS] + ([(-1.5, 1.)] if m.alpha is None else [])).T,
                              x_scale='jac', ftol=1e-10, xtol=1e-10, gtol=1e-10)
        delta = float(trial.fun @ trial.fun)-b['chi2']
        results['constant_temperature_comparison'][name] = dict(
            shared_temperature_K=float(trial.x[0]), chi2=float(trial.fun @ trial.fun),
            delta_chi2_for_two_extra_temperature_parameters=delta,
            nominal_p_value=float(chi2.sf(delta, 2)), converged=bool(trial.success),
            caveat='Nominal likelihood comparison only: both models have poor absolute goodness of fit; '
                   'covariance and model error are not included. Not a robust variability detection.')
    results['numerical_checks'] = dict(max_band_quadrature_fractional_change=fractional_quadrature_change,
                                       noiseless_temperature_recovery_max_error_K=recovery_error)
    results['redshift'] = Z
    results['input_sha256'] = hashlib.sha256(CSV_PATHS['P9694'].read_bytes()).hexdigest()
    results['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    results['amplitude_order'] = ['shared_host_at_rest_1.6um',
                                 'visit1_disc_at_rest_1um', 'visit1_dust_at_rest_2um',
                                 'visit2_disc_at_rest_1um', 'visit2_dust_at_rest_2um',
                                 'visit3_disc_at_rest_1um', 'visit3_dust_at_rest_2um']
    results['assumptions'] = dict(
        error_model='Supplied total errors, diagonal covariance; no error inflation',
        response='Uniform wavelength top-hat from exported centre +/- half-width; 32-point quadrature',
        extinction='No foreground or intrinsic extinction correction applied',
        normalization='Observer-frame mJy at the named rest wavelengths; no luminosity-distance factor needed',
        temperature_bounds_K=list(TBOUNDS),
        model='Fnu = shared host + per-visit A_disc * lambda_rest^(-alpha) + '
              'per-visit A_dust * Bnu(lambda_rest,T) * (lambda_rest/2um)^(-beta), '
              'with dust normalized at rest 2um',
        use_for_proposal_temperature_labels=False)
    (OUT / 'results.json').write_text(json.dumps(results, indent=2)+'\n')
    write_report(results, model.groups)
    with (OUT / 'measurements.csv').open('w') as f:
        writer = csv.writer(f)
        writer.writerow(['input_row_zero_based', 'epoch', 'wavelength_um', 'flux_mjy',
                         'error_mjy', 'used_in_fit', 'mask_reason', 'baseline_model_mjy',
                         'residual_sigma'])
        for i in range(len(pred)):
            d = model.data
            writer.writerow([i, model.groups[model.epoch[i]]['label'], d['wavelength_um'][i],
                             d['flux_mjy'][i], d['flux_err_mjy'][i], bool(model.keep[i]),
                             ';'.join(model.reasons[i]), pred[i],
                             (d['flux_mjy'][i]-pred[i])/d['flux_err_mjy'][i]])
    plot_fit(model, best, 'continuum_fits')
    plot_fit(models['free_disc_slope'], results['free_disc_slope'], 'free_slope_fits')
    print('Wrote', OUT, flush=True)


def write_report(results, groups):
    names = [('baseline', 'Ell5 host, fixed disc slope 1/3'),
             ('free_disc_slope', 'Ell5 host, free common power-law slope'),
             ('young_host', 'Ell2 host, fixed disc slope 1/3'),
             ('S0_host', 'S0 host, fixed disc slope 1/3'),
             ('modified_blackbody_beta1', 'Ell5 host, modified blackbody β=1'),
             ('wider_line_masks', 'Ell5 host, wider line masks')]
    table = '\n'.join('| ' + label + ' | ' + ' | '.join(f'{t:.0f}' for t in results[key]['parameters'][:3])
                      + f" | {results[key]['reduced_chi2']:.2f} |" for key, label in names)
    base = results['baseline']
    intervals = '\n'.join(f"| {p['epoch']} | {p['temperature_K']:.0f} | "
                          f"{p['conditional_68pct_K'][0]:.0f}–{p['conditional_68pct_K'][1]:.0f} |"
                          for p in results['profiles'])
    comparisons = '\n'.join(
        f"- {key}: shared T = {r['shared_temperature_K']:.0f} K; allowing separate visit "
        f"temperatures improves χ² by {r['delta_chi2_for_two_extra_temperature_parameters']:.2f} "
        f"for two extra parameters (nominal p = {r['nominal_p_value']:.3g})."
        for key, r in results['constant_temperature_comparison'].items())
    report = f"""# P9694 exploratory hot-dust continuum fits

The fits give an approximate hot-dust colour temperature, but **do not support
precise temperature labels in the proposal or a robust temperature-variability
claim**. The fixed-disc baseline gives about 1160–1180 K; a freely fitted red
power law gives about 1330–1390 K. Both leave structured residuals larger than
the supplied errors. Figure 1 and the proposal PDF were therefore not changed.

| Model | Jul 2025 T (K) | Dec 2025–Jan 2026 T (K) | Jun–Jul 2026 T (K) | χ² / nominal dof |
| --- | ---: | ---: | ---: | ---: |
{table}

These alternatives are sensitivity checks, **not a confidence interval**. In
particular β=1 changes the assumed emissivity as well as the inferred temperature
and provides a worse fit. The free power-law index is α =
{results['free_disc_slope']['alpha']:.3f}, where Fν ∝ ν^α. This red slope is a
phenomenological continuum; it is not evidence for a physical thin-disc spectrum.
The shared host amplitude falls from {base['amplitudes_mjy'][0]:.2f} mJy in the
baseline to {results['free_disc_slope']['amplitudes_mjy'][0]:.2f} mJy with the free
slope (both at rest 1.6 µm), demonstrating a substantial decomposition degeneracy.
Dust amplitudes and any extrapolated dust luminosities are correspondingly
model dependent; no physical bolometric luminosity is reported.

## Inputs and method

- Source P9694, J004607.98+090720.9, z = {Z}; 289 native, unbinned measurements.
  Rest wavelength coverage is approximately 0.60–4.03 µm.
- Simultaneous fit to three visits: constant host normalization and common disc
  slope; independent disc amplitude, dust amplitude and dust temperature per visit.
  All component amplitudes are nonnegative. Temperatures are bounded at 600–2400 K.
- Baseline host: SWIRE Ell5; Ell2 and S0 are sensitivity alternatives. The library
  Fλ spectra are converted to Fν before fitting. Ell13 was downloaded but not used.
- Each predicted Fν is averaged uniformly in wavelength over the exported centre
  ± half-width, using 32-point Gauss–Legendre quadrature. This is a top-hat
  approximation: the full instrument response curves were not supplied.
- Fixed masks exclude any band touching ±2% in wavelength around Hα, Paδ,
  He I/Paγ, Paβ, Paα, Brγ, Brβ or Brα, and rest 3.20–3.35 µm around PAH 3.3.
  The wider-mask check uses ±3%. No residual-based clipping is applied.
  The baseline uses {base['n_used']} points ({', '.join(map(str, base['n_per_visit']))}
  per visit); all supplied points remain visible in the diagnostic plots.
- The supplied total flux errors are used without rescaling. Their 2% calibration
  term is retained, but covariance is unavailable and the fit treats errors as
  independent. Nominal χ² probabilities and profile intervals are conditional on
  this approximation; active zero-amplitude boundaries further limit their interpretation.
- No foreground or intrinsic extinction correction is applied. Each visit spans
  roughly two to three weeks; within-visit source variability is not modelled.
  Noncontemporaneous optical/NGPS spectra are not used as rigid flux constraints.
- Model components are evaluated in the rest frame, with amplitudes in observer
  mJy; the common redshift/distance normalization is absorbed in the amplitudes.

## Formal errors and visit comparison

For reproducibility, the baseline profile-likelihood intervals (Δχ²=1, all other
parameters refitted) are recorded below. **These are not credible total errors**:
model dependence and structured residuals dominate this apparent precision.

| Visit | Baseline T (K) | Conditional 68% interval (K) |
| --- | ---: | ---: |
{intervals}

{comparisons}

These nominal comparisons assume an adequate spectral model and independent
errors, which the residuals do not establish. They must not be presented as a
robust detection of heating or cooling. An extinction-aware continuum model,
better host constraints, instrumental response/covariance and a check of
within-visit variability are needed before making that claim.

## Files and verification

- [Baseline components and residuals](continuum_fits.pdf)
- [Free-slope components and residuals](free_slope_fits.pdf)
- [Results, profiles, assumptions and input/script hashes](results.json)
- [Every measurement, mask decision and baseline residual](measurements.csv)

Reproduce from the project root with
`python jwst_proposal/fit_spherex_hot_dust.py`.
Noiseless synthetic data recover their input temperatures with maximum error
{results['numerical_checks']['noiseless_temperature_recovery_max_error_K']:.3g} K.
Increasing the quadrature to 64 nodes changes predicted band fluxes by at most
{100*results['numerical_checks']['max_band_quadrature_fractional_change']:.3f}%.
Multistart optimization was used for all six model variants. The diagnostics
show that the main limitation is model adequacy, not numerical convergence.

Host templates and their units: [SWIRE library, Polletta et al. (2007)](https://www.iasf-milano.inaf.it/~polletta/templates/swire_templates.html).
The use of near-IR spectral decomposition to infer dust colour temperature, and
its dependence on emissivity assumptions, has precedent in
[Landt et al. (2019)](https://arxiv.org/abs/1908.01627).
Template download provenance and checksums are in `../../inputs/host_templates/provenance.json`.
"""
    (OUT / 'README.md').write_text(report)


def plot_fit(model, best, filename):
    pred = model.design(best['parameters']) @ best['amplitudes_mjy']
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.2), sharex='col',
                             gridspec_kw={'height_ratios': [3, 1]}, layout='constrained')
    wave = np.geomspace(.72, 5.05, 600)
    rest = wave/(1+Z)
    a = best['amplitudes_mjy']
    for e, group in enumerate(model.groups):
        ax, residual_ax = axes[:, e]
        idx = group['indices']
        used, masked = idx[model.keep[idx]], idx[~model.keep[idx]]
        for rows, color, alpha in [(used, group['color'], .8), (masked, '.7', .45)]:
            ax.errorbar(model.data['wavelength_um'][rows], model.data['flux_mjy'][rows],
                        yerr=model.data['flux_err_mjy'][rows], fmt='.', ms=4, lw=.6,
                        color=color, alpha=alpha)
        host = a[0]*model.host_at(rest)
        disc = a[1+2*e]*rest**(-best['alpha'])
        dust = a[2+2*e]*blackbody(rest, best['parameters'][e])
        for flux, color, label, style in [(host, '.5', 'Host', ':'),
                                         (disc, '#bc772e', 'Disc', '--'),
                                         (dust, '#b3493e', 'Hot dust', '--'),
                                         (host+disc+dust, '.15', 'Total', '-')]:
            ax.plot(wave, flux, color=color, lw=1.3, ls=style, label=label)
        ax.set(yscale='log', ylim=(.15, 20), xlim=(.72, 5.05),
               title=f"{group['label']}\nConditional colour T = {best['parameters'][e]:.0f} K")
        residual_ax.axhline(0, color='.4', lw=.6)
        residual_ax.axhspan(-2, 2, color='.8', alpha=.3)
        residual_ax.plot(model.data['wavelength_um'][used],
                         (model.data['flux_mjy'][used]-pred[used])/model.data['flux_err_mjy'][used],
                         '.', color=group['color'], ms=4)
        residual_ax.set_xlabel('Observed wavelength (µm)')
    axes[0, 0].set_ylabel('Flux density (mJy)')
    axes[1, 0].set_ylabel('Residual / σ')
    axes[0, 0].legend(fontsize=8, ncol=2)
    slope_label = 'Fixed disc slope 1/3' if model.alpha is not None else f"Free power-law slope {best['alpha']:.2f}"
    fig.suptitle('P9694: exploratory host + power law + single-temperature blackbody\n'
                 f"{slope_label}; shared Ell5 host; χ² / dof = {best['chi2']:.1f} / {best['dof']}"
                 '   (grey points excluded from fit)', fontsize=11)
    for suffix in ('png', 'pdf'):
        fig.savefig(OUT / f'{filename}.{suffix}', dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
