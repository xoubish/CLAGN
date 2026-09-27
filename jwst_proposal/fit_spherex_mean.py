"""Fit one error-weighted continuum to all P9694 SPHEREx visits.

Retains the earlier line masks and bandpass approximation. The curve is a
year-averaged illustration, not an instantaneous spectrum or a variability fit.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, nnls
from fit_spherex_hot_dust import Fit, Z, TBOUNDS, blackbody
from spherex_data import CSV_PATHS

HERE = Path(__file__).resolve().parent
OUT = HERE / 'review' / 'p9694_hot_dust'


class MeanFit(Fit):
    def design(self, p):
        expanded = [p[0]]*3 + ([p[1]] if self.alpha is None else [])
        separate = super().design(expanded)
        return np.column_stack([separate[:, 0], separate[:, [1, 3, 5]].sum(axis=1),
                                separate[:, [2, 4, 6]].sum(axis=1)])

    def optimize_mean(self, y=None):
        bounds = [TBOUNDS] + ([(-1.5, 1.)] if self.alpha is None else [])
        trials = []
        for t in (850., 1400., 2000.):
            for alpha in ((-.8, 1/3) if self.alpha is None else (1/3,)):
                start = [t] + ([alpha] if self.alpha is None else [])
                fit = least_squares(lambda p: self.solve_linear(p, y=y)[0], start,
                                    bounds=np.array(bounds).T, x_scale='jac',
                                    ftol=1e-10, xtol=1e-10, gtol=1e-10)
                trials.append(fit)
        fit = min(trials, key=lambda f: f.fun @ f.fun)
        residual, amps = self.solve_linear(fit.x, y=y)
        dof = len(self.y)-len(bounds)-3
        return dict(temperature_K=float(fit.x[0]), parameters=fit.x.tolist(),
                    alpha=float(fit.x[1] if self.alpha is None else self.alpha),
                    beta=self.beta, host_template=self.template,
                    amplitudes_mjy=amps.tolist(), amplitude_order=['host_rest1.6um', 'disc_rest1um', 'dust_rest2um'],
                    chi2=float(residual @ residual), dof=dof,
                    reduced_chi2=float(residual @ residual)/dof,
                    n_used=len(self.y), converged=bool(fit.success))


def main():
    variants = [('baseline', {}), ('free_disc_slope', {'alpha': None}),
                ('young_host', {'template': 'Ell2'}), ('S0_host', {'template': 'S0'}),
                ('modified_blackbody_beta1', {'beta': 1}),
                ('wider_line_masks', {'mask_velocity': .03})]
    results = {}
    for key, kwargs in variants:
        fit = MeanFit(**kwargs)
        results[key] = fit.optimize_mean()
        assert results[key]['converged']
        print(key, json.dumps(results[key]), flush=True)
    fit = MeanFit()
    best = results['baseline']
    amps = np.array(best['amplitudes_mjy'])
    prediction = fit.design(best['parameters']) @ amps
    fine = MeanFit(quadrature=64)
    fine_prediction = fine.design(best['parameters']) @ amps
    quadrature_error = float(np.max(abs(prediction-fine_prediction)/fine_prediction))
    synthetic = fit.design([1450.])[fit.keep] @ np.array([1.5, .8, 3.])
    recovery = fit.optimize_mean(y=synthetic)
    assert abs(recovery['temperature_K']-1450.) < 1
    assert quadrature_error < .005
    previous = json.loads((OUT / 'results.json').read_text())
    assert previous['input_sha256'] == hashlib.sha256(CSV_PATHS['P9694'].read_bytes()).hexdigest()
    for key in ('baseline', 'free_disc_slope'):
        results[key]['separate_visit_chi2'] = previous[key]['chi2']
        results[key]['delta_chi2_vs_separate_visits'] = results[key]['chi2']-previous[key]['chi2']
        results[key]['extra_parameters_in_separate_visits'] = 6
    observed = np.linspace(.70, 5.05, 1200)
    rest = observed/(1+Z)
    host = amps[0]*fit.host_at(rest)
    disc = amps[1]*rest**(-best['alpha'])
    dust = amps[2]*blackbody(rest, best['temperature_K'])
    curve_path = HERE/'inputs'/'P9694_spherex_mean_fit.csv'
    np.savetxt(curve_path, np.column_stack([observed, host+disc+dust, host, disc, dust]),
               delimiter=',', header='wavelength_observed_um,total_mjy,host_mjy,disc_mjy,dust_mjy', comments='')
    results['provenance'] = dict(
        input_sha256=hashlib.sha256(CSV_PATHS['P9694'].read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        model_script_sha256=hashlib.sha256((HERE/'fit_spherex_hot_dust.py').read_bytes()).hexdigest(),
        curve_path=str(curve_path.relative_to(HERE)),
        curve_sha256=hashlib.sha256(curve_path.read_bytes()).hexdigest(),
        redshift=Z, total_measurements=len(fit.keep), fit_row_indices=np.flatnonzero(fit.keep).tolist(),
        visits=[{k:g[k] for k in ('label', 'start_mjd', 'end_mjd')} for g in fit.groups],
        interpretation='Single inverse-variance-weighted continuum across July 2025–July 2026; '
                       'all amplitudes, slope and temperature shared. Line bands excluded; '
                       'native points remain unbinned and unscaled. Conditional colour temperature, '
                       'not a physical single-temperature dust distribution or instantaneous spectrum.',
        plot_model='baseline', plot_temperature_label_K=round(best['temperature_K']/100)*100,
        numerical_checks=dict(max_band_quadrature_fractional_change=quadrature_error,
                              noiseless_temperature_error_K=abs(recovery['temperature_K']-1450.)))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'mean_fit.json').write_text(json.dumps(results, indent=2)+'\n')
    table = '\n'.join(f"| {k} | {results[k]['temperature_K']:.0f} | {results[k]['reduced_chi2']:.2f} |"
                      for k, _ in variants)
    (OUT/'mean_fit.md').write_text(f"""# P9694 pooled SPHEREx continuum

The figure now shows one dashed **total host + disc + blackbody** curve fitted
to the continuum measurements from all three visits, with a rounded,
model-dependent hot-dust colour temperature of **1200 K**. The exact baseline
optimum is {best['temperature_K']:.1f} K. All three component amplitudes and the
temperature are shared; this is not an average of the earlier fitted temperatures.
The fixed-disc baseline uses Fν ∝ ν^(1/3) and the SWIRE Ell5 host template.

The visits span **July 2025–July 2026**, approximately one year. Pooling gives a
useful descriptive average; it does not establish simultaneity or remove
variability, host/disc degeneracy, extinction or calibration systematics.
All 289 measurements are displayed unchanged; 209 continuum points enter the
fit with their supplied errors and band widths. The earlier flux-independent
emission-line/PAH masks are retained. No binning or inter-visit rescaling is used.

| Pooled model | Colour T (K) | χ² / nominal dof |
| --- | ---: | ---: |
{table}

The baseline has χ² = {best['chi2']:.1f} for {best['dof']} nominal degrees of freedom,
compared with {best['separate_visit_chi2']:.1f} for the earlier model allowing
separate visit disc amplitudes, dust amplitudes and temperatures (six additional
parameters). These residuals exceed the supplied errors. The figure's 1200 K
label is consequently an illustrative, conditional colour temperature without
a formal precision claim. The alternative models in this table are sensitivity
checks, not a confidence interval. Pooling does not resolve the earlier model
dependence, so no temperature-change or bolometric luminosity claim is made.

Reproduce with `python jwst_proposal/fit_spherex_mean.py`, then regenerate Figure 1.
[Machine-readable fit and curve provenance](mean_fit.json).
The [earlier fit report](README.md) documents shared assumptions, line masks,
template sources, uncertainty limitations and the original separate-visit fits.
Numerical checks recover a known synthetic temperature to better than 1 K;
64-node quadrature changes predicted band fluxes by at most {100*quadrature_error:.3f}%.
""")


if __name__ == '__main__':
    main()
