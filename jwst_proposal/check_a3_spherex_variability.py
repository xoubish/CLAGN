"""Exploratory continuum variability check for a supplied SPHEREx CSV.

Fit each epoch independently at its own wavelength samples, then integrate
the fitted continuum over identical observed-wavelength intervals. This avoids
comparing unequal raw wavelength averages. No input data are modified.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import chi2
from spherex_line_labels import CSV_PATHS, load_spherex, retained_measurements

HERE = Path(__file__).resolve().parent
WINDOWS = [(2.6, 3.5), (3.5, 4.0), (4.0, 4.5), (4.5, 5.0), (4.0, 5.0)]


def solve(design, flux, covariance):
    weighted = cho_solve(cho_factor(covariance), design)
    parameter_covariance = np.linalg.inv(design.T @ weighted)
    return parameter_covariance @ weighted.T @ flux, parameter_covariance


def summarize(flux, covariance):
    contrasts = []
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        change = flux[j] - flux[i]
        error = np.sqrt(covariance[i, i] + covariance[j, j] - 2 * covariance[i, j])
        contrasts.append(dict(epochs=[i, j], change_mjy=float(change),
                              error_mjy=float(error),
                              change_percent=float(100 * change / flux[i]),
                              signed_sigma=float(change / error)))
    weight = np.linalg.inv(covariance)
    common = np.sum(weight @ flux) / np.sum(weight)
    residual = flux - common
    statistic = float(residual @ weight @ residual)
    return dict(mean_flux_mjy=flux.tolist(), covariance_mjy2=covariance.tolist(),
                contrasts=contrasts, constant_flux_chi2=statistic,
                constant_flux_dof=2, constant_flux_p=float(chi2.sf(statistic, 2)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target', choices=['A3', 'P2190', 'P1823'], default='A3')
    parser.add_argument('--include-excluded', action='store_true',
                        help='Sensitivity check including the user-excluded P1823 point.')
    args = parser.parse_args()
    prefix = args.target
    output_name = prefix + ('_including_outlier' if args.include_excluded else '')
    source = CSV_PATHS[prefix]
    data, groups, rows = load_spherex(prefix)
    assert len(groups) == 3, 'This comparison requires exactly three observing periods.'
    wave, flux, errors = (data[k] for k in ['wavelength_um', 'flux_mjy', 'flux_err_mjy'])
    retained = np.ones(len(rows), bool) if args.include_excluded else retained_measurements(prefix, data)
    epoch = np.empty(len(rows), dtype=int)
    for i, group in enumerate(groups):
        epoch[group['indices']] = i
    measurement_var = np.array([float(r['variance_measurement_mjy2']) for r in rows])
    calibration_var = np.array([float(r['variance_calibration_mjy2']) for r in rows])
    # CSV errors/variances are rounded independently to 5-6 significant digits.
    np.testing.assert_allclose(measurement_var + calibration_var, errors**2, rtol=2e-4)
    cal_fraction = np.median(np.sqrt(calibration_var) / np.abs(flux))
    assert np.isclose(cal_fraction, .02, atol=1e-5)
    results = []
    windows = ([(.85, 1.), (1.15, 1.5), (2.2, 2.6)] + WINDOWS
               if prefix == 'P1823' else WINDOWS)
    for lower, upper in windows:
        selected = (wave >= lower) & (wave <= upper) & retained
        x = (wave[selected] - (lower + upper) / 2) / (upper - lower)
        y, err, ids = flux[selected], errors[selected], epoch[selected]
        for degree in (1, 2, 3):
            design = np.column_stack([(ids == i) * x**k
                                      for i in range(3) for k in range(degree + 1)])
            # Integral over x=-0.5..0.5, divided by the interval width (1).
            integral = np.zeros((3, design.shape[1]))
            for i in range(3):
                for k in range(degree + 1):
                    integral[i, i * (degree + 1) + k] = 0 if k % 2 else 1 / (2**k * (k + 1))
            b, cov = solve(design, y, np.diag(err**2))
            predicted = design @ b
            independent = summarize(integral @ b, integral @ cov @ integral.T)
            independent['fit_chi2'] = float(np.sum(((y - predicted) / err)**2))
            independent['fit_dof'] = len(y) - len(b)
            independent['fit_p'] = float(chi2.sf(independent['fit_chi2'], independent['fit_dof']))
            # Sensitivity scenario, not a measured covariance: move the CSV's
            # existing 2% term out of the diagonal into one shared scale per
            # epoch. Independent epochs; no double counting the 2% term.
            covariance = np.diag(measurement_var[selected])
            for i in range(3):
                calibration = cal_fraction * predicted * (ids == i)
                covariance += np.outer(calibration, calibration)
            bc, cc = solve(design, y, covariance)
            correlated = summarize(integral @ bc, integral @ cc @ integral.T)
            residual = y - design @ bc
            correlated['fit_chi2'] = float(residual @ cho_solve(cho_factor(covariance), residual))
            correlated['fit_dof'] = len(y) - len(b)
            # Sensitivity to excess point-to-point scatter; only inflate the
            # measurement variance, leaving the existing shared 2% term intact.
            inflation = max(1., correlated['fit_chi2'] / correlated['fit_dof'])
            inflated_covariance = covariance + (inflation - 1.) * np.diag(measurement_var[selected])
            bi, ci = solve(design, y, inflated_covariance)
            inflated = summarize(integral @ bi, integral @ ci @ integral.T)
            inflated['measurement_variance_inflation'] = inflation
            results.append(dict(observed_window_um=[lower, upper], degree=degree,
                                counts=np.bincount(ids, minlength=3).tolist(),
                                independent_errors=independent,
                                epoch_correlated_calibration=correlated,
                                excess_scatter_and_correlated_calibration=inflated))
    summary = dict(
        input=source.name, sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        epochs=[g['label'] for g in groups],
        excluded_csv_row_indices=np.flatnonzero(~retained).tolist(),
        method='Independent polynomial per epoch; uniform wavelength average of F_nu over common observed intervals. Quadratic adopted, linear/cubic sensitivity checks.',
        calibration_fraction_in_csv=float(cal_fraction),
        covariance_scenario='Independent measurement variances plus a 2% multiplicative term fully correlated within each epoch and independent between epochs; replaces the original diagonal calibration variance.',
        limitations=[
            'Exploratory, uncorrected local significances; aggregate 4-5 um overlaps two narrower windows and is not independent of them.',
            'Smooth fits use wavelength centres, not full spectral response curves. Band means are model estimates, not direct synthetic photometry.',
            'CSV does not specify cross-measurement or cross-epoch calibration covariance; the correlated scenario is illustrative, not established.',
            'Background, aperture, processing-version and other systematics have not been independently validated with control sources.',
            'Nominal significances assume the fitted smooth continuum is adequate; inspect fit chi-square and degree sensitivity.',
        ], results=results)
    (HERE / f'spherex_{output_name}_variability.json').write_text(json.dumps(summary, indent=2) + '\n')
    lines = [
        f'# {prefix} SPHEREx continuum variability check', '',
        f'Regenerate with `python check_a3_spherex_variability.py --target {prefix}`.', '',
        'The three epochs are ' + ', '.join(g['label'] for g in groups) + '. ' +
        'Fit an independent quadratic continuum to each epoch at its actual wavelengths, '
        'then average the fitted F_nu over identical wavelength intervals. '
        'This accounts for unequal wavelength sampling without imposing the same continuum shape.', '',
        '| Observed interval (µm) | Mean fluxes in chronological order (mJy) | First-to-last change | Nominal significance | Shared 2% scenario | Shared 2% + excess scatter |',
        '|---|---|---|---|---|---|',
    ]
    for result in results:
        if result['degree'] != 2:
            continue
        nominal = result['independent_errors']
        comparison = nominal['contrasts'][1]
        shared = result['epoch_correlated_calibration']['contrasts'][1]
        inflated = result['excess_scatter_and_correlated_calibration']['contrasts'][1]
        means = ', '.join(f'{f:.3f}' for f in nominal['mean_flux_mjy'])
        lo, hi = result['observed_window_um']
        lines.append(f'| {lo:g}–{hi:g} | {means} | {comparison["change_percent"]:+.1f}% | '
                     f'{comparison["signed_sigma"]:.1f}σ | {shared["signed_sigma"]:.1f}σ | {inflated["signed_sigma"]:.1f}σ |')
    broad_results = [r for r in results if r['observed_window_um'] == [4.0, 5.0]]
    changes = [r['independent_errors']['contrasts'][1]['change_percent'] for r in broad_results]
    interpretation = (
        'P1823 interpretations must use the global constant-flux tests, all epoch pairs, '
        'continuum-model checks and fit residuals recorded in the JSON. The isolated '
        '4.06996 micron point is excluded in the primary calculation and restored '
        'with --include-excluded for sensitivity analysis.'
        if prefix == 'P1823' else
        'P2190 shows a candidate long-wavelength brightening concentrated between '
        'summer 2025 and early 2026; the last two epochs are consistent with each other. '
        'There is no comparable significant change at 3.5–4 µm. '
        'Proposed wording: “Three SPHEREx epochs reveal a rising near-infrared '
        'continuum with tentative long-wavelength brightening, providing the '
        'short-wavelength context for MIRI spectroscopy.”'
        if prefix == 'P2190' else
        'The 4.5–5 µm epoch ordering is not monotonic. These results support a candidate '
        'few-percent continuum change, not a secure claim of monotonic brightening. '
        'Proposed wording: “Three SPHEREx epochs reveal a rising near-infrared '
        'continuum with indications of continued evolution, providing the '
        'short-wavelength context for the proposed MIRI spectroscopy.”')
    lines += ['',
        'The nominal column treats the supplied total errors as independent. The CSV includes '
        'a 2% calibration term. In the sensitivity scenario this same term is shared within '
        'each epoch, independently between epochs, while retaining the supplied measurement '
        'variance; it is not added twice. A scale error common to all epochs would instead '
        'largely cancel. The CSV alone cannot determine the correct covariance.', '',
        f'Linear, quadratic and cubic fits give first-to-last changes of '
        f'{min(changes):.1f}–{max(changes):.1f}% for 4–5 µm. ' + interpretation, '',
        '## Limits', '',
    ] + ['- ' + s for s in summary['limitations']]
    lines += ['',
        'Instrument context: [IRSA spectral calibration products](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/spherex_spectral_calibrations.html) '
        'provide the spectral response curves needed for a full response-based comparison. '
        'The numerical findings above come from the supplied CSV, not from that documentation.', '']
    (HERE / f'spherex_{output_name}_variability.md').write_text('\n'.join(lines))
    print('\n'.join(lines[:12]))


if __name__ == '__main__':
    main()
