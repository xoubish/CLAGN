"""Independent quadrature checks and compact plot from saved fitted parameters."""
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss
from analyse_forecasts import Audit

HERE = Path(__file__).resolve().parent


def main():
    path = HERE / 'adequacy_and_forecasts.json'
    results = json.loads(path.read_text())
    audit = Audit()
    nodes, weights = leggauss(64)
    audit.delay_nodes = (nodes + 1) / 2
    audit.delay_weights = weights / 2
    wave = np.load(HERE / 'conditional_curves.npz')['rest_um']
    observed = (wave * (1 + .2377224))[:, None]
    band_weights = np.ones_like(observed)
    times = np.full(len(wave), 2028.)
    fixed = results['sensitivity_5pct_extra_scatter']['0.0']
    for item in results['whole_spectrum_counterexamples']:
        eta = item['evolution']
        evolving = results['sensitivity_5pct_extra_scatter'][str(eta)]
        fm, fo = audit.design(fixed['parameters'], 0., observed, band_weights, times, 1.)
        em, eo = audit.design(evolving['parameters'], eta, observed, band_weights, times, 1.)
        residual = (fm @ item['fixed_amplitudes'] + fo) / (em @ evolving['amplitudes'] + eo) - 1
        item['independent_64_node_rms_fractional_mismatch'] = float(np.sqrt(np.mean(residual**2)))
        item['independent_64_node_maximum_fractional_mismatch'] = float(np.max(abs(residual)))
        print('64-node counterexample', eta, 'maximum difference', np.max(abs(residual)))
    original_driver = audit.driver
    audit.driver = lambda t, p, future=1.: np.ones_like(np.asarray(t), dtype=float)
    baseline, offset = audit.design(fixed['parameters'], 0.)
    for eta in [-.5, -1.]:
        matrix, other_offset = audit.design(fixed['parameters'], eta)
        np.testing.assert_allclose(matrix, baseline, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(other_offset, offset, rtol=1e-12, atol=1e-12)
    audit.driver = original_driver
    results['constant_unit_illumination_fixed_equals_evolving'] = True
    path.write_text(json.dumps(results, indent=2) + '\n')

    data = np.genfromtxt(HERE / 'counterexample_spectra.csv', delimiter=',', names=True)
    fig, axes = plt.subplots(2, 1, figsize=(7, 5), sharex=True, layout='constrained',
                             gridspec_kw={'height_ratios': [3, 1]})
    axes[0].plot(data['rest_um'], data['evolving_observed_mjy'], color='#be513b',
                 label='Evolving inner dust mass (η = −0.5)')
    axes[0].plot(data['rest_um'], data['fixed_observed_mjy'], '--', color='#24658a',
                 label='Refitted fixed dust distribution')
    axes[0].set(yscale='log', ylabel='Predicted flux density (mJy)',
                title='Different dust responses, nearly identical MIRI spectra')
    axes[0].legend(fontsize=9)
    axes[1].plot(data['rest_um'], 100 * data['fixed_over_evolving_minus1'], color='#24658a')
    axes[1].axhline(0, color='.5', lw=.6)
    axes[1].set(xlabel='Rest wavelength (µm)', ylabel='Fixed / evolving\n− 1 (%)')
    for axis in axes:
        for wavelength in [9.7, 18.]:
            axis.axvline(wavelength, color='.7', ls=':', lw=.6)
    fig.suptitle('Same date and future factor; additional 5% historical scatter assumed', fontsize=9)
    for extension in ['pdf', 'png']:
        fig.savefig(HERE / f'counterexample.{extension}', dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
