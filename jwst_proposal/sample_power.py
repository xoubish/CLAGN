"""Conditional sample-size calculation for the primary rising/fading contrast.

No scatter or effect size is inferred from the illustrative Figure 2.
Run with /opt/anaconda3/bin/python jwst_proposal/sample_power.py.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.stats import nct, t

HERE = Path(__file__).resolve().parent


def power(n_per_group, scatter_dex, offset_dex, alpha=.05):
    df = 2*n_per_group-2
    se = scatter_dex*np.sqrt(2/n_per_group)
    critical = t.ppf(1-alpha/2, df)
    nc = offset_dex/se
    return float(nct.sf(critical, df, nc)+nct.cdf(-critical, df, nc))


def main():
    baseline = dict(n_per_group=12, scatter_dex=.20, offset_dex=.25, alpha=.05)
    scenarios = []
    for scatter in [.15, .20, .25, .30]:
        threshold = brentq(lambda offset: power(12, scatter, offset)-.8, 0., 2.)
        scenarios.append(dict(scatter_dex=scatter,
                              mean_difference_se_dex=float(scatter*np.sqrt(2/12)),
                              offset_dex_for_80_percent_power=float(threshold),
                              ratio_for_80_percent_power=float(10**threshold),
                              power_for_025_dex=power(12, scatter, .25)))
    sample_sizes = [dict(n_per_group=n, total=2*n, power=power(n, .2, .25))
                    for n in [8, 10, 11, 12, 15, 18]]
    # Independent simulation checks the t-statistic implementation and both tails.
    rng = np.random.default_rng(260926)
    trials = 200000
    a = rng.normal(0., .20, size=(trials, 12))
    b = rng.normal(0., .20, size=(trials, 12))
    denominator = np.sqrt((a.var(axis=1, ddof=1)+b.var(axis=1, ddof=1))/12)
    mean_difference = a.mean(axis=1)-b.mean(axis=1)
    critical = t.ppf(.975, 22)
    false_positive = float(np.mean(np.abs(mean_difference/denominator) > critical))
    simulated_power = float(np.mean(np.abs((mean_difference+.25)/denominator) > critical))
    exact = power(**baseline)
    assert abs(false_positive-.05) < 5*np.sqrt(.05*.95/trials)
    assert abs(simulated_power-exact) < 5*np.sqrt(exact*(1-exact)/trials)
    result = dict(
        status='Conditional planning sensitivity, not a measured scatter or predicted population effect.',
        diagnostic='Difference of group means in delta_warm = log10(observed/predicted nuclear 8--13 micron continuum luminosity).',
        baseline=baseline, baseline_power=exact,
        baseline_mean_difference_se_dex=float(.2*np.sqrt(2/12)),
        baseline_offset_ratio=float(10**.25),
        test='Independent equal-variance Gaussian residuals; pooled two-sample t-test; two-sided alpha=0.05; df=22.',
        assumptions=[
            '0.20 dex is an explicit planning assumption, not estimated from the sample or Figure 2.',
            'Scatter includes intrinsic dispersion and independent per-source measurement/model errors.',
            '0.25 dex is a design effect-size benchmark, not a dust-model population prediction.',
            'The calculation includes no shared/group-dependent systematic errors or covariate-estimation penalty.',
            'A 5% two-sided false-positive threshold is not a three-sigma detection criterion.',
            'Final inference must account for unequal errors, covariates, and correlated uncertainties.',
        ],
        scatter_scenarios=scenarios, sample_size_scenarios=sample_sizes,
        stricter_threshold_check=dict(alpha=.0027, power=power(12, .2, .25, alpha=.0027)),
        monte_carlo=dict(seed=260926, trials=trials, false_positive_rate=false_positive,
                         empirical_power=simulated_power, analytic_power=exact),
        references=[
            'https://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm',
            'https://www.itl.nist.gov/div898/handbook/prc/section2/prc222.htm',
        ])
    (HERE/'inputs/sample_power.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(baseline_power=exact, scenarios=scenarios,
                          monte_carlo=result['monte_carlo']), indent=2))


if __name__ == '__main__':
    main()
