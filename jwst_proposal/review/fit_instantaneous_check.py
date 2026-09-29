"""Test a current-power shell model against the same R01 historical data.

Finite multistart profile search, not a posterior or likelihood-ratio test.
"""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
from pathlib import Path
import sys, json
import numpy as np
from scipy.optimize import minimize
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent), str(HERE/'p9694_predictivity')]
from analyse_forecasts import Audit, Constrained, BOUNDS

def main():
    a = Audit(); c = Constrained(a)
    archive = json.loads((HERE/'p9694_predictivity/adequacy_and_forecasts.json').read_text())
    records = []
    for seed in archive['sensitivity_5pct_extra_scatter'].values():
        p = np.array(seed['parameters']); p[-1] = 0
        def expand(x): return np.r_[x, 0.]
        opt = minimize(lambda x: c.fit(expand(x), 0.)[0], p[:-1], method='Powell',
                       bounds=BOUNDS[:-1], options={'maxfev':1100, 'ftol':.0002, 'xtol':.002})
        p = expand(opt.x)
        stat, amps, valid = c.fit(p, 0.)
        r = dict(parameters=p.tolist(), amplitudes=amps.tolist(), chi2_extra5pct=stat,
                 primary_chi2_at_same_parameters=a.fit(p, 0., floor=False)[0],
                 success=bool(opt.success), valid=bool(valid), nfev=int(opt.nfev))
        records.append(r); print(json.dumps(r), flush=True)
        (HERE/'instantaneous_check.json').write_text(json.dumps(dict(
            scope=__doc__, records=records,
            threshold=archive['adequacy_chi2_threshold_99pct']), indent=2)+'\n')

if __name__ == '__main__': main()
