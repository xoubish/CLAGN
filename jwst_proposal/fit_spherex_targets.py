"""Host + disc + blackbody continuum fit to a SPHEREx spectrum, for any target.

Same model as fit_spherex_hot_dust.py (Ell5 host template normalised at rest 1.6 um,
disc F_nu ~ lambda^-1/3, blackbody with free temperature, NNLS amplitudes, channel
bandpass integration, emission-line and 3.3-um PAH windows masked), applied to one
pooled spectrum per target. Writes inputs/spherex_fits/<id>.json and <id>_curve.csv.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import least_squares, nnls
from fit_spherex_hot_dust import blackbody, LINES, TBOUNDS

HERE = Path(__file__).resolve().parent; OUT = HERE/'inputs/spherex_fits'; OUT.mkdir(exist_ok=True)
SPECTRA = {'R01': HERE/'inputs/spaxel_scryer_P9694_ra11p5333_dec9p1225_cleaned.csv'}
ALPHA, TEMPLATE, MASK_V = 1/3, 'Ell5', .02


def load(tid):
    path = SPECTRA.get(tid, HERE/f'inputs/spherex_extracted/{tid}.csv')
    d = pd.read_csv(path)
    if 'spectrum_type' in d: d = d[d.spectrum_type.eq('cleaned')]
    d = d[np.isfinite(d.flux_mjy) & (d.flux_err_mjy > 0)]
    return d.reset_index(drop=True), str(path.relative_to(HERE))


class TargetFit:
    def __init__(self, tid, z, quadrature=32):
        d, self.source = load(tid); self.z = z
        lo = (d.wavelength_um-d.wavelength_half_width_um).values/(1+z)
        hi = (d.wavelength_um+d.wavelength_half_width_um).values/(1+z)
        keep = np.ones(len(lo), bool)
        for l in LINES: keep &= ~((lo <= l*(1+MASK_V)) & (hi >= l*(1-MASK_V)))
        keep &= ~((lo <= 3.35) & (hi >= 3.20))
        nodes, w = np.polynomial.legendre.leggauss(quadrature)
        self.waves = (lo[:, None]+hi[:, None])/2+(hi-lo)[:, None]*nodes/2; self.w = w/2
        tab = np.loadtxt(HERE/'inputs/host_templates'/f'{TEMPLATE}_template_norm.sed')
        hw = tab[:, 0]/1e4; hf = tab[:, 1]*hw**2; hf /= np.interp(1.6, hw, hf)
        self.host = np.interp(self.waves, hw, hf) @ self.w
        self.disc = self.waves**(-ALPHA) @ self.w
        self.keep = keep; self.y = d.flux_mjy.values[keep]; self.err = d.flux_err_mjy.values[keep]
        self.n_total = len(d)

    def design(self, T):
        return np.column_stack([self.host, self.disc, blackbody(self.waves, T) @ self.w])

    def solve(self, T):
        a = self.design(T)[self.keep]
        amps, _ = nnls(a/self.err[:, None], self.y/self.err)
        return (a @ amps-self.y)/self.err, amps

    def fit(self):
        best = None
        for t0 in (850., 1200., 1600., 2000.):
            r = least_squares(lambda p: self.solve(p[0])[0], [t0], bounds=([TBOUNDS[0]], [TBOUNDS[1]]), x_scale='jac')
            if best is None or r.fun @ r.fun < best.fun @ best.fun: best = r
        T = float(best.x[0]); res, amps = self.solve(T)
        # 1-sigma from delta chi2 = 1 scan
        grid = np.linspace(max(T-400, TBOUNDS[0]), min(T+400, TBOUNDS[1]), 161)
        chi = np.array([self.solve(t)[0] @ self.solve(t)[0] for t in grid]); c0 = res @ res
        ok = grid[chi <= c0+1]
        return dict(temperature_K=T, temperature_1sigma_K=[float(ok.min()), float(ok.max())] if len(ok) else None,
                    amplitudes_mjy=amps.tolist(), amplitude_order=['host_rest1.6um', 'disc_rest1um', 'dust_rest2um'],
                    chi2=float(c0), n_used=int(self.keep.sum()), n_total=int(self.n_total), dof=int(self.keep.sum()-4),
                    alpha=ALPHA, host_template=TEMPLATE, mask_velocity=MASK_V, source=self.source, z=self.z)

    def curve(self, T, amps, obs_wave):
        rest = obs_wave/(1+self.z)
        tab = np.loadtxt(HERE/'inputs/host_templates'/f'{TEMPLATE}_template_norm.sed')
        hw = tab[:, 0]/1e4; hf = tab[:, 1]*hw**2; hf /= np.interp(1.6, hw, hf)
        host = amps[0]*np.interp(rest, hw, hf); disc = amps[1]*rest**(-ALPHA); dust = amps[2]*blackbody(rest, T)
        return pd.DataFrame(dict(wavelength_observed_um=obs_wave, total_mjy=host+disc+dust, host_mjy=host, disc_mjy=disc, dust_mjy=dust))


def main(ids):
    sample = pd.read_csv(HERE/'inputs/jwst_sample_cycle6.csv').set_index('id')
    for tid in ids:
        f = TargetFit(tid, float(sample.loc[tid, 'z'])); res = f.fit()
        (OUT/f'{tid}.json').write_text(json.dumps(res, indent=2))
        f.curve(res['temperature_K'], res['amplitudes_mjy'], np.geomspace(.72, 5.1, 400)).to_csv(OUT/f'{tid}_curve.csv', index=False)
        print(f"{tid}: T = {res['temperature_K']:.0f} K (1σ {res['temperature_1sigma_K']}), amps host/disc/dust = "
              f"{[round(a, 2) for a in res['amplitudes_mjy']]} mJy, chi2/dof = {res['chi2']:.0f}/{res['dof']}, n = {res['n_used']}/{res['n_total']}")


if __name__ == '__main__':
    main(sys.argv[1:] or ['R01', 'F04'])
