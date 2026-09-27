"""Conditional P9694 forecast audit, not a torus posterior or observing simulation.

Fits static nonnegative dust-shell masses to dated SPHEREx and WISE data under
gridded optical driving histories. An alternative allows causal inner-dust mass
response. Full radiative transfer, extinction and calibration covariance remain
outside this deliberately transparent feasibility calculation.
"""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
from pathlib import Path
import json
import hashlib
import itertools

import numpy as np
import pandas as pd
from astropy.time import Time
from scipy.linalg import solve_triangular
from scipy.optimize import nnls, minimize
from scipy.stats import qmc, chi2
from numpy.polynomial.legendre import leggauss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from fit_spherex_hot_dust import Fit, Z
from make_fig2_diagnostics import GRAINS, efficiency, planck_nu

HERE = Path(__file__).resolve().parent
OUT = HERE/'review'/'p9694_predictivity'
PARAM_NAMES = ['optical_host_fraction', 'pre2010_relative_flux',
               '2014_relative_flux', 'log10_hot_radius_lightyears',
               'disc_alpha', 'NGPS_relative_calibration', 'delay_extent']
BOUNDS = np.array([[0., .75], [.5, 1.5], [.5, 1.3], [-1.3, .4],
                   [-1.5, 1/3], [.8, 1.2], [1., 2.]])


def year(mjd):
    return Time(mjd, format='mjd').decimalyear


class Audit:
    def __init__(self):
        self.target = json.loads((HERE/'inputs/P9694_figure_target.json').read_text())
        self.spx = Fit()
        self.optical = np.array(self.target['ztf']['g'], dtype=float)
        self.optical[:, 1] = 3631e3*10**(-.4*self.optical[:, 1])
        # The exported ZTF errors are magnitudes; fluxes above are AB mJy.
        self.optical[:, 2] *= self.optical[:, 1]*np.log(10)/2.5
        self.optical[:, 0] = year(self.optical[:, 0])
        # 60-day medians suppress intra-season fluctuations without inventing
        # unobserved monitoring. Interpolation is a scenario, not a recovered UV curve.
        bins = np.floor(self.optical[:, 0]*365.25/60)
        self.opt_t = np.array([np.median(self.optical[bins == b, 0]) for b in np.unique(bins)])
        self.opt_f = np.array([np.median(self.optical[bins == b, 1]) for b in np.unique(bins)])
        self.ref_flux = float(np.median(self.optical[self.optical[:, 0] < 2021, 1]))
        spec = []
        for e in self.target['spec']['epochs']:
            w, f = np.asarray(e['wave'], float), np.asarray(e['flux'], float)
            f = f*1e-17*w*w/2.99792458e18/1e-26
            val = float(np.nanmedian(f[(w > 8600) & (w < 9000)]))
            spec.append((float(year(e['mjd'])), val))
        self.spec = np.array(spec)
        self.ngps_ratio = self.spec[-1, 1]/self.spec[1, 1]
        self.ngps_g_proxy = np.interp(self.spec[1, 0], self.opt_t, self.opt_f)*self.ngps_ratio
        nodes, weights = leggauss(32)
        self.records, waves, ww = [], [], []
        d = self.spx.data
        for i in np.flatnonzero(self.spx.keep):
            waves.append(d['wavelength_um'][i]+d['wavelength_half_width_um'][i]*nodes)
            ww.append(weights/2)
            self.records.append(dict(kind='SPHEREx', time=float(year(d['mjds'][i])),
                                     flux=d['flux_mjy'][i], error=d['flux_err_mjy'][i],
                                     wave=d['wavelength_um'][i], input_row=int(i)))
        raw = pd.read_csv(HERE/'inputs/figure_neowise_exposures.csv', dtype={'cc_flags': str})
        raw = raw[raw['name'].eq('P9694')].copy()
        flags = raw['cc_flags'].str.strip().str.zfill(4)
        good = ((raw['qual_frame'] > 0) & flags.str[:2].eq('00') & (raw['dist_x'] < 3)
                & (raw['w1sigmpro'] > 0) & (raw['w2sigmpro'] > 0))
        raw = raw.loc[good].copy()
        raw['visit'] = np.round(raw['mjd']/180)
        self.band = {}
        for band, zero, iso in [('W1', 309.540, 3.3526), ('W2', 171.787, 4.6028),
                                ('W3', 31.674, 11.5608), ('W4', 8.363, 22.0883)]:
            response = np.loadtxt(HERE/f'inputs/wise_response/{band}.txt')
            use = response[:, 1] > .0001*response[:, 1].max()
            lo, hi = response[use, 0].min(), response[use, 0].max()
            wave = (lo+hi)/2+(hi-lo)*nodes/2
            weight = weights*np.interp(wave, response[:, 0], response[:, 1])/wave
            # WISE convention has zero colour correction for constant F_lambda:
            # Fnu_ref/Fnu_iso = (lambda/lambda_iso)^2.
            weight /= np.sum(weight*(wave/iso)**2)
            self.band[band] = (wave, weight)
            if band not in ('W1', 'W2'):
                continue
            col = band.lower()+'mpro'
            for _, group in raw.groupby('visit'):
                if len(group) < 3:
                    continue
                flux = zero*1e3*10**(-.4*float(group[col].median()))
                magerr = 1.2533*group[col].std(ddof=1)/np.sqrt(len(group))
                error = flux*np.log(10)/2.5*magerr
                waves.append(wave); ww.append(weight)
                self.records.append(dict(kind=band, time=float(year(group['mjd'].median())),
                                         flux=flux, error=error, wave=iso))
        # Only the two pre-NEOWISE W1 visits; no duplicate archival/NEOWISE rows.
        for mjd, flux, error in self.target['wise']['W1']:
            if year(mjd) >= 2011.5:
                continue
            waves.append(self.band['W1'][0]); ww.append(self.band['W1'][1])
            self.records.append(dict(kind='W1', time=float(year(mjd)), flux=flux, error=error, wave=3.3526))
        # Optical continuum snapshots have no supplied absolute calibration errors.
        # A stated 15% calibration scenario is used, not a measured uncertainty.
        for t, flux in self.spec:
            waves.append(np.full(32, .88)); ww.append(np.ones(32)/32)
            self.records.append(dict(kind='optical_15pct_assumed', time=t, flux=flux,
                                     error=.15*flux, wave=.88))
        self.waves, self.weights = np.array(waves), np.array(ww)
        self.times = np.array([r['time'] for r in self.records])
        self.y = np.array([r['flux'] for r in self.records])
        self.err = np.array([r['error'] for r in self.records])
        self.kind = np.array([r['kind'] for r in self.records])
        self.cov = np.diag(self.err**2)
        for name, fraction in [('W1', .024), ('W2', .028)]:
            vector = fraction*self.y*(self.kind == name)
            self.cov += np.outer(vector, vector)
        self.chol = np.linalg.cholesky(self.cov)
        self.chol5 = np.linalg.cholesky(self.cov+np.diag((.05*self.y)**2))
        self.shells = []
        reference = next(g for g in GRAINS if g['name'].startswith('Gra'))
        reference_heat = np.exp(np.interp(np.log(1500.), reference['log_t'], reference['log_cool']))/reference['q_heat']
        for g in GRAINS:
            temperatures = ([200., 350., 550., 750., 1000., 1200.] if g['name'].startswith('Sil')
                            else [250., 450., 700., 1000., 1300., 1600.])
            for temp in temperatures:
                cooling = np.exp(np.interp(np.log(temp), g['log_t'], g['log_cool']))
                radius = np.sqrt(reference_heat*g['q_heat']/cooling)
                normwave = np.geomspace(.2, 100, 1000)
                norm = np.max(efficiency(g['name'], normwave)*planck_nu(normwave, temp))
                self.shells.append(dict(grain=g, temperature=temp, cooling=cooling,
                                        radius=radius, norm=norm,
                                        inner=(temp >= (750 if g['name'].startswith('Sil') else 1000))))
        nodes, weights = leggauss(8)
        self.delay_nodes = (nodes+1)/2
        self.delay_weights = weights/2

    def driver(self, time, p, future=1.):
        host, early, middle, loglag, alpha, ngps, extent = p
        # Before the optical coverage, explicitly vary unmeasured illumination.
        ts = np.r_[1900., 2010., 2014., self.opt_t, self.spec[-1, 0], 2028.5, 2035.]
        fs = np.r_[early, early, middle, self.opt_f/self.ref_flux,
                   self.ngps_g_proxy/self.ref_flux*ngps,
                   self.ngps_g_proxy/self.ref_flux*ngps*future,
                   self.ngps_g_proxy/self.ref_flux*ngps*future]
        # Host fraction refers to pre-2021 reference g flux; optical proxy stays positive.
        nuclear = np.maximum(fs-host, .05)/(1-host)
        return np.interp(time, ts, nuclear)

    def design(self, p, evolution=0., waves=None, weights=None, times=None, future=1.):
        waves = self.waves if waves is None else np.asarray(waves)
        weights = self.weights if weights is None else np.asarray(weights)
        times = self.times if times is None else np.asarray(times)
        rest = waves/(1+Z)
        hostfraction, _, _, loglag, alpha, _, extent = p
        hostnorm = hostfraction*self.ref_flux/self.spx.host_at(.4722/(1+Z))
        discnorm = (1-hostfraction)*self.ref_flux/(.4722/(1+Z))**(-alpha)
        host = hostnorm*np.sum(self.spx.host_at(rest)*weights, axis=1)
        disc = discnorm*self.driver(times, p, future)*np.sum(rest**(-alpha)*weights, axis=1)
        matrix = np.zeros((len(times), len(self.shells)))
        for j, shell in enumerate(self.shells):
            g = shell['grain']
            delays = extent*10**loglag*shell['radius']*(1+Z)*self.delay_nodes
            illumination = self.driver(times[:, None]-delays, p, future)
            temp = np.exp(np.interp(np.log(shell['cooling']*illumination), g['log_cool'], g['log_t']))
            # Causal inner mass response: eta=0 fixed; eta<0 destruction on
            # brightening / replenishment on fading, with no extra response delay.
            mass_factor = illumination**evolution if shell['inner'] else np.ones_like(illumination)
            q = efficiency(g['name'], rest)
            spectrum = q[:, None, :]*planck_nu(rest[:, None, :], temp)/shell['norm']
            matrix[:, j] = np.sum(np.sum(spectrum*weights[:, None, :], axis=2)
                                   *mass_factor*self.delay_weights, axis=1)
        return matrix, host+disc

    def fit(self, p, evolution=0., floor=False):
        matrix, offset = self.design(p, evolution)
        chol = self.chol5 if floor else self.chol
        aw = solve_triangular(chol, matrix, lower=True, check_finite=False)
        yw = solve_triangular(chol, self.y-offset, lower=True, check_finite=False)
        amps, _ = nnls(aw, yw, maxiter=300)
        residual = aw@amps-yw
        return float(residual@residual), amps, matrix@amps+offset

    def warm_bounds(self, p, evolution):
        waves = np.array([self.band[b][0] for b in ('W3','W4')])
        weights = np.array([self.band[b][1] for b in ('W3','W4')])
        return self.design(p, evolution, waves, weights, np.array([2010.5,2010.5]))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    audit = Audit()
    samples = qmc.Sobol(len(BOUNDS), scramble=True, seed=9694).random_base2(7)
    parameters = qmc.scale(samples, BOUNDS[:,0], BOUNDS[:,1])
    records = []
    for evolution in (0., -.5, -1.):
        for i,p in enumerate(parameters):
            stat, amps, prediction = audit.fit(p, evolution)
            records.append(dict(evolution=evolution, parameters=p.tolist(), chi2=stat,
                                amplitudes=amps.tolist()))
        best = min([r for r in records if r['evolution']==evolution], key=lambda r:r['chi2'])
        print('GRID', evolution, best['chi2'], flush=True)
        # Refine the best starts; this is a profile fit, not posterior sampling.
        starts = sorted([r for r in records if r['evolution']==evolution],key=lambda r:r['chi2'])[:2]
        for start in starts:
            opt = minimize(lambda p:audit.fit(p,evolution)[0],start['parameters'],
                           method='Powell',bounds=BOUNDS,options={'maxfev':800,'ftol':1e-5,'xtol':.002})
            stat,amps,prediction = audit.fit(opt.x,evolution)
            records.append(dict(evolution=evolution,parameters=opt.x.tolist(),chi2=stat,
                                amplitudes=amps.tolist(),optimizer_success=bool(opt.success)))
        best = min([r for r in records if r['evolution']==evolution],key=lambda r:r['chi2'])
        print('REFINED',evolution,best['chi2'],best['parameters'],flush=True)
    payload=dict(n_measurements=len(audit.y),n_shell_coefficients=len(audit.shells),
                 parameter_names=PARAM_NAMES,parameter_bounds=BOUNDS.tolist(),
                 records=records,
                 data=audit.records,
                 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'profile_fits.json').write_text(json.dumps(payload,indent=2)+'\n')
    print('Saved profile_fits.json',flush=True)


if __name__=='__main__':
    main()
