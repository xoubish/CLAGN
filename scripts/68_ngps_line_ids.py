"""Identify stellar absorption lines in a reduced NGPS spectrum and measure its velocity.

For targets whose spectrum turns out to be a star (H5MUSES-327 on October 26). Each line is a
Gaussian absorption on a linear continuum within +-25 A of its vacuum rest wavelength. The adopted
velocity is the inverse-variance mean of the clean lines; its error is the larger of the formal error
and the line-to-line scatter / sqrt(N), because U/G wavelengths have no sky-line check.
Writes <products>/line_ids/NAME.json and <products>/plots/NAME_lines.png; observed_ngps.py
puts both on the observer page.

Usage: python scripts/68_ngps_line_ids.py H5MUSES-327 [--products oct26_data/reduction_20261027/products_feige34]
"""
import argparse
import json
from pathlib import Path
import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
C = 299792.458
# Vacuum rest wavelengths. Groups share one label on the plots.
# Label row 0-2 (bottom to top) keeps neighbouring labels apart.
GROUPS = [('Ca II K, H', [3934.777, 3969.591], 2), ('Hδ', [4102.892], 1), ('G band', [4305.6], 0), ('Hγ', [4341.684], 1),
          ('Hβ', [4862.683], 0), ('Mg b', [5168.76, 5174.125, 5185.048], 1), ('Na D', [5891.583, 5897.558], 0),
          ('Hα (emission)', [6564.61], 1), ('Ca II triplet', [8500.36, 8544.44, 8664.52], 0)]
# Lines used for the velocity: isolated and well fitted by one Gaussian. Na D (blended doublet),
# the Mg b blend, the G band (molecular) and the weak Hα are marked but not measured.
MEASURE = {'Ca II K': 3934.777, 'Ca II H': 3969.591, 'Hδ': 4102.892, 'Hγ': 4341.684, 'Hβ': 4862.683,
           'Mg I 5185': 5185.048, 'Ca II 8500': 8500.36, 'Ca II 8544': 8544.44, 'Ca II 8664': 8664.52}


def measure(w, f, e, m):
    rows = []
    for name, rest in MEASURE.items():
        q = m & (abs(w-rest) < 25) & np.isfinite(f)
        x, y, s = w[q], f[q], e[q]
        model = lambda x, a, mu, sg, b0, b1: b0+b1*(x-rest)-a*np.exp(-.5*((x-mu)/sg)**2)
        p, cov = curve_fit(model, x, y, p0=[.3*np.median(y), rest, 2, np.median(y), 0], sigma=s,
                           bounds=([0, rest-15, .5, -np.inf, -np.inf], [np.inf, rest+15, 10, np.inf, np.inf]))
        rows.append(dict(line=name, rest_vac_A=rest, observed_A=float(p[1]),
                         velocity_kms=float((p[1]-rest)/rest*C), error_kms=float(np.sqrt(cov[1, 1])/rest*C),
                         ew_A=float(p[0]*p[2]*np.sqrt(2*np.pi)/p[3]), depth_over_noise=float(p[0]/np.median(s))))
    good = [r for r in rows if r['depth_over_noise'] > 5]
    # H-alpha is checked separately: in H5MUSES-327 it is a narrow emission line on the stellar trace.
    q = m & (abs(w-6564.61) < 30)
    model = lambda x, a, mu, sg, b0, b1: b0+b1*(x-6564.61)+a*np.exp(-.5*((x-mu)/sg)**2)
    p, cov = curve_fit(model, w[q], f[q], p0=[.1*np.median(f[q]), 6563, 1.5, np.median(f[q]), 0], sigma=e[q])
    emission = dict(line='Hα emission', rest_vac_A=6564.61, observed_A=float(p[1]),
                    velocity_kms=float((p[1]-6564.61)/6564.61*C), error_kms=float(np.sqrt(cov[1, 1])/6564.61*C),
                    ew_A=float(-p[0]*abs(p[2])*np.sqrt(2*np.pi)/p[3]), fwhm_kms=float(2.3548*abs(p[2])/6564.61*C),
                    flux_1e17=float(p[0]*abs(p[2])*np.sqrt(2*np.pi)), snr=float(p[0]/np.sqrt(cov[0, 0])))
    v = np.array([r['velocity_kms'] for r in good]); err = np.array([r['error_kms'] for r in good])
    mean = float(np.average(v, weights=1/err**2))
    error = float(max(1/np.sqrt(np.sum(1/err**2)), np.std(v)/np.sqrt(len(v))))
    return rows, mean, error, emission


def plot(name, w, f, m, velocity, out):
    shift = 1+velocity/C
    fig = plt.figure(figsize=(13, 7.5))
    top = fig.add_axes([.06, .45, .92, .5]); zooms = [fig.add_axes([.06+i*.315, .07, .29, .28]) for i in range(3)]
    q = m & np.isfinite(f)
    top.plot(w[q], f[q], lw=.5, color='#333')
    lo, hi = np.percentile(f[q], [1, 99.5]); top.set_ylim(min(0, lo), hi*1.3)
    for label, rests, row in GROUPS:
        for rest in rests:
            top.axvline(rest*shift, color='#d94c36', lw=.6, alpha=.7)
        top.text(np.mean(rests)*shift, hi*(1.07+.06*row), label, ha='center', fontsize=8, color='#d94c36')
    top.set(xlim=(3240, 10235), xlabel='Vacuum heliocentric wavelength (Å)', ylabel='Flux (10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹)',
            title=f'{name}: stellar absorption lines at z ≈ 0 (cz = {velocity:.0f} km/s); markers at that velocity')
    for ax, (label, (a, b)) in zip(zooms, [('Ca II K, H', (3900, 4000)), ('Mg b · Na D', (5140, 5920)), ('Ca II triplet', (8470, 8690))]):
        s = q & (w > a) & (w < b)
        ax.plot(w[s], f[s], lw=.7, color='#333')
        for _, rests, _ in GROUPS:
            for rest in rests:
                if a < rest*shift < b:
                    ax.axvline(rest*shift, color='#d94c36', lw=.7, alpha=.7)
        ax.set_title(label, fontsize=9); ax.set_xlim(a, b); ax.tick_params(labelsize=8)
    fig.savefig(out, dpi=150); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('name')
    p.add_argument('--products', type=Path, default=ROOT/'oct26_data/reduction_20261027/products_feige34')
    args = p.parse_args()
    t = Table.read(args.products/f'spectra/{args.name}.fits')
    w = np.asarray(t['WAVE_VAC_HELIO_A']); f = np.asarray(t['FLUX']); e = np.asarray(t['FLUX_ERR']); m = np.asarray(t['MASK'])
    rows, velocity, error, emission = measure(w, f, e, m)
    (args.products/'line_ids').mkdir(exist_ok=True)
    plot(args.name, w, f, m, velocity, args.products/f'plots/{args.name}_lines.png')
    result = dict(name=args.name, classification='star', velocity_kms=round(velocity, 1), velocity_error_kms=round(error, 1),
                  z=round(velocity/C, 6), frame='heliocentric, vacuum wavelengths',
                  marks=[dict(label=label, rest=rests, row=row) for label, rests, row in GROUPS], lines=rows,
                  halpha_emission=emission if emission['snr'] > 5 and emission['ew_A'] < 0 else None,
                  note='Absorption-line velocity of a star; not a cosmological redshift.')
    (args.products/f'line_ids/{args.name}.json').write_text(json.dumps(result, indent=2))
    for r in rows:
        print(f"{r['line']:11s} {r['velocity_kms']:7.1f} ± {r['error_kms']:5.1f} km/s  EW {r['ew_A']:5.2f} Å")
    print(f'adopted cz = {velocity:.1f} ± {error:.1f} km/s, z = {velocity/C:.6f}')
    print('H-alpha:', {k: round(v, 2) if isinstance(v, float) else v for k, v in emission.items()})


if __name__ == '__main__':
    main()
