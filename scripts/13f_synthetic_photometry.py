"""
13f_synthetic_photometry.py  --  synthetic ZTF g and r AB magnitudes from every cached archival (DR16/DR17 spec-lite) spectrum,
so the archival epoch can be compared directly with ZTF photometry (the CLAS+ test of Nakazono+2026) without the PSF-vs-fibre
systematic.  Photon-weighted AB magnitude through the SVO Palomar/ZTF.g and ZTF.r curves (data/external/ZTF_*.dat).
Also returns the fraction of the filter covered by good pixels (ivar > 0) and the median S/N in the band.
Usage: 13f_synthetic_photometry.py [--names csv]   (default: every P* directory in data/spectra_cache with a spec-*.fits that is NOT
an SDSS-V file, i.e. the DR16 spectrum)  ->  data/synthetic_ztf_dr16.csv
"""
import os, sys, glob, time, argparse, numpy as np, pandas as pd
from astropy.io import fits
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(HERE, 'data'); CACHE = os.path.join(DATA, 'spectra_cache')
C_AA = 2.99792458e18                     # speed of light in Angstrom/s
FILT = {b: np.loadtxt(os.path.join(DATA, 'external', f'ZTF_{b}.dat')) for b in ('g', 'r')}


def synth(lam, flam, ivar, band):
    """AB mag of a spectrum (flam in 1e-17 erg/s/cm2/A) through filter `band`; photon-weighted mean f_nu."""
    fl, ft = FILT[band][:, 0], FILT[band][:, 1]
    T = np.interp(lam, fl, ft, left=0.0, right=0.0); good = (ivar > 0) & np.isfinite(flam)
    cover = np.trapz(T * good, lam) / np.trapz(T, lam)
    if cover < 0.9:
        return np.nan, cover, np.nan
    w = T * good; num = np.trapz(w * flam * 1e-17 * lam, lam); den = np.trapz(w * C_AA / lam, lam)
    fnu = num / den
    snr = np.nanmedian((flam * np.sqrt(ivar))[(T > 0.1 * ft.max()) & good])
    return (-2.5 * np.log10(fnu) - 48.6) if fnu > 0 else np.nan, cover, snr


def one(name):
    d = os.path.join(CACHE, name)
    files = [f for f in glob.glob(os.path.join(d, 'spec-*.fits')) if 'v6_' not in f and '/daily/' not in f]
    # the DR16 file is the one whose name is spec-PLATE-MJD-FIBER with a 4-digit fibre and MJD < 59000 (SDSS-I..IV)
    files = [f for f in files if int(os.path.basename(f).split('-')[2]) < 59000] or files
    if not files:
        return dict(name=name)
    out = dict(name=name, dr16_file=os.path.basename(files[0]))
    try:
        h = fits.open(files[0]); t = h[1].data; lam = 10 ** t['loglam']; flam = t['flux'].astype(float); ivar = t['ivar'].astype(float)
        out['mjd_spec'] = int(os.path.basename(files[0]).split('-')[2])
        for b in ('g', 'r'):
            m, cov, snr = synth(lam, flam, ivar, b); out[f'syn_{b}'] = m; out[f'cover_{b}'] = cov; out[f'snr_{b}'] = snr
        # SDSS pipeline spectrophotometric flux in the same file, for a cross-check (nanomaggies)
        for hdu in h[1:]:
            cols = getattr(hdu, 'columns', None)
            if cols is not None and 'SPECTROFLUX' in cols.names:
                sf = np.atleast_2d(hdu.data['SPECTROFLUX'])[0]; out['pipe_g'] = 22.5 - 2.5 * np.log10(sf[1]) if sf[1] > 0 else np.nan; out['pipe_r'] = 22.5 - 2.5 * np.log10(sf[2]) if sf[2] > 0 else np.nan; break
        h.close()
    except Exception as e:
        out['error'] = str(e)[:80]
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--names', default=None); ap.add_argument('--out', default='synthetic_ztf_dr16.csv'); a = ap.parse_args()
    names = pd.read_csv(a.names).name.astype(str).tolist() if a.names else sorted(n for n in os.listdir(CACHE) if n.startswith('P') and os.path.isdir(os.path.join(CACHE, n)))
    t0 = time.time(); rows = []
    with ThreadPoolExecutor(max_workers=6) as ex:
        for k, r in enumerate(ex.map(one, names)):
            rows.append(r)
            if (k + 1) % 2000 == 0: print(f'[{time.time()-t0:5.0f}s] {k+1}/{len(names)}', flush=True)
    out = pd.DataFrame(rows); out.to_csv(os.path.join(DATA, a.out), index=False)
    ok = out.syn_r.notna() if 'syn_r' in out else pd.Series(False, index=out.index)
    print(f'done in {time.time()-t0:.0f}s: {int(ok.sum())} of {len(out)} objects with synthetic r (g: {int(out.syn_g.notna().sum()) if "syn_g" in out else 0}); '
          f'median syn_r - pipeline r = {(out.syn_r - out.pipe_r).median():+.3f} mag' if 'pipe_r' in out else '', flush=True)
