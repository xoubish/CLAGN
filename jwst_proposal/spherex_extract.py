"""Extract a SPHEREx spectrum from QR2 level-2b cutouts (IRSA IBE) at a source position.

Per image: PSF-weighted photometry of the background-subtracted IMAGE (MJy/sr) at the
source pixel, using the oversampled PSF nearest to the detector position; wavelength and
bandwidth from the WCS-WAVE lookup table at the source's active-pixel coordinates.
Output CSV columns match the spaxel_scryer files used by the figures:
wavelength_um, wavelength_half_width_um, flux_mjy, flux_err_mjy, mjds, plus detector/obsid.
"""
import sys, glob, warnings
from pathlib import Path
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import shift as ndshift
warnings.filterwarnings('ignore')
from astropy import log as _alog; _alog.setLevel('ERROR')
import logging; logging.getLogger('astropy').setLevel(logging.ERROR)

PIX_ARCSEC = 6.15
SR_PER_PIX = (PIX_ARCSEC/206264.806)**2
MJY_PER_MJYSR = SR_PER_PIX*1e9
# FLAGS bits that make a pixel unusable (MP_SOURCE=21 marks known sources and is expected at the target)
BAD_BITS = [0, 1, 2, 4, 5, 6, 7, 9, 10, 11, 15, 17, 19]
BAD_MASK = sum(1 << b for b in BAD_BITS)
SOURCE_BIT = 1 << 21


def psf_pixel_weights(psfcube, xa, ya, xc, yc, shape, oversamp=10):
    """PSF fraction per cutout pixel for a point source at (xc, yc) (0-based pixel coords)."""
    n = psfcube.shape[0]; side = int(round(np.sqrt(n)))
    gx = np.linspace(0, 2040, side); ix = int(np.argmin(abs(gx-xa))); iy = int(np.argmin(abs(gx-ya)))
    psf = psfcube[iy*side+ix].astype(float); psf /= psf.sum()
    # shift the oversampled PSF so its peak sits at the fractional pixel position, then bin
    c = (psf.shape[0]-1)/2.
    fx, fy = (xc-np.floor(xc))*oversamp, (yc-np.floor(yc))*oversamp
    sh = ndshift(psf, (fy, fx), order=1, mode='constant', cval=0.)
    # bin to detector pixels: the PSF centre (c) maps to pixel floor(xc)
    npix = int(np.ceil(psf.shape[0]/oversamp))+2
    w = np.zeros(shape)
    for j in range(shape[0]):
        for i in range(shape[1]):
            y0 = (j-np.floor(yc))*oversamp+c-oversamp/2.+.5; x0 = (i-np.floor(xc))*oversamp+c-oversamp/2.+.5
            ys, xs = slice(int(max(y0, 0)), int(max(min(y0+oversamp, psf.shape[0]), 0))), slice(int(max(x0, 0)), int(max(min(x0+oversamp, psf.shape[1]), 0)))
            if ys.stop > ys.start and xs.stop > xs.start: w[j, i] = sh[ys, xs].sum()
    return w


def extract_one(path, ra, dec):
    with fits.open(path) as h:
        img = h['IMAGE'].data.astype(float); hdr = h['IMAGE'].header
        flags = h['FLAGS'].data; var = h['VARIANCE'].data.astype(float)
        wcs = WCS(hdr); x, y = [float(v) for v in wcs.all_world2pix(ra, dec, 0)]
        if not (2 <= x <= img.shape[1]-3 and 2 <= y <= img.shape[0]-3): return None
        xa, ya = [float(v) for v in WCS(hdr, fobj=h, key='A').all_pix2world(x, y, 0)]
        tab = h['WCS-WAVE'].data; X, Y, V = tab['X'][0], tab['Y'][0], tab['VALUES'][0]
        lam, bw = RegularGridInterpolator((Y, X), V, bounds_error=False, fill_value=None)((ya, xa))
        w = psf_pixel_weights(h['PSF'].data, xa, ya, x, y, img.shape)
        mjd = hdr.get('MJD-AVG', hdr.get('MJD-OBS')); det = hdr.get('DETECTOR'); obsid = hdr.get('OBSID')
    yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]; rr = np.hypot(xx-x, yy-y)
    good = np.isfinite(img) & ((flags & BAD_MASK) == 0) & np.isfinite(var) & (var > 0)
    bkg_px = good & ((flags & SOURCE_BIT) == 0) & (rr > 2.5) & (rr < 7)
    if bkg_px.sum() < 15 or not good[int(round(y)), int(round(x))]: return None
    bkg = np.median(img[bkg_px]); bkg_var = (1.2533*np.std(img[bkg_px])/np.sqrt(bkg_px.sum()))**2
    ap = good & (rr <= 3.2)
    wa = w*ap; norm = (wa**2/var).sum()
    if norm <= 0 or wa.sum() < .5: return None
    flux = ((img-bkg)*wa/var).sum()/norm            # PSF-weighted, variance-weighted estimate of total source flux
    err = np.sqrt(1./norm+bkg_var*(wa/var).sum()**2/norm**2)
    return dict(wavelength_um=float(lam), wavelength_half_width_um=float(bw)/2., flux_mjy=float(flux*MJY_PER_MJYSR),
                flux_err_mjy=float(err*MJY_PER_MJYSR), mjds=float(mjd), detector=int(det), obsid=obsid,
                psf_fraction=float(wa.sum()), x=x, y=y)


def main(folder, ra, dec, out):
    rows = []
    for p in sorted(glob.glob(f'{folder}/*.fits')):
        try:
            r = extract_one(p, ra, dec)
        except Exception as e:
            r = None
        if r: rows.append(r)
    d = pd.DataFrame(rows).sort_values('wavelength_um'); d.to_csv(out, index=False)
    print(f'{out}: {len(d)} measurements from {len(glob.glob(f"{folder}/*.fits"))} cutouts; wavelength {d.wavelength_um.min():.2f}-{d.wavelength_um.max():.2f} um')


if __name__ == '__main__':
    main(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4])
