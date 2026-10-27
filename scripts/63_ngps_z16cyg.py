"""16 Cyg A/B/C (Z16CygABC, night 2026-10-26): all-slice reduction and bright-star products.

The regular chain reduces one object in the central slice. 16 Cyg needs all three
slices and several objects per slice, so it gets its own PypeIt run in work/<night>/z16
using the night's calibrations (U: the helium-anchored solution). Products: the
unsaturated central-slice star from the two 1 s frames, fluxed with the nearest
standard, plus spatial-profile and outer/central diagnostics for the companion C.

Run with NGPS_NIGHT_CONFIG pointing at the October night.json:
    63_ngps_z16cyg.py reduce u      # PypeIt run(s) for the listed channels
    63_ngps_z16cyg.py products [--standard FEIGE34]   # default: nearest standard in time
"""
import importlib
import json
import shutil
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
env = importlib.import_module('54_reduce_ngps_sep23')
from pypeit import specobjs, coadd1d, fluxcalibrate
from pypeit.spectrographs.util import load_spectrograph
import matplotlib.pyplot as plt

Z16 = env.NIGHT/'z16'
OUT = env.BASE/'products_16cyg'
SHORT = ['0051', '0053']          # 1 s exposures; the 10 s frame saturates R/I/U
# Unsaturated star near the left end of the central slice (U has its own pixel layout).
STAR_SPAT = dict(u=(490, 520), g=(170, 205), r=(170, 205), i=(170, 205))
LIMITS = dict(u=(3241, 4340), g=(4340, 5790), r=(5790, 7740), i=(7740, 10234))
UNIT = '1e-17 erg s-1 cm-2 Angstrom-1'


def reduce(ch):
    directory = Z16/f'setup_{ch}'
    if not (directory/'Calibrations').exists():
        shutil.copytree(env.pypeit_file(ch).parent/'Calibrations', directory/'Calibrations')
    env.run([env.BIN/'run_pypeit', f'z16_{ch}.pypeit'], f'z16_run_{ch}.log', cwd=directory)


def star(ch, frame):
    path = next((Z16/f'setup_{ch}/Science').glob(f'spec1d_ngps_*_{frame}-*.fits'))
    sobjs = specobjs.SpecObjs.from_fitsfile(str(path), chk_version=False)
    central = [o for o in sobjs if o.SLITID == env.CONFIG['slits'][ch]
               and STAR_SPAT[ch][0] < o.SPAT_PIXPOS < STAR_SPAT[ch][1]]
    assert len(central) <= 1, (path, [o.NAME for o in central])
    # Object finding keeps at most three objects per slice and can miss the star.
    return path, sobjs, central[0] if central else None


def products(standard=None):
    for directory in ['channels', 'work']:
        (OUT/directory).mkdir(parents=True, exist_ok=True)
    report = dict(frames=SHORT, star_spat_window=STAR_SPAT, channels={})
    segments = []
    for ch in 'ugri':
        files = []
        for frame in SHORT:
            path, sobjs, obj = star(ch, frame)
            if obj is None:
                report['channels'].setdefault(ch, []).append(dict(file=path.name, star_found=False))
                continue
            one = specobjs.SpecObjs(specobjs=np.array([obj]), header=sobjs.header)
            out = OUT/'work'/f'{ch}_{path.name}'
            one.write_to_fits(sobjs.header, str(out), overwrite=True)
            files.append(out)
            report['channels'].setdefault(ch, []).append(dict(file=path.name,
                spat=float(obj.SPAT_PIXPOS), masked_box_pixels=int((~obj.BOX_MASK).sum())))
        mjd = np.mean([fits.getheader(p)['MJD'] for p in files])
        sens = []
        for p in sorted((env.BASE/'audit/standards'/ch).glob('sens_spec1d*fits')):
            h = fits.getheader(env.BASE/'audit/standards'/ch/p.name.removeprefix('sens_'))
            if standard is None or env.canonical(h['TARGET']) == standard:
                sens.append((abs(float(h['MJD'])-mjd), p, env.canonical(h['TARGET'])))
        _, sensfile, used = min(sens, key=lambda row: row[0])
        fluxed = [OUT/'work'/('fluxed_'+p.name) for p in files]
        fluxcalibrate.flux_calibrate([str(p) for p in files], [str(sensfile)]*len(files),
                                     outfiles=[str(p) for p in fluxed])
        spec = load_spectrograph(f'p200_ngps_{ch}')
        par = spec.default_pypeit_par()['coadd1d']
        par['scale_method'] = 'none'; par['ex_value'] = 'OPT'; par['flux_value'] = True
        names = [fits.getheader(p, 1)['NAME'] for p in fluxed]
        c = coadd1d.CoAdd1D.get_instance([str(p) for p in fluxed], names, spectrograph=spec, par=par)
        c.run()
        coadd = OUT/'work'/f'16CYG_STAR_{ch.upper()}.fits'
        c.save(str(coadd))
        with fits.open(coadd) as h:
            d = h[1].data
            w = np.asarray(d['wave']); f = np.asarray(d['flux']); e = np.asarray(d['sigma'])
            m = np.asarray(d['mask']).astype(bool) & (w > 100)
        tab = Table(dict(WAVE_VAC_HELIO_A=w[m], FLUX=f[m], FLUX_ERR=e[m],
                         CHANNEL=np.full(m.sum(), ch.upper())))
        tab.meta = dict(TARGET='16 Cyg (central-slice star)', FLUXSTD=used, BUNIT=UNIT,
                        NEXP=len(files), TEXPTIME=float(len(files)), SLICE='CENTRAL', EXTRACT='OPT',
                        NOTE='Not telluric-corrected; slit losses for a 6th-mag star off the slice centre')
        tab.write(OUT/f'channels/16CYG_STAR_{ch.upper()}.fits', overwrite=True)
        report['channels'][ch].append(dict(flux_standard=used, sensfile=sensfile.name))
        lo, hi = LIMITS[ch]
        segments.append(tab[(tab['WAVE_VAC_HELIO_A'] >= lo) & (tab['WAVE_VAC_HELIO_A'] < hi)])
    merged = vstack(segments, metadata_conflicts='silent')
    merged.sort('WAVE_VAC_HELIO_A')
    merged.write(OUT/'16CYG_STAR.fits', overwrite=True)
    merged.write(OUT/'16CYG_STAR.csv', format='ascii.csv', overwrite=True)
    fig, ax = plt.subplots(figsize=(13, 5))
    for ch, color in zip('UGRI', ['#7851a9', '#16843b', '#d94c36', '#8a451e']):
        q = merged['CHANNEL'] == ch
        ax.plot(merged['WAVE_VAC_HELIO_A'][q], merged['FLUX'][q], lw=.5, color=color, label=ch)
    for lo, hi in [(6850, 6940), (7160, 7330), (7570, 7700), (8120, 8350), (8900, 9700)]:
        ax.axvspan(lo, hi, color='gray', alpha=.08)
    ax.set(xlabel='Vacuum heliocentric wavelength (Å)', ylabel='Flux density (10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹)',
           title='16 Cyg: unsaturated central-slice star, 1 s frames (not telluric-corrected; slit losses not included)')
    ax.legend(ncol=4)
    fig.tight_layout(); fig.savefig(OUT/'16CYG_STAR.png', dpi=150); plt.close(fig)
    for name in ['spatial_profiles.png', 'outer_over_central_ratio.png']:
        if (env.BASE/'audit/z16'/name).exists():
            shutil.copy2(env.BASE/'audit/z16'/name, OUT/name)
    (OUT/'z16_report.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    if sys.argv[1] == 'reduce':
        for ch in sys.argv[2:]:
            reduce(ch)
    else:
        products(sys.argv[sys.argv.index('--standard')+1] if '--standard' in sys.argv else None)
