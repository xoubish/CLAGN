"""
03e_blending.py  --  is the W1 photometry of each candidate contaminated by neighbours inside the WISE beam?

WISE W1 has a 6.1" FWHM PSF on 2.75" unWISE pixels, so anything within ~5-8" blends into the light curve. Two
independent measures per target:
  * unWISE DR1 object catalogue (Schlafly+2019) via NOIRLab Data Lab TAP: fracflux_w1 / fracflux_w2 = fraction of the
    flux at the target's position that the deblender attributes to the target itself (1 = clean); flags; number of
    other unWISE sources within 10".
  * Pan-STARRS1 DR2 mean catalogue via MAST: neighbours within 8" with nDetections >= 3; separation and z-band
    magnitude difference to the target; a rough W1 contamination estimate assuming neighbours are 1.5 mag bluer in
    z - W1 (AB) than a quasar.
Usage: /opt/anaconda3/bin/python 03e_blending.py <in.csv> <tag> [ra_col dec_col name_col]
Output: data/blending_<tag>.csv with name, fracflux_w1, fracflux_w2, flags_unwise_w1, n_unwise_10as, n_ps1_8as,
        nbr_min_sep_as, nbr_min_dz, w1_contam_est, blend_flag
"""
import io, os, sys, time
import numpy as np
import pandas as pd
import requests
import pyvo
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, 'data')
TAP = pyvo.dal.TAPService('https://datalab.noirlab.edu/tap')
PS1 = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv'
BOX = 0.0028          # deg half-width (~10") for the unWISE box


def unwise(targets):
    out = []
    for i in range(0, len(targets), 20):
        sub = targets.iloc[i:i + 20]
        cond = ' OR '.join(f"(ra BETWEEN {r.ra - BOX / np.cos(np.radians(r.dec)):.6f} AND {r.ra + BOX / np.cos(np.radians(r.dec)):.6f} AND dec BETWEEN {r.dec - BOX:.6f} AND {r.dec + BOX:.6f})"
                           for r in sub.itertuples())
        # note: 'primary' is a reserved word for the ADQL parser, so it cannot be used in the WHERE clause; duplicates
        # from overlapping coadd tiles are handled by taking the nearest source below
        q = f"SELECT ra, dec, flux_w1, fracflux_w1, fracflux_w2, flags_unwise_w1, flags_info_w1, nm_w1 FROM unwise_dr1.object WHERE {cond}"
        for attempt in range(3):
            try:
                t = TAP.search(q).to_table().to_pandas(); break
            except Exception as e:
                print(f'   unWISE batch {i//20}: {str(e)[:80]} (attempt {attempt+1})', flush=True); time.sleep(10); t = None
        if t is None or not len(t):
            continue
        for r in sub.itertuples():
            d = np.hypot((t.ra - r.ra) * np.cos(np.radians(r.dec)), t.dec - r.dec) * 3600
            near = t[d < 3.0]
            n10 = int(((d >= 3.0) & (d < 10.0)).sum())
            if len(near):
                j = near.index[np.argmin(d[near.index])]
                out.append(dict(name=r.name, fracflux_w1=float(t.at[j, 'fracflux_w1']), fracflux_w2=float(t.at[j, 'fracflux_w2']),
                                flags_unwise_w1=int(t.at[j, 'flags_unwise_w1']), nm_w1=int(t.at[j, 'nm_w1']), n_unwise_10as=n10))
            else:
                out.append(dict(name=r.name, fracflux_w1=np.nan, fracflux_w2=np.nan, flags_unwise_w1=-1, nm_w1=0, n_unwise_10as=n10))
        print(f'   unWISE {min(i+20, len(targets))}/{len(targets)}', flush=True)
    return pd.DataFrame(out)


def ps1_neighbours(row):
    name, ra, dec = row['name'], float(row['ra']), float(row['dec'])
    res = dict(name=name, n_ps1_8as=np.nan, nbr_min_sep_as=np.nan, nbr_min_dz=np.nan, w1_contam_est=np.nan, target_z_ps1=np.nan)
    for attempt in range(3):
        try:
            r = requests.get(PS1, params={'ra': ra, 'dec': dec, 'radius': 8 / 3600, 'nDetections.gte': 3,
                                          'columns': '[raMean,decMean,zMeanPSFMag,rMeanPSFMag,nDetections]'}, timeout=120)
            r.raise_for_status(); break
        except Exception:
            time.sleep(5); r = None
    if r is None or not r.text.strip():
        return res
    t = pd.read_csv(io.StringIO(r.text))
    if not len(t):
        res['n_ps1_8as'] = 0; return res
    t['sep'] = np.hypot((t.raMean - ra) * np.cos(np.radians(dec)), t.decMean - dec) * 3600
    t.loc[(t.zMeanPSFMag <= 0) | (t.zMeanPSFMag > 30), 'zMeanPSFMag'] = np.nan     # PS1 uses -999 for missing magnitudes
    tgt = t[t.sep < 1.0].sort_values('sep').head(1)
    nbr = t[t.sep >= 1.0]
    tz = float(tgt.zMeanPSFMag.iloc[0]) if len(tgt) and pd.notna(tgt.zMeanPSFMag.iloc[0]) else np.nan
    res['target_z_ps1'] = tz; res['n_ps1_8as'] = int(len(nbr))
    if len(nbr):
        res['nbr_min_sep_as'] = float(nbr.sep.min())
        dz = (nbr.zMeanPSFMag - tz).values.astype(float)
        dz = dz[np.isfinite(dz)]
        if len(dz):
            res['nbr_min_dz'] = float(dz.min())
            # neighbour W1 flux relative to the quasar: 10^(-0.4 dz) x 10^(-0.4 x 1.5) (stars/galaxies ~1.5 mag bluer in z-W1)
            frac = np.sum(10 ** (-0.4 * (dz + 1.5)))
            res['w1_contam_est'] = float(frac / (1 + frac))
    return res


def main():
    inp, tag = sys.argv[1], sys.argv[2]
    ra_col, dec_col, name_col = (sys.argv[3:6] if len(sys.argv) >= 6 else ('ra', 'dec', 'name'))
    t = pd.read_csv(inp).drop_duplicates(name_col)
    targets = pd.DataFrame({'name': t[name_col].astype(str), 'ra': t[ra_col].astype(float), 'dec': t[dec_col].astype(float)}).reset_index(drop=True)
    t0 = time.time()
    u = unwise(targets); print(f'[{time.time()-t0:.0f}s] unWISE done: {len(u)} rows', flush=True)
    with ThreadPoolExecutor(max_workers=6) as ex:
        p = pd.DataFrame(list(ex.map(ps1_neighbours, [r for _, r in targets.iterrows()])))
    print(f'[{time.time()-t0:.0f}s] PS1 done', flush=True)
    out = targets.merge(u, on='name', how='left').merge(p, on='name', how='left')
    # two kinds of "blend": a separate PS1 neighbour inside the beam (another object's light, possibly variable) versus
    # a low unWISE fraction with no PS1 neighbour, which is the target's own extended host being split by the deblender
    # (constant light: dilutes the fractional variability but cannot fake a change)
    major = (out.n_ps1_8as > 0) & ((out.w1_contam_est > 0.15) | (out.fracflux_w1 < 0.8))     # neighbour matters for W1
    minor = (out.n_ps1_8as > 0) & (out.nbr_min_dz < 1.5) & ~major                          # neighbour present, < 15% of W1
    out['blend_kind'] = np.select([major, out.fracflux_w1 < 0.8, minor], ['neighbour', 'extended host', 'minor neighbour'], 'clean')
    out['blend_flag'] = out.blend_kind == 'neighbour'
    out.to_csv(os.path.join(DATA, f'blending_{tag}.csv'), index=False)
    print(f'wrote data/blending_{tag}.csv: {len(out)} targets | unWISE fracflux_w1 < 0.8: {int((out.fracflux_w1 < 0.8).sum())} | < 0.5: {int((out.fracflux_w1 < 0.5).sum())} '
          f'| PS1 neighbour within 8": {int((out.n_ps1_8as > 0).sum())} | neighbour blends (flagged): {int(out.blend_flag.sum())} | extended-host only: {int((out.blend_kind == "extended host").sum())}')


if __name__ == '__main__':
    main()
