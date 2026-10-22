"""Archive acquisition and reproducible review of the requested H5MUSES-327.

Run with /opt/anaconda3/bin/python; fetch is read-only against public archives.
The saved W1 model is transformed, never refitted. This is a catalogue addition,
not an automatic change to either observing queue.
"""
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import base64
import hashlib
import importlib
import io
import json
import os
import sys

import numpy as np
import pandas as pd
import requests
from astropy.io import ascii
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'data/additional_targets/H5MUSES-327'
NAME = 'H5MUSES-327'
RA, DEC = 261.137, 59.44636111111111
CACHE = ROOT / 'data/reselection_2026-09-20'
os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-additional-matplotlib')


def save_json(path, value):
    # NumPy scalar values and absent measurements must remain valid JSON.
    native = importlib.import_module('16_candidate_webpage').native
    path.write_text(json.dumps(native(value), indent=2, allow_nan=False)+'\n')


def date(mjd):
    return Time(float(mjd), format='mjd').strftime('%Y-%m-%d')


def fetch():
    DEST.mkdir(parents=True, exist_ok=True)
    base = dict(spatial='Cone', objstr=f'{RA} {DEC}', radius=3,
                radunits='arcsec', outfmt=1)
    gator = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query'
    sql = f'''SELECT p.objID,p.ra,p.dec,p.type,p.psfMag_g,p.psfMag_r,
        p.psfMag_i,p.modelMag_g,p.modelMag_r,p.modelMag_i,
        s.specObjID,s.plate,s.mjd,s.fiberID,s.run2d,s.z,s.zWarning,s.class,s.subClass
        FROM dbo.fGetNearbyObjEq({RA},{DEC},0.05) n
        JOIN PhotoObjAll p ON p.objID=n.objID
        LEFT JOIN SpecObjAll s ON s.bestObjID=p.objID'''
    jobs = {
        'ztf_objects.tbl': (gator, dict(base, catalog='ztf_objects_dr24',
            selcols='oid,ra,dec,filtercode,ngoodobs,medianmag,magrms,minmag,maxmag')),
        'neowise.tbl': (gator, dict(base, catalog='neowiser_p1bs_psd',
            selcols='ra,dec,mjd,w1mpro,w1sigmpro,w2mpro,w2sigmpro,qual_frame,cc_flags')),
        'allwise.tbl': (gator, dict(base, catalog='allwise_p3as_psd')),
        'ztf_snad.json': ('https://db.ztf.snad.space/api/v3/data/dr24/circle/full/json',
            dict(ra=RA, dec=DEC, radius_arcsec=3)),
        'sdss.json': ('https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch',
            dict(cmd=sql, format='json')),
        'sdss_spectra_cone.json': ('https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch',
            dict(cmd=f'SELECT s.specObjID,s.ra,s.dec,s.plate,s.mjd,s.fiberID,s.run2d,s.z,s.zWarning,s.class,s.subClass FROM SpecObjAll s JOIN dbo.fGetNearbySpecObjEq({RA},{DEC},0.05) n ON s.specObjID=n.specObjID', format='json')),
        'gaia.tbl': (gator, dict(base, catalog='gaia_dr3_source', radius=5,
            selcols='source_id,ra,dec,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,phot_g_mean_mag,bp_rp,ruwe')),
        'spitzer_ch0.tbl': ('https://irsa.ipac.caltech.edu/data/SPITZER/5MUSES/spectra/r24183552b_ch0.tbl', {}),
        'spitzer_ch2.tbl': ('https://irsa.ipac.caltech.edu/data/SPITZER/5MUSES/spectra/r24183552b_ch2.tbl', {}),
    }
    def one(item):
        filename, (url, params) = item
        path = DEST / filename
        if path.exists():
            return filename, dict(status='cached', url=requests.Request('GET',url,params=params).prepare().url,
                                  bytes=path.stat().st_size)
        try:
            r = requests.get(url, params=params, timeout=(20, 120))
            r.raise_for_status()
            if filename.endswith('.json'):
                r.json()
            elif 'stat="ERROR"' in r.text[:2000] or '<html' in r.text[:300].lower():
                raise ValueError(r.text[:300])
            path.write_bytes(r.content)
            return filename, dict(status='downloaded', url=r.url, bytes=len(r.content),
                queried_utc=datetime.now(timezone.utc).isoformat())
        except Exception as e:
            return filename, dict(status='failed', error=str(e))
    manifest_path = DEST / 'archive_queries.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    with ThreadPoolExecutor(max_workers=5) as ex:
        for future in as_completed([ex.submit(one, item) for item in jobs.items()]):
            name, result = future.result()
            if result['status'] != 'cached' or name not in manifest: manifest[name] = result
            manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
            print(name, result, flush=True)
    images = importlib.import_module('18_fetch_review_images')
    from types import SimpleNamespace
    for field, pixels in [(40, 256), (240, 384)]:
        print('SDSS cutout', images.one(SimpleNamespace(name=NAME, ra=RA, dec=DEC), field, pixels), flush=True)
    pd.DataFrame([dict(name=NAME, ra=RA, dec=DEC, r_planning=15.755)]).to_csv(DEST/'target.csv', index=False)


def geometry(ra=RA, dec=DEC, dest=DEST, visit_minutes=14):
    """20-second full-visit checks, including the existing standard reservations."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord, AltAz, get_body
    from astropy.utils import iers
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    obs = importlib.import_module('05_observability')
    coord = SkyCoord(ra*u.deg, dec*u.deg)
    standards = pd.read_csv(ROOT/'observing/oct26_27_2026/standards.csv')
    nights, windows = [], {}
    for night in ['oct26', 'oct27']:
        dusk, dawn, _, _ = obs.night_window(obs.NIGHTS[night][0], 'full')
        grid = np.arange(np.ceil(dusk.unix/20)*20, dawn.unix, 20)
        times = Time(grid, format='unix')
        frame = AltAz(obstime=times, location=obs.PALOMAR.location, pressure=0*u.hPa)
        aa = coord.transform_to(frame)
        x = np.where(aa.alt.deg > 0, aa.secz.value, np.inf)
        moon = aa.separation(get_body('moon', times, obs.PALOMAR.location).transform_to(frame)).deg
        duration = visit_minutes*60
        count = int(duration/20)+1  # endpoints inclusive
        xmax = np.lib.stride_tricks.sliding_window_view(x, count).max(axis=1)
        mmin = np.lib.stride_tricks.sliding_window_view(moon, count).min(axis=1)
        starts = grid[:len(xmax)]
        free = np.ones(len(starts), dtype=bool)
        for r in standards[standards.NIGHT == night].itertuples():
            free &= ~((starts < pd.Timestamp(r.END_PDT).timestamp()) &
                      (starts+duration > pd.Timestamp(r.START_PDT).timestamp()))
        def local(seconds):
            return pd.Timestamp(float(seconds), unit='s', tz='UTC').tz_convert('America/Los_Angeles')
        tier_ranges, audit = {}, {}
        for tier, limit in [('preferred', 1.5), ('extended', 1.8), ('fallback', 2.)]:
            valid = free & (xmax <= limit-.001) & (mmin >= 40.02)
            changes = np.diff(np.r_[False, valid, False].astype(int))
            ranges, checks = [], []
            for a, b in zip(np.where(changes == 1)[0], np.where(changes == -1)[0]):
                first, last = starts[a], starts[b-1]
                ranges.append(dict(start_utc=Time(first, format='unix').isot,
                    end_utc=Time(last+duration, format='unix').isot,
                    start=local(first).strftime('%H:%M'), end=local(last+duration).strftime('%H:%M'),
                    minutes=float((last+duration-first)/60)))
                checks.append(dict(first_start_pdt=local(first).isoformat(),
                    last_start_pdt=local(last).isoformat(), last_end_pdt=local(last+duration).isoformat(),
                    min_moon_deg=float(mmin[a:b].min()), max_airmass=float(xmax[a:b].max())))
            tier_ranges[tier], audit[str(limit)] = ranges, checks
        good = (x >= 1) & (x <= 1.8) & (moon >= 40)
        nights.append(dict(night=night, tier_ranges=tier_ranges, ranges=tier_ranges['extended'],
            longest_minutes=max((r['minutes'] for r in tier_ranges['extended']), default=0),
            min_airmass=float(x[good].min()) if good.any() else None,
            moon_min=float(moon[good].min()) if good.any() else None,
            moon_max=float(moon[good].max()) if good.any() else None, minutes_airmass_le1p3=float(((x<=1.3)&good).sum()/3),
            curve=[[float(times[j].mjd), float(x[j]) if x[j]<4 else None, float(moon[j])]
                   for j in range(0, len(times), 15)]))
        windows[night] = audit
    save_json(dest/'visibility.json', windows)
    return nights


def analyze():
    """Use one ZTF object per band; retain all other matches as separate archives."""
    raw = json.loads((DEST/'ztf_snad.json').read_text())
    summary, ztf, weekly, selected_frames = {}, {}, {}, []
    for band in ['g', 'r', 'i']:
        matches = [(oid, r) for oid, r in raw.items() if r['meta']['filter'] == 'z'+band]
        oid, record = max(matches, key=lambda item: len(item[1]['lc']))
        frame = pd.DataFrame(record['lc'])
        frame = frame[np.isfinite(frame.mag) & np.isfinite(frame.mjd) & (frame.magerr>0)]
        frame = frame.drop_duplicates(['mjd', 'mag', 'magerr']).sort_values('mjd')
        frame['filtercode'], frame['oid'] = 'z'+band, oid
        selected_frames.append(frame)
        def bins(days):
            return frame.groupby(np.floor(frame.mjd/days)).agg(
                mjd=('mjd', 'median'), mag=('mag', 'median'), err=('magerr', 'median'), n=('mag', 'size'))
        w = bins(7)
        weekly[band] = w
        ztf[band] = w[['mjd', 'mag', 'err']].round(5).values.tolist()
        seasonal = bins(90)
        seasonal = seasonal[seasonal.n>=5]
        seasonal.to_csv(DEST/f'ztf_{band}_90day.csv', index=False)
        first = frame[frame.mjd <= frame.mjd.min()+365].mag.median()
        last = frame[frame.mjd >= frame.mjd.max()-365].mag.median()
        summary[band] = dict(oid=oid, n=len(frame), first_date=date(frame.mjd.min()),
            last_date=date(frame.mjd.max()), median_mag=frame.mag.median(),
            latest180_mag=frame[frame.mjd>=frame.mjd.max()-180].mag.median(),
            scatter_mag=frame.mag.std(), median_error_mag=frame.magerr.median(),
            p95_minus_p5_mag=frame.mag.quantile(.95)-frame.mag.quantile(.05),
            seasonal_range_mag=seasonal.mag.max()-seasonal.mag.min(),
            last_minus_first_year_mag=last-first)
    pd.concat(selected_frames).to_csv(DEST/'ztf_selected.csv', index=False)
    neo = ascii.read(DEST/'neowise.tbl', format='ipac').to_pandas().sort_values('mjd')
    # Cluster physical visits by >50-day gaps, rather than splitting visits at
    # arbitrary calendar-bin edges. Assess quality independently in each band.
    visits, ir = {}, {}
    for b in [1, 2]:
        flags = neo.cc_flags.astype(str).str.strip().str.zfill(4)
        f = neo[(neo.qual_frame>0) & flags.str[b-1].eq('0') &
                np.isfinite(neo[f'w{b}mpro']) & np.isfinite(neo[f'w{b}sigmpro']) &
                neo[f'w{b}sigmpro'].gt(0)].copy()
        f['visit'] = (f.mjd.diff().fillna(0)>50).cumsum()
        v = f.groupby('visit').agg(mjd=('mjd','median'), mag=(f'w{b}mpro','median'),
            err=(f'w{b}mpro', lambda s: 1.2533*s.std(ddof=1)/np.sqrt(len(s))), n=('mjd','size'))
        v = v[v.n>=3]
        v.to_csv(DEST/f'neowise_w{b}_visits.csv', index=False)
        visits[b] = v
        ir[f'w{b}'] = dict(n_frames=len(f), n_visits=len(v), first_date=date(v.mjd.min()),
            last_date=date(v.mjd.max()), visit_range_mag=v.mag.max()-v.mag.min(),
            p90_minus_p10_mag=v.mag.quantile(.9)-v.mag.quantile(.1),
            last_minus_first_visit_mag=v.mag.iloc[-1]-v.mag.iloc[0],
            typical_visit_stat_error_mag=v.err.median())
    projection = pd.read_csv(CACHE/'wise_h5muses327/projected.csv').iloc[0].to_dict()
    unwise = pd.concat([pd.read_parquet(p) for p in (CACHE/'wise_h5muses327').glob('pixel_*.parquet')])
    unwise = unwise.drop_duplicates(['time','flux','err']).sort_values('time')
    unwise.to_csv(DEST/'unwise_w1.csv', index=False)
    projection['largest_gap_days'] = float(unwise.time.diff().max())
    projection['model_sha256'] = hashlib.sha256((ROOT/'data/umap_w1_model.pkl').read_bytes()).hexdigest()
    # Same labelled reference cells and nearest-neighbour definition as the pool.
    from scipy.spatial import cKDTree
    reference = pd.read_csv(ROOT/'data/sampleA_embedding_objectid.csv')
    _, xe, ye = np.histogram2d(reference.umap_x, reference.umap_y, bins=10)
    px, py = [int(np.clip(np.searchsorted(edges[1:], projection[col]),0,9))
              for edges,col in [(xe,'umap_x'),(ye,'umap_y')]]
    ix = np.clip(np.searchsorted(xe[1:],reference.umap_x),0,9)
    iy = np.clip(np.searchsorted(ye[1:],reference.umap_y),0,9)
    regions=[]
    for col,label in [('in_region_clagn','literature'),('in_region_zeltyn','Zeltyn')]:
        cells=set(zip(ix[reference[col].fillna(False)],iy[reference[col].fillna(False)]))
        if (px,py) in cells: regions.append(label)
    dist, idx = cKDTree(reference[['umap_x','umap_y']]).query(
        [projection['umap_x'],projection['umap_y']], k=51)
    idx = idx[1:] if dist[0]<1e-6 else idx[:50]
    bits = reference.label_bits.to_numpy().astype(int)[idx]
    projection.update(region=' + '.join(regions) if regions else 'outside selected',
        on_fraction=float(((bits&16)>0).mean()), off_fraction=float(((bits&32)>0).mean()))
    gaia = ascii.read(DEST/'gaia.tbl',format='ipac').to_pandas().iloc[0].to_dict()
    gaia['source_id'] = str(int(gaia['source_id']))
    gaia['parallax_snr'] = gaia['parallax']/gaia['parallax_error']
    result = dict(name=NAME, ra=RA, dec=DEC, ztf=summary, neowise=ir, manifold=projection,
        gaia=gaia, interpretation='Likely foreground stellar optical counterpart; AGN identity unconfirmed.',
        ztf_quality='SNAD DR24 returned light curves; use longest single object per filter. Other field IDs retained separately; no inter-field offsets mixed.',
        variability_caution='Observed changes, not a confirmed intrinsic-variability classification. WISE visit errors are statistical; calibration, scan-angle and unresolved-blend effects remain possible.',
        manifold_caution='Original saved-model projection. A 2010–2014 gap allows the zero-mean GP to revert toward zero; GP-based Fvar is not the observed variability amplitude. Location and neighbour fractions do not establish AGN identity.')
    save_json(DEST/'summary.json', result)
    nights = geometry()
    build_payload(result, ztf, visits, unwise, nights)
    plot(result, weekly, visits, unwise, reference)
    report(result)
    print(json.dumps(result, default=str, indent=2), flush=True)


def build_payload(s, ztf, visits, unwise, nights):
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    coord = SkyCoord(RA*u.deg, DEC*u.deg)
    g, r = s['ztf']['g'], s['ztf']['r']
    m = s['manifold']
    science = dict(reference_private=False, ztf_r_latest180_mag=r['latest180_mag'],
        ztf_g_latest180_mag=g['latest180_mag'], ztf_last_date=r['last_date'],
        state_classification='Likely foreground stellar counterpart; AGN identity unconfirmed',
        science_question='Establish optical source identity and assess the infrared variability; no optical redshift or SDSS spectrum was returned by the archive queries.',
        field_notes='Gaia DR3 parallax 0.4935 ± 0.0277 mas; significant proper motion and RUWE 0.914. Projection is for comparison only; the GP interpolates across a long W1 gap. The unWISE blending-diagnostic query timed out; an unresolved blend is not excluded.')
    v = visits[1]; flux = 309.54e3*10**(-.4*v.mag.to_numpy())
    screen_path=CACHE/'neighbour_cache'/f'{NAME}_screen.json'
    screen=json.loads(screen_path.read_text()) if screen_path.exists() else dict(status='pending', flags=[])
    screen['flags'] = list(screen.get('flags',[]))+['Likely stellar optical counterpart; verify identity before AGN follow-up']
    screen['status'] = 'review'
    manifest = CACHE/'neighbour_unwise_manifest.csv'
    if manifest.exists():
        uw=pd.read_csv(manifest);uw=uw[uw.name==NAME]
        if len(uw):screen['unwise']=uw.iloc[0].to_dict()
    image_path=ROOT/'data/cutouts'/f'{NAME}_sdss_wide.jpg'
    image_meta=json.loads(image_path.with_suffix('.json').read_text())
    target=dict(name=NAME,jname='SDSS J172432.88+592646.9',aliases=['5MUSES-327','5MUSES 327'],
        ra=RA,dec=DEC,sex=coord.to_string('hmsdms',sep=':',precision=2),z=None,
        rmag=r['median_mag'],r_source='ZTF DR24 median; not a current measurement',
        ux=m['umap_x'],uy=m['umap_y'],region=m['region'],on_fraction=m['on_fraction'],off_fraction=m['off_fraction'],
        status='unmatched',known=[],n_spec=0,n_spec_available_dates=0,epochs=[],nights=nights,
        ztf=ztf,ztf_status=dict(status='available',collection='ZTF DR24 via SNAD',
            first_date=r['first_date'],last_date=r['last_date'],n_g=g['n'],n_r=r['n']),
        wise={'W1':unwise[['time','flux','err']].values.tolist()},
        neo=np.column_stack([v.mjd,flux,flux*.921034*v.err]).tolist(),
        pool_role='requested',prepared_nights=['oct26','oct27'],field_status='review',neighbour_screen=screen,
        spec=None,lines=[],science=science,exposure_plans={},image_status='available',
        cut=dict(source='SDSS DR18',field_arcsec=image_meta['field_arcsec'],
            image='data:image/jpeg;base64,'+base64.b64encode(image_path.read_bytes()).decode()),
        archive_review=dict(summary=f"ZTF {r['first_date']}–{r['last_date']}: g/r scatter {g['scatter_mag']:.3f}/{r['scatter_mag']:.3f} mag. NEOWISE W1/W2 visit ranges {s['neowise']['w1']['visit_range_mag']:.3f}/{s['neowise']['w2']['visit_range_mag']:.3f} mag; intrinsic origin not established. Gaia strongly favors a stellar optical counterpart. No SDSS DR18 spectrum within 3 arcsec. Spitzer IRS SL/LL spectra are available (2008-01-14).",
            note=s['manifold_caution'],sources=[
                dict(label='ZTF r light curve',url='https://ztf.snad.space/dr24/view/825202200007543'),
                dict(label='Spitzer IRS short-low spectrum',url='https://irsa.ipac.caltech.edu/data/SPITZER/5MUSES/spectra/r24183552b_ch0.tbl'),
                dict(label='Spitzer IRS long-low spectrum',url='https://irsa.ipac.caltech.edu/data/SPITZER/5MUSES/spectra/r24183552b_ch2.tbl')]))
    save_json(DEST/'target_payload.json',target)


def plot(s, weekly, visits, unwise, reference):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), constrained_layout=True)
    for band,color in [('g','#268b54'),('r','#c34e3e'),('i','#7854a2')]:
        f=weekly[band]
        axes[0].errorbar(Time(f.mjd.to_numpy(),format='mjd').decimalyear,f.mag,yerr=f.err,
                        fmt='.',ms=3,alpha=.65,color=color,label=f'ZTF {band}, 7-day medians')
    axes[0].invert_yaxis();axes[0].set_ylabel('Optical magnitude (AB)');axes[0].legend(ncol=3,fontsize=8)
    for b,color in [(1,'#286fa1'),(2,'#c77a21')]:
        v=visits[b]
        axes[1].errorbar(Time(v.mjd.to_numpy(),format='mjd').decimalyear,v.mag,yerr=v.err,
                        fmt='o-',ms=3,lw=.7,color=color,label=f'NEOWISE W{b} visits')
    wmag=23.9-2.699-2.5*np.log10(unwise.flux.to_numpy()*1000)
    axes[1].plot(Time(unwise.time.to_numpy(),format='mjd').decimalyear,wmag,'.',color='#222222',label='unWISE W1; manifold input')
    axes[1].invert_yaxis();axes[1].set_ylabel('Infrared magnitude (Vega)');axes[1].legend(ncol=3,fontsize=8)
    for ax in axes[:2]:ax.set_xlabel('Year');ax.grid(alpha=.2)
    axes[2].scatter(reference.umap_x,reference.umap_y,s=2,c='#999999',alpha=.4)
    m=s['manifold'];axes[2].scatter([m['umap_x']],[m['umap_y']],s=80,c='#d64631',marker='*',label=NAME)
    axes[2].set(xlabel='W1 UMAP 1',ylabel='W1 UMAP 2');axes[2].legend()
    fig.suptitle(f'{NAME}: requested comparison target; likely stellar optical counterpart\nObserved photometric changes do not establish an AGN transition',fontsize=12)
    for suffix in ['png','pdf']:fig.savefig(DEST/f'archive_review.{suffix}',dpi=170)
    plt.close(fig)


def report(s):
    lines=[f'# {NAME}', '',
        'Requested catalogue/manifold addition, October 22, 2026. Scheduling is maintained separately in observing/oct26_27_2026/all_80_targets.csv by october_additions.py.', '',
        'RA 17:24:32.88, Dec +59:26:46.9 (J2000); 261.137°, +59.4463611111°.', '',
        '**Identity:** Gaia DR3 1434108585365868160 matches within 0.105 arcsec. '
        'Parallax 0.4935 ± 0.0277 mas (17.8 sigma), proper motion '
        '(+3.763 ± 0.039, −1.294 ± 0.038) mas/yr, RUWE 0.914. '
        'This strongly favors a stellar optical counterpart; a blended infrared association is not excluded. '
        'SDSS is point-like and returns no DR18 spectrum within 3 arcsec. Redshift is unknown.', '',
        '## Variability', '', '| Band | Dates | Measurements | Scatter / visit range |',
        '| --- | --- | ---: | ---: |']
    for band,r in s['ztf'].items():
        lines.append(f"| ZTF {band} | {r['first_date']}–{r['last_date']} | {r['n']} | {r['scatter_mag']:.3f} mag RMS |")
    for band,r in s['neowise'].items():
        lines.append(f"| NEOWISE {band.upper()} | {r['first_date']}–{r['last_date']} | {r['n_visits']} visits | {r['visit_range_mag']:.3f} mag peak-to-peak |")
    lines += ['', 'The optical changes are small, with no large sustained brightening or fading. '
        'The IR visit medians show modest changes. Statistical errors do not include calibration, '
        'scan-angle or unresolved-blend systematics; intrinsic variability remains unconfirmed. '
        'The unWISE blending-diagnostic query timed out. SDSS/Gaia neighbour queries returned no '
        'nearby-source flag, but that does not exclude an unresolved blend.', '',
        'The ZTF analysis uses one field/object ID per band, selected by returned coverage. '
        'The SNAD mirror returns fewer measurements than the IRSA object-catalogue counts; '
        'the analysis uses the actual returned curves, not those catalogue counts.', '',
        '## Observability', '', 'Palomar, 14-minute visits, Moon ≥40°, airmass ≤1.8, '
        '20-second checks, existing opening/closing standards excluded. Times are PDT. '
        'These are feasible starts, not free queue slots; consult the current queue for the assigned visit. '
        'Geometry uses the project’s Palomar location, unrefracted topocentric coordinates, '
        'the built-in ephemeris and bundled IERS tables.', '',
        '| Night | First start | Latest start | Latest finish |', '| --- | --- | --- | --- |']
    windows=json.loads((DEST/'visibility.json').read_text())
    for night,tiers in windows.items():
        for w in tiers['1.8']:
            lines.append(f"| {night} | {w['first_start_pdt'][11:19]} | {w['last_start_pdt'][11:19]} | {w['last_end_pdt'][11:19]} |")
    m=s['manifold']
    lines += ['', '## Manifold', '',
        f"Saved-model transform: ({m['umap_x']:.6f}, {m['umap_y']:.6f}); {m['region']} region. "
        f"{m['n_w1']} unWISE W1 measurements; {m['gp_extrapolated_fraction']:.1%} of the 4000-day input grid lies beyond the measured span.", '',
        s['manifold_caution'], '', '## Files and provenance', '',
        '- `archive_queries.json`: public query URLs and timestamps.',
        '- `summary.json`, `visibility.json`: derived quantities and geometry.',
        '- `ztf_selected.csv`, `ztf_*_90day.csv`, `neowise_w*_visits.csv`, `unwise_w1.csv`: analysis inputs and summaries.',
        '- `spitzer_ch0.tbl`, `spitzer_ch2.tbl`: public IRS spectra, observed January 14, 2008.',
        '- `archive_review.png` / `.pdf`: light curves and manifold location.',
        '- `target_payload.json`: persistent addition loaded by the observer-page generator.', '',
        'Rebuild: run `63_5muses327.py fetch`; run `15c_project_expanded_pool.py fetch '
        '--targets data/additional_targets/H5MUSES-327/target.csv --batch wise_h5muses327`; '
        'then `15c_project_expanded_pool.py project --batch wise_h5muses327`, '
        '`63_5muses327.py analyze` and `51_observer_page.py`. '
        'The model requires NumPy <2.3 with the installed Numba; use an isolated compatible environment. '
        'No refitting, commit, upload or publication is performed.']
    (DEST/'README.md').write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['fetch', 'analyze'])
    args = parser.parse_args()
    fetch() if args.stage == 'fetch' else analyze()
