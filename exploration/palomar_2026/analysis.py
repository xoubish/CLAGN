"""Read-only census and reproducible screening of the completed Palomar campaign.

All outputs are local. Line indices are pseudo-continuum diagnostics, never fitted
broad-line fluxes. No missing errors are imputed and no CL classifications are made.
"""
from pathlib import Path
import hashlib
import json
import sys
from datetime import datetime
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.io import fits
import astropy.units as u

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'outputs'
sys.path.insert(0, str(ROOT / 'scripts'))
from spectral_utils import centered_record, quality

NIGHTS = {'sep23': 'sep23_p330e', 'oct26': 'oct26_feige34', 'oct27': 'oct27_bd284211'}
LINES = {
    'Hb': dict(center=4862.68, sides=[(4750, 4790), (5100, 5140)],
               wings=[(4800, 4848), (4878, 4930)], panel=(4680, 5160)),
    'Ha': dict(center=6564.61, sides=[(6400, 6450), (6750, 6800)],
               wings=[(6460, 6520), (6610, 6690)], panel=(6340, 6840)),
    'MgII': dict(center=2798.75, sides=[(2650, 2700), (2930, 2970)],
                 wings=[(2720, 2770), (2825, 2880)], panel=(2600, 3010)),
}
VARIANTS = {'nominal': 0., 'blue_sides': -15., 'red_sides': 15.}
SOURCES = {}


def track(path):
    path = Path(path)
    key = str(path.relative_to(ROOT))
    if key not in SOURCES:
        SOURCES[key] = dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(), bytes=path.stat().st_size)
    return path


def read_csv(path, **kwargs):
    return pd.read_csv(track(path), **kwargs)


def read_json(path):
    return json.loads(track(path).read_text())


def canonical(name, jwst_names):
    if name.startswith('J') and name[1:] in jwst_names:
        return name[1:]
    return {'P8828': 'P8823', 'J085337.26014303.6': 'J085337.26+014303.6'}.get(name, name)


def census():
    jwst = read_csv(ROOT/'jwst_proposal/inputs/jwst_sample_cycle6.csv')
    jnames = set(jwst.internal_id)
    payload = read_json(ROOT/'data/reselection_2026-09-20/candidate_payload_local.json')
    targets = {t['name']: t for t in payload['targets']}
    for p in sorted((ROOT/'data/additional_targets').glob('*/target_payload.json')):
        t = read_json(p)
        extra = p.with_name('archive_payload.json')
        if extra.exists():
            t.update(read_json(extra))
        targets[t['name']] = t
    schedule = read_csv(ROOT/'docs/october/all_80_targets.csv')
    for r in schedule.to_dict('records'):
        name = canonical(r['NAME'], jnames)
        targets.setdefault(name, {}).update(schedule_reason=r['SCIENCE_QUESTION'])
        if 'z' not in targets[name]:
            targets[name]['z'] = r['Z']
    for r in jwst.to_dict('records'):
        targets.setdefault(r['internal_id'], {}).update(z=r['z'], ra=r['ra'], dec=r['dec'],
            jwst_id=r['id'], jwst_family=r['family'], jwst_transition=r['transition'])
    parent = read_csv(ROOT/'data/master_list_scored.csv', low_memory=False).dropna(subset=['ra', 'dec'])
    pc = SkyCoord(parent.ra.to_numpy()*u.deg, parent.dec.to_numpy()*u.deg)
    visits = []
    objects = {}
    for night, folder in NIGHTS.items():
        path = ROOT/'docs/observed'/folder/'manifest.json'
        coordfile = 'observed_targets.csv' if night == 'sep23' else f'observed_{night}_targets.csv'
        coordinates = read_csv(ROOT/'docs/targets'/coordfile)
        coords = {canonical(r['NAME'], jnames): r for r in coordinates.to_dict('records')}
        for v in read_json(path):
            name = v['name']
            t = targets.get(name, {})
            co = coords[name]
            sky = SkyCoord(co['RA'], co['DECL'], unit=(u.hourangle, u.deg))
            idx, sep, _ = sky.match_to_catalog_sky(pc)
            aliases = [name]
            if name in jnames:
                aliases.append(jwst.set_index('internal_id').loc[name, 'target'])
            if sep.arcsec.item() < 1.5:
                aliases.append(str(parent.iloc[int(idx)]['name']))
            is_star = name in ['Z16CygABC', 'H5MUSES-327']
            z = t.get('z', np.nan)
            if is_star:
                z = 0.
            objects[name] = dict(name=name, z=z, ra=sky.ra.deg, dec=sky.dec.deg,
                metadata_available=bool(t), aliases=';'.join(dict.fromkeys(aliases)),
                population='stellar/special' if is_star else ('nearby-galaxy special' if name.startswith('H') else 'AGN sample'),
                region=t.get('region', 'unassigned'), pool_role=t.get('pool_role', 'unassigned'),
                umap_x=t.get('ux', np.nan), umap_y=t.get('uy', np.nan),
                known_catalog_matches=len(t.get('known') or []),
                jwst_id=t.get('jwst_id', ''), jwst_family=t.get('jwst_family', ''),
                prior_question=t.get('schedule_reason') or (t.get('science') or {}).get('science_question', ''),
                association_warning=(t.get('ztf_status') or {}).get('association_warning', False))
            visits.append(dict(name=name, night=night, mjd=v['start_mjd'], end_mjd=Time(datetime.fromisoformat(v['end_utc'])).mjd,
                exposures=v['exposures'], exposure_seconds=v['exposure_seconds'], standard=v['flux_standard'],
                spectrum_path=str((path.parent/v['csv']).relative_to(ROOT)),
                **{f'snr_{ch}': (v.get('snr') or {}).get(ch, np.nan) for ch in 'UGRI'}))
    return pd.DataFrame(objects.values()).sort_values('name'), pd.DataFrame(visits), targets, jwst


def archived_spectra(obj, native_index):
    """Keep one reduction per survey/day; retain explicit rejections in audit.

    Bad SDSS-V field flags propagate across duplicate releases of that field/day.
    Short DESI coadds retain their start/end; long/unknown intervals are excluded.
    """
    raw = []
    audit = []
    paths = set()
    for name in obj.aliases.split(';'):
        for folder in ['data/spectra_dl', 'data/reselection_2026-09-20/sdssv_spectra']:
            p = ROOT/folder/f'{name}.json'
            if p.exists() and p not in paths:
                paths.add(p)
                for i, r in enumerate(read_json(p)):
                    raw.append((centered_record(r), str(p.relative_to(ROOT)), i))
    bad = {(str(r.get('source', 'SDSS')), int(np.floor(r['mjd']))) for r, _, _ in raw
           if r.get('mjd') and (r.get('fieldquality') == 'bad' or r.get('metadata_quality_ok') is False)}
    selected = {}
    for r, path, index in raw:
        source = r.get('source', 'SDSS')
        mjd = r.get('mjd')
        w = np.asarray(r.get('wave', []), float)
        f = np.asarray(r.get('flux', []), float)
        reason = ''
        if mjd is None or not np.isfinite(mjd) or not 40000 < mjd < 100000:
            reason = 'missing_or_invalid_date'
        elif (source, int(np.floor(mjd))) in bad:
            reason = 'bad_metadata_or_duplicate_bad_field'
        elif len(w) != len(f) or len(w) < 20 or not np.all(np.diff(w) > 0):
            reason = 'invalid_arrays'
        elif r.get('meta', {}).get('zwarning', 0) not in (None, 0, 0.):
            reason = 'redshift_warning'
        lo = hi = mjd
        if r.get('coadd'):
            lo, hi = r.get('min_mjd'), r.get('max_mjd')
            if source != 'DESI' or not r.get('date_verified') or lo is None or hi is None:
                reason = reason or 'coadd_date_unverified'
            elif not (np.isfinite(lo) and np.isfinite(hi) and 0 <= hi-lo <= 30):
                reason = reason or 'coadd_interval_long_or_invalid'
        record = dict(name=obj.name, path=path, record_index=index, source=source, mjd=mjd,
                      status='rejected' if reason else 'eligible', reason=reason)
        audit.append(record)
        if reason:
            continue
        key = (source, int(np.floor(hi)))
        if key in selected:
            old = selected[key]
            if quality(r) <= quality(old[0]):
                record.update(status='duplicate', reason='lower_ranked_same_survey_day')
                continue
            old[3].update(status='duplicate', reason='lower_ranked_same_survey_day')
        selected[key] = (r, path, index, record, lo, hi)
    spectra = []
    for key, (r, path, index, rec, lo, hi) in sorted(selected.items()):
        w, f = np.asarray(r['wave'], float), np.asarray(r['flux'], float)
        # Legacy `ok` is rounded n_good_in_bin / n_all_spectrum_pixels, often
        # identically zero. It is NOT a boolean mask. Flux nulls carry the mask.
        err = np.full(len(f), np.nan)
        native_path = find_native(obj, r) if key[0] == 'SDSS' else None
        native = native_index[(native_index.internal_id == obj.name) & (abs(native_index.mjd-float(r['mjd'])) < .01)]
        if native_path is not None:
            with fits.open(track(native_path), memmap=False) as hdus:
                d = hdus[1].data
                w, f = 10**np.asarray(d['loglam'], float), np.asarray(d['flux'], float)
                iv = np.asarray(d['ivar'], float)
                err = np.where(iv > 0, 1/np.sqrt(np.maximum(iv, 1e-300)), np.nan)
                f[~np.isfinite(err)] = np.nan
            path = str(native_path.relative_to(ROOT))
        elif len(native) and key[0] == 'SDSS':
            p = ROOT/'jwst_proposal'/native.iloc[0]['file']
            if p.exists():
                d = read_csv(p)
                w, f = d.wave_A.to_numpy(float), d.flux.to_numpy(float)
                iv = d.ivar.to_numpy(float)
                err = np.where(iv > 0, 1/np.sqrt(np.maximum(iv, 1e-300)), np.nan)
                f[~np.isfinite(err)] = np.nan
                path = str(p.relative_to(ROOT))
        record_z = r.get('meta', {}).get('z')
        z_bad = record_z is not None and np.isfinite(record_z) and abs(record_z-obj.z)/(1+obj.z) > .02
        if z_bad:
            rec.update(status='rejected', reason='catalog_redshift_disagrees')
            continue
        rec['status'] = 'accepted_screen'
        spectra.append(dict(name=obj.name, epoch_id=f'{obj.name}_{key[0]}_{key[1]}', source=key[0],
            night='archive', mjd=(lo+hi)/2, min_mjd=lo, max_mjd=hi, z=obj.z, path=path,
            proprietary=bool(r.get('proprietary')), error_kind='native diagonal' if np.isfinite(err).any() else 'missing',
            metadata_status='not_rejected; visual QA still required', wave=w, flux=f, error=err,
            quality=np.zeros(len(w), int), channel=np.full(len(w), 'archive'), alternate=None, uncorrected=None))
    return spectra, audit


def find_native(obj, record):
    """Only recover the exact URL cache or a unique field/day/version product."""
    url = record.get('url')
    for alias in obj.aliases.split(';'):
        if url:
            directory = ROOT/'data/spectra_cache'/alias
            p = directory/hashlib.sha256(url.encode()).hexdigest()[:20]/url.split('/')[-1]
            if p.exists():
                return p
            p = directory/url.split('/')[-1]
            if record.get('cache_identity_version') != 2 and p.exists():
                return p
        if record.get('proprietary'):
            version = record.get('archive_version', record.get('run2d', 'v6_2_1'))
            field = record.get('meta', {}).get('field')
            if field is not None:
                directory = ROOT/'data/reselection_2026-09-20/spectra_files'/version/alias
                matches = list(directory.glob(f"spec-{int(field):06d}-{int(record['mjd'])}-*.fits"))
                if len(matches) == 1:
                    return matches[0]
    return None


def load_spectra(objects, visits):
    native = read_csv(ROOT/'jwst_proposal/inputs/archival_spectra/index.csv')
    spectra, audit = [], []
    for obj in objects.itertuples(index=False):
        old, history = archived_spectra(obj, native)
        spectra.extend(old)
        audit.extend(history)
        for v in visits[visits.name == obj.name].itertuples(index=False):
            d = read_csv(ROOT/v.spectrum_path)
            # The special stellar spectrum uses the same standard columns.
            w, f = d.WAVE_VAC_HELIO_A.to_numpy(float), d.FLUX.to_numpy(float)
            err = d.FLUX_ERR.to_numpy(float)
            mask = d.MASK.astype(str).str.lower().isin(['true', '1']).to_numpy() if 'MASK' in d else np.ones(len(d), bool)
            good = mask & np.isfinite(err) & (err > 0)
            f[~good] = np.nan
            alternate = d.FLUX_ALTSTD.to_numpy(float) if 'FLUX_ALTSTD' in d else None
            unc = d.FLUX_UNCORR.to_numpy(float) if 'FLUX_UNCORR' in d else None
            for arr in [alternate, unc]:
                if arr is not None:
                    arr[~good] = np.nan
            spectra.append(dict(name=obj.name, epoch_id=f'{obj.name}_NGPS_{v.night}', source='NGPS', night=v.night,
                mjd=(v.mjd+v.end_mjd)/2, min_mjd=v.mjd, max_mjd=v.end_mjd, z=obj.z,
                path=v.spectrum_path, proprietary=False, error_kind='native diagonal; correlated pixels',
                metadata_status='delivered', wave=w, flux=f, error=err,
                quality=d.QUALITY.to_numpy(int) if 'QUALITY' in d else np.zeros(len(w), int),
                channel=d.CHANNEL.to_numpy() if 'CHANNEL' in d else np.full(len(w), 'stellar'),
                alternate=alternate, uncorrected=unc))
    return spectra, pd.DataFrame(audit)


def integration_weights(wave, good, windows):
    """Exact piecewise-linear integral on common endpoints, never over masked gaps.

    Returns weights on input samples plus fractional wavelength coverage. A
    missing run or a gap >2.5 times median sampling is never interpolated across.
    """
    w = np.asarray(wave, float)
    weights = np.zeros(len(w))
    valid = np.asarray(good, bool)
    step = np.median(np.diff(w))
    left, right = w[:-1], w[1:]
    pair = valid[:-1] & valid[1:] & (right-left <= 2.5*step)
    covered = 0.
    for a, b in windows:
        l, r = np.maximum(left, a), np.minimum(right, b)
        sel = pair & (r > l)
        idx = np.flatnonzero(sel)
        dx = right[idx]-left[idx]
        wr = ((r[idx]-left[idx])**2-(l[idx]-left[idx])**2)/(2*dx)
        wl = r[idx]-l[idx]-wr
        np.add.at(weights, idx, wl)
        np.add.at(weights, idx+1, wr)
        covered += (r[idx]-l[idx]).sum()
    return weights, covered/sum(b-a for a, b in windows)


def line_index(s, line, shift=0., flux_key='flux'):
    """Linear sideband fit and signed wing integral; errors propagate sidebands.

    Rest EW is integral(f/continuum - 1) d(lambda_rest); flux uses the observed
    wavelength differential. Host, FeII and narrow blends are not decomposed.
    """
    cfg = LINES[line]
    w = s['wave']/(1+s['z'])
    f = s[flux_key]
    empty = dict(coverage=0., valid=False)
    if f is None or not np.isfinite(s['z']):
        return empty
    good = np.isfinite(f)
    sides = [(a+shift, b+shift) for a, b in cfg['sides']]
    masks = [(w >= a) & (w <= b) & good for a, b in sides]
    if any(m.sum() < 4 for m in masks):
        return empty
    if any(integration_weights(w, good, [side])[1] < .85 for side in sides):
        return empty
    side = masks[0] | masks[1]
    x = (w-cfg['center'])/1000
    design = np.column_stack([np.ones(len(w)), x])
    inverse = np.linalg.pinv(design[side])
    coeff = inverse @ f[side]
    cont = design @ coeff
    weights, cov = integration_weights(w, good, cfg['wings'])
    if coeff[0] <= 0 or np.any(cont[weights != 0] <= 0):
        return dict(coverage=cov, valid=False)
    h = weights.copy()
    h[side] -= (weights @ design) @ inverse
    area_rest = h @ np.nan_to_num(f)
    # Exact EW relative to the fitted local slope, no constant-continuum shortcut.
    ew = weights @ np.where(good, np.nan_to_num(f)/np.where(cont != 0, cont, np.nan)-1, 0)
    grad = np.divide(weights, cont, out=np.zeros(len(w)), where=cont != 0)
    grad[side] -= ((weights*np.nan_to_num(f)/np.maximum(cont**2, 1e-100)) @ design) @ inverse
    err = s['error']
    supported = (h != 0) | (grad != 0)
    have_error = np.isfinite(err[supported]).all() and (err[supported] > 0).all()
    area_err = np.sqrt(np.sum((h[supported]*err[supported])**2))*(1+s['z']) if have_error else np.nan
    ew_err = np.sqrt(np.sum((grad[supported]*err[supported])**2)) if have_error else np.nan
    relevant = side | (weights != 0)
    q = s['quality'][relevant]
    telluric_fraction = np.mean((q & 4) != 0) if len(q) else 0.
    join = len(set(s['channel'][relevant])) > 1 and s['source'] == 'NGPS'
    result = dict(coverage=cov, valid=cov >= .98, wing_ew=ew, wing_flux=area_rest*(1+s['z']),
        wing_ew_err_diagonal=ew_err, wing_flux_err_diagonal=area_err, continuum=coeff[0],
        continuum_slope_per_1000A=coeff[1], telluric_fraction=telluric_fraction,
        cross_channel=join, outside_wavelength_validation=bool(np.any((q & 16) != 0)))
    if line == 'Hb':
        ow, oc = integration_weights(w, good, [(4990, 5025)])
        oh = ow.copy()
        oh[side] -= (ow @ design) @ inverse
        oarea = oh @ np.nan_to_num(f)
        result.update(oiii_coverage=oc, oiii_flux=oarea*(1+s['z']),
            wing_over_oiii=area_rest/oarea if oc >= .98 and oarea > 0 else np.nan)
    return result


def continuum_color(s):
    """Observed total pseudo-continuum slope; deliberately not an AGN-only slope."""
    w = s['wave']/(1+s['z'])
    x, y = [], []
    for a, b in [(4200, 4230), (5080, 5120), (5600, 5650)]:
        good = np.isfinite(s['flux'])
        sel = good & (w >= a) & (w <= b)
        if sel.sum() >= 4 and integration_weights(w, good, [(a, b)])[1] >= .9:
            level = np.median(s['flux'][sel])
            if level > 0:
                x.append((a+b)/2); y.append(level)
    out = dict(color_n_windows=len(x), color_windows=';'.join(map(str, x)), alpha_lambda_total=np.nan, c5100=np.nan, c4200=np.nan)
    if len(x) >= 2:
        out['alpha_lambda_total'] = np.polyfit(np.log(x), np.log(y), 1)[0]
    for xx, yy in zip(x, y):
        if xx == 5100:
            out['c5100'] = yy
        if xx == 4215:
            out['c4200'] = yy
    return out


def common_support_pair(a, b, observed_step=3.):
    """Sensitivity diagnostic on equal support, with no interpolation over gaps.

    No line-spread matching is implied. Resampling induces covariance, so errors
    are explicitly unavailable here; original native errors remain in the catalog.
    """
    grid = np.arange(max(a['wave'][0], b['wave'][0]), min(a['wave'][-1], b['wave'][-1]), observed_step)
    result = []
    for s in [a, b]:
        w, f = s['wave'], s['flux']
        idx = np.flatnonzero(np.isfinite(f))
        gap = (np.diff(idx) > 1) | (np.diff(w[idx]) > 2.5*np.median(np.diff(w)))
        y = np.full(len(grid), np.nan)
        for run in np.split(idx, np.flatnonzero(gap)+1):
            if len(run) < 2:
                continue
            mask = (grid >= w[run[0]]) & (grid <= w[run[-1]])
            y[mask] = np.interp(grid[mask], w[run], f[run])
        result.append(dict(s, wave=grid, flux=y, error=np.full(len(grid), np.nan),
            quality=np.zeros(len(grid), int), channel=np.full(len(grid), 'common-grid')))
    common = np.isfinite(result[0]['flux']) & np.isfinite(result[1]['flux'])
    for s in result:
        s['flux'][~common] = np.nan
    return result


def measure(spectra, objects):
    rows, epochs = [], []
    pops = objects.set_index('name').population.to_dict()
    for s in spectra:
        meta = {k: s[k] for k in ['name', 'epoch_id', 'source', 'night', 'mjd', 'min_mjd', 'max_mjd', 'z', 'path', 'proprietary', 'error_kind']}
        epochs.append(dict(**meta, good_pixels=int(np.isfinite(s['flux']).sum()), **continuum_color(s)))
        if pops[s['name']] == 'stellar/special':
            continue
        for line in LINES:
            for variant, shift in VARIANTS.items():
                metric = line_index(s, line, shift)
                rows.append(dict(**meta, line=line, variant=variant, **metric))
            if s['source'] == 'NGPS':
                for key in ['uncorrected', 'alternate']:
                    rows.append(dict(**meta, line=line, variant=key, **line_index(s, line, flux_key=key)))
    return pd.DataFrame(rows), pd.DataFrame(epochs)


def comparisons(metrics, epochs):
    """Latest previous accepted archive per line, plus first archive and NGPS repeat."""
    out = []
    for (name, line), group in metrics.groupby(['name', 'line']):
        nom = group[group.variant == 'nominal'].sort_values('mjd')
        for new in nom[nom.source == 'NGPS'].itertuples(index=False):
            if not new.valid:
                continue
            prev = nom[(nom.max_mjd < new.min_mjd) & nom.valid]
            archive = prev[prev.source != 'NGPS']
            pairs = []
            if len(archive):
                pairs.append(('latest_archive', archive.iloc[-1]))
                if len(archive) > 1:
                    pairs.append(('first_archive', archive.iloc[0]))
            repeats = prev[prev.source == 'NGPS']
            if len(repeats):
                pairs.append(('previous_NGPS', repeats.iloc[-1]))
            for baseline, old in pairs:
                ratios, deltas = [], []
                for variant in VARIANTS:
                    a = group[(group.epoch_id == old.epoch_id) & (group.variant == variant)]
                    b = group[(group.epoch_id == new.epoch_id) & (group.variant == variant)]
                    if len(a) and len(b) and a.iloc[0].valid and b.iloc[0].valid:
                        av, bv = a.iloc[0].wing_ew, b.iloc[0].wing_ew
                        deltas.append(bv-av)
                        if av > 2 and bv > 2:
                            ratios.append(bv/av)
                ratio = new.wing_ew/old.wing_ew if old.wing_ew > 2 and new.wing_ew > 2 else np.nan
                amin = min(ratios) if len(ratios) == 3 else np.nan
                amax = max(ratios) if len(ratios) == 3 else np.nan
                robust_direction = len(deltas) == 3 and (min(deltas) > 0 or max(deltas) < 0)
                row = dict(name=name, line=line, baseline=baseline, old_epoch=old.epoch_id, new_epoch=new.epoch_id,
                    night=new.night, old_mjd=old.mjd, new_mjd=new.mjd,
                    baseline_rest_years=(new.mjd-old.mjd)/365.25/(1+new.z),
                    old_ew=old.wing_ew, new_ew=new.wing_ew, delta_ew=new.wing_ew-old.wing_ew,
                    ew_ratio=ratio, ratio_min=amin, ratio_max=amax, continuum_variants=len(deltas),
                    direction_stable=robust_direction, telluric_fraction=new.telluric_fraction,
                    cross_channel=new.cross_channel, old_error_kind=old.error_kind,
                    new_error_kind=new.error_kind,
                    review_large_change=bool(robust_direction and np.isfinite(amin) and (amin > 1.5 or amax < 2/3)),
                    delta_ew_err_diagonal=np.hypot(new.wing_ew_err_diagonal, old.wing_ew_err_diagonal))
                if line == 'Hb':
                    a, b = old.get('wing_over_oiii', np.nan), getattr(new, 'wing_over_oiii', np.nan)
                    row['oiii_normalized_ratio'] = b/a if a > 0 and b > 0 else np.nan
                variants = group[(group.epoch_id == new.epoch_id) & group.variant.isin(['uncorrected', 'alternate'])]
                for v in variants.itertuples(index=False):
                    row[v.variant+'_ew_ratio_to_primary'] = v.wing_ew/new.wing_ew if v.valid and abs(new.wing_ew) > 2 else np.nan
                out.append(row)
    return pd.DataFrame(out)


def photometry(objects, targets):
    """Separate survey, corrected-alert, and IR streams; one ZTF object per band."""
    rows = []
    history = read_csv(ROOT/'jwst_proposal/review/history_audit/binned_measurements.csv')
    for obj in objects.itertuples(index=False):
        survey_found = False
        candidates = []
        if obj.jwst_id:
            candidates.append(ROOT/f'jwst_proposal/review/history_audit/raw_ztf/{obj.jwst_id}.csv')
        for alias in obj.aliases.split(';'):
            for tag in ['review_dr24_20260920', 'pool', 'zeltyn', 'calib', 'v2']:
                candidates.append(ROOT/f'data/ztf_cache/{tag}/{alias}.csv')
        for p in candidates:
            if not p.exists() or p.stat().st_size < 20:
                continue
            try:
                d = read_csv(p, dtype={'oid': str}, low_memory=False)
            except pd.errors.EmptyDataError:
                continue
            if not {'mjd', 'mag', 'magerr', 'filtercode', 'oid'}.issubset(d.columns):
                continue
            for col in ['mjd', 'mag', 'magerr', 'catflags', 'ra', 'dec']:
                if col in d:
                    d[col] = pd.to_numeric(d[col], errors='coerce')
            good = np.isfinite(d.mjd) & np.isfinite(d.mag) & np.isfinite(d.magerr) & (d.magerr > 0) & (d.magerr < .3)
            if 'catflags' in d:
                good &= d.catflags.fillna(32768).astype(int).eq(0)
            if {'ra', 'dec'}.issubset(d.columns):
                pos = SkyCoord(d.ra.to_numpy()*u.deg, d.dec.to_numpy()*u.deg)
                good &= pos.separation(SkyCoord(obj.ra*u.deg, obj.dec*u.deg)).arcsec < 1.5
            d = d[good].copy()
            for band, b in d.groupby('filtercode'):
                if band not in ['zg', 'zr', 'zi']:
                    continue
                oid = b.groupby('oid').size().idxmax()
                b = b[b.oid == oid].sort_values('magerr').drop_duplicates('mjd')
                for r in b.itertuples(index=False):
                    rows.append(dict(name=obj.name, band=band[1:], mjd=r.mjd, value=r.mag, error=r.magerr,
                        units='AB mag', stream='ZTF survey', source=str(p.relative_to(ROOT)), oid=oid, bin_days=0))
                survey_found = True
            if survey_found:
                break
        t = targets.get(obj.name, {})
        if not survey_found:
            h = history[(history.internal_id == obj.name) & history.band.isin(['g', 'r'])]
            for r in h.itertuples(index=False):
                rows.append(dict(name=obj.name, band=r.band, mjd=r.mjd, value=r.mag, error=r.mag_se,
                    units='AB mag', stream='ZTF 90-day cache', source='jwst_proposal/review/history_audit/binned_measurements.csv', oid=r.oid, bin_days=90))
            if h.empty:
                for band, points in (t.get('ztf') or {}).items():
                    for mjd, mag, err in points:
                        rows.append(dict(name=obj.name, band=band, mjd=mjd, value=mag, error=err,
                            units='AB mag', stream='ZTF display cache', source='candidate/target payload', oid='', bin_days=7))
        alert = ROOT/f'observing/science_review/alerts/{obj.name}_corrected.csv'
        if alert.exists():
            a = read_csv(alert)
            required = {'mjd', 'fid', 'magpsf_corr', 'sigmapsf_corr_ext', 'dubious', 'corrected'}
            if required.issubset(a.columns):
                ok = a.corrected.astype(str).str.lower().eq('true') & a.dubious.astype(str).str.lower().eq('false')
                ok &= np.isfinite(a.magpsf_corr) & a.sigmapsf_corr_ext.between(0, .3, inclusive='neither')
                if 'match_arcsec' in a:
                    ok &= a.match_arcsec < 1.5
                a = a[ok].sort_values('sigmapsf_corr_ext').drop_duplicates(['mjd', 'fid'])
                for r in a.itertuples(index=False):
                    rows.append(dict(name=obj.name, band={1: 'g', 2: 'r', 3: 'i'}.get(r.fid, '?'),
                        mjd=r.mjd, value=r.magpsf_corr, error=r.sigmapsf_corr_ext, units='AB mag',
                        stream='ZTF corrected alerts', source=str(alert.relative_to(ROOT)), oid=r.oid, bin_days=0))
        h = history[(history.internal_id == obj.name) & history.band.isin(['W1', 'W2'])]
        for r in h.itertuples(index=False):
            rows.append(dict(name=obj.name, band=r.band, mjd=r.mjd, value=r.mag, error=r.mag_se,
                units='Vega mag', stream='NEOWISE visits', source='jwst_proposal/review/history_audit/binned_measurements.csv', oid='', bin_days=10))
        if h.empty:
            for band, points in (t.get('wise') or {}).items():
                for mjd, flux, err in points:
                    if flux > 0 and band in ['W1', 'W2']:
                        zp = 309540. if band == 'W1' else 171787.
                        rows.append(dict(name=obj.name, band=band, mjd=mjd, value=-2.5*np.log10(flux/zp),
                            error=2.5/np.log(10)*err/flux, units='Vega mag', stream='unWISE visits',
                            source='candidate/target payload; mJy converted with WISE Vega zero point', oid='', bin_days=10))
            # NEOWISE and unWISE stay distinct to avoid a false offset or double counting.
            for mjd, flux, err in t.get('neo', []) or []:
                if flux > 0:
                    rows.append(dict(name=obj.name, band='W1', mjd=mjd, value=-2.5*np.log10(flux/309540.),
                        error=2.5/np.log(10)*err/flux, units='Vega mag', stream='NEOWISE visits',
                        source='candidate/target payload', oid='', bin_days=10))
    return pd.DataFrame(rows).drop_duplicates(['name', 'band', 'stream', 'mjd']).sort_values(['name', 'stream', 'band', 'mjd'])


def photometric_summary(photo, visits):
    rows = []
    for (name, band, stream), g in photo.groupby(['name', 'band', 'stream']):
        g = g.sort_values('mjd')
        # Equal time bins prevent the most densely sampled season dominating slope/amplitude.
        width = 90 if band.startswith('W') else 30
        b = g.groupby(np.floor(g.mjd/width)).agg(mjd=('mjd', 'median'), value=('value', 'median'))
        last = b[b.mjd >= b.mjd.max()-3*365.25]
        slope = np.polyfit((last.mjd-last.mjd.min())/365.25, last.value, 1)[0] if len(last) >= 4 and last.mjd.max()-last.mjd.min() > 180 else np.nan
        rows.append(dict(name=name, band=band, stream=stream, n_points=len(g), n_bins=len(b),
            first_mjd=g.mjd.min(), last_mjd=g.mjd.max(), p95_p05_mag=np.ptp(np.percentile(b.value, [5, 95])),
            recent_slope_mag_per_observed_year=slope, last_minus_first_mag=b.value.iloc[-1]-b.value.iloc[0],
            gap_to_latest_NGPS_days=visits[visits.name == name].mjd.max()-g.mjd.max()))
    return pd.DataFrame(rows)


def match_photometry(photo, epochs, window=30.):
    rows = []
    for e in epochs.itertuples(index=False):
        for band in ['g', 'r']:
            p = photo[(photo.name == e.name) & (photo.band == band) & photo.stream.isin(['ZTF survey', 'ZTF corrected alerts'])]
            for stream, g in p.groupby('stream'):
                nearest = float(np.min(abs(g.mjd-e.mjd)))
                m = g[abs(g.mjd-e.mjd) <= window]
                rows.append(dict(name=e.name, epoch_id=e.epoch_id, source=e.source, band=band, stream=stream,
                    spectral_mjd=e.mjd, nearest_days=nearest, window_days=window, n_points=len(m),
                    matched_mag=m.value.median() if len(m) else np.nan,
                    matched_scatter=m.value.std() if len(m) > 1 else np.nan,
                    coadd_span_days=e.max_mjd-e.min_mjd))
    return pd.DataFrame(rows)


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    objects, visits, targets, jwst = census()
    evidence = read_csv(ROOT/'observing/science_review/all_candidate_evidence.csv')
    groupcols = ['name', 'experiment_group', 'known_state_status', 'post_spectrum_optical_trigger', 'post_spectrum_ir_flag']
    objects = objects.merge(evidence[groupcols].drop_duplicates('name'), on='name', how='left')
    objects['experiment_group'] = objects.experiment_group.fillna('outside pre-run review')
    objects['known_state_status'] = objects.known_state_status.fillna('not audited in pre-run review')
    print(f'Census: {len(objects)} objects, {len(visits)} visits', flush=True)
    spectra, audit = load_spectra(objects, visits)
    print(f'Spectral screen: {len(spectra)} retained epochs', flush=True)
    metrics, epochs = measure(spectra, objects)
    pairs = comparisons(metrics, epochs)
    photo = photometry(objects, targets)
    summary = photometric_summary(photo, visits)
    matches = match_photometry(photo, epochs)
    coverage = epochs.groupby('name').agg(n_epochs=('epoch_id', 'size'), first_spectrum_mjd=('mjd', 'min'), last_spectrum_mjd=('mjd', 'max'))
    objects = objects.merge(coverage, on='name', how='left')
    objects['n_NGPS'] = objects.name.map(visits.groupby('name').size())
    objects['n_archive'] = objects.n_epochs-objects.n_NGPS
    objects['n_native_error_epochs'] = objects.name.map(epochs[epochs.error_kind != 'missing'].groupby('name').size()).fillna(0).astype(int)
    objects['spectral_baseline_rest_years'] = (objects.last_spectrum_mjd-objects.first_spectrum_mjd)/365.25/(1+objects.z)
    objects['has_optical'] = objects.name.isin(photo.loc[photo.band.isin(['g', 'r']), 'name'])
    objects['has_IR'] = objects.name.isin(photo.loc[photo.band.isin(['W1', 'W2']), 'name'])
    readiness = jwst.merge(objects, left_on='internal_id', right_on='name', how='left', suffixes=('_proposal', ''))
    tables = dict(objects=objects, visits=visits, epochs=epochs, archive_audit=audit, line_indices=metrics,
        comparisons=pairs, photometry=photo, photometry_summary=summary, photometry_matches=matches, jwst_coverage=readiness)
    for name, table in tables.items():
        table.to_csv(OUT/f'{name}.csv', index=False)
    track(ROOT/'jwst_proposal/proposal.tex')
    track(ROOT/'observing/science_review/README.md')
    manifest = dict(snapshot_date='2026-11-03', inputs=SOURCES,
        counts={k: len(v) for k, v in tables.items()},
        definitions=dict(line_indices='signed wing pseudo-EW; no host/FeII decomposition',
            error_budget='diagonal statistical only where actual pixel errors exist; no imputed archival errors',
            baseline='latest accepted covered archive before each NGPS visit; first and repeat separately',
            coadds='DESI date-verified intervals <=30 days; midpoint used and endpoints retained'),
        software=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'provenance.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(manifest['counts'], indent=2), flush=True)
    return tables, spectra


if __name__ == '__main__':
    run()
