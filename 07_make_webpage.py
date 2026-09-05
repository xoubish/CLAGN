"""
07_make_webpage.py  --  build the self-contained CLAGN Night Sheet (web/clagn_night_sheet.html).

Reads the pipeline outputs in data/ (targets_<night>.csv, master_list_scored.csv, spectra epochs, ZTF and unWISE
caches, Sample A / Zeltyn embeddings) and writes one HTML file with the data embedded as JSON and the charts drawn
client-side (no external libraries). Re-run after every pipeline update; drop observed NGPS spectra as
data/ngps_spectra/<name>.csv (columns wave_A, flux) and they appear in each target's spectrum slot.

Usage: /opt/anaconda3/bin/python 07_make_webpage.py
"""
import os, json, glob, datetime
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
OUT = os.path.join(HERE, 'web', 'clagn_night_sheet.html')
NIGHT_META = {
    'sep23': dict(label='Sep 23', date='2026-09-23', part='first half', window='20:05 – 00:39 PDT', lst='19.5h – 0.1h',
                  moon='93 %', moon_pos='RA 22.2h  Dec −12°', moon_note='up all window', slots=14, mjd=61307.2),
    'oct26': dict(label='Oct 26', date='2026-10-26', part='full night', window='19:24 – 05:38 PST', lst='21.0h – 7.2h',
                  moon='98 %', moon_pos='RA 3.0h  Dec +22°', moon_note='rises ~3.5 h after twilight', slots=34, mjd=61340.3),
    'oct27': dict(label='Oct 27', date='2026-10-27', part='full night', window='19:23 – 05:39 PST', lst='21.0h – 7.3h',
                  moon='94 %', moon_pos='RA 4.1h  Dec +26°', moon_note='rises ~3.5 h after twilight', slots=34, mjd=61341.3),
}
TIER_LABEL = {'T1': 'Tier 1 · manifold-selected', 'T2': 'Tier 2 · EVQ completing transition',
              'T3': 'Tier 3 · confirmed CLAGN, revisit', 'T4': 'Tier 4 · control'}
LINES = [('Mg II', 2798.0), ('Hβ', 4861.0), ('[O III]', 5007.0), ('Hα', 6563.0), ('Ca II', 8600.0), ('[S III]', 9531.0)]


def mjd_to_date(mjd):
    return Time(mjd, format='mjd').to_datetime().strftime('%Y-%m-%d') if np.isfinite(mjd) else ''


def load_targets():
    rows = []
    for night in NIGHT_META:
        p = os.path.join(DATA, f'targets_{night}.csv')
        if os.path.exists(p):
            t = pd.read_csv(p); t['night'] = night; rows.append(t)
    t = pd.concat(rows, ignore_index=True)
    m = pd.read_csv(os.path.join(DATA, 'master_list_scored.csv'))
    return t, m


def ztf_series(name, tags=('zeltyn', 'pool')):
    for tag in tags:
        p = os.path.join(DATA, 'ztf_cache', tag, f'{name}.csv')
        if os.path.exists(p) and os.path.getsize(p) > 5:
            df = pd.read_csv(p)
            if 'mjd' not in df or df.mjd.isna().all():
                return {}
            out = {}
            for fc, key in [('zg', 'g'), ('zr', 'r')]:
                b = df[(df.filtercode == fc) & np.isfinite(df.mag)].copy()
                if not len(b):
                    continue
                b['night'] = np.floor(b.mjd)
                g = b.groupby('night').agg(mjd=('mjd', 'median'), mag=('mag', 'median'), err=('magerr', 'median'), n=('mag', 'size'))
                out[key] = [[round(r.mjd, 2), round(r.mag, 3), round(float(r.err), 3)] for r in g.itertuples()]
            return out
    return {}


def load_wise(cache_dir, name_by_objectid=None):
    files = sorted(glob.glob(os.path.join(cache_dir, '*.parquet')))
    if not files:
        return {}
    lc = pd.concat([pd.read_parquet(f) for f in files]).reset_index()
    lc['band'] = lc['band'].str.replace('WISE_', '', regex=False)
    out = {}
    for (oid, band), g in lc.groupby(['objectid', 'band']):
        key = name_by_objectid.get(oid) if name_by_objectid else f'P{oid}'
        if key is None:
            continue
        g = g.sort_values('time')
        out.setdefault(key, {})[band] = [[round(t, 1), round(f, 4), round(e, 4)] for t, f, e in zip(g.time, g.flux, g.err)]
    return out


def epochs_for(name, ep_tables):
    rows = []
    for ep in ep_tables:
        e = ep[ep.name == name]
        for r in e.sort_values('mjd').itertuples():
            if getattr(r, 'is_coadd', False) is True:
                continue
            src = getattr(r, 'source', 'SDSS')
            prog = str(getattr(r, 'programname', '') if src == 'SDSS' else f"{getattr(r, 'survey', '')}/{getattr(r, 'program', '')}")
            rows.append(dict(mjd=round(float(r.mjd), 1), date=mjd_to_date(float(r.mjd)), src=src, prog=prog.replace('nan', ''),
                             cls=str(getattr(r, 'class', '') or '').replace('nan', ''), sub=str(getattr(r, 'subclass', '') or '').replace('nan', ''),
                             z=(round(float(r.z), 4) if pd.notna(getattr(r, 'z', np.nan)) else None)))
    return rows


def main():
    t, m = load_targets()
    names = t.name.unique().tolist()
    mm = m.set_index('name')
    ep_tables = []
    for tag in ['zeltyn', 'pool']:
        p = os.path.join(DATA, f'spectra_epochs_{tag}.csv')
        if os.path.exists(p):
            ep_tables.append(pd.read_csv(p, low_memory=False))
    # unWISE: pool cache keyed by poolid; Zeltyn cache keyed by 1..N in zeltyn_coords order
    zc = pd.read_csv(os.path.join(DATA, 'zeltyn_coords.csv'))
    wise = load_wise(os.path.join(DATA, 'wise_cache', 'pool'))
    wise.update(load_wise(os.path.join(DATA, 'wise_cache', 'zeltyn'), {i + 1: n for i, n in enumerate(zc.name)}))
    # NEOWISE-R visits to 2024 (03c_neowise_now.py), converted to mJy so they share the W1 panel
    neo = {}
    for tag in ['zeltyn', 'pool']:
        p = os.path.join(DATA, f'neowise_visits_{tag}.csv')
        if os.path.exists(p):
            v = pd.read_csv(p)
            for nm, g in v.groupby('name'):
                g = g.sort_values('mjd')
                fl = 309.54e3 * 10 ** (-0.4 * g.w1.values); er = fl * 0.921 * g.w1err.fillna(0.05).values
                neo[nm] = [[round(a, 1), round(b, 4), round(c, 4)] for a, b, c in zip(g.mjd, fl, er)]
    # manifold background
    A = pd.read_csv(os.path.join(DATA, 'sampleA_embedding_objectid.csv'))
    Z = pd.read_csv(os.path.join(DATA, 'zeltyn_embedding_full.csv'))
    Zc = Z[Z.projected & Z.is_clagn & (Z.z < 1)]
    NB = 10
    hist, xe, ye = np.histogram2d(A.umap_x, A.umap_y, bins=NB)
    ix = np.clip(np.searchsorted(xe[1:], A.umap_x), 0, NB - 1); iy = np.clip(np.searchsorted(ye[1:], A.umap_y), 0, NB - 1)
    rects = []
    for col, key in [('in_region_zeltyn', 'zeltyn'), ('in_region_clagn', 'clagn')]:
        flagged = set(zip(ix[A[col].values], iy[A[col].values]))
        rects += [dict(kind=key, x0=float(xe[i]), x1=float(xe[i + 1]), y0=float(ye[j]), y1=float(ye[j + 1])) for i, j in flagged]
    Ze = Z[Z.projected & ~Z.is_clagn & (Z.z < 1)]
    code = ((A.label_bits & 16) > 0).astype(int) + 2 * ((A.label_bits & 32) > 0).astype(int)   # 0 none, 1 turn-on, 2 turn-off, 3 both
    manifold = dict(
        A=[[round(x, 2), round(y, 2), int(c)] for x, y, c in zip(A.umap_x, A.umap_y, code)],
        Z=[[round(x, 2), round(y, 2)] for x, y in zip(Zc.umap_x, Zc.umap_y)],
        E=[[round(x, 2), round(y, 2)] for x, y in zip(Ze.umap_x, Ze.umap_y)],
        rects=rects, xlim=[float(xe[0]), float(xe[-1])], ylim=[float(ye[0]), float(ye[-1])])

    targets = []
    for name in names:
        r = mm.loc[name]
        if isinstance(r, pd.DataFrame):
            r = r.iloc[0]
        sc = SkyCoord(r.ra * u.deg, r.dec * u.deg)
        nights = []
        for tr in t[t.name == name].itertuples():
            nights.append(dict(night=tr.night, rank=int(tr.rank), prio=float(tr.priority_night) if pd.notna(tr.priority_night) else None,
                               hrs=round(float(getattr(tr, f'hrs_{tr.night}', np.nan)), 1), moonsep=float(getattr(tr, f'moonsep_{tr.night}', np.nan)),
                               minx=float(getattr(tr, f'minX_{tr.night}', np.nan)),
                               texp=float(getattr(tr, 't_exp_min', np.nan)) if pd.notna(getattr(tr, 't_exp_min', np.nan)) else None,
                               plan=str(getattr(tr, 'exp_plan', '')) if pd.notna(getattr(tr, 'exp_plan', np.nan)) else '',
                               pph=float(getattr(tr, 'prio_per_hour', np.nan)) if pd.notna(getattr(tr, 'prio_per_hour', np.nan)) else None))
        z = float(r.z) if pd.notna(r.z) else None
        lines = [dict(name=n, obs=round(w * (1 + z), 0), inrange=bool(3200 <= w * (1 + z) <= 10400)) for n, w in LINES] if z is not None else []
        cut = {}
        for kind, ext, mime in [('sdss', 'jpg', 'image/jpeg'), ('ps1_g', 'png', 'image/png'), ('ps1_r', 'png', 'image/png'),
                                ('ztf_g', 'png', 'image/png'), ('ztf_r', 'png', 'image/png')]:
            cp = os.path.join(DATA, 'cutouts', f'{name}_{kind}.{ext}')
            if os.path.exists(cp) and os.path.getsize(cp) > 100:
                import base64
                cut[kind] = f'data:{mime};base64,' + base64.b64encode(open(cp, 'rb').read()).decode('ascii')
        # archival spectra (03d_fetch_spectra.py): up to four epochs for the overlay + the full EW history
        spec = None
        sp_p = os.path.join(DATA, 'spectra_dl', f'{name}.json')
        if os.path.exists(sp_p):
            recs = json.load(open(sp_p))
            recs = [r for r in recs if r.get('flux') and any(v is not None for v in r['flux'])]
            if recs:
                sdss = [r for r in recs if r.get('source', 'SDSS') == 'SDSS' and not r.get('coadd')]
                coadd = [r for r in recs if r.get('source', 'SDSS') == 'SDSS' and r.get('coadd')]
                desi = [r for r in recs if r.get('source') == 'DESI']
                pick = []
                if sdss:
                    pick.append(sdss[0])
                    if len(sdss) > 1 and sdss[-1] is not sdss[0]:
                        pick.append(sdss[-1])
                if coadd:
                    pick.append(coadd[-1])
                if desi:
                    pick.append(desi[0])
                def rebin12(fl):
                    a = np.array([np.nan if v is None else v for v in fl], float)
                    a = a[:len(a) // 2 * 2].reshape(-1, 2)
                    with np.errstate(all='ignore'):
                        m = np.nanmean(a, axis=1)
                    return [None if np.isnan(v) else round(float(v), 2) for v in m]
                w12 = [round(float(w), 1) for w in np.array(recs[0]['wave'])[:len(recs[0]['wave']) // 2 * 2].reshape(-1, 2).mean(axis=1)]
                def lab(r):
                    d = mjd_to_date(r['mjd']) if r.get('mjd') and np.isfinite(r['mjd']) else ''
                    src = 'DESI' if r.get('source') == 'DESI' else ('SDSS-V coadd' if r.get('coadd') else f"SDSS {str(r.get('program', '')).strip() or ''}".strip())
                    return f'{d[:7]} {src}'.strip()
                spec = dict(wave=w12, epochs=[dict(label=lab(r), flux=rebin12(r['flux']), cls=str(r['meta'].get('class', '') or ''),
                                                    z=r['meta'].get('z'), ew_hb=r.get('ew', {}).get('Hb'), ew_ha=r.get('ew', {}).get('Ha')) for r in pick],
                            history=[dict(date=mjd_to_date(r['mjd']) if r.get('mjd') and np.isfinite(r['mjd']) else '', src=('DESI' if r.get('source') == 'DESI' else 'SDSS'),
                                          prog=str(r.get('program', '')).strip(), coadd=bool(r.get('coadd')), cls=str(r['meta'].get('class', '') or ''),
                                          ew_hb=(None if r.get('ew', {}).get('Hb') is None or not np.isfinite(r['ew']['Hb']) else round(r['ew']['Hb'], 1)),
                                          ew_ha=(None if r.get('ew', {}).get('Ha') is None or not np.isfinite(r['ew']['Ha']) else round(r['ew']['Ha'], 1)),
                                          sn=r['meta'].get('sn_median_all')) for r in recs])
        spec_p = os.path.join(DATA, 'ngps_spectra', f'{name}.csv')
        ngps = None
        if os.path.exists(spec_p):
            s = pd.read_csv(spec_p); ngps = [[round(a, 1), round(b, 4)] for a, b in zip(s.iloc[:, 0], s.iloc[:, 1])]
        def f(col, nd=3):
            v = r.get(col, np.nan)
            return (round(float(v), nd) if pd.notna(v) and np.isfinite(float(v)) else None) if not isinstance(v, str) else v
        jname = None
        if str(name).startswith('P'):        # pool objects: show an SDSS-style J-name, keep the pool id as the key
            jname = 'J' + sc.ra.to_string(u.hour, sep='', precision=2, pad=True) + sc.dec.to_string(u.deg, sep='', precision=1, alwayssign=True, pad=True)
        targets.append(dict(
            name=name, jname=jname, tier=r.tier, source=r.source_catalog, ra=round(float(r.ra), 6), dec=round(float(r.dec), 6),
            sex=sc.to_string('hmsdms', sep=':', precision=1), z=z, rmag=f('r_mag', 2), priority=f('priority'),
            M=f('M'), P=f('P'), S=f('S'), B=f('B'), trend=str(r.get('trend', '')), notes=str(r.get('notes', '')),
            clagn_score=f('clagn_score'), density=f('zeltyn_density_ratio', 2), m_comb=f('M_combined', 2), in_zeltyn=bool(r.get('in_region_zeltyn', False)),
            in_clagn=bool(r.get('in_region_clagn', False)), ux=f('umap_x'), uy=f('umap_y'),
            n_spec=f('n_spec', 0), yrs=f('years_since_last_spec', 1), last_class=str(r.get('last_class', '')),
            mjd_last=f('mjd_last_spec', 1), w1_ratio=f('w1_ratio_now_over_spec', 2), dr_ref=f('dr_since_ref', 2),
            r_last=f('r_last', 2), g_last=f('g_last', 2), mjd_last_ztf=f('mjd_last_ztf', 1),
            lines=lines, nights=nights, epochs=epochs_for(name, ep_tables), ztf=ztf_series(name), wise=wise.get(name, {}), ngps=ngps,
            cut=cut, neo=neo.get(name, []), spec=spec))

    # night summary counts
    for n, meta in NIGHT_META.items():
        sub = t[(t.night == n) & (t['rank'] > 0)]
        meta['counts'] = sub.tier.value_counts().to_dict(); meta['n_primary'] = int(len(sub)); meta['n_backup'] = int(((t.night == n) & (t['rank'] == 0)).sum())

    # program numbers for the About tab
    def nrows(fn):
        p = os.path.join(DATA, fn); return int(sum(1 for _ in open(p)) - 1) if os.path.exists(p) else None
    pool_scored = pd.read_csv(os.path.join(DATA, 'parent_pool_scored.csv'), usecols=['projected', 'psfmag_r']) if os.path.exists(os.path.join(DATA, 'parent_pool_scored.csv')) else None
    stats = dict(n_sampleA=len(A), n_sampleA_clagn=int(A.is_known_clagn.sum()), n_zeltyn=len(zc), n_zeltyn_clagn=int(zc['class_zeltyn'].str.startswith('CL-AGN').sum()),
                 n_pool=nrows('parent_pool_dr16qso.csv'), n_pool_projected=int(pool_scored.projected.fillna(False).sum()) if pool_scored is not None else None,
                 n_T1=int((m.tier == 'T1').sum()), n_T4=int((m.tier == 'T4').sum()), n_lit=nrows('literature_clagn.csv'),
                 n_targets=len(targets), n_primary={k: v['n_primary'] for k, v in NIGHT_META.items()})
    payload = dict(generated=datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), nights=NIGHT_META, tiers=TIER_LABEL,
                   manifold=manifold, targets=targets, stats=stats)
    html = TEMPLATE.replace('__DATA__', json.dumps(payload, separators=(',', ':'), allow_nan=False))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w') as fh:
        fh.write(html)
    print(f'wrote {OUT}: {len(targets)} targets, {os.path.getsize(OUT)/1e6:.1f} MB')
    # standalone copy for GitHub Pages (docs/index.html): the artifact host wraps the fragment in a document, GitHub does not
    head_end = html.find('<header>')
    standalone = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width,initial-scale=1">\n'
                  + html[:head_end] + '</head>\n<body>\n' + html[head_end:] + '\n</body>\n</html>\n')
    os.makedirs(os.path.join(HERE, 'docs'), exist_ok=True)
    with open(os.path.join(HERE, 'docs', 'index.html'), 'w') as fh:
        fh.write(standalone)
    print(f'wrote docs/index.html (standalone, {os.path.getsize(os.path.join(HERE, "docs", "index.html"))/1e6:.1f} MB)')


TEMPLATE = r'''<title>CLAGN Night Sheet</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Barlow+Condensed:wght@500;600;700&family=Source+Sans+3:wght@400;600&family=JetBrains+Mono:wght@400;500&display=swap">
<style>
:root{
  color-scheme:dark;
  --ground:#0e1420; --surface:#151d2c; --raised:#1b2536; --hair:#26314a; --hair2:#31405d;
  --ink:#e6e9ef; --ink2:#aab2c2; --ink3:#727c91; --accent:#d9a441; --accent-soft:rgba(217,164,65,.14);
  --g:#199e70; --r:#d95926; --w1:#9085e9;
  --t1:#3987e5; --t2:#d95926; --t3:#199e70; --t4:#727c91;
  --sdss:#aab2c2; --desi:#d9a441;
}
@media (prefers-color-scheme: light){
  :root:not([data-theme="dark"]){
    color-scheme:light;
    --ground:#f4f6fa; --surface:#ffffff; --raised:#eaeef5; --hair:#d6dce8; --hair2:#c2cad9;
    --ink:#171c26; --ink2:#4b5568; --ink3:#7d869a; --accent:#a8761a; --accent-soft:rgba(168,118,26,.14);
    --g:#158f63; --r:#d6552a; --w1:#4a3aa7;
    --t1:#2a78d6; --t2:#d6552a; --t3:#158f63; --t4:#7d869a;
    --sdss:#4b5568; --desi:#a8761a;
  }
}
:root[data-theme="light"]{
  color-scheme:light;
  --ground:#f4f6fa; --surface:#ffffff; --raised:#eaeef5; --hair:#d6dce8; --hair2:#c2cad9;
  --ink:#171c26; --ink2:#4b5568; --ink3:#7d869a; --accent:#a8761a; --accent-soft:rgba(168,118,26,.14);
  --g:#158f63; --r:#d6552a; --w1:#4a3aa7;
  --t1:#2a78d6; --t2:#d6552a; --t3:#158f63; --t4:#7d869a;
  --sdss:#4b5568; --desi:#a8761a;
}
*{box-sizing:border-box}
body{margin:0;background:var(--ground);color:var(--ink);font:15px/1.45 "Source Sans 3",system-ui,sans-serif}
.mono{font-family:"JetBrains Mono",ui-monospace,Menlo,monospace;font-variant-numeric:tabular-nums}
.disp{font-family:"Barlow Condensed","Arial Narrow",sans-serif}
a{color:var(--accent)}
header{position:sticky;top:0;z-index:5;background:var(--ground);border-bottom:1px solid var(--hair)}
.bar{display:flex;flex-wrap:wrap;align-items:center;gap:14px 22px;padding:12px 22px;max-width:1500px;margin:0 auto}
.bar h1{margin:0;font:600 26px/1 "Barlow Condensed",sans-serif;letter-spacing:.01em}
.bar h1 small{display:block;font:400 12px/1.3 "Source Sans 3",sans-serif;color:var(--ink3);letter-spacing:.08em;text-transform:uppercase;margin-top:4px}
.tabs{display:flex;gap:4px}
.tab{background:none;border:1px solid var(--hair);color:var(--ink2);padding:6px 12px;border-radius:3px;font:600 15px "Barlow Condensed",sans-serif;letter-spacing:.03em;cursor:pointer}
.tab[aria-pressed="true"]{border-color:var(--accent);color:var(--ink);background:var(--accent-soft)}
.tab:focus-visible,.chip:focus-visible,input:focus-visible,.toggle:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.chips{display:flex;gap:6px;flex-wrap:wrap}
.chip{background:none;border:1px solid var(--hair);color:var(--ink2);padding:4px 9px;border-radius:999px;font:600 12px "Barlow Condensed",sans-serif;letter-spacing:.06em;cursor:pointer;display:inline-flex;align-items:center;gap:6px}
.chip i{width:8px;height:8px;border-radius:50%;display:inline-block;background:var(--c)}
.chip[aria-pressed="false"]{opacity:.45}
input[type=search]{background:var(--surface);border:1px solid var(--hair);color:var(--ink);padding:6px 10px;border-radius:3px;font:14px "Source Sans 3",sans-serif;min-width:200px}
.night{max-width:1500px;margin:0 auto;padding:14px 22px 6px;display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:10px 26px;color:var(--ink2);font-size:13.5px}
.night b{display:block;color:var(--ink3);font:600 11px "Barlow Condensed",sans-serif;letter-spacing:.1em;text-transform:uppercase;margin-bottom:2px}
.night .v{color:var(--ink)}
main{max-width:1500px;margin:0 auto;padding:10px 22px 60px}
.section{font:600 13px "Barlow Condensed",sans-serif;letter-spacing:.12em;text-transform:uppercase;color:var(--ink3);margin:22px 0 8px;display:flex;align-items:center;gap:12px}
.section::after{content:"";flex:1;height:1px;background:var(--hair)}
.card{background:var(--surface);border:1px solid var(--hair);border-radius:4px;margin:0 0 12px;display:grid;grid-template-columns:230px minmax(360px,1fr) 300px;gap:0}
.card > div{padding:14px 16px;min-width:0}
.card .id{border-right:1px solid var(--hair)}
.card .side{border-left:1px solid var(--hair);display:flex;flex-direction:column;gap:12px}
.card .foot{grid-column:1 / -1;border-top:1px solid var(--hair);display:grid;grid-template-columns:minmax(300px,1fr) minmax(300px,1.4fr);gap:18px}
.rank{font:700 30px/1 "Barlow Condensed",sans-serif;color:var(--ink);display:flex;align-items:baseline;gap:8px}
.rank small{font:600 12px "Barlow Condensed",sans-serif;letter-spacing:.1em;color:var(--ink3);text-transform:uppercase}
.tier{display:inline-flex;align-items:center;gap:6px;margin-top:8px;font:600 12px "Barlow Condensed",sans-serif;letter-spacing:.08em;text-transform:uppercase;color:var(--ink2)}
.tier i{width:10px;height:10px;border-radius:2px;background:var(--c)}
h2{margin:10px 0 2px;font:600 21px/1.1 "Barlow Condensed",sans-serif;letter-spacing:.01em;word-break:break-all}
.coords{color:var(--ink2);font-size:12.5px;margin-bottom:10px}
.kv{display:grid;grid-template-columns:auto 1fr;gap:3px 10px;font-size:13px;color:var(--ink2)}
.kv b{font-weight:600;color:var(--ink3);font-family:"Barlow Condensed",sans-serif;letter-spacing:.06em;text-transform:uppercase;font-size:11.5px;padding-top:2px}
.kv .v{color:var(--ink)}
.meters{margin-top:12px;display:grid;gap:5px}
.meter{display:grid;grid-template-columns:22px 1fr 34px;align-items:center;gap:8px;font-size:12px;color:var(--ink2)}
.meter b{font:600 12px "Barlow Condensed",sans-serif;color:var(--ink3);letter-spacing:.06em}
.meter .track{height:5px;background:var(--raised);border-radius:3px;overflow:hidden}
.meter .fill{height:100%;background:var(--accent);border-radius:3px}
.meter .fill.dim{background:var(--hair2)}
.legend{display:flex;gap:14px;flex-wrap:wrap;font-size:12px;color:var(--ink2);margin-bottom:4px}
.legend span{display:inline-flex;align-items:center;gap:6px}
.legend i{width:10px;height:10px;border-radius:50%;background:var(--c);display:inline-block}
.legend i.tri{width:0;height:0;border-radius:0;background:none;border-left:5px solid transparent;border-right:5px solid transparent;border-top:9px solid var(--c)}
.legend i.band{width:14px;height:10px;border-radius:2px;background:var(--accent-soft);border:1px solid var(--accent)}
svg text{font-family:"JetBrains Mono",ui-monospace,monospace;font-size:10.5px;fill:var(--ink3)}
svg .axis{stroke:var(--hair);stroke-width:1}
svg .grid{stroke:var(--hair);stroke-width:1}
.chart{position:relative}
.tip{position:absolute;pointer-events:none;background:var(--raised);border:1px solid var(--hair2);color:var(--ink);font:12px "JetBrains Mono",monospace;padding:5px 8px;border-radius:3px;white-space:nowrap;display:none;z-index:3}
.empty{color:var(--ink3);font-size:13px;padding:8px 0}
.mini h4,.foot h4{margin:0 0 6px;font:600 12px "Barlow Condensed",sans-serif;letter-spacing:.12em;text-transform:uppercase;color:var(--ink3)}
table{border-collapse:collapse;width:100%;font-size:12.5px}
td,th{padding:3px 5px;text-align:left;border-bottom:1px solid var(--hair);vertical-align:top}
td:first-child{white-space:nowrap}
th{font:600 11px "Barlow Condensed",sans-serif;letter-spacing:.1em;text-transform:uppercase;color:var(--ink3)}
td.n{text-align:right}
.toggle{background:none;border:none;color:var(--accent);font:inherit;font-size:12.5px;cursor:pointer;padding:4px 0}
.badge{display:inline-block;padding:1px 7px;border-radius:999px;border:1px solid var(--hair2);font:600 11px "Barlow Condensed",sans-serif;letter-spacing:.08em;text-transform:uppercase;color:var(--ink2);margin-left:8px;vertical-align:middle}
.about{display:grid;grid-template-columns:minmax(0,68ch) minmax(280px,380px);gap:40px;padding:18px 4px 30px}
.about .col{min-width:0}
.about h2{font:600 34px/1.05 "Barlow Condensed",sans-serif;margin:6px 0 10px;text-wrap:balance}
.about .lede{font-size:18px;line-height:1.45;color:var(--ink2);margin:0 0 18px}
.about h3{font:600 13px "Barlow Condensed",sans-serif;letter-spacing:.14em;text-transform:uppercase;color:var(--ink3);margin:26px 0 8px}
.about p,.about li{font-size:16px;line-height:1.55;color:var(--ink)}
.about ul{padding-left:20px;margin:6px 0}
.about li{margin:4px 0}
.about .fine{font-size:13px;color:var(--ink3);margin-top:28px}
.about table.tiers td{padding:8px 8px 10px 0;font-size:14px;line-height:1.45;color:var(--ink2)}
.about table.tiers td:first-child{padding-right:12px;white-space:nowrap}
.about table.tiers b{color:var(--ink);display:block;margin-bottom:2px}
.side-col h3:first-child{margin-top:8px}
@media (max-width:1000px){.about{grid-template-columns:1fr}}
.card .foot .specrow{grid-column:1 / -1}
.ewhist{margin-top:6px;max-width:640px}
.ewhist td,.ewhist th{padding:2px 6px}
.over{background:var(--surface);border:1px solid var(--hair);border-radius:4px;padding:6px 16px 12px;margin:8px 0 18px}
.over svg{display:block;max-height:420px}
.tgt:hover{stroke:var(--ink)}
.cuts{display:flex;gap:8px}
.cuts figure{margin:0;flex:1;min-width:0}
.cuts img{width:100%;height:auto;aspect-ratio:1;display:block;border:1px solid var(--hair);border-radius:2px;background:var(--raised);image-rendering:auto}
.cuts figcaption{font:600 10.5px "Barlow Condensed",sans-serif;letter-spacing:.08em;text-transform:uppercase;color:var(--ink3);margin-top:3px;text-align:center}
.spec-slot{border:1px dashed var(--hair2);border-radius:4px;min-height:90px;display:flex;align-items:center;justify-content:center;color:var(--ink3);font-size:13px;text-align:center;padding:10px}
footer{max-width:1500px;margin:20px auto 0;padding:14px 22px;color:var(--ink3);font-size:12.5px;border-top:1px solid var(--hair)}
@media (max-width:1100px){.card{grid-template-columns:200px 1fr}.card .side{grid-column:1 / -1;border-left:none;border-top:1px solid var(--hair);flex-direction:row;flex-wrap:wrap}.card .foot{grid-template-columns:1fr}}
@media (max-width:700px){.card{grid-template-columns:1fr}.card .id{border-right:none;border-bottom:1px solid var(--hair)}}
</style>

<header>
  <div class="bar">
    <h1>CLAGN Night Sheet<small>Palomar 200-inch · NGPS · 2026B</small></h1>
    <div class="tabs" id="tabs" role="group" aria-label="Night"></div>
    <div class="chips" id="chips" role="group" aria-label="Tiers"></div>
    <input type="search" id="q" placeholder="Search name…" aria-label="Search by name">
  </div>
  <div class="night" id="nightinfo"></div>
</header>
<main id="main"></main>
<footer id="foot"></footer>

<script type="application/json" id="data">__DATA__</script>
<script>
const D = JSON.parse(document.getElementById('data').textContent);
const TIERC = {T1:'var(--t1)',T2:'var(--t2)',T3:'var(--t3)',T4:'var(--t4)'};
const state = {night:'sep23', tiers:new Set(['T1','T2','T3','T4']), q:''};
try{ const s=JSON.parse(localStorage.getItem('clagn_ns')||'{}'); if(s.night) state.night=s.night; }catch(e){}
function save(){ try{ localStorage.setItem('clagn_ns', JSON.stringify({night:state.night})); }catch(e){} }

/* ---------- header controls ---------- */
const tabs = document.getElementById('tabs');
[...Object.entries(D.nights).map(([k,v])=>[k,v.label]), ['all','All nights'], ['about','About']].forEach(([k,lab])=>{
  const b=document.createElement('button'); b.className='tab'; b.textContent=lab; b.dataset.k=k;
  b.setAttribute('aria-pressed', String(state.night===k));
  b.onclick=()=>{state.night=k; save(); render();}; tabs.appendChild(b);
});
const chips=document.getElementById('chips');
Object.entries(D.tiers).forEach(([k,lab])=>{
  const b=document.createElement('button'); b.className='chip'; b.style.setProperty('--c',TIERC[k]); b.innerHTML=`<i></i>${k}`; b.title=lab;
  b.setAttribute('aria-pressed','true'); b.onclick=()=>{ state.tiers.has(k)?state.tiers.delete(k):state.tiers.add(k); b.setAttribute('aria-pressed',String(state.tiers.has(k))); render(); };
  chips.appendChild(b);
});
document.getElementById('q').addEventListener('input', e=>{state.q=e.target.value.trim().toLowerCase(); render();});

/* ---------- helpers ---------- */
const fmt=(v,d=2)=> (v==null||Number.isNaN(v))?'—':Number(v).toFixed(d);
const esc=s=>String(s??'').replace(/[&<>"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
function mjdToYear(m){ return 2000 + (m-51544.5)/365.25; }
function yearToMjd(y){ return 51544.5 + (y-2000)*365.25; }
function lin(d0,d1,r0,r1){ return v=> r0 + (v-d0)*(r1-r0)/(d1-d0); }
function niceTicks(a,b,n){ const span=b-a, step0=span/n, p=Math.pow(10,Math.floor(Math.log10(step0))), c=step0/p, step=(c<1.5?1:c<3.5?2:c<7.5?5:10)*p; const out=[]; for(let v=Math.ceil(a/step)*step; v<=b+1e-9; v+=step) out.push(+v.toFixed(6)); return out; }

/* ---------- light curve chart (two panels, shared time axis) ---------- */
function lightCurve(t){
  const W=720, padL=54, padR=14, hTop=150, hBot=110, gap=26, top=8;
  const runs=Object.values(D.nights).map(n=>n.mjd);
  let all=[];
  ['g','r'].forEach(k=>(t.ztf[k]||[]).forEach(p=>all.push(p[0])));
  Object.values(t.wise).forEach(s=>s.forEach(p=>all.push(p[0])));
  (t.neo||[]).forEach(p=>all.push(p[0]));
  t.epochs.forEach(e=>all.push(e.mjd));
  if(!all.length) return '<div class="empty">No light curves fetched for this target yet.</div>';
  let x0=Math.min(...all)-120, x1=Math.max(...all, Math.max(...runs))+120;
  const X=lin(x0,x1,padL,W-padR);
  const H=top+hTop+gap+hBot+34;
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="Light curves of ${esc(t.name)}">`;
  // run bands
  runs.forEach(m=>{ s+=`<rect x="${X(m)-2}" y="${top}" width="4" height="${hTop+gap+hBot}" fill="var(--accent)" opacity=".28"/>`; });
  // epochs markers
  t.epochs.forEach(e=>{ const x=X(e.mjd), c=e.src==='DESI'?'var(--desi)':'var(--sdss)'; s+=`<line x1="${x}" y1="${top}" x2="${x}" y2="${top+hTop+gap+hBot}" stroke="${c}" stroke-width="1" opacity=".55"/><polygon points="${x-4},${top} ${x+4},${top} ${x},${top+7}" fill="${c}"/>`; });
  // top panel: ZTF mags (inverted)
  const zpts=[...(t.ztf.g||[]).map(p=>({k:'g',m:p[0],v:p[1],e:p[2]})),...(t.ztf.r||[]).map(p=>({k:'r',m:p[0],v:p[1],e:p[2]}))];
  if(zpts.length){
    const vals=zpts.map(p=>p.v).sort((a,b)=>a-b); const lo=vals[Math.floor(vals.length*.02)], hi=vals[Math.ceil(vals.length*.98)-1];
    const y0=hi+0.15, y1=lo-0.15; const Y=lin(y1,y0,top,top+hTop); // brighter (smaller mag) at top
    niceTicks(y1,y0,4).forEach(v=>{ s+=`<line class="grid" x1="${padL}" x2="${W-padR}" y1="${Y(v)}" y2="${Y(v)}"/><text x="${padL-6}" y="${Y(v)+3.5}" text-anchor="end">${v.toFixed(1)}</text>`; });
    s+=`<text x="${padL-6}" y="${top-1}" text-anchor="end" style="fill:var(--ink2)">mag</text>`;
    zpts.forEach(p=>{ if(p.v<y1||p.v>y0) return; s+=`<circle cx="${X(p.m)}" cy="${Y(p.v)}" r="2.6" fill="var(--${p.k})" opacity=".85"><title>ZTF ${p.k}  MJD ${p.m.toFixed(1)}  ${p.v.toFixed(2)} ± ${p.e.toFixed(2)}</title></circle>`; });
    if(t.r_last!=null){ /* reference: r at last spectrum, if computed */ }
  } else {
    s+=`<text x="${padL+8}" y="${top+hTop/2}" style="fill:var(--ink3)">no ZTF light curve</text>`;
  }
  // bottom panel: unWISE flux
  const w1=t.wise.W1||[], w2=t.wise.W2||[], neo=t.neo||[];
  const wb=top+hTop+gap;
  if(w1.length||neo.length){
    const fl=[...w1.map(p=>p[1]),...w2.map(p=>p[1]),...neo.map(p=>p[1])]; const lo=Math.min(...fl)*0.9, hi=Math.max(...fl)*1.08;
    const Y=lin(lo,hi,wb+hBot,wb);
    niceTicks(lo,hi,3).forEach(v=>{ s+=`<line class="grid" x1="${padL}" x2="${W-padR}" y1="${Y(v)}" y2="${Y(v)}"/><text x="${padL-6}" y="${Y(v)+3.5}" text-anchor="end">${v.toFixed(2)}</text>`; });
    s+=`<text x="${padL-6}" y="${wb-3}" text-anchor="end" style="fill:var(--ink2)">mJy</text>`;
    const path=arr=>arr.map((p,i)=>`${i?'L':'M'}${X(p[0]).toFixed(1)},${Y(p[1]).toFixed(1)}`).join('');
    if(w2.length) s+=`<path d="${path(w2)}" fill="none" stroke="var(--w1)" stroke-width="1.5" opacity=".45"/>`;
    s+=`<path d="${path(w1)}" fill="none" stroke="var(--w1)" stroke-width="2" stroke-linejoin="round"/>`;
    w1.forEach(p=>{ s+=`<circle cx="${X(p[0])}" cy="${Y(p[1])}" r="3.5" fill="var(--w1)" stroke="var(--surface)" stroke-width="2"><title>unWISE W1  MJD ${p[0].toFixed(0)}  ${p[1].toFixed(3)} ± ${p[2].toFixed(3)} mJy</title></circle>`; });
    if(neo.length){ s+=`<path d="${path(neo)}" fill="none" stroke="var(--w1)" stroke-width="1.5" stroke-dasharray="3 3" opacity=".8"/>`;
      neo.forEach(p=>{ s+=`<circle cx="${X(p[0])}" cy="${Y(p[1])}" r="3.5" fill="var(--surface)" stroke="var(--w1)" stroke-width="2"><title>NEOWISE-R visit  MJD ${p[0].toFixed(0)}  ${p[1].toFixed(3)} ± ${p[2].toFixed(3)} mJy</title></circle>`; }); }
  } else {
    s+=`<text x="${padL+8}" y="${wb+hBot/2}" style="fill:var(--ink3)">no unWISE light curve</text>`;
  }
  // time axis (years)
  const yA=top+hTop+gap+hBot;
  s+=`<line class="axis" x1="${padL}" x2="${W-padR}" y1="${yA}" y2="${yA}"/>`;
  niceTicks(mjdToYear(x0),mjdToYear(x1),8).forEach(y=>{ const x=X(yearToMjd(y)); if(x<padL||x>W-padR) return; s+=`<line class="axis" x1="${x}" x2="${x}" y1="${yA}" y2="${yA+5}"/><text x="${x}" y="${yA+17}" text-anchor="middle">${y}</text>`; });
  s+=`<text x="${W-padR}" y="${yA+31}" text-anchor="end">year</text>`;
  s+='</svg>';
  return s;
}

/* ---------- manifold thumbnail ---------- */
function manifold(t){
  const M=D.manifold, W=260, H=200, p=6;
  const X=lin(M.xlim[0],M.xlim[1],p,W-p), Y=lin(M.ylim[0],M.ylim[1],H-p,p);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="Position of ${esc(t.name)} on the W1 manifold">`;
  M.rects.forEach(r=>{ s+=`<rect x="${X(r.x0)}" y="${Y(r.y1)}" width="${X(r.x1)-X(r.x0)}" height="${Y(r.y0)-Y(r.y1)}" fill="${r.kind==='zeltyn'?'var(--accent)':'var(--t3)'}" opacity=".13"/>`; });
  s+=`<g fill="var(--ink3)" opacity=".45">`+M.A.map(a=>`<circle cx="${X(a[0]).toFixed(1)}" cy="${Y(a[1]).toFixed(1)}" r="${a[2]?1.9:1.3}" ${a[2]?'fill="var(--ink2)"':''}/>`).join('')+`</g>`;
  s+=`<g fill="none" stroke="var(--desi)" stroke-width="1" opacity=".7">`+M.E.map(a=>`<circle cx="${X(a[0]).toFixed(1)}" cy="${Y(a[1]).toFixed(1)}" r="2.2"/>`).join('')+`</g>`;
  s+=`<g fill="var(--desi)" opacity=".8">`+M.Z.map(a=>`<circle cx="${X(a[0]).toFixed(1)}" cy="${Y(a[1]).toFixed(1)}" r="1.8"/>`).join('')+`</g>`;
  if(t.ux!=null){ const x=X(t.ux), y=Y(t.uy); s+=`<circle cx="${x}" cy="${y}" r="7" fill="${TIERC[t.tier]}" stroke="var(--surface)" stroke-width="2.5"/>`; }
  else s+=`<text x="${W/2}" y="${H/2}" text-anchor="middle" style="fill:var(--ink3)">not projected</text>`;
  s+='</svg>';
  return s;
}

/* ---------- overview manifold (top of page) ---------- */
function overview(list, nk){
  const M=D.manifold, W=980, H=400, p=14;
  const X=lin(M.xlim[0],M.xlim[1],p,W-p), Y=lin(M.ylim[0],M.ylim[1],H-p,p);
  const tri=(x,y,up,c)=>{ x=+x; y=+y; return `<polygon points="${x},${(up?y-4.5:y+4.5).toFixed(1)} ${(x-4).toFixed(1)},${(up?y+3:y-3).toFixed(1)} ${(x+4).toFixed(1)},${(up?y+3:y-3).toFixed(1)}" fill="${c}"/>`; };
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="W1 variability manifold with changing-look AGN and tonight's targets">`;
  M.rects.forEach(r=>{ s+=`<rect x="${X(r.x0)}" y="${Y(r.y1)}" width="${X(r.x1)-X(r.x0)}" height="${Y(r.y0)-Y(r.y1)}" fill="${r.kind==='zeltyn'?'var(--accent)':'var(--t3)'}" opacity=".12"/>`; });
  s+=`<g fill="var(--ink3)" opacity=".4">`+M.A.filter(a=>!a[2]).map(a=>`<circle cx="${X(a[0]).toFixed(1)}" cy="${Y(a[1]).toFixed(1)}" r="2.2"/>`).join('')+`</g>`;
  s+=`<g opacity=".9">`+M.A.filter(a=>a[2]).map(a=>tri(X(a[0]).toFixed(1),Y(a[1]).toFixed(1),a[2]!==2,'var(--ink2)')).join('')+`</g>`;
  s+=`<g fill="none" stroke="var(--desi)" stroke-width="1.2" opacity=".85">`+M.E.map(a=>`<circle cx="${X(a[0]).toFixed(1)}" cy="${Y(a[1]).toFixed(1)}" r="3.2"/>`).join('')+`</g>`;
  s+=`<g fill="var(--desi)">`+M.Z.map(a=>`<circle cx="${X(a[0]).toFixed(1)}" cy="${Y(a[1]).toFixed(1)}" r="3"/>`).join('')+`</g>`;
  list.filter(t=>t.ux!=null).forEach(t=>{ const n=nk==='all'?null:t.nights.find(x=>x.night===nk); const primary=nk==='all'||(n&&n.rank>0);
    s+=`<circle class="tgt" data-name="${esc(t.name)}" cx="${X(t.ux)}" cy="${Y(t.uy)}" r="${primary?6.5:4.5}" fill="${TIERC[t.tier]}" stroke="var(--surface)" stroke-width="2" opacity="${primary?1:.6}" style="cursor:pointer"><title>${esc(t.jname||t.name)} · ${t.tier}${n?` · rank ${n.rank||'backup'}`:''} · z ${fmt(t.z,2)} · r ${fmt(t.rmag,1)}</title></circle>`; });
  s+='</svg>';
  return `<section class="over"><div class="section" style="margin-top:6px">W1 manifold · where the changing-look AGN sit</div>${s}
    <div class="legend" style="margin-top:6px"><span><i style="--c:var(--ink3);opacity:.6"></i>Sample A (Hemmati+2026)</span><span><i class="tri" style="--c:var(--ink2)"></i>literature turn-on ▲ / turn-off ▼</span><span><i style="--c:var(--desi)"></i>Zeltyn+2024 CL-AGN</span><span><i style="--c:transparent;border:1.5px solid var(--desi)"></i>Zeltyn+2024 EVQ</span><span><i class="band"></i>Zeltyn-enriched bins</span><span><i class="band" style="background:rgba(25,158,112,.18);border-color:var(--t3)"></i>turn-on/off-enriched bins</span><span style="margin-left:auto">large dots: ${nk==='all'?'all':'this night’s'} targets by tier (click to jump)</span></div></section>`;
}

/* ---------- About tab ---------- */
function about(){
  const S=D.stats||{}, N=D.nights, n=v=>v==null?'—':Number(v).toLocaleString();
  const tierRows=[
    ['T1','Manifold-selected discovery targets',`Bright SDSS quasars (z &lt; 0.8) in the run's RA windows whose WISE W1 light curve places them in the changing-look-enriched part of the manifold and which have never been reported to change. Each has an archival SDSS spectrum from 2000–2018 as the baseline. ${n(S.n_T1)} qualify; the nights take the best per hour of telescope time.`,'var(--t1)'],
    ['T2','EVQs completing a transition',`Extremely variable quasars from Zeltyn et al. (2024) whose broad lines had dimmed strongly by 2020–21 without disappearing, and which sit in the CLAGN region. A 2026 spectrum tests whether the transition completed.`,'var(--t2)'],
    ['T3','Confirmed CLAGNs, revisited',`Spectroscopically confirmed SDSS-V changing-look AGN (${n(S.n_zeltyn_clagn)} in the sample). A few of the brightest per night, to test for state reversal and to check that the manifold region really predicts change.`,'var(--t3)'],
    ['T4','Controls',`Bright, photometrically quiet quasars from the low-variability side of the manifold, outside both enriched regions, where no spectral change is expected. Needed to claim that the region predicts change rather than that all quasars change.`,'var(--t4)']];
  return `<section class="about">
  <div class="col">
    <h2>Catching AGN before they turn off</h2>
    <p class="lede">A spectroscopic follow-up of changing-look AGN candidates with the Palomar 200-inch and NGPS in 2026B, selected from the shape of their mid-infrared light curves rather than from a prior spectral change.</p>
    <h3>Why</h3>
    <p>Changing-look AGN (CLAGN) gain or lose their broad emission lines on timescales of months to years, together with large continuum changes. They are rare, roughly 0.4–1.25 % of re-observed quasars, and almost all have been found by chance in repeat spectroscopy. A purely photometric way to pick them out would let surveys like LSST target them deliberately.</p>
    <p>Hemmati et al. (2026, ApJ 998, 130) built a low-dimensional map, a UMAP manifold, of the WISE/NEOWISE W1 light curves of ~${n(S.n_sampleA)} AGN at z &lt; 1 (Sample A), without using any labels. Known turn-on and turn-off CLAGNs from the literature occupy distinct parts of that map. Projecting the ${n(S.n_zeltyn)} SDSS-V CLAGNs and extremely variable quasars of Zeltyn et al. (2024), which were not used in training, shows them concentrating in a compact region too: the Zeltyn EVQs land there at 6.3× the rate of the general sample. Interestingly the SDSS-V CLAGNs and the literature CLAGNs sit in <em>different</em> regions. Their mean W1 curves show why: the literature objects changed before or around the start of WISE (2010–13), the SDSS-V ones faded slowly through the whole decade. Both regions are used here.</p>
    <h3>The runs</h3>
    <p>${Object.values(N).map(v=>`<b>${v.date}</b>, ${v.part} (${v.window}; moon ${v.moon}, ${v.moon_note})`).join('; ')}. All three are bright time, so the ranking favours bright targets and large moon separation, and every exposure estimate includes the moon.</p>
    <h3>How targets are ranked</h3>
    <p>Every candidate gets <span class="mono">priority = B · S · (M + P)</span>, plus 0.5 if the pipeline class of its archival spectra changed between epochs.</p>
    <ul>
      <li><b>M</b>, manifold likeness: the larger of the local density of Zeltyn CLAGNs around the object on the map and the fraction of literature CLAGNs among its 50 nearest neighbours, each normalised so 1 is the threshold; capped at 2.</li>
      <li><b>P</b>, photometric change since the last spectrum: the largest of the W1 change to 2020 (unWISE), to 2024 (NEOWISE-R) and the ZTF r change to 2025, in units where 1 means a factor 1.5; capped at 2. This is the term that says <em>something is happening now</em>.</li>
      <li><b>S</b>, staleness: 1 if the last spectrum is at least three years old, else 0.3.</li>
      <li><b>B</b>, brightness: 1.0, 0.8, 0.5, 0.25 for r brighter than 18.5, 19, 19.5 and fainter. There is no magnitude cut.</li>
    </ul>
    <p>Per night, the priority is multiplied by a moon-distance weight (1 beyond 60°, 0.7 at 40–60°, 0.4 at 30–40°; closer is excluded) and divided by the estimated exposure time, so the list is ordered by <em>science return per hour</em>. The exposure model scales the proposal's NGPS ETC point (r = 18.5, dark: 9 min at S/N 10) to S/N ≈ 7 on the diagnostic broad line: Hα for z ≤ 0.55, else Hβ, with the bright moon adding about 1 mag of sky in the red and 2 in the blue-green. Roughly 4 targets per hour at r = 18, 2 per hour at r = 19 with Hα, 1 per hour if Hβ is needed. Each night's usable hours are filled in that order, with floors so that every tier is represented and caps on revisits and controls. Objects within 2″ of a published CLAGN or of a Zeltyn object are excluded from Tier 1.</p>
    <h3>What a card shows</h3>
    <ul>
      <li><b>Light curves:</b> ZTF g and r nightly medians (2018–2025), unWISE W1 and W2 (2010–2020) and NEOWISE-R W1 visits (2014–2024) in mJy; every archival spectral epoch as a marker (grey SDSS, gold DESI) and the three run dates as gold bands.</li>
      <li><b>Manifold:</b> the object's position among Sample A (grey), the literature turn-on/off objects (▲ ▼), the Zeltyn CLAGNs and EVQs (gold) and the two enriched regions (shaded).</li>
      <li><b>Archival spectra:</b> every SDSS epoch from DR19, including SDSS-V through 2022–23, and DESI DR1, with the pipeline class and redshift where available.</li>
      <li><b>What NGPS sees:</b> the lines that fall inside 3200–10400 Å at the object's redshift, against the narrower SDSS range.</li>
      <li><b>Imaging:</b> 64″ cutouts from SDSS (gri) and the ZTF g and r reference images, north up and east left.</li>
      <li><b>NGPS spectrum:</b> empty until observed; drop <span class="mono">data/ngps_spectra/&lt;name&gt;.csv</span> in the repository and regenerate.</li>
    </ul>
    <h3>Status and caveats</h3>
    <p>Lists are current as of ${D.generated}. Candidate pool: ${n(S.n_pool)} DR16 quasars in the RA windows, ${n(S.n_pool_projected)} with WISE light curves and manifold positions; ${n(S.n_lit)} published CLAGN positions used for exclusion. The exposure model is a scaling and should be checked against the NGPS ETC with the moon phase set. SDSS-V spectra taken after the DR19 cutoff (2023 onward) are not public; a check against the collaboration's internal spAll will de-prioritise anything recently re-observed. The unWISE record ends in December 2020 and NEOWISE in early 2024, so the most recent mid-infrared behaviour is unknown; ZTF covers the optical to 2025.</p>
    <p class="fine">Pipeline: github repository CLAGN (scripts 00–08), data from IRSA (unWISE, NEOWISE-R, ZTF), SDSS DR16/DR19, DESI DR1, NOIRLab Data Lab, MAST/PS1. Contact: S. Hemmati (Caltech/IPAC).</p>
  </div>
  <aside class="col side-col">
    <h3>Tiers</h3>
    <table class="tiers">${tierRows.map(([k,t,d,c])=>`<tr><td><span class="tier" style="--c:${c};margin:0"><i></i>${k}</span></td><td><b>${t}</b><div>${d}</div></td></tr>`).join('')}</table>
    <h3>Tonight's counts</h3>
    <table><thead><tr><th>night</th><th>primary</th><th>backups</th><th>by tier</th></tr></thead><tbody>${Object.values(N).map(v=>`<tr><td>${v.label}</td><td class="n mono">${v.n_primary}</td><td class="n mono">${v.n_backup}</td><td>${Object.entries(v.counts||{}).sort().map(([k,c])=>`${k} ${c}`).join(' · ')}</td></tr>`).join('')}</tbody></table>
    <h3>Proposal in one line</h3>
    <p>Are objects in the CLAGN-enriched part of a photometric variability manifold more likely to show changing-look behaviour? A spectrum now, against a 2000–2018 SDSS baseline, answers it per object; the controls answer it for the method.</p>
  </aside>
</section>`;
}

/* ---------- archival spectra overlay + EW history ---------- */
const SPECC=['var(--ink3)','var(--t1)','var(--accent)','var(--t2)'];
function specPanel(t){
  const S=t.spec; if(!S||!S.epochs.length) return `<div class="empty">No archival spectrum file fetched yet.</div>`;
  const W=900,H=240,pL=46,pR=12,pT=30,pB=34; const w=S.wave; const X=lin(w[0],w[w.length-1],pL,W-pR);
  const vals=[]; S.epochs.forEach(e=>e.flux.forEach(v=>{ if(v!=null) vals.push(v); })); vals.sort((a,b)=>a-b);
  const lo=Math.min(0,vals[Math.floor(vals.length*.01)]), hi=vals[Math.ceil(vals.length*.995)-1]*1.05; const Y=lin(lo,hi,H-pB,pT);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="Archival spectra of ${esc(t.name)}">`;
  niceTicks(lo,hi,4).forEach(v=>{ s+=`<line class="grid" x1="${pL}" x2="${W-pR}" y1="${Y(v)}" y2="${Y(v)}"/><text x="${pL-5}" y="${Y(v)+3.5}" text-anchor="end">${v.toFixed(0)}</text>`; });
  s+=`<text x="${pL+4}" y="${pT+10}" style="fill:var(--ink2)">F<tspan baseline-shift="sub" font-size="8">λ</tspan> · 10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹</text>`;
  // rest-frame line markers, labels on two rows so Hβ and [O III] never collide
  if(t.z!=null){ [['Mg II',2798,0],['Hβ',4861,1],['[O III]',5007,0],['Hα',6563,1]].forEach(([nm,l0,row])=>{ const lo_=l0*(1+t.z); if(lo_<w[0]||lo_>w[w.length-1]) return; const x=X(lo_); s+=`<line x1="${x}" x2="${x}" y1="${pT-2}" y2="${H-pB}" stroke="var(--hair2)" stroke-width="1"/><text x="${x}" y="${row?pT-4:pT-16}" text-anchor="middle" style="fill:var(--ink3);font-size:10px">${nm}</text>`; }); }
  S.epochs.forEach((e,i)=>{ let d='',pen=false; e.flux.forEach((v,k)=>{ if(v==null){pen=false;return;} const x=X(w[k]).toFixed(1), y=Y(Math.min(v,hi)).toFixed(1); d+=(pen?'L':'M')+x+','+y; pen=true; }); s+=`<path d="${d}" fill="none" stroke="${SPECC[i%SPECC.length]}" stroke-width="${i===S.epochs.length-1?1.6:1.1}" opacity="${i===0?.9:.95}"/>`; });
  s+=`<line class="axis" x1="${pL}" x2="${W-pR}" y1="${H-pB}" y2="${H-pB}"/>`; niceTicks(w[0],w[w.length-1],7).forEach(v=>{ const x=X(v); if(x<pL||x>W-pR) return; s+=`<text x="${x}" y="${H-pB+14}" text-anchor="middle">${v}</text>`; }); s+=`<text x="${W-pR}" y="${H-3}" text-anchor="end">observed wavelength, Å</text>`;
  s+='</svg>';
  const leg=`<div class="legend" style="margin-top:4px">${S.epochs.map((e,i)=>`<span><i style="--c:${SPECC[i%SPECC.length]};border-radius:0;height:3px;width:16px"></i>${esc(e.label)} <span style="color:var(--ink3)">${esc(e.cls)}${e.ew_hb!=null?` · EW(Hβ) ${e.ew_hb.toFixed(0)} Å`:''}${e.ew_ha!=null?` · EW(Hα) ${e.ew_ha.toFixed(0)} Å`:''}</span></span>`).join('')}</div>`;
  const hist=S.history.length>1?`<table class="ewhist"><thead><tr><th>epoch</th><th>source</th><th>class</th><th>EW Hβ</th><th>EW Hα</th></tr></thead><tbody>${S.history.map(h=>`<tr><td class="mono">${h.date}</td><td>${h.src}${h.coadd?' coadd':''} <span style="color:var(--ink3)">${esc(h.prog)}</span></td><td>${esc(h.cls)}</td><td class="n mono">${h.ew_hb==null?'—':h.ew_hb}</td><td class="n mono">${h.ew_ha==null?'—':h.ew_ha}</td></tr>`).join('')}</tbody></table>`:'';
  return s+leg+hist;
}

/* ---------- NGPS wavelength ruler ---------- */
function ruler(t){
  const W=560, H=84, p=28, a=3200, b=10400; const X=lin(a,b,p,W-p);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="Emission lines in the NGPS range for ${esc(t.name)}">`;
  // title row, then the NGPS bar, then the thinner SDSS/BOSS bar for comparison, then tick labels
  s+=`<text x="${p}" y="11" style="fill:var(--ink2)">NGPS 3200–10400 Å</text>`;
  s+=`<text x="${W-p}" y="11" text-anchor="end">thin bar: SDSS / BOSS 3800–9200 Å</text>`;
  s+=`<rect x="${X(a)}" y="30" width="${X(b)-X(a)}" height="14" rx="2" fill="var(--raised)"/>`;
  s+=`<rect x="${X(3800)}" y="48" width="${X(9200)-X(3800)}" height="3" rx="1" fill="var(--hair2)"/>`;
  [3200,4000,5000,6000,7000,8000,9000,10400].forEach(w=>{ s+=`<text x="${X(w)}" y="76" text-anchor="middle">${w}</text>`; });
  // line labels alternate above / below the bar so neighbours (Hβ, [O III]) never collide
  t.lines.forEach((l,i)=>{ if(l.obs<2600||l.obs>11400) return; const x=Math.min(Math.max(X(l.obs),p-14),W-p+14); const ok=l.inrange, above=(i%2===0);
    s+=`<line x1="${x}" x2="${x}" y1="27" y2="47" stroke="${ok?'var(--accent)':'var(--ink3)'}" stroke-width="2"/><text x="${x}" y="${above?24:62}" text-anchor="middle" style="fill:${ok?'var(--ink)':'var(--ink3)'};font-size:10px">${esc(l.name)}</text>`; });
  s+='</svg>';
  return s;
}

/* ---------- card ---------- */
function meter(label,val,max,title){ const f=val==null?0:Math.min(1,val/max); return `<div class="meter" title="${esc(title)}"><b>${label}</b><div class="track"><div class="fill ${val==null?'dim':''}" style="width:${(f*100).toFixed(0)}%"></div></div><span class="mono">${fmt(val,2)}</span></div>`; }
function card(t, nightKey){
  const nk = nightKey==='all' ? (t.nights.slice().sort((a,b)=>(b.prio??0)-(a.prio??0))[0]?.night) : nightKey;
  const n = t.nights.find(x=>x.night===nk) || {};
  const isBackup = n.rank===0;
  const epochs = t.epochs;
  const epRows = epochs.map(e=>`<tr><td class="mono">${e.date}</td><td>${esc(e.src)} <span style="color:var(--ink3)">${esc(e.prog)}</span></td><td>${esc(e.cls)}${e.sub?` <span style="color:var(--ink3)">${esc(e.sub)}</span>`:''}</td><td class="n mono">${e.z==null?'—':e.z.toFixed(3)}</td></tr>`);
  const shown = epRows.slice(0,6).join(''), hidden = epRows.slice(6).join('');
  const change = t.dr_ref!=null ? `Δr since last spectrum ${t.dr_ref>0?'+':''}${fmt(t.dr_ref)} mag` : (t.w1_ratio!=null ? `W1 now / at spectrum = ${fmt(t.w1_ratio)}` : 'no change measure yet');
  return `<article class="card" id="c-${esc(t.name)}">
    <div class="id">
      <div class="rank">${isBackup?'<small>backup</small>':`${n.rank??'—'}<small>${D.nights[nk]?.label??''}</small>`}</div>
      <div class="tier" style="--c:${TIERC[t.tier]}"><i></i>${esc(D.tiers[t.tier]||t.tier)}</div>
      <h2>${esc(t.jname||t.name)}</h2>
      <div class="coords mono">${esc(t.sex)}${t.jname?` <span style="color:var(--ink3)">· ${esc(t.name)}</span>`:''}</div>
      <div class="kv">
        <b>z</b><span class="v mono">${fmt(t.z,3)}</span>
        <b>r</b><span class="v mono">${fmt(t.rmag,1)}</span>
        <b>hours</b><span class="v mono">${fmt(n.hrs,1)} <span style="color:var(--ink3)">best airmass ${fmt(n.minx,2)}</span></span>
        <b>moon</b><span class="v mono">${fmt(n.moonsep,0)}°</span>
        <b>exposure</b><span class="v mono">${n.plan||(n.texp!=null?`~${fmt(n.texp,0)} min`:'—')} <span style="color:var(--ink3)">${n.pph!=null?`${fmt(n.pph,1)} /h`:''}</span></span>
        <b>last spec</b><span class="v">${t.mjd_last!=null?`${fmt(t.yrs,1)} yr ago`:'none'} <span style="color:var(--ink3)">${esc(t.last_class)}</span></span>
        <b>trend</b><span class="v">${esc(t.trend)}</span>
      </div>
      <div class="meters">
        ${meter('M',t.M,2,'manifold CLAGN-likeness (1 = threshold)')}
        ${meter('P',t.P,2,'photometric change since last spectrum (1 = factor 1.5)')}
        ${meter('S',t.S,1,'staleness of last spectrum')}
        ${meter('B',t.B,1,'brightness weight')}
        <div class="meter"><b>Σ</b><span></span><span class="mono" style="color:var(--ink)">${fmt(n.prio??t.priority,2)}</span></div>
      </div>
    </div>
    <div class="chart">
      <div class="legend"><span><i style="--c:var(--g)"></i>ZTF g</span><span><i style="--c:var(--r)"></i>ZTF r</span><span><i style="--c:var(--w1)"></i>unWISE W1 (W2 faint)</span><span><i style="--c:transparent;border:2px solid var(--w1)"></i>NEOWISE-R visits to 2024</span><span><i class="tri" style="--c:var(--sdss)"></i>SDSS epoch</span><span><i class="tri" style="--c:var(--desi)"></i>DESI epoch</span><span><i class="band"></i>NGPS run</span></div>
      ${lightCurve(t)}
      <div style="color:var(--ink3);font-size:12.5px;margin-top:2px">${esc(change)}${t.mjd_last_ztf?` · ZTF through ${mjdToYear(t.mjd_last_ztf).toFixed(1)}`:''}</div>
    </div>
    <div class="side">
      ${Object.keys(t.cut||{}).length?`<div class="mini"><h4>Imaging · 40″ · N up, E left</h4><div class="cuts">${(t.cut.ps1_g||t.cut.ps1_r?[['sdss','SDSS gri'],['ps1_g','PS1 g'],['ps1_r','PS1 r']]:[['sdss','SDSS gri'],['ztf_g','ZTF g ref'],['ztf_r','ZTF r ref']]).filter(([k])=>t.cut[k]).map(([k,lab])=>`<figure><img src="${t.cut[k]}" alt="${lab} cutout of ${esc(t.name)}" width="88" height="88"><figcaption>${lab}</figcaption></figure>`).join('')}</div></div>`:''}
      <div class="mini"><h4>W1 manifold</h4>${manifold(t)}
        <div style="font-size:12px;color:var(--ink3)">kNN score ${fmt(t.clagn_score,2)}${t.density!=null?` · Zeltyn density ${fmt(t.density,1)}×`:''}${t.in_zeltyn?' · in Zeltyn region':''}${t.in_clagn?' · in turn-on/off region':''}${t.m_comb!=null?`<br>ZTF+W1 manifold score ${fmt(t.m_comb,2)} ${t.m_comb>=1?'(agrees)':t.m_comb<0.5?'(disagrees)':''}`:''}</div></div>
      <div class="mini"><h4>Archival spectra <span class="badge">${epochs.length}</span></h4>
        ${epochs.length?`<table><thead><tr><th>date</th><th>survey</th><th>class</th><th>z</th></tr></thead><tbody>${shown}</tbody>${hidden?`<tbody hidden>${hidden}</tbody>`:''}</table>${hidden?`<button class="toggle" onclick="const b=this.previousElementSibling.querySelector('tbody[hidden],tbody[data-open]');if(b.hidden){b.hidden=false;b.dataset.open=1;this.textContent='fewer epochs'}else{b.hidden=true;delete b.dataset.open;this.textContent='all ${epochs.length} epochs'}">all ${epochs.length} epochs</button>`:''}`:'<div class="empty">No epoch inventory yet.</div>'}
      </div>
    </div>
    <div class="foot">
      <div class="specrow"><h4>Archival spectra <span style="color:var(--ink3);font-weight:400;text-transform:none;letter-spacing:0">rest-frame EW from fixed windows, positive = emission</span></h4>${specPanel(t)}</div>
      <div><h4>What NGPS sees at z = ${fmt(t.z,3)}</h4>${ruler(t)}<div style="font-size:12.5px;color:var(--ink2)">${esc(t.notes)}</div></div>
      <div><h4>NGPS spectrum</h4>${t.ngps?spectrum(t):`<div class="spec-slot"><div>Not yet observed · planned ${esc(D.nights[nk]?.label??'')}${isBackup?' (backup)':''}<br><span style="font-size:12px">drop data/ngps_spectra/${esc(t.name)}.csv and regenerate</span></div></div>`}</div>
    </div>
  </article>`;
}
function spectrum(t){
  const W=560,H=150,pL=44,pR=10,pT=8,pB=26; const w=t.ngps.map(p=>p[0]), f=t.ngps.map(p=>p[1]);
  const X=lin(Math.min(...w),Math.max(...w),pL,W-pR), lo=Math.min(...f), hi=Math.max(...f), Y=lin(lo,hi,H-pB,pT);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="NGPS spectrum"><path d="${t.ngps.map((p,i)=>`${i?'L':'M'}${X(p[0]).toFixed(1)},${Y(p[1]).toFixed(1)}`).join('')}" fill="none" stroke="var(--accent)" stroke-width="1.2"/>`;
  s+=`<line class="axis" x1="${pL}" x2="${W-pR}" y1="${H-pB}" y2="${H-pB}"/>`; niceTicks(Math.min(...w),Math.max(...w),6).forEach(v=>{ s+=`<text x="${X(v)}" y="${H-8}" text-anchor="middle">${v}</text>`; }); s+='</svg>'; return s;
}

/* ---------- render ---------- */
function render(){
  document.querySelectorAll('.tab').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.k===state.night)));
  const nk=state.night, main=document.getElementById('main'), info=document.getElementById('nightinfo');
  if(nk==='about'){ info.innerHTML=`<div><b>About</b><span class="v">what this program is, how targets are chosen, and how to read a card</span></div>`; main.innerHTML=about(); return; }
  if(nk!=='all'){ const n=D.nights[nk]; const counts=Object.entries(n.counts||{}).sort().map(([k,v])=>`${k} ${v}`).join(' · ');
    info.innerHTML=`<div><b>Night</b><span class="v">${n.date} · ${n.part}</span></div><div><b>Window</b><span class="v mono">${n.window}</span></div><div><b>LST</b><span class="v mono">${n.lst}</span></div><div><b>Moon</b><span class="v">${n.moon} · ${n.moon_pos}</span><br>${n.moon_note}</div><div><b>Slots</b><span class="v">${n.n_primary} primary · ${n.n_backup} backups</span><br>${counts}</div>`;
  } else info.innerHTML=`<div><b>View</b><span class="v">all three nights, best priority first</span></div>`;
  let list=D.targets.filter(t=>state.tiers.has(t.tier) && (!state.q || t.name.toLowerCase().includes(state.q)));
  if(nk!=='all') list=list.filter(t=>t.nights.some(x=>x.night===nk));
  const key=t=>{ const n=nk==='all'?null:t.nights.find(x=>x.night===nk); return n?(n.prio??0):Math.max(...t.nights.map(x=>x.prio??0)); };
  const prim=list.filter(t=>nk==='all'||t.nights.find(x=>x.night===nk).rank>0).sort((a,b)=>{ if(nk!=='all'){return a.nights.find(x=>x.night===nk).rank-b.nights.find(x=>x.night===nk).rank;} return key(b)-key(a); });
  const back=nk==='all'?[]:list.filter(t=>t.nights.find(x=>x.night===nk).rank===0).sort((a,b)=>key(b)-key(a));
  let html=overview([...prim,...back], nk);
  if(prim.length) html+=`<div class="section">${nk==='all'?'All targets':'Primary list'} · ${prim.length}</div>`+prim.map(t=>card(t,nk)).join('');
  if(back.length) html+=`<div class="section">Backups · ${back.length}</div>`+back.map(t=>card(t,nk)).join('');
  if(!html) html='<div class="empty">Nothing matches the current filters.</div>';
  main.innerHTML=html;
  main.querySelectorAll('.tgt').forEach(c=>c.addEventListener('click',()=>{ const el=document.getElementById('c-'+c.dataset.name); if(el){ el.scrollIntoView({behavior:'smooth',block:'start'}); el.style.outline='2px solid var(--accent)'; setTimeout(()=>el.style.outline='',1800);} }));
}
document.getElementById('foot').innerHTML=`Generated ${D.generated} from the CLAGN pipeline (Hemmati et al. 2026 W1 manifold; Zeltyn et al. 2024 SDSS-V CLAGNs; SDSS DR19 + DESI DR1 spectra; ZTF + unWISE light curves). Priority = B·S·(M+P) weighted by moon distance per night; nights are filled by priority per hour of telescope time using an exposure model scaled from the NGPS ETC (S/N ~7 on the diagnostic line under a bright moon; the "/h" figure on each card). Grey dots: Sample A; gold dots: Zeltyn CL-AGNs; shaded: enriched regions.`;
render();
</script>
'''

if __name__ == '__main__':
    main()
