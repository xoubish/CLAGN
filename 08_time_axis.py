"""
08_time_axis.py  --  Is the offset between the literature-CLAGN region and the Zeltyn region on the W1 manifold a
*time* axis?  All W1 light curves cover the same 2010-2020 window in absolute time, so the epoch of the largest
change in each curve can be compared with the object's position along the axis joining the two regions.

For every Sample A Turn-on/Turn-off object (parquet W1) and every projected Zeltyn CL-AGN / EVQ (unWISE cache):
  t_change = epoch that maximises |median(flux after) - median(flux before)|  (single change-point estimate)
  s        = projection of (umap_x, umap_y) onto the unit vector from the Turn-on/off-region centroid to the
             Zeltyn-region centroid (0 at the former, 1 at the latter)
Writes data/manifold_time_axis.csv and data/manifold_time_axis.png; prints correlations.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.time import Time

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, 'data')


def change_epoch(t, f, min_side=3):
    """MJD of the largest step between the medians before/after an epoch; also the step sign and relative size."""
    o = np.argsort(t); t, f = np.asarray(t)[o], np.asarray(f)[o]
    if len(t) < 2 * min_side + 1:
        return np.nan, np.nan, np.nan
    best = (0.0, np.nan, 0.0)
    for i in range(min_side, len(t) - min_side + 1):
        a, b = np.median(f[:i]), np.median(f[i:])
        step = (b - a) / max(np.median(f), 1e-9)
        if abs(step) > abs(best[0]):
            best = (step, 0.5 * (t[i - 1] + t[i]), abs(step))
    return best[1], np.sign(best[0]), best[2]


A = pd.read_csv(os.path.join(DATA, 'sampleA_embedding_objectid.csv'))
Z = pd.read_csv(os.path.join(DATA, 'zeltyn_embedding_full.csv'))
zc = pd.read_csv(os.path.join(DATA, 'zeltyn_coords.csv'))

# --- axis between the two regions
c_lit = A.loc[A.in_region_clagn, ['umap_x', 'umap_y']].mean().values
c_zel = A.loc[A.in_region_zeltyn, ['umap_x', 'umap_y']].mean().values
u = (c_zel - c_lit); L = np.linalg.norm(u); u = u / L
def s_of(x, y):
    return ((np.c_[x, y] - c_lit) @ u) / L
print(f'axis: literature-region centroid {c_lit.round(2)} -> Zeltyn-region centroid {c_zel.round(2)}, length {L:.2f} UMAP units')

# --- Sample A CLAGNs: W1 from the parquet
df = pd.read_parquet(os.path.join(DATA, 'df_lc_020724.parquet.gzip'))
w1 = df[df.index.get_level_values('band') == 'W1'].reset_index()[['objectid', 'time', 'flux']]
rows = []
cl = A[A.is_known_clagn]
for r in cl.itertuples():
    g = w1[w1.objectid == r.objectid]
    tc, sgn, size = change_epoch(g.time.values, g.flux.values)
    kind = 'Turn-on' if (r.label_bits & 16) else 'Turn-off'
    rows.append(dict(group='SampleA ' + kind, name=str(r.objectid), umap_x=r.umap_x, umap_y=r.umap_y, t_change=tc, sign=sgn, step=size))
# also all non-CLAGN Sample A (background) for the correlation over the whole manifold
bg = A[~A.is_known_clagn].sample(n=min(600, (~A.is_known_clagn).sum()), random_state=0)
for r in bg.itertuples():
    g = w1[w1.objectid == r.objectid]
    tc, sgn, size = change_epoch(g.time.values, g.flux.values)
    rows.append(dict(group='SampleA other', name=str(r.objectid), umap_x=r.umap_x, umap_y=r.umap_y, t_change=tc, sign=sgn, step=size))

# --- Zeltyn: W1 from the unWISE cache (objectid 1..N in zeltyn_coords order)
import glob
lc = pd.concat([pd.read_parquet(f) for f in glob.glob(os.path.join(DATA, 'wise_cache', 'zeltyn', '*.parquet'))]).reset_index()
lc = lc[lc.band == 'WISE_W1']
for r in Z[Z.projected].itertuples():
    oid = int(zc.index[zc.name == r.Name][0]) + 1
    g = lc[lc.objectid == oid]
    tc, sgn, size = change_epoch(g.time.values, g.flux.values)
    rows.append(dict(group='Zeltyn CL-AGN' if r.is_clagn else 'Zeltyn EVQ', name=r.Name, umap_x=r.umap_x, umap_y=r.umap_y, t_change=tc, sign=sgn, step=size))

R = pd.DataFrame(rows)
R['s'] = s_of(R.umap_x.values, R.umap_y.values)
R['year_change'] = 2000 + (R.t_change - 51544.5) / 365.25
R.to_csv(os.path.join(DATA, 'manifold_time_axis.csv'), index=False)

print('\nmedian change epoch (year) and axis coordinate s by group:')
print(R.groupby('group').agg(n=('s', 'size'), s_median=('s', 'median'), year_change_median=('year_change', 'median'),
                             step_median=('step', 'median')).round(2).to_string())
for grp, m in [('all CLAGN (SampleA on/off + Zeltyn CL-AGN)', R.group.isin(['SampleA Turn-on', 'SampleA Turn-off', 'Zeltyn CL-AGN'])),
               ('CLAGN + EVQ', R.group.str.contains('Turn|Zeltyn')), ('Sample A other (background)', R.group == 'SampleA other'), ('everything', R.s.notna())]:
    sub = R[m & R.t_change.notna()]
    rho, p = spearmanr(sub.s, sub.year_change)
    print(f'Spearman rho(s, year of largest W1 change) {grp}: {rho:+.2f} (p={p:.1e}, n={len(sub)})')
# large-step objects only (a real transition, not noise)
sub = R[R.group.str.contains('Turn|Zeltyn') & (R.step > 0.3)]
rho, p = spearmanr(sub.s, sub.year_change); print(f'  ... CLAGN/EVQ with |step| > 30%: rho {rho:+.2f} (p={p:.1e}, n={len(sub)})')

# --- figure
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
cols = {'SampleA other': '#9aa3b5', 'SampleA Turn-on': '#3987e5', 'SampleA Turn-off': '#1b6fc2', 'Zeltyn CL-AGN': '#d9a441', 'Zeltyn EVQ': '#e8c67a'}
for grp, g in R.groupby('group'):
    ax[0].scatter(g.s, g.year_change, s=14 if grp == 'SampleA other' else 26, c=cols[grp], alpha=.35 if grp == 'SampleA other' else .85, label=f'{grp} (n={len(g)})', edgecolor='none')
ax[0].axvline(0, color='#3987e5', lw=1, ls=':'); ax[0].axvline(1, color='#d9a441', lw=1, ls=':')
ax[0].set_xlabel('s: position along literature-region (0) → Zeltyn-region (1) axis'); ax[0].set_ylabel('year of largest W1 step'); ax[0].legend(fontsize=8, loc='lower right')
ax[0].set_title('Is the region offset a time axis?')
for grp, g in R.groupby('group'):
    if grp == 'SampleA other': continue
    ax[1].hist(g.year_change.dropna(), bins=np.arange(2010, 2021.5, 1), histtype='step', lw=1.6, color=cols[grp], label=grp, density=True)
ax[1].set_xlabel('year of largest W1 step'); ax[1].set_ylabel('density'); ax[1].legend(fontsize=8); ax[1].set_title('When did the W1 curves change?')
fig.tight_layout(); fig.savefig(os.path.join(DATA, 'manifold_time_axis.png'), dpi=130)
print('wrote data/manifold_time_axis.csv and .png')
