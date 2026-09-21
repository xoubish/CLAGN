"""Window manifold prototype: UMAP (DTW) on same-length NEOWISE W1 windows so an object's motion between windows is meaningful.
Windows (6.0 yr each): A 2014.0-2020.0, B 2016.5-2022.5, C 2018.5-2024.5.  GP (RBF 200 d, fixed, normalize_y) on a 60-point grid,
normalised by the window's clipped max (shape only), like the paper.  Fit on the calibration objects' windows, transform the rest.
Usage: window_manifold.py <tag> <visits csv> [<visits csv> ...]"""
import os, sys
os.environ.setdefault('NUMBA_THREADING_LAYER', 'workqueue')
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numba, umap
from sklearn.gaussian_process import GaussianProcessRegressor; from sklearn.gaussian_process.kernels import RBF
from scipy.spatial import cKDTree; from sklearn.metrics import roc_auc_score
from scipy import stats
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
SCR = '/private/tmp/claude-502/-Users-shemmati-Desktop-CLAGN/21a23a78-c7ac-423b-bad2-e7c1affbd39c/scratchpad'
tag = sys.argv[1]; files = sys.argv[2:]
ZP_MJY = 309.54e3; L = 2190.0; NG = 60; MINV = 8
WIN = {'A': 56658.0, 'B': 57570.0, 'C': 58300.0}           # window starts: 2014.0, 2016.5, 2018.5
v = pd.concat([pd.read_csv(f) for f in files]).drop_duplicates(['name', 'mjd']); v = v[np.isfinite(v.w1)]
v['flux'] = ZP_MJY * 10 ** (-0.4 * v.w1); v['err'] = 0.921 * v.flux * v.w1err.fillna(0.05).clip(lower=0.005)
grid = np.linspace(0, L, NG).reshape(-1, 1)
def rep(t, f, e):
    gp = GaussianProcessRegressor(kernel=RBF(length_scale=200.0), alpha=e ** 2, optimizer=None, normalize_y=True).fit(t.reshape(-1, 1), f)
    y = gp.predict(grid); cl, _, _ = stats.sigmaclip(y, 5, 5); mx = cl.max() if len(cl) else y.max()
    return y / mx if mx > 0 else None
rows = []; feats = []
for name, g in v.groupby('name'):
    for w, t0 in WIN.items():
        s = g[(g.mjd >= t0) & (g.mjd < t0 + L)]
        if len(s) < MINV or (s.mjd.max() - s.mjd.min()) < 0.8 * L:
            continue
        r = rep(s.mjd.values - t0, s.flux.values, s.err.values)
        if r is None or not np.all(np.isfinite(r)):
            continue
        rows.append((name, w, len(s), -2.5 * np.log10(s.flux.iloc[-3:].median() / s.flux.iloc[:3].median()))); feats.append(r)
W = pd.DataFrame(rows, columns=['name', 'win', 'nvis', 'fade_win']); X = np.array(feats)
print(f'{len(W)} windows for {W.name.nunique()} objects: ' + ', '.join(f'{w}: {n}' for w, n in W.win.value_counts().sort_index().items()))
@numba.njit(fastmath=True)
def dtw(a, b):
    n = a.shape[0]; m = b.shape[0]; E = np.empty((n, m))
    E[0, 0] = (a[0] - b[0]) ** 2
    for i in range(1, n): E[i, 0] = E[i - 1, 0] + (a[i] - b[0]) ** 2
    for j in range(1, m): E[0, j] = E[0, j - 1] + (a[0] - b[j]) ** 2
    for i in range(1, n):
        for j in range(1, m):
            v = (a[i] - b[j]) ** 2; v1 = E[i - 1, j]; v2 = E[i - 1, j - 1]; v3 = E[i, j - 1]
            E[i, j] = v + min(v1, min(v2, v3))
    return np.sqrt(E[n - 1, m - 1])
c = pd.read_csv('data/sdssv_calibration_pool.csv').set_index('name')
W['is_cal'] = W.name.isin(c.index)
fit_idx = np.where(W.is_cal.values)[0] if W.is_cal.sum() >= 500 else np.arange(len(W))
mapp = umap.UMAP(n_neighbors=50, min_dist=0.3, metric=dtw, random_state=3).fit(X[fit_idx])
emb = np.full((len(W), 2), np.nan); emb[fit_idx] = mapp.embedding_
rest = np.setdiff1d(np.arange(len(W)), fit_idx)
if len(rest): emb[rest] = mapp.transform(X[rest])
W['x'], W['y'] = emb[:, 0], emb[:, 1]
W.to_csv(f'{SCR}/window_manifold_{tag}.csv', index=False)
# per-object table: positions per window, displacement A->C
P = W.pivot(index='name', columns='win', values=['x', 'y', 'fade_win']); P.columns = [f'{a}_{b}' for a, b in P.columns]
P = P.join(c[['dimmed_any', 'bright07', 'dW1_full', 'M', 'z', 'dr', 'sv_mjd']], how='left')
P['dx'] = P.x_C - P.x_A; P['dy'] = P.y_C - P.y_A; P['disp'] = np.hypot(P.dx, P.dy)
neo = pd.concat([pd.read_csv(f.replace('visits', 'now')) for f in files]).drop_duplicates('name').set_index('name'); P = P.join(neo[['dw1_neowise']], how='left')
# within-SDSS-V recent change, if available
wsv = f'{SCR}/sdssv_within_change.csv'
if os.path.exists(wsv):
    s = pd.read_csv(wsv).set_index('name'); P = P.join(s[['dr_sv', 'span_yr', 'mjd_first', 'mjd_last']].rename(columns={'mjd_first': 'sv_first', 'mjd_last': 'sv_last'}), how='left')
cal = P[P.dimmed_any.notna() & P.x_A.notna() & P.x_C.notna()].copy(); cal['dimmed'] = cal.dimmed_any.astype(bool)
print(f'calibration objects with windows A and C: {len(cal)}, dimmed {cal.dimmed.sum()}')
# local dimming-rate fields (leave-one-out KNN, k=100) in each window's own positions
def loo_rate(XY, y, k=100):
    _, i = cKDTree(XY).query(XY, k=min(k + 1, len(XY))); return y[i[:, 1:]].mean(axis=1)
for w in 'ABC':
    ok = cal[f'x_{w}'].notna(); cal.loc[ok, f'rate_{w}'] = loo_rate(cal.loc[ok, [f'x_{w}', f'y_{w}']].values, cal.loc[ok, 'dimmed'].values.astype(float))
print('\nAUC for the DR16 -> SDSS-V dimming outcome:')
for nm, col in [('local rate, window A (2014-20)', 'rate_A'), ('local rate, window B (2016.5-22.5)', 'rate_B'), ('local rate, window C (2018.5-24.5)', 'rate_C'), ('rate change C - A (motion)', None),
                ('NEOWISE fade 2014-24 (scalar)', 'dw1_neowise'), ('unWISE fade 2010-20 (scalar)', 'dW1_full'), ('fade inside window C', 'fade_win_C'), ('old manifold M', 'M')]:
    s = cal.dropna(subset=[col] if col else ['rate_A', 'rate_C']); val = s[col] if col else (s.rate_C - s.rate_A)
    print(f'   {nm:36s} AUC {roc_auc_score(s.dimmed, val.fillna(0)):.3f}  (N={len(s)})')
if 'dr_sv' in P:
    r = P[P.dr_sv.notna() & P.x_C.notna()].copy(); r['rdim'] = r.dr_sv > 0.5
    print(f'\nwithin-SDSS-V outcome (r change between first and last SDSS-V epoch, >= 1.5 yr apart): N={len(r)}, dimmed > 0.5 mag: {r.rdim.sum()}')
    if r.rdim.sum() >= 8:
        for w in 'ABC':
            ok = r[f'x_{w}'].notna(); r.loc[ok, f'rrate_{w}'] = loo_rate(r.loc[ok, [f'x_{w}', f'y_{w}']].values, r.loc[ok, 'rdim'].values.astype(float))
        for nm, val in [('local recent-dimming rate, window C', r.rrate_C), ('window A', r.rrate_A), ('motion (rate C - rate A)', r.rrate_C - r.rrate_A), ('NEOWISE fade 2014-24', r.dw1_neowise), ('fade inside window C', r.fade_win_C), ('unWISE fade 2010-20', r.dW1_full), ('old manifold M', r.M)]:
            ok = val.notna(); print(f'   {nm:36s} AUC {roc_auc_score(r.rdim[ok], val[ok]):.3f}  (N={ok.sum()})')
for lab, s in [('dimmed', cal[cal.dimmed]), ('unchanged', cal[~cal.dimmed & ~cal.bright07.astype(bool)])]:
    print(f'   {lab:10s} N={len(s):4d} median displacement A->C {s.disp.median():.2f}, rate A {100*s.rate_A.mean():.1f}% -> C {100*s.rate_C.mean():.1f}%')
P.to_csv(f'{SCR}/window_manifold_{tag}_objects.csv')
# figure
fig, axs = plt.subplots(2, 3, figsize=(19, 11)); axs = axs.ravel()
allW = W.dropna(subset=['x'])
for k, w in enumerate('ABC'):
    ax = axs[k]; s = allW[allW.win == w]; ax.scatter(s.x, s.y, s=3, c='lightgrey')
    sc = ax.scatter(s.x, s.y, s=6, c=s.fade_win, cmap='RdBu_r', vmin=-0.4, vmax=0.4); plt.colorbar(sc, ax=ax, fraction=0.04, label='W1 change inside window (mag, + = faded)')
    cc = cal[cal[f'x_{w}'].notna()]; d = cc[cc.dimmed]; ax.scatter(d[f'x_{w}'], d[f'y_{w}'], s=26, c='red', ec='k', lw=0.4, label=f'dimmed at SDSS-V ({len(d)})')
    b = cc[cc.bright07.astype(bool)]; ax.scatter(b[f'x_{w}'], b[f'y_{w}'], s=26, c='blue', ec='k', lw=0.4, label=f'brightened ({len(b)})'); ax.legend(fontsize=8)
    ax.set_title(f'window {w}: {2014.0 + [0, 2.5, 4.5][k]:.1f}-{2020.0 + [0, 2.5, 4.5][k]:.1f}, {len(s)} objects')
ax = axs[3]; u = cal[~cal.dimmed]; ax.quiver(u.x_A, u.y_A, u.dx, u.dy, angles='xy', scale_units='xy', scale=1, color='grey', alpha=0.35, width=0.002)
d = cal[cal.dimmed]; ax.quiver(d.x_A, d.y_A, d.dx, d.dy, angles='xy', scale_units='xy', scale=1, color='red', width=0.004); ax.set_title('motion A -> C: dimmed (red) vs unchanged (grey)')
ax = axs[4]
if 'dr_sv' in P and P.dr_sv.notna().sum() > 20:
    r2 = P[P.dr_sv.notna() & P.x_A.notna() & P.x_C.notna()]; q = ax.quiver(r2.x_A, r2.y_A, r2.dx, r2.dy, r2.dr_sv.clip(-1, 1), angles='xy', scale_units='xy', scale=1, cmap='RdBu_r', clim=(-1, 1), width=0.003)
    plt.colorbar(q, ax=ax, fraction=0.04, label='r change between SDSS-V epochs (mag, + = fainter)'); ax.set_title(f'motion A -> C coloured by the change WITHIN SDSS-V ({len(r2)} objects)')
ax = axs[5]; cc = cal.dropna(subset=['rate_C'])
for w, col in [('A', 'grey'), ('C', 'k')]:
    ok = cc[f'rate_{w}'].notna(); b = pd.qcut(cc.loc[ok, f'rate_{w}'], 5, duplicates='drop'); g = cc[ok].groupby(b, observed=False).dimmed.agg(['mean', 'size'])
    ax.errorbar(range(len(g)), 100 * g['mean'], yerr=100 * np.sqrt(g['mean'] * (1 - g['mean']) / g['size']), fmt='o-', color=col, label=f'window {w}', capsize=3)
ax.set_xlabel('quintile of local dimming rate (leave-one-out)'); ax.set_ylabel('% dimmed at SDSS-V'); ax.legend(); ax.set_title('calibration: does position in each window predict dimming?'); ax.grid(alpha=0.3)
for ax in axs[:5]: ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
plt.suptitle(f'Window manifold ({tag}): UMAP+DTW on 6-yr NEOWISE W1 windows; motion = same object in successive windows'); plt.tight_layout()
plt.savefig(f'{SCR}/window_manifold_{tag}.png', dpi=100); print(f'wrote window_manifold_{tag}.png')
