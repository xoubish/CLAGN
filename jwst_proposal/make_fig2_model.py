"""Figure 2 (v0.3): schematic fixed-dust vs evolving-dust prediction for anchor A1 (P11530).

Toy model, stated as such in the caption: optically thin shells of dust in radiative
equilibrium with retarded illumination. Illumination history from the data:
ZTF g (host g = 18.8 mag assumed and subtracted) for 2018.5-2025.8, NEOWISE W1 minus a
1.0 mJy floor before that (scaled to join), constant after 2025.8.
Opacity: lambda^-1.5 continuum with Lorentzian 9.7 and 18 um silicate features.
Fixed dust: inner edge at r_sub(L_peak), unchanged. Evolving dust: inner edge follows
r_sub(L) after the fade (reformation complete by 2027), same density law.
Prints the diagnostic ratios used in the text."""
from pathlib import Path
import json
import numpy as np
from astropy.time import Time
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SRC = HERE / 'inputs/figure_targets.json'
data = json.loads(SRC.read_text(encoding='utf-8'))
t = {x['name']: x for x in data['targets']}['P11530']; z = t['z']
yr = lambda m: Time(np.asarray(m, float), format='mjd').decimalyear
BLUE, ORANGE, INK, INK2, GRID, AXIS, MUTED = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e', '#e1e0d9', '#c3c2b7', '#898781'
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Times New Roman', 'Times'], 'font.size': 8.5, 'axes.labelsize': 8.5,
                     'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 7.5, 'pdf.fonttype': 42,
                     'axes.edgecolor': AXIS, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .6})

# ---- illumination history L(t), relative to the 2019 peak ----
g = np.asarray(t['ztf']['g'], float); tg = yr(g[:, 0]); fg = 10**(-0.4*g[:, 1]) - 10**(-0.4*19.3)
o = np.argsort(tg); tg, fg = tg[o], np.maximum(fg[o], 0.02*np.nanmax(fg))
neo = np.asarray(t['neo'], float); tn = yr(neo[:, 0]); fn = np.maximum(neo[:, 1] - 1.0, 0.05)
wise = np.asarray(t['wise']['W1'], float); tw = yr(wise[:, 0]); fw = np.maximum(wise[:, 1] - 1.0, 0.05)
grid = np.arange(2009.0, 2028.6, 0.02)
Lz = np.interp(grid, tg, fg); scale = np.interp(2018.6, tg, fg) / np.interp(2018.6, tn, fn)
Lpre = np.interp(grid, np.r_[tw[tw < 2013.5], tn], np.r_[fw[tw < 2013.5], fn]) * scale
L = np.where(grid < tg.min(), Lpre, Lz); L[grid > tg.max()] = np.mean(fg[tg > tg.max()-1])
L = L / L.max()
Lpeak, tpeak = 1.0, grid[L.argmax()]

# ---- dust model ----
lam = np.geomspace(1.0, 35.0, 900)
lor = lambda l0, gam: gam**2/((lam-l0)**2+gam**2)
kappa = lam**-1.5 * (1 + 2.2*lor(9.7, 1.2) + 0.9*lor(18.0, 2.6))
beta, Tsub = 1.5, 1500.0
def T_of(Lr, r_over_rsub):  # radiative equilibrium: T^(4+beta) ∝ L/r^2, T=Tsub at r_sub(Lpeak) when L=Lpeak
    return Tsub * (Lr/Lpeak)**(1/(4+beta)) * r_over_rsub**(-2/(4+beta))
c_ltyr = 1.0  # radii in light-years; delays in years
rsub_peak = 0.25  # light-years at the 2019 peak: 0.4 pc (L/1e45)^0.5 for L ~ 4e43 erg/s
def Leff(r, tobs):  # thin-shell response: mean of L over the last 2r/c
    m = (grid <= tobs) & (grid >= tobs - 2*r/c_ltyr); return L[m].mean() if m.any() else L[grid <= tobs][-1]
def B_nu(T):  # in arbitrary units ∝ lam^-3 / (exp(hc/lam k T)-1)
    x = 14387.77/(lam*T); return lam**-3/np.expm1(np.clip(x, 1e-6, 500))
def spectrum(r_in, r_out, tobs, p=1.5):
    r = np.geomspace(r_in, r_out, 220); dM = r**(2-p) * np.gradient(r)
    F = np.zeros_like(lam)
    for rk, mk in zip(r, dM):
        F += mk * kappa * B_nu(T_of(Leff(rk, tobs), rk/rsub_peak))
    return F
tobs = 2027.8
r_out = 300*rsub_peak
Lnow = Leff(rsub_peak, tobs); r_in_new = rsub_peak*np.sqrt(Lnow/Lpeak)
cases = {}
for p in (1.0, 1.5):
    Ff = spectrum(rsub_peak, r_out, tobs, p=p); Fe = spectrum(r_in_new, r_out, tobs, p=p); Fp = spectrum(rsub_peak, r_out, tpeak+0.3, p=p)
    cases[p] = (Ff, Fe, Fp)
lobs = lam*(1+z)
norm = 3.5 / np.interp(12.0, lobs, cases[1.0][0])  # p=1 fixed model = 3.5 mJy at observed 12 um (0.7 x AllWISE W3)
F_fixed, F_evol, F_peak = [x*norm for x in cases[1.0]]
ratio = {p: cases[p][1]/cases[p][0] for p in cases}
# statistical S/N per R=100 bin from the JDox MRS sensitivity curve (S/N 10 in 10 ks per pixel), 60 groups x 4 dithers per setting
slam = np.array([5,7,10,12,14,15.5,17,18,19,20,21,22,23,24,25,26,27,27.5]); slim = np.array([.08,.07,.09,.10,.13,.17,.22,.30,.5,.8,1,1.5,2,2.7,4,6,15,30])
Flim = np.exp(np.interp(lobs, slam, np.log(slim))); Rnat = np.interp(lobs, [5,10,15,20,25,28], [3500,3000,2500,2000,1700,1500])
texp = 4*60*2.775
sn_pix = np.where(F_fixed > Flim, 10*np.sqrt(F_fixed/Flim), 10*F_fixed/Flim) * np.sqrt(texp/1e4)
sn = sn_pix*np.sqrt(2*Rnat/100); frac = 1/np.maximum(sn, 1e-3)
mrs = (lobs >= 4.9) & (lobs <= 27.9)
print('L_now/L_peak =', round(Lnow, 3), ' r_in,new/r_in,peak =', round(r_in_new/rsub_peak, 2), ' T_in fixed now =', round(T_of(Lnow, 1.0)), 'K')
for l in [2, 3, 4, 5, 6, 8, 10.8, 12, 15, 20, 24]:
    i = np.argmin(abs(lobs-l)); print(f'  obs {l:5.1f} um: re-formed/fixed p=1.0 {ratio[1.0][i]:5.2f}  p=1.5 {ratio[1.5][i]:5.2f}   MRS 1-sigma stat = {100*frac[i]:4.1f}%')
print('fade at obs 12 um (p=1, 2027.8 / 2019.8):', round(np.interp(12, lobs, F_fixed)/np.interp(12, lobs, F_peak), 2))

fig = plt.figure(figsize=(6.5, 2.55), facecolor='white')
gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.3, 1.1], wspace=.45, left=.07, right=.985, top=.86, bottom=.19)
ax = fig.add_subplot(gs[0]); ax.plot(grid, L, color=INK, lw=1.2)
ax.axvspan(2027.5, 2028.5, color=GRID, alpha=.7, lw=0); ax.axvline(tobs, color=ORANGE, lw=.8, ls='--')
ax.set_xlim(2009.5, 2028.7); ax.set_ylim(0, 1.1); ax.set_xticks([2010, 2016, 2022, 2028]); ax.set_yticks([0, .5, 1])
ax.set_xlabel('Year', labelpad=1); ax.set_ylabel(r'$L / L_{\rm peak}$'); ax.set_title('(a) Illumination history of A1', loc='left', fontsize=8.5, pad=3)
ax.text(2026.4, .93, 'JWST', color=ORANGE, fontsize=7.5, ha='right'); ax.grid(axis='y', color=GRID, lw=.5); ax.set_axisbelow(True)
ax = fig.add_subplot(gs[1])
ax.axvspan(0.75, 5.0, color='#f3e7d3', lw=0); ax.axvspan(4.9, 27.9, color=GRID, alpha=.6, lw=0)
ax.plot(lobs, F_peak, color=MUTED, lw=.8, ls=':', label='same dust near the 2019 peak')
ax.plot(lobs, F_fixed, color=BLUE, lw=1.4, label='fixed dust, 2027.8')
ax.plot(lobs, F_evol, color=ORANGE, lw=1.4, label='dust re-formed inward, 2027.8')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1, 32); ax.set_ylim(0.1, 40)
ax.set_xticks([1, 2, 5, 10, 20]); ax.set_xticklabels(['1', '2', '5', '10', '20'])
ax.set_xlabel('Observed wavelength (μm)', labelpad=1); ax.set_ylabel(r'$F_\nu$ (mJy)'); ax.set_title('(b) Predicted dust spectrum', loc='left', fontsize=8.5, pad=3)
ax.legend(frameon=False, loc='lower right', handlelength=1.6, labelcolor=INK2, fontsize=7)
ax.text(1.9, 28, 'SPHEREx', color=INK2, fontsize=7, ha='center'); ax.text(11.5, 28, 'MIRI/MRS', color=INK2, fontsize=7, ha='center')
ax = fig.add_subplot(gs[2])
ax.axvspan(0.75, 5.0, color='#f3e7d3', lw=0); ax.axvspan(4.9, 27.9, color=GRID, alpha=.6, lw=0)
lo = np.minimum(ratio[1.0], ratio[1.5]); hi = np.maximum(ratio[1.0], ratio[1.5])
ax.fill_between(lobs, lo, hi, color=ORANGE, alpha=.35, lw=0, label='re-formed / fixed')
ax.plot(lobs, ratio[1.0], color=ORANGE, lw=1.0)
ax.fill_between(lobs[mrs], 1-frac[mrs], 1+frac[mrs], color=BLUE, alpha=.3, lw=0, label='MRS ±1σ stat., 33 min')
ax.axhline(1, color=AXIS, lw=.6); ax.axhline(1.056, color=MUTED, lw=.6, ls='--'); ax.text(26, 1.075, '5.6% abs. cal.', fontsize=6.5, color=MUTED, ha='right')
ax.set_xscale('log'); ax.set_xlim(1, 32); ax.set_xticks([1, 2, 5, 10, 20]); ax.set_xticklabels(['1', '2', '5', '10', '20'])
ax.set_ylim(0.9, 2.6); ax.set_xlabel('Observed wavelength (μm)', labelpad=1); ax.set_ylabel('flux ratio')
ax.set_title('(c) Discriminant and precision', loc='left', fontsize=8.5, pad=3); ax.legend(frameon=True, framealpha=.9, edgecolor='none', loc='center right', handlelength=1.6, labelcolor=INK2, fontsize=6.8)
for l0, lab in [(9.7*(1+z), '9.7'), (18*(1+z), '18')]:
    ax.axvline(l0, color=AXIS, lw=.5, zorder=0); ax.text(l0, .03, lab, transform=ax.get_xaxis_transform(), ha='center', fontsize=7, color=INK2)
fig.savefig(HERE/'fig2_model.pdf'); fig.savefig(HERE/'fig2_model.png', dpi=220); print('ok')
