"""Illustrative dust response, not a target fit or an ETC prediction.

A factor-six step fade motivates two limiting responses at A1's redshift.
Optically thin shells integrate the emitted spectra over the exact top-hat
light-travel response of each spherical shell. The fixed model retains the
pre-fade inner edge. The re-formed model adds dust down to the post-fade
sublimation radius, with the same radial density law and outer distribution.
Both spectra are normalized to their own mean 13-14 micron rest continuum.
No inferred heating history is constructed from an infrared light curve.
"""
from pathlib import Path
import json
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

HERE = Path(__file__).resolve().parent
snapshot = json.loads((HERE / 'inputs/figure_targets.json').read_text())
z = next(t['z'] for t in snapshot['targets'] if t['name'] == 'P11530')
fade, t_fade, t_obs = 1 / 6, 2019.5, 2027.5
elapsed = (t_obs - t_fade) / (1 + z)  # rest years
r_old = .25  # light-years; illustrative radius, not fitted
r_new, r_outer = r_old * np.sqrt(fade), 300 * r_old
beta, t_sub = 1.5, 1500.
lam = np.geomspace(1., 35., 1200)  # rest microns
lorentz = lambda center, width: width**2 / ((lam-center)**2 + width**2)
kappa = lam**-beta * (1 + 2.2*lorentz(9.7, 1.2) + .9*lorentz(18., 2.6))


def planck_nu(temperature):
    return lam**-3 / np.expm1(np.clip(14387.77 / (lam * temperature), 1e-6, 500))


def shells(rin, rout, density_index, formed=False, peak=False):
    radii = np.geomspace(rin, rout, 600)
    weight = radii**(2-density_index) * np.gradient(radii)
    result = np.zeros_like(lam)
    for radius, mass in zip(radii, weight):
        hot = t_sub * (radius/r_old)**(-2/(4+beta))
        cool = hot * fade**(1/(4+beta))
        low_fraction = np.clip(elapsed/(2*radius), 0, 1)
        if peak:
            emission = planck_nu(hot)
        elif formed:
            # Newly formed dust has no pre-fade emission in the empty cavity.
            emission = low_fraction * planck_nu(cool)
        else:
            emission = ((1-low_fraction)*planck_nu(hot)
                        + low_fraction*planck_nu(cool))
        result += mass * kappa * emission
    return result


def mean_band(flux, lo, hi):
    grid = np.linspace(lo, hi, 401)
    return np.trapezoid(np.interp(grid, lam, flux), grid)/(hi-lo)


def normalized(flux):
    return flux/mean_band(flux, 13., 14.)


cases, diagnostics = {}, []
for p in (1., 1.5):
    fixed = shells(r_old, r_outer, p)
    reformed = fixed + shells(r_new, r_old, p, formed=True)
    peak = shells(r_old, r_outer, p, peak=True)
    cases[p] = tuple(normalized(f) for f in (fixed, reformed, peak))
    hot_fixed = mean_band(fixed, 4.5, 5.5)/mean_band(fixed, 13., 14.)
    hot_reformed = mean_band(reformed, 4.5, 5.5)/mean_band(reformed, 13., 14.)
    diagnostics.append(dict(density_index=p, hot_colour_fixed=float(hot_fixed),
        hot_colour_reformed=float(hot_reformed),
        hot_colour_excess_percent=float(100*(hot_reformed/hot_fixed-1))))
summary = dict(status='Illustrative limiting cases; neither a fit to A1 nor an ETC simulation',
    redshift=z, luminosity_ratio=fade, fade_observed_year=t_fade,
    observing_year=t_obs, elapsed_rest_years=elapsed, old_inner_radius_light_years=r_old,
    new_inner_radius_light_years=float(r_new), emissivity_index=beta,
    sublimation_temperature_K=t_sub, fixed_inner_temperature_K=float(t_sub*fade**(1/(4+beta))),
    reference_rest_um=[13.,14.], hot_rest_um=[4.5,5.5], cases=diagnostics)
(HERE/'inputs/illustrative_model_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))

BLUE, ORANGE, GREY = '#244c78', '#d25a2f', '#73808b'
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','DejaVu Serif'],
    'font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,
    'legend.fontsize':6.5,'axes.spines.top':False,'axes.spines.right':False,
    'axes.linewidth':.6,'axes.edgecolor':'#9ba5ad','pdf.fonttype':42,'ps.fonttype':42})
fig, axes = plt.subplots(1,3,figsize=(6.5,2.25),
    gridspec_kw=dict(width_ratios=[.85,1.2,1.2],left=.065,right=.98,bottom=.21,top=.86,wspace=.48))
a=axes[0]
a.plot([2014,t_fade,t_fade,2029],[1,1,fade,fade],color=GREY,lw=1.4)
a.axvline(t_obs,color=ORANGE,ls='--',lw=.9)
a.set(xlim=(2014,2029),ylim=(0,1.12),xticks=[2015,2020,2025],yticks=[0,.5,1],
      xlabel='Year',ylabel=r'$L/L_{\rm initial}$')
a.set_title('(a) Illustrative fade',loc='left',fontsize=8,pad=5)
a.text(2028,.92,'MIRI',ha='right',fontsize=7,color=ORANGE)
lobs=lam*(1+z)
for a in axes[1:]:
    a.axvspan(.75,5,facecolor='#f3e7d3',zorder=0)
    a.axvspan(4.9,27.9,facecolor='#e9edf1',alpha=.7,zorder=0)
    a.set_xscale('log');a.set_xlim(1,30)
    a.xaxis.set_major_locator(FixedLocator([1,2,5,10,20]))
    a.xaxis.set_major_formatter(FuncFormatter(lambda v,_:f'{v:g}'))
    a.xaxis.set_minor_formatter(NullFormatter())
    a.set_xlabel('Observed wavelength (µm)',labelpad=2)
a=axes[1]
fixed,reformed,peak=cases[1.5]
a.plot(lobs,peak,color=GREY,ls=':',lw=.8,label='Before fade')
a.plot(lobs,fixed,color=BLUE,lw=1.2,label='Fixed dust')
a.plot(lobs,reformed,color=ORANGE,lw=1.2,label='Re-formed dust')
a.set_yscale('log');a.set_ylim(.004,4)
a.set_ylabel(r'$F_\nu / \langle F_\nu\rangle_{13-14\,\mu m}$',labelpad=2)
a.set_title('(b) Dust spectra',loc='left',fontsize=8,pad=5)
a.legend(loc='lower right',frameon=False,handlelength=1.3)
a.text(1.7,2.3,'SPHEREx',fontsize=6,ha='center',color=GREY)
a.text(11.5,2.3,'MIRI',fontsize=6,ha='center',color=GREY)
a=axes[2]
ratios=[v[1]/v[0] for v in cases.values()]
a.fill_between(lobs,np.minimum(*ratios),np.maximum(*ratios),color=ORANGE,alpha=.3)
a.plot(lobs,ratios[0],color=ORANGE,lw=1)
a.axhline(1,color=GREY,lw=.7)
a.set_ylim(.94,2.7);a.set_ylabel('Re-formed / fixed',labelpad=2)
a.set_title('(c) Spectral contrast',loc='left',fontsize=8,pad=5)
for rest,label in [(9.7,'9.7'),(18.,'18')]:
    a.axvline(rest*(1+z),color=GREY,lw=.5,alpha=.6)
    a.text(rest*(1+z),.94,label,transform=a.get_xaxis_transform(),
           fontsize=6.5,ha='center',va='top',color=GREY)
for a in axes:
    a.grid(axis='y',color='#e3e7eb',lw=.4);a.set_axisbelow(True)
for ext in ['pdf','png']:
    fig.savefig(HERE/f'fig2_model.{ext}',dpi=240)
plt.close(fig)
