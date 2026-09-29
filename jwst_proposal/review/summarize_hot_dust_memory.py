"""Render the completed conditional memory tests and optical extensions."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/clagn-matplotlib')
import json,html
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.time import Time

HERE=Path(__file__).resolve().parent
OUT=HERE/'hot_dust_memory_test'


def main():
    results=[]
    for sub,label in [(OUT,'catalogue + extrapolation'),(OUT/'with_saved_alerts','saved alerts + extrapolation')]:
        summary=json.loads((sub/'summary.json').read_text())
        for key,rec in summary['reports'].items():
            if 'fractional_rmse_reduction' not in rec:continue
            band,host,continuation,fraction,family,grid,cohort,gain=key.split('|')
            results.append(dict(input=label,band=band,template=host,continuation=continuation,optical_host_fraction=float(fraction),
                family=family,grid=grid,cohort=cohort,gain=gain,n=rec['n'],
                current_rmse=rec['rmse_instantaneous'],delayed_rmse=rec['rmse_nested_heldout'],
                improvement_percent=100*rec['fractional_rmse_reduction'],best_delay_parameter=rec['best_grid_delay_parameter']))
    table=pd.DataFrame(results);table.to_csv(OUT/'all_comparisons.csv',index=False)
    current=json.loads((OUT/'with_saved_alerts/summary.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(12,10),constrained_layout=True)
    for row,(family,grid,gain,label) in enumerate([('fixed_years','0to1yr','unit_gain','One common delay'),
                                                ('luminosity_scaled','scaled','free_gain','Luminosity-scaled delays + fitted response')]):
        for col,band in enumerate(['g','r']):
            ax=axes[row,col];key=f'{band}|Ell5|linear|0.0|{family}|{grid}|screened|{gain}'
            rec=current['reports'][key]
            values=rec['per_target']
            for r in values:
                ax.plot([r['instantaneous'],r['heldout_prediction']],[r['observed']]*2,color='#bbb',lw=.7)
                ax.scatter(r['instantaneous'],r['observed'],marker='x',color='#888',s=25)
                ax.errorbar(r['heldout_prediction'],r['observed'],xerr=r['optical_error'],yerr=r['dust_error'],
                            fmt='o',color='#b66d35',ecolor='#dbc0a8',ms=3.5,lw=.6)
                ax.annotate(r['id'],(r['heldout_prediction'],r['observed']),xytext=(3,3),textcoords='offset points',fontsize=6)
            lim=max(.2,max(abs(r[k]) for r in values for k in ['observed','instantaneous','heldout_prediction']))*1.2
            ax.plot([-lim,lim],[-lim,lim],':',color='#999');ax.set(xlim=(-lim,lim),ylim=(-lim,lim),
                xlabel='Predicted ln(last / first dust flux)',ylabel='Fitted ln(last / first dust flux)',
                title=f'{label}\n{band}: n={rec["n"]}; prediction error {abs(100*rec["fractional_rmse_reduction"]):.1f}% '+
                      ('lower' if rec['fractional_rmse_reduction']>0 else 'higher'))
    fig.suptitle('Conditional echo test · linear optical extensions where needed\nGrey ×: current-state model · orange: delayed model · provisional decomposition errors',fontsize=12)
    fig.savefig(OUT/'model_comparison.png',dpi=150);fig.savefig(OUT/'model_comparison.pdf');plt.close(fig)
    # All tested optical extensions are reviewable, including sparse targets.
    sample=pd.read_csv(HERE.parent/'inputs/jwst_sample_cycle6.csv')
    sub=OUT/'with_saved_alerts';slopes=pd.read_csv(sub/'optical_extrapolations.csv')
    spectra=pd.read_csv(sub/'spectral_fits.csv')
    for page in range(4):
        fig,axes=plt.subplots(3,2,figsize=(12,10),constrained_layout=True)
        for ax,target in zip(axes.flat,list(sample.itertuples())[6*page:6*(page+1)]):
            for band,color in [('g','#267967'),('r','#a64f5e')]:
                rec=slopes[slopes.id.eq(target.id)&slopes.band.eq(band)]
                if rec.empty:continue
                r=rec.iloc[0];d=pd.read_csv(sub/f'{target.id}_{band}_optical_bins.csv')
                norm=d.flux_mjy.iloc[-1]
                ax.errorbar(Time(d.mjd,format='mjd').decimalyear,d.flux_mjy/norm,yerr=d.error_mjy/norm,
                            fmt='.',color=color,ms=2,lw=.4,alpha=.65,label=band)
                stop=max(61500,r.last_mjd+30);time=np.linspace(r.last_mjd,stop,80);dt=(time-r.last_mjd)/365.25
                extrap=1+r.slope_mjy_per_observed_year*dt/norm
                err=np.sqrt(d.error_mjy.iloc[-1]**2+(r.slope_error*dt)**2)/norm
                ax.plot(Time(time,format='mjd').decimalyear,extrap,'--',color=color,lw=1)
                ax.fill_between(Time(time,format='mjd').decimalyear,extrap-err,extrap+err,color=color,alpha=.10)
            for t in spectra[spectra.id.eq(target.id)&spectra.template.eq('Ell5')].mjd:
                ax.axvline(Time(t,format='mjd').decimalyear,color='#aaa',lw=.7)
            ax.axhline(0,color='grey',lw=.5)
            ax.set(xlim=(2023,2027),title=target.id,xlabel='Observed year',ylabel='Flux / last measured bin')
            ax.legend(fontsize=7,loc='best')
        fig.suptitle(f'Optical history and estimated continuation · {page+1}/4\nDashed: last-year flux slope; shading: approximate trend error; grey: SPHEREx dates',fontsize=12)
        fig.savefig(OUT/f'optical_extensions_{page+1}.png',dpi=130)
        fig.savefig(OUT/f'optical_extensions_{page+1}.pdf');plt.close(fig)
    notes='''# Conditional hot-dust memory test

The requested last-year linear ZTF extrapolation has been tested. Some delayed
models improve held-out prediction, but the result is not robust enough to
claim a detected hot-dust memory signal or a measured response time.

The initial catalogue-only pass and a separate pass incorporating available
saved alert photometry are both retained. The latter adds measured points for
F02, F04, F06, F07 and R01, after matching their magnitude zero points to the
catalogue light curves. Matching offsets require at least five near-simultaneous
points. This avoids extrapolating over saved measurements; no new alert query
was made. Alert-reference calibration uncertainties are not fully modeled.

## What was tested

- The final 365.25 observed days of each optical series were fit in **flux**,
  using 14-day bins. At least three bins spanning 90 days were required. The
  slope is extended continuously from the final measured bin. Flat continuation
  is a sensitivity test. Negative predicted illumination is excluded, not clipped.
- Spectra were fit jointly per target with a constant host, variable disc, and
  one blackbody per sufficiently covered visit. Ell5 and S0 host templates were
  tested. The thermal quantity is a model-dependent bolometric blackbody proxy.
- Observed changes between first and last adequately covered SPHEREx visits
  were compared with changes in the instantaneous or delayed optical driver.
  A constant source-specific reprocessing normalization cancels in these ratios.
- The response is a uniform delay distribution from zero to twice its mean.
  Fixed mean lags of 0, 0.1, 0.25, 0.5, 1, 2 and 3 rest years were explored.
  The 0–1 year and 0–3 year ranges are reported separately because requiring
  prehistory coverage for the longer grid removes additional targets.
- A second family uses luminosity-scaled delays. The reference follows
  log10(tau/day) = -2.11 - 0.2 M_V (Koshida et al. 2014), with approximate V-band
  luminosities inferred from the median optical flux, a fixed Fnu ∝ nu^(1/3)
  continuum, and Planck18 distances. AB/Vega differences, extinction and host
  luminosity make this an approximate prior, not a measured lag. Factors
  0, 0.5, 1, 2 and 4 times that reference were tested.
- Both unit response and one globally fitted response exponent in [0,3] were
  evaluated. Each target was left out when selecting the lag and, where free,
  fitting that exponent. The current-state comparison has the same freedom.
  Targets within any one comparison are identical across all tested lags.
- Optical host sensitivity subtracts either nothing or 50% of the minimum
  observed optical flux. The latter is an assumption, not a measured host.

## Interpretation

With saved alerts included, the simple common-lag/unit-response model worsens
held-out prediction by about 12% in g and 6% in r for the baseline host template.
The luminosity-scaled/fitted-response family improves it by about 7% in g and
18% in r. These families use different common-coverage cohorts, so compare
current and delayed models **within** each row, not raw RMSE between rows.

The improvements depend on optical host subtraction, response freedom,
continuation and host template. Some gains reverse under the host-subtraction
sensitivity case. The trial families were explored in this session, including
the shorter lag range after the initial broad-grid pass; this is exploratory,
not a preregistered significance test. All comparisons are in
`all_comparisons.csv`, including less favorable results.

The main limitation is spectral model adequacy: **none of these simple joint
spectral fits passes a p > 0.01 goodness-of-fit check with the supplied independent
errors**. Reduced chi-squared ranges from about 1.6 to 62. Scaling local parameter
covariances by reduced chi-squared does not fix those physical/model residuals.
Thus even the more promising predictive improvements are conditional on
unvalidated dust decompositions. The `screened` cohort only removes boundary,
rank, extreme-uncertainty and predominantly nonphysical-driver cases; it is
not a statistically adequate-fit cohort. The `adequate` cohort is empty.

Errors propagate local spectral covariance (including shared host), an assumed
3% independent scale error per visit, and Monte Carlo optical/trend errors.
The model-selection score uses conditional mean predictions, not a fully
marginalized optical-history likelihood. It is equal-target log-change RMSE.
Bootstrap intervals resample targets only. They are not detection significances.

Saved spectra/optical points, exclusions, slope errors, extrapolation horizons,
all parameter-grid predictions and full comparison reports are retained.
SPHEREx visits are represented by median times, and unmodeled intra-visit
evolution, extinction, complex dust SEDs and extraction/calibration effects remain.
The experiment tests an ordinary echo. It does not measure a dust formation or
destruction timescale and does not rule out structural memory.

## Reproduce

From the project root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/anaconda3/bin/python jwst_proposal/review/test_hot_dust_memory.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/anaconda3/bin/python jwst_proposal/review/test_hot_dust_memory.py --use-saved-alerts
/opt/anaconda3/bin/python jwst_proposal/review/summarize_hot_dust_memory.py
```

No proposal source or proposal PDF is modified.

Reference: https://arxiv.org/abs/1406.2078
'''
    (OUT/'README.md').write_text(notes)
    chosen=table[(table.input=='saved alerts + extrapolation')&(table.continuation=='linear')&
        (table.template=='Ell5')&(table.optical_host_fraction==0)&(table.cohort=='screened')&
        (((table.family=='fixed_years')&(table.grid=='0to1yr')&(table.gain=='unit_gain'))|
         ((table.family=='luminosity_scaled')&(table.gain=='free_gain')))]
    parts=['<!doctype html><meta charset="utf-8"><title>Hot-dust memory test</title>',
           '<style>body{font:16px system-ui;max-width:1200px;margin:40px auto;color:#26333c;padding:0 20px}img{width:100%}td,th{padding:7px;border-bottom:1px solid #ddd}pre{white-space:pre-wrap}</style>',
           '<h1>Conditional hot-dust memory test</h1>',
           '<p>Some delayed models improve prediction; the result depends on assumptions and is not a validated detection. None of the simple spectral decompositions passes the supplied-error adequacy test.</p>',
           chosen[['band','family','n','current_rmse','delayed_rmse','improvement_percent']].round(3).to_html(index=False),
           '<p>Positive improvement means lower error. Each row compares identical targets; different rows can have different coverage selections.</p>',
           '<img src="model_comparison.png" alt="Current versus delayed predictions">',
           '<h2>Review the optical extrapolations</h2>']
    for i in range(1,5):parts.append(f'<img src="optical_extensions_{i}.png" alt="Optical extensions page {i}">')
    parts+=['<h2>Full methods and limitations</h2><pre>'+html.escape(notes)+'</pre>',
            '<a href="all_comparisons.csv">All comparisons CSV</a>']
    (OUT/'index.html').write_text('\n'.join(parts))
    print(chosen[['band','family','n','improvement_percent']].to_string(index=False))


if __name__=='__main__':main()
