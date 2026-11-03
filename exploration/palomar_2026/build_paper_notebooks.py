"""Build only the two Palomar paper pilot notebooks (06 and 07)."""
from pathlib import Path
import nbformat as nbf

HERE=Path(__file__).resolve().parent


def md(s):return nbf.v4.new_markdown_cell(s.strip())
def code(s):return nbf.v4.new_code_cell(s.strip())


SETUP=code('''
from pathlib import Path
import sys, json, inspect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, Image, Markdown
HERE = Path.cwd().parent if Path.cwd().name == 'notebooks' else Path.cwd()
sys.path.insert(0, str(HERE))
from paper_protocol import OUT
pd.set_option('display.max_colwidth', 80)
sample = pd.read_csv(OUT/'analysis_sample.csv')
pilot = pd.read_csv(OUT/'pilot_sample.csv')
epochs = pd.read_csv(OUT/'pilot_epochs.csv')
''')


def save(name,cells):
    book=nbf.v4.new_notebook(cells=cells,metadata=dict(kernelspec=dict(display_name='Python 3',language='python',name='python3')))
    nbf.write(book,HERE/'notebooks'/name)
    print(name)


def main():
    save('06_palomar_paper_protocol_and_pilot.ipynb',[
        md('''# Palomar paper: sample, timing, and pilot

**Question:** does the pre-run light-curve manifold identify broad-line evolution beyond simple variability information?
The likely first paper is a systematic Palomar spectral-evolution and selection-validation study, supported by a few well-validated individual cases. A successful forecast is a hypothesis to test, not the present conclusion. Dust physics is deferred to JWST.

This protocol was written **after** the initial spectral inspection. It is not a preregistration. Eleven science cases were chosen using that first look; four comparison sources were matched using pre-run covariates. Consequently, the pilot cannot estimate selection yield, significance, or predictive performance.

[Pilot fits and results → notebook 07](07_palomar_paper_fits_and_results.ipynb) · [Paper plan](../PALOMAR_PAPER.md)'''),SETUP,
        md('''## What is fixed for this pilot

Keep every Palomar visit for the selected sources, the earliest accepted archive, and the latest accepted archive preceding each visit. The reference is chosen by date and existing quality criteria, not by the largest line change. A reference lacking usable native errors stays unavailable; an older spectrum is not silently substituted.

The primary fitting experiment concerns broad Hβ. Total-continuum EW and broad-Hβ/[O III] are complementary diagnostics: host dilution affects the former; extended narrow emission and aperture differences affect the latter. Neither alone is an intrinsic luminosity measurement.'''),
        code('''
print(f'{len(sample)} campaign sources; {len(pilot)} pilot sources; {len(epochs)} selected epochs')
display(pilot[['name','pilot_role','experiment_group','known_state_status','z','r_pre_run','comparison_baseline_years']])
display(epochs.groupby(['source','night'], dropna=False).size().rename('epochs').to_frame())
'''),
        md('''## Comparison matching

Candidate comparisons come from pre-run group D and exclude catalog-confirmed CLAGN and the eleven science cases. Optimize distinct matches using redshift, historical optical brightness and archival baseline. Calipers: |Δz| ≤ 0.15, |Δr| ≤ 1.5 mag, |log₁₀ baseline ratio| ≤ 1, with a 0.1-year offset in baseline ratios. Retain four low-cost matches.

These comparisons are not assumed nonvariable. P15270 is a relatively weak match to P9694 in baseline; the matching table makes that visible. Matching a deliberately selected pilot does not remove its outcome-selection bias.'''),
        code("display(pd.read_csv(OUT/'pilot_matches.csv'))"),
        md('''## Is the spectral endpoint genuinely prospective?

For original DR16-parent sources, the original W1 pool cache links `objectid` to `poolid` (see `scripts/02_parent_pool.py`). Recover the last valid input measurement and compare it with the accepted archival spectrum. Do **not** equate pool IDs with Sample A training IDs; they are different namespaces.

An archival spectrum after all linked W1 measurements is a timing candidate for a forecast test. The historical projection snapshot and training-identity overlap still need auditing. If the archival baseline predates some input photometry, part of the eventual spectral change may already have been visible to the selection method. Unresolved timing is retained explicitly.

The on/off neighbor fractions describe labeled training neighborhoods. They are not calibrated probabilities or proven directional forecasts.'''),
        code('''
display(sample.timing_interpretation.value_counts().rename('sources').to_frame())
display(pilot[['name','manifold_direction','predicted_direction','latest_archive_mjd','w1_input_last_mjd','timing_interpretation']])
fig,ax=plt.subplots(figsize=(9,3.5))
g=sample[sample.origin.eq('old_DR16_parent') & sample.w1_input_last_mjd.notna()].copy()
delta=(g.latest_archive_mjd-g.w1_input_last_mjd)/365.25
ax.hist(delta.dropna(),bins=16,color='#346b73');ax.axvline(0,color='#b53c4a',ls='--')
ax.set(xlabel='Latest archival spectrum − last linked W1 input [observed years]',ylabel='Sources',title='Negative values allow retrospective information; positive values still require provenance checks')
plt.show()
'''),
        md('''## Native data and instrumental assumptions

Three of four selected DESI epochs were recovered from cached native files by exact time interval and a coordinate match within 2 arcsec. The remaining P10381 DESI reference has no uniquely identified native file and is not fit. SDSS supplies native inverse variance and dispersion. NGPS uses delivered diagonal errors and an assumed R=1650. DESI resolution matrices are approximated by their second moments; the final analysis should forward-convolve through the actual instrumental response.

This pilot fits 42 of 43 epochs with three host/Fe II variants, then refits all 20 NGPS visits under the alternative standard with all three variants. Alternative standards are **stress tests**: notably, October 27 Feige 34 was observed through cloud and is not equally credible with the adopted BD+28 4211 calibration.'''),
        code('''
audit=pd.read_csv(OUT/'native_data_audit.csv')
display(audit.groupby(['source','status']).size().rename('epochs').to_frame())
display(audit[['epoch_id','error_kind','resolution_status','status']])
'''),
        md('''## Inspect and reproduce the selection code

The snapshot records hashes of the actual input tables and W1 cache files. The builder refuses to replace it when those inputs change. Any revised input set should get a versioned snapshot. This preserves a post-first-look audit trail; it does not prove that today's files are byte-identical to their pre-observation versions.'''),
        code('''
import paper_protocol
print(inspect.getsource(paper_protocol.build))
snapshot=json.loads((OUT/'protocol_snapshot.json').read_text())
display(pd.Series({k:v for k,v in snapshot.items() if k not in ['inputs','restrictions','science_cases','covariates']}))
print('Hashed input files:',len(snapshot['inputs']))
'''),
        md('''## Decision gates before a population claim

1. Finish the historical feature and training-overlap audit; separate known-source follow-up from new discoveries.
2. Validate broad-line changes with Hα where covered, shared host constraints, exposure checks and supported calibration choices.
3. Extend the frozen outcome rules to the full Hβ-accessible observed sample, preserving failures and upper limits.
4. Measure per-source detectability with injections into actual spectra. The present toy simulations are not that completeness experiment.
5. Compare manifold information with simple pre-run variability at the source level; account for time baseline, brightness and detectability. The tiny C group and small comparison pool limit model complexity.

Until these gates pass, describe a targeted spectral-evolution sample. Do not claim an AGN incidence rate, a measured enrichment factor, or successful prospective prediction.''')
    ])
    save('07_palomar_paper_fits_and_results.ipynb',[
        md('''# Palomar paper pilot: decompositions, stress tests, and priorities

**P21422 is the cleanest current weakening lead. P7947 is a strengthening lead with greater host-model sensitivity. P9506 deserves promotion as a repeat-confirmed long-baseline weakening case.** These are pilot conclusions about line variability, not new confirmed changing-look classifications.

The first-look wing indices and the fitted full broad-Hβ measurements below are different quantities. In particular, the P24930 and P17881 interpretations change after decomposition. Dust modeling remains deferred.

[Protocol and sample → notebook 06](06_palomar_paper_protocol_and_pilot.ipynb) · [All 15 fit pages, PDF](../outputs/palomar_paper/pilot_fit_atlas.pdf) · [Paper plan](../PALOMAR_PAPER.md)'''),SETUP,
        code('''
fits=pd.read_csv(OUT/'fit_results.csv')
pairs=pd.read_csv(OUT/'fit_comparisons.csv')
envelope=pd.read_csv(OUT/'model_sensitivity_envelope.csv')
stress=pd.read_csv(OUT/'calibration_model_envelope.csv')
display(fits.groupby(['status','variant']).size().rename('fits').to_frame())
'''),
        md('''## The actual fitted model

Fit rest 4200–5500 Å independently at each epoch with nonnegative coefficients for an AGN power law, two BC03 stellar templates, optical Fe II, narrow and broad Balmer components, [O III], and broad He II. Broad Hβ is a two-Gaussian sum. Narrow-line kinematics are shared within each epoch. Fix [O III] λ5008/λ4960 integrated flux to 2.98 and include a fixed-shape, free-amplitude blue wing. Pixel-integrated Gaussian line profiles include approximate instrumental broadening.

Nonlinear parameters vary continuum slope and narrow/broad kinematics; NNLS solves component amplitudes at each step. Two starting points reduce sensitivity to local minima. Three model families vary the young stellar component or Fe II smoothing. Fixed host dispersion and limited Fe II shapes are assumptions, not measured physical parameters.

Masks retain finite positive native errors and exclude selected sky/telluric intervals. Report coverage explicitly. Reduced χ² assumes diagonal errors and is only a residual diagnostic. **No statistical detection significance, posterior, black-hole mass, or disappearance label is inferred.** Galactic reddening and calibration/aperture covariance are not yet in a joint model; continuum slopes/host fractions remain conditional nuisance estimates.

Templates come from the pinned [PyQSOFit repository](https://github.com/legolason/PyQSOFit), with checksums and its license retained locally. This small local fitter does not run PyQSOFit itself. [Ren et al. (2024)](https://arxiv.org/abs/2406.17598) motivates the need to test host decomposition rather than assume a unique solution.'''),
        code('''
from pilot_fits import Basis, fit_arrays, one, VARIANTS
print(inspect.getsource(Basis))
print(inspect.getsource(fit_arrays))
display(VARIANTS)
'''),
        md('''## Refit an epoch here

Change `SOURCE` and `VARIANT` to inspect another object or model. This cell loads read-only inputs and recomputes one fit without replacing saved fit products.'''),
        code('''
from pilot_data import load
spectra=load()
SOURCE='P21422'; VARIANT='baseline'
chosen=[s for s in spectra if s['name']==SOURCE and s['source']=='NGPS'][-1]
display(pd.Series(one(chosen,VARIANT,save=False)))
'''),
        md('''## P21422: strongest current lead

The broad line stays present. Across the primary-standard host/Fe variants, its EW is 0.66–0.69 of the 2022 DESI reference and broad Hβ/[O III] is 0.56–0.61. Including the alternative-standard stress tests expands the ranges to roughly 0.59–0.69 and 0.52–0.61. These ranges are **not confidence intervals**.

The two independently extracted Palomar exposures agree in wing EW to about 2.4%, supporting a repeatable feature in the delivered spectrum. This does not eliminate aperture or common calibration effects. A local literature association points to Graham 2019, but its confirmation status remains unresolved; this is not yet a new-source discovery claim. The older [Graham et al. study](https://doi.org/10.1093/mnras/stz3244) already studies changing-state quasars, so novelty must be established from the source's history and new behavior.'''),
        code("display(Image(filename=str(OUT/'figures/P21422_decomposition.png')))"),
        md('''## Contrasting strengthening case and repeat-confirmed fading case

P7947's full broad-Hβ EW rises by 1.70–2.11× under the adopted calibration, and 1.51–2.11× across the tested standards. Its host decomposition and broad centroid are less secure; the strict [O III] window support fails because of the conservative sky mask. Treat the [O III]-normalized change as a model-interpolation diagnostic pending reinspection.

P9506 shows long-baseline weakening in both September and October. October/reference EW spans 0.67–0.77 across the tested models/standards, and broad Hβ/[O III] spans 0.53–0.55. Individual-exposure wing indices are similar across both nights. This supports persistent weakening; it does not establish an abrupt month-scale transition.'''),
        code('''
for name in ['P7947','P9506']:
    display(Image(filename=str(OUT/'figures'/f'{name}_decomposition.png')))
'''),
        md('''## Full pilot, including checks that weaken the first interpretation

P24930's fitted full-line EW ratio is about 0.87–0.88 under the adopted standard, considerably less dramatic than its wing-index change; missing wing coverage also matters. P17881 has EW growth but nearly unchanged broad-Hβ/[O III] relative to February 2025, so continuum dilution/calibration must be separated before claiming new broad-flux growth. Some comparison sources have opposite EW and narrow-normalized trends, reinforcing the aperture/calibration concern.

The plotted bars span model choices, not measurement errors. Orange marks inadequate window support. Missing fits are not non-detections.'''),
        code('''
display(Image(filename=str(OUT/'figures/latest_archive_model_envelope.png')))
cols=['name','epoch_id','hb_ew_total_cont_ratio_min','hb_ew_total_cont_ratio_max','hb_oiii_ratio_min','hb_oiii_ratio_max','both_hb_supported','both_oiii_supported']
display(stress[stress.comparison.eq('latest_archive')][cols])
display(pd.read_csv(OUT/'unavailable_fit_comparisons.csv'))
'''),
        md('''## Repeat visits and alternative standards

P9694's apparent September-to-October full-line growth is not stable across standards: the EW stress range includes unity. P8544 and P10381 likewise do not give a clean EW direction across all tests. P7837's apparent month-scale fading weakens to almost no change under one combination and has missing spectral support. None is promoted here as a validated rapid transition.

The alternative calibration on October 27 comes from a cloud-attenuated standard; it is intentionally a severe diagnostic, not an equally probable solution. Independent-night alternatives are allowed separately when forming repeat ratios. Narrow-line invariance is also an assumption to test across slit/fiber apertures.'''),
        code('''
display(stress[stress.comparison.eq('previous_ngps')][cols])
alt=pd.read_csv(OUT/'alternate_standard_fits.csv')
display(alt[alt.variant.eq('baseline')][['epoch_id','hb_flux_alternate_over_primary','hb_ew_total_cont_alternate_over_primary','hb_oiii_alternate_over_primary']])
'''),
        md('''## Independent exposures

These checks use native R-channel exposures under the adopted standard. Their Hβ and sideband windows avoid the strong telluric bands. They are uncoadded wing diagnostics, not replacements for the full-line decomposition. P7947 has greater exposure-to-exposure scatter than the two weakening cases.'''),
        code('''
display(Image(filename=str(OUT/'figures/independent_exposures.png')))
exp=pd.read_csv(OUT/'independent_exposure_indices.csv')
display(exp[exp.variant.eq('nominal')][['name','night','mjd','wing_ew','wing_ew_err_diagonal','wing_over_oiii','valid']])
'''),
        md('''## Recovery checks and their limits

Unit checks test grey-calibration invariance, host dilution without spurious broad-flux growth, known line changes across resolutions, flux units/redshift, and missing-data handling. Seventy-two two-epoch synthetic trials use 12 noise seeds, three true changes, and two residual scenarios. The extra-curvature/correlated-residual challenge exposes model bias: in the fading case the median recovered ratio is about 14% low.

These trials share much of the fitting model by design. They are a limited validation and a warning against taking the small model-family spread as the full uncertainty. They do **not** calibrate real-source completeness, false positives, or CLAGN discovery rates.'''),
        code('''
inj=pd.read_csv(OUT/'injection_recovery.csv')
display(inj.groupby(['challenge','true_hb_oiii_ratio']).fractional_error.agg(['median','std','min','max']))
fig,ax=plt.subplots(figsize=(8,4))
for (challenge,truth),g in inj.groupby(['challenge','true_hb_oiii_ratio']):
    offset=-.025 if challenge=='matched_model' else .025
    ax.scatter(np.full(len(g),truth+offset),g.recovered_hb_oiii_ratio,s=18,alpha=.6,
               color='#227981' if challenge=='matched_model' else '#bf7134',label=challenge if truth==.5 else None)
ax.plot([.35,2.4],[.35,2.4],'k--');ax.set(xlabel='True broad Hβ / [O III] change factor',ylabel='Recovered change factor')
ax.legend(fontsize=8);plt.show()
'''),
        md('''## Refreshed optical data

Public ALeRCE/ZTF corrected-alert photometry was queried separately from all frozen selection inputs. Require reliable correction flags, drb ≥ 0.8, error < 0.2 mag and a close reference association. The local timing criterion is at least three distinct nights from one object ID within ±30 days. It is a coverage diagnostic, not a precise flux-calibration prior.

There are now five pilot visits with at least one band meeting that criterion, including three October visits. P21422 has six r-band nights, but the nearest is ~19 days before spectroscopy; temporal uncertainty remains. Alert photometry is detection-selected and reference dependent. Lack of alerts is not evidence for an unchanged or undetected source. Do not mix these points with forced or survey photometry without addressing their systematics.'''),
        code('''
opt=pd.read_csv(OUT/'optical_visit_matches.csv')
display(opt[opt.n_nights_30d.gt(0)])
print('Visits with a qualifying band:',opt.loc[opt.usable_local_anchor,'epoch_id'].nunique())
display(Image(filename=str(OUT/'figures/refreshed_optical.png')))
'''),
        md('''## What to push next

**First:** jointly fit Hβ and Hα for P21422, P9506 and P7947, retaining all usable archival epochs, with shared host information and explicit aperture/calibration terms. Inspect wavelength-dependent residuals and compare supported photometry. Establish literature identity/status before using “new.”

**Then:** apply a validated measurement and detectability procedure to the full Hβ-accessible sample, including comparison and null-result sources. Test whether the manifold adds information beyond simple pre-run variability without reusing new outcomes in feature selection. The observed sample may support spectral-evolution associations more strongly than an independent forecast claim.

The strongest current paper path is the complete, carefully calibrated Palomar sample with a few convincing case studies. A P21422-only discovery paper and a successful manifold-prediction paper are both premature. The executed pilot has now shown which assumptions matter enough to change those decisions.''')
    ])


if __name__=='__main__':main()
