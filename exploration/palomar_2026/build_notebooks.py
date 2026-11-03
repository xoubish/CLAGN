"""Create reviewable notebooks. Generated/executed copies stay local and ignored."""
from pathlib import Path
import sys
import nbformat as nb

HERE = Path(__file__).resolve().parent
DEST = HERE/'notebooks'
SETUP = '''from pathlib import Path
import sys, os
ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p / 'scripts/observed_ngps.py').exists())
HERE = ROOT / 'exploration/palomar_2026'
os.environ['MPLCONFIGDIR'] = str(HERE / '.cache/matplotlib')
sys.path.insert(0, str(HERE))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, Markdown, Image
import analysis as A
from plots import load_tables, configure, dossier, year
configure()
T = load_tables()
objects, visits, epochs = T['objects'], T['visits'], T['epochs']
metrics, pairs, photo = T['line_indices'], T['comparisons'], T['photometry']
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 110)
'''
SPECTRA = '''# Recover the selected native spectra with their actual error arrays.
base_objects, base_visits, target_payloads, jwst = A.census()
spectra, archive_audit = A.load_spectra(base_objects, base_visits)
by_epoch = {s['epoch_id']: s for s in spectra}
'''


def md(text):
    return nb.v4.new_markdown_cell(text.strip())


def code(text):
    return nb.v4.new_code_cell(text.strip())


def save(name, cells):
    if sys.argv[1:] and name[:2] not in sys.argv[1:]:
        return
    book = nb.v4.new_notebook(cells=cells, metadata={
        'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
        'language_info': {'name': 'python', 'version': '3.12'}})
    nb.validate(book)
    nb.write(book, DEST/f'{name}.ipynb')


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    save('00_start_here_census_and_quality', [
        md('''# Palomar 2026 — start here
**Snapshot: 2026-11-03.** This is an executed discovery workbook for the completed
September 23 / October 26 / October 27 campaign. It uses the Palomar science notes
and current JWST proposal as requested. The outputs contain actual local data;
no simulated source measurements are included.

Read in order: **00 census → 01 spectra → 02 photometry → 03 source atlas →
04 population tests → 05 proposal and work plan**. The same notebooks are exported
as HTML under `outputs/html/`. The full source atlas is `outputs/all_sources_atlas.pdf`.
Code shared between notebooks is in `analysis.py` and `plots.py` one directory up.

These are screening measurements, not confirmed changing-look classifications.
The dated input manifest includes SHA-256 hashes; archives and original reductions
are read only. Executed notebooks and outputs remain local because some spectra
are collaboration data.'''), code(SETUP),
        code('''summary = visits.groupby('night', sort=False).agg(visits=('name','size'),
    distinct_targets=('name','nunique'), exposures=('exposures','sum'),
    total_seconds=('exposure_seconds','sum'), median_snr_G=('snr_G','median'))
display(summary)
display(pd.Series({'distinct targets': len(objects), 'delivered visits': len(visits),
    'archival epochs after screening': (epochs.source != 'NGPS').sum(),
    'targets with archives': (objects.n_archive > 0).sum(),
    'targets with repeated Palomar visits': (objects.n_NGPS > 1).sum(),
    'JWST targets observed': objects.jwst_id.notna().sum(),
    'targets with optical cache': objects.has_optical.sum(), 'targets with IR cache': objects.has_IR.sum()}))
display(objects[objects.n_NGPS > 1][['name','z','n_NGPS','n_archive','jwst_id']])'''),
        md('''## What was actually observed?
Use delivered visit manifests, not the planned lists, as the denominator. Names
are canonicalized: the telescope header P8828 is P8823, and prefixed JWST names
are linked to their proposal IDs. Z16CygABC and the stellar optical counterpart
of H5MUSES-327 remain in the inventory and atlas but are excluded from AGN line
screening. HIC-1953 and HIRAS-F23290-1252 are separate nearby-galaxy targets.'''),
        code('''fig, ax = plt.subplots(1,3,figsize=(15,4))
for night,g in visits.groupby('night',sort=False):
    ax[0].scatter(g.mjd, g.snr_G, label=night, s=24)
ax[0].set(xlabel='MJD (UTC)',ylabel='G continuum S/N per native pixel'); ax[0].legend()
ax[1].hist(objects.n_archive, bins=np.arange(0,objects.n_archive.max()+2)-.5)
ax[1].set(xlabel='Retained archival epochs / source',ylabel='Sources')
ax[2].scatter(objects.z, objects.spectral_baseline_rest_years, c=objects.n_archive, cmap='viridis')
ax[2].set(xlabel='Redshift',ylabel='Spectral baseline (rest years)')
plt.show()
display(objects[objects.n_archive == 0][['name','population','jwst_id','prior_question']])'''),
        md('''## Quality gates before discovery claims
* **September:** P330E and BD+28 disagree by roughly 41–53% in flux scale. This is
  an ambiguity between calibrations, not a Gaussian error to add in quadrature.
* **October 26:** Feige 34 is the delivered scale. Early transparency variations,
  slit throughput and the P8823 single accepted exposure need source-level checks.
* **October 27:** partly cloudy; observed absolute flux can be suppressed by grey
  losses. Comparisons in flux alone are unsafe. EW cancels a grey factor but still
  depends on host light, aperture, color calibration and blends.
* **All nights:** resampled-pixel errors omit covariance and calibration/host
  systematics. Telluric correction, arm joins, narrow-line spatial extent and
  spectral resolution are not interchangeable across instruments.

Bad metadata propagates across duplicate reductions of the same survey/day.
SDSS coadds without verified intervals are excluded. Verified DESI coadds up to
30 days are retained with both endpoints and cannot imply a one-night transition.
Native SDSS FITS are preferred. The historical JSON `ok` array is a rounded
occupancy fraction, **not** a pixel mask; its zero values must not delete spectra.'''),
        code('''display(T['archive_audit'].groupby(['status','reason'],dropna=False).size().rename('records').to_frame())
display(epochs.groupby(['source','error_kind']).size().rename('epochs').to_frame())
nominal = metrics[metrics.variant == 'nominal']
display(nominal.groupby(['source','line']).agg(attempts=('valid','size'),
    full_coverage=('valid','sum')))
# A rejected/partial index is missing evidence, never a non-detection.
display(nominal[(nominal.source=='NGPS') & ~nominal.valid][['name','night','line','coverage']].head(25))'''),
        md('''## Reproduce or inspect
From the repository root:
```sh
/opt/anaconda3/bin/python exploration/palomar_2026/analysis.py
/opt/anaconda3/bin/python exploration/palomar_2026/build_notebooks.py
/opt/anaconda3/bin/python exploration/palomar_2026/execute_notebooks.py
```
The notebook build regenerates notebook sources; save personal exploratory edits
under a different name before rebuilding. `outputs/provenance.json` gives the
input snapshot and software versions. `README.md` summarizes the findings and
`SCIENCE_PLAN.md` specifies the next measurements.''')])

    save('01_emission_lines_and_continuum', [
        md('''# Emission lines and continuum through time
**Question:** which sources changed their broad-profile wings, continuum color,
or line-to-narrow-line contrast? Three continuum choices and the corrected,
uncorrected, and alternative-standard NGPS spectra expose several systematics.

We measure signed wing indices of Hβ, Hα and Mg II. Narrow cores are excluded,
but Fe II, host absorption and other blends remain. These are **not decomposed
broad-line measurements**; widths, black-hole masses and Eddington ratios require
a later joint spectral fit.'''), code(SETUP), code(SPECTRA),
        md(r'''## Definition and inspectable implementation
For a local linear pseudo-continuum C(λ), the rest-frame wing index is
\( W=\int_{\mathrm{wings}}[f_\lambda/C(\lambda)-1]d\lambda_{\rm rest}\).
The signed line-flux proxy uses \(\int(f_\lambda-C)d\lambda_{\rm obs}\), hence the
extra factor of (1+z). Integrals use common interval endpoints and never cross
masked or large wavelength gaps. At least 98% wing coverage and both continuum
sidebands are required. Statistical errors propagate sideband uncertainty but
assume independent pixels; they are not total errors.

Ratios are only reported if both indices exceed 2 Å, an explicit exploratory
denominator floor rather than a detection threshold. Signed differences are
retained for weak/absorbing profiles. A 1.5× or 2/3× ratio stable across the three
continuum choices triggers review, not a discovery declaration.'''),
        code('''import inspect
display(pd.DataFrame(A.LINES).T)
print(inspect.getsource(A.line_index))'''),
        code('''# Change NAME and LINE and rerun this cell to inspect any source.
NAME, LINE = 'P17881', 'Hb'
trial = []
for s in spectra:
    if s['name'] == NAME:
        for variant,shift in A.VARIANTS.items():
            trial.append(dict(epoch=s['epoch_id'],variant=variant,**A.line_index(s,LINE,shift)))
display(pd.DataFrame(trial))
fig = dossier(NAME,spectra,T); plt.show()'''),
        code('''latest = pairs[pairs.baseline=='latest_archive'].copy()
leads = latest[latest.review_large_change].merge(objects[['name','known_state_status','jwst_id']],on='name')
display(leads[['name','night','line','old_mjd','ew_ratio','ratio_min','ratio_max',
    'oiii_normalized_ratio','telluric_fraction','cross_channel','known_state_status']].sort_values(['line','ew_ratio']))
fig,axes = plt.subplots(1,3,figsize=(15,4))
for ax,line in zip(axes,['Hb','Ha','MgII']):
    d = latest[(latest.line==line)&np.isfinite(latest.ew_ratio)]
    ax.scatter(d.old_ew,d.new_ew,c=d.telluric_fraction,cmap='magma',vmin=0,vmax=.5)
    lim=max(d.old_ew.max(),d.new_ew.max())*1.05
    ax.plot([0,lim],[0,lim],color='grey',ls='--')
    ax.set(xlabel='Previous wing index (Å)',ylabel='NGPS wing index (Å)',title=line)
plt.show()'''),
        md('''## The month-scale repeat experiment
Each source is its own comparison across September–October. Examine both EW and
wing/[O III] diagnostics. Disagreement between those normalizations is informative
about continuum, narrow-line aperture or calibration effects. A missing index
does not mean no change; inspect the coverage and native-channel products.'''),
        code('''repeat = pairs[pairs.baseline=='previous_NGPS']
display(repeat[['name','line','baseline_rest_years','ew_ratio','ratio_min','ratio_max',
    'oiii_normalized_ratio','uncorrected_ew_ratio_to_primary','delta_ew_err_diagonal']])
fig, ax = plt.subplots(figsize=(9,4))
r = repeat[(repeat.line=='Hb')&repeat.ew_ratio.notna()]
ax.scatter(r.name,r.ew_ratio,label='Hβ wing EW ratio',s=50)
ax.scatter(r.name,r.oiii_normalized_ratio,marker='s',label='wing/[O III] ratio',s=45)
ax.axhline(1,color='grey',ls='--'); ax.set_ylabel('October / September'); ax.legend(); plt.show()'''),
        md('''## Sensitivity to missing pixels in the five repeats
P10381 fails the strict September Hβ coverage cut. Here both NGPS epochs are
interpolated only within valid segments onto a common 3 Å observed-frame grid,
then restricted to their **shared** valid pixels. The relaxed 90% threshold is
only a labeled sensitivity experiment; these results do not replace the primary
98%-coverage catalog. No errors are invented after resampling and no
instrumental-resolution match is implied.'''),
        code('''shared=[]
for name in objects.loc[objects.n_NGPS>1,'name']:
    s=sorted([s for s in spectra if s['name']==name and s['source']=='NGPS'],key=lambda s:s['mjd'])
    a,b=A.common_support_pair(s[0],s[-1])
    for variant,shift in A.VARIANTS.items():
        x,y=A.line_index(a,'Hb',shift),A.line_index(b,'Hb',shift)
        coverage=min(x.get('coverage',0),y.get('coverage',0))
        ratio=y.get('wing_ew',np.nan)/x.get('wing_ew',np.nan) if x.get('wing_ew',0)>2 and coverage>=.9 else np.nan
        shared.append(dict(name=name,variant=variant,shared_coverage=coverage,ew_ratio=ratio))
shared=pd.DataFrame(shared)
display(shared)
shared.to_csv(A.OUT/'repeat_common_support_sensitivity.csv',index=False)'''),
        md('''## Continuum slope is a measurement, not automatically a disk slope
The slope below fits log fλ against log λ using positive median flux in available
4200–4230, 5080–5120 and 5600–5650 Å rest windows. Host light and Fe II remain;
windows that fall across different NGPS channels are particularly sensitive to
color calibration. Compare only identical window sets. Galactic reddening is
constant within a source but matters for cross-source absolute slopes.'''),
        code('''hb = latest[latest.line=='Hb'].merge(
    epochs[['epoch_id','alpha_lambda_total','color_windows']],left_on='old_epoch',right_on='epoch_id').rename(
    columns={'alpha_lambda_total':'alpha_old','color_windows':'windows_old'})
hb = hb.merge(epochs[['epoch_id','alpha_lambda_total','color_windows']],left_on='new_epoch',right_on='epoch_id',suffixes=('','_new'))
hb = hb[(hb.windows_old==hb.color_windows)&np.isfinite(hb.alpha_old)&np.isfinite(hb.alpha_lambda_total)]
hb['delta_alpha'] = hb.alpha_lambda_total-hb.alpha_old
fig,ax=plt.subplots(figsize=(7,4))
ax.scatter(hb.delta_alpha,hb.delta_ew)
ax.axvline(0,color='grey',ls=':'); ax.axhline(0,color='grey',ls=':')
ax.set(xlabel='Δ total pseudo-continuum αλ',ylabel='Δ Hβ wing pseudo-EW (Å)')
plt.show()
display(hb[['name','night','delta_alpha','delta_ew','windows_old']].sort_values('delta_alpha').head(12))'''),
        md('''This scatter is descriptive: EW and slope share continuum measurements,
so their errors are coupled. A line-responsivity claim must use joint posterior
draws, host subtraction, matched resolution and calibration nuisance parameters.
Do not assign a formal changing-look significance from diagonal errors alone.''')])

    save('02_photometry_timing_and_cross_band', [
        md('''# Photometry, timing and optical/IR context
**Question:** do line changes track optical variations, changing colors or delayed
IR emission? First establish actual overlap in time. Survey PSF photometry uses
one object ID per filter, clean flags and a 1.5″ association cut. Reference-corrected
alerts remain a separate stream. ZTF is AB; WISE is Vega; neither is an AGN-only flux.
No bins from different pipelines are merged to create an artificial trend.'''),code(SETUP),
        code('''summary = T['photometry_summary']
display(summary.groupby(['stream','band']).agg(sources=('name','nunique'),points=('n_points','sum')))
display(summary.groupby(['stream','band']).gap_to_latest_NGPS_days.median().rename('median gap days'))
fig,axes=plt.subplots(1,2,figsize=(13,4))
for band in ['g','r','W1','W2']:
    d=summary[(summary.band==band)&summary.stream.isin(['ZTF survey','NEOWISE visits'])]
    axes[0].hist(d.gap_to_latest_NGPS_days,bins=18,alpha=.45,label=band)
axes[0].set(xlabel='NGPS minus latest cached photometry (days)',ylabel='Sources'); axes[0].legend()
opt=summary[(summary.band=='g')&(summary.stream=='ZTF survey')]
ir=summary[(summary.band=='W1')&(summary.stream=='NEOWISE visits')]
joined=opt.merge(ir,on='name',suffixes=('_g','_W1'))
axes[1].scatter(joined.p95_p05_mag_g,joined.p95_p05_mag_W1)
axes[1].set(xlabel='g 95–5% amplitude (mag)',ylabel='W1 95–5% amplitude (mag)')
plt.show()'''),
        md('''The amplitude plot compares the cached baselines, which can differ in
duration and cadence. It is useful for identifying outliers, not for estimating
a transfer function or a population optical-to-IR response ratio.'''),
        code('''matches = T['photometry_matches']
ngps=matches[matches.source=='NGPS']
display(ngps[ngps.n_points>0])
# Explicit sensitivity to timing tolerance; no interpolation over seasonal gaps.
timing=[]
for window in [7,14,30,60,90]:
    matched=A.match_photometry(photo,epochs[epochs.source=='NGPS'],window)
    ok=matched[matched.n_points>0]
    timing.append(dict(window_days=window,matched_sources=ok.name.nunique(),matched_epochs=ok.epoch_id.nunique()))
display(pd.DataFrame(timing))
nearest=ngps.groupby('name').nearest_days.min().sort_values()
display(nearest.head(20).rename('nearest optical detection days'))'''),
        md('''## Optical color evolution with measured simultaneity
Pair g and r within one day, within the same stream. The following cell uses
one median per night and one match per g night. No cross-stream zero point is
assumed. Coloring by date can reveal a loop, but a loop alone does not establish
hysteresis: seasonal sampling and broad-line contamination require checks.'''),
        code('''NAME='P9694'
d=photo[(photo.name==NAME)&(photo.stream=='ZTF survey')].copy()
d['night']=np.floor(d.mjd).astype(int)
nightly=d.groupby(['night','band']).agg(mjd=('mjd','median'),mag=('value','median')).reset_index()
g=nightly[nightly.band=='g'][['mjd','mag']].sort_values('mjd').rename(columns={'mag':'g'})
r=nightly[nightly.band=='r'][['mjd','mag']].sort_values('mjd').rename(columns={'mag':'r','mjd':'r_mjd'})
matched=pd.merge_asof(g,r,left_on='mjd',right_on='r_mjd',direction='nearest',tolerance=1).dropna()
# Retain one g point per matched r night to avoid duplicating a noisy datum.
matched['dt']=abs(matched.mjd-matched.r_mjd)
matched=matched.sort_values('dt').drop_duplicates('r_mjd').sort_values('mjd')
fig,axes=plt.subplots(1,2,figsize=(12,4))
axes[0].scatter(year(matched.mjd),matched.g-matched.r,s=10)
axes[0].set(xlabel='Year',ylabel='g − r (AB)',title=NAME)
sc=axes[1].scatter(matched.g,matched.g-matched.r,c=year(matched.mjd),s=12,cmap='viridis')
axes[1].invert_xaxis(); axes[1].set(xlabel='g (brighter →)',ylabel='g − r (AB)')
fig.colorbar(sc,ax=axes[1],label='Year');plt.show()
display(matched[['mjd','r_mjd','g','r','dt']].head())'''),
        md('''## What is justified next?
1. Retrieve total-flux optical photometry around each Palomar visit, with
   association/field/reference checks. Preserve upper limits and signed fluxes.
2. Compute true synthetic g/r fluxes using throughput curves on fully covered
   spectra; compare to date-matched photometry. A pivot-wavelength continuum is
   not a synthetic magnitude. Fit grey/color calibration factors jointly with
   host light; never rescale to stale photometry.
3. Measure line/continuum response and color–luminosity tracks per source using
   joint uncertainties, repeated spectra and actual time intervals.
4. Test optical–IR delays only where optical driver and IR coverage overlap.
   Use red-noise and cadence-matched null simulations, lag-recovery injections,
   and held-out seasons; a maximum correlation in sparse WISE sampling is not
   a lag detection. Optical–IR luminosity-state mismatch can be reported earlier.

The cached data are rich in history but sparse at the new spectral endpoints.
The timing table above is the retrieval specification, not evidence of stability.''')])

    save('03_individual_sources_and_atlas', [
        md('''# Individual-source discovery atlas
Every observed source has a six-panel dossier and a row in the audit tables.
The atlas includes stars and ancillary galaxies so that unexpected objects are
not erased by an AGN-only selection. The spectra show all retained epochs,
with early/latest archives and every NGPS visit emphasized.

The shortlist below mixes new screening leads with deliberately retained
time-domain cases. No source is called a new CLAGN on the basis of an index.
Novelty requires a source-specific literature and catalog audit.'''),code(SETUP),code(SPECTRA),
        code('''# Rebuild all 99 figures and a single PDF. Set False to inspect an existing atlas.
REBUILD_ATLAS=True
if REBUILD_ATLAS:
    from plots import export_atlas
    atlas=export_atlas(spectra,T)
else:
    atlas=pd.read_csv(A.OUT/'dossier_index.csv')
display(Markdown('[Open all-source PDF](../outputs/all_sources_atlas.pdf) · '
                 '[Browse individual sources](../outputs/source_index.html)'))
display(atlas)'''),
        md('''## Sources to inspect first
**P17881:** compare both Balmer lines with the recent archive; establish whether
strengthening is nuclear, host dilution, aperture, or calibration. It already
belongs to the JWST sample, so novelty is not assumed.

**P21422, P24930, P7947:** Hβ weakening/strengthening review leads. Check archival
baseline age, broad/narrow separation, narrow-line behavior and optical history.

**P9694:** the September strengthening lead now has an October repeat. Distinguish
the persistence of a stronger wing index from a new September–October change.

**P8544:** test the earlier weakening interpretation against the October repeat,
including continuum-window and telluric sensitivity.

**P11530 and P1823:** preserve the rich long-term histories even if the latest
Palomar spectrum does not establish a new transition.

**Mg II leads:** Fe II and pseudo-continuum placement make these valuable fitting
targets; the index itself cannot establish Mg II broad-line variability.

**H5MUSES-327:** existing work identifies a stellar optical spectrum at the IR
position. Dusty-star versus blended-background-source interpretations require
spatial/spectral separation; this is a distinct discovery track.'''),
        code('''for name in ['P17881','P21422','P24930','P7947','P9694','P8544','P11530','P1823']:
    display(Markdown('### '+name))
    display(Image(filename=str(A.OUT/'dossiers'/f'{name}.png')))
    display(pairs[(pairs.name==name)&(pairs.line=='Hb')][['baseline','night','old_mjd','new_mjd',
        'ew_ratio','ratio_min','ratio_max','oiii_normalized_ratio']])'''),
        md('''## Any source, editable in one place
Use this cell to regenerate a dossier from its native spectra. For a line near
an NGPS arm boundary, return to the individual channel files documented in the
reduction README; the merged spectrum has not been renormalized across joins.'''),
        code('''NAME='P10252'  # replace with any name in objects.name
fig=dossier(NAME,spectra,T);plt.show()
display(epochs[epochs.name==NAME][['epoch_id','mjd','min_mjd','max_mjd','source','path','error_kind']])''')])

    save('04_population_and_selection_experiments', [
        md('''# Population and selection experiments
**Primary question from the Palomar notes:** does the time-domain manifold enrich
the yield of spectroscopic state changes beyond brightness, redshift, variability
amplitude, coverage and observational selection?

The run is a targeted heterogeneous follow-up sample. These 99 sources are not a
random sample of all AGN. Planned experimental labels were attached before the
new spectra; use actual observed sources for denominators and preserve known
CLAGN follow-up as a separate stratum. Missing data are not null outcomes.'''),code(SETUP),
        code('''latest=pairs[(pairs.baseline=='latest_archive')&(pairs.line=='Hb')].sort_values('new_mjd').drop_duplicates('name',keep='last')
sample=objects.merge(latest,on='name',how='left',suffixes=('','_pair'))
sample['has_Hb_pair']=sample.new_epoch.notna()
sample['review_flag']=sample.review_large_change.fillna(False).astype(bool)
display(sample.groupby('experiment_group').agg(observed=('name','size'),
    archived=('n_archive',lambda x:int((x>0).sum())),with_Hb_pair=('has_Hb_pair','sum'),
    Hbeta_review_flags=('review_flag','sum')))
display(sample.groupby('known_state_status').agg(observed=('name','size'),with_Hb_pair=('has_Hb_pair','sum')))
fig,axes=plt.subplots(1,2,figsize=(13,5))
for label,g in sample.groupby('experiment_group'):
    axes[0].scatter(g.z,g.spectral_baseline_rest_years,s=24,label=label)
axes[0].set(xlabel='Redshift',ylabel='Available spectral baseline (rest yr)')
axes[0].legend(fontsize=7)
ok=sample.ew_ratio.gt(0)&sample.umap_x.notna()
sc=axes[1].scatter(sample.loc[ok,'umap_x'],sample.loc[ok,'umap_y'],
    c=np.log10(sample.loc[ok,'ew_ratio']),cmap='coolwarm',vmin=-.5,vmax=.5,s=38)
axes[1].set(xlabel='Frozen UMAP x',ylabel='Frozen UMAP y',title='Latest Hβ pair per source')
fig.colorbar(sc,ax=axes[1],label='log10 NGPS/archive wing EW ratio');plt.show()'''),
        md('''## Source-level uncertainty for a descriptive statistic
This bootstrap resamples sources, not individual spectral pairs. Its interval
describes sampling variability of the median log index ratio in the observed,
measurable sample. It omits measurement/calibration uncertainty and cannot be
interpreted as a population-changing-look rate or evidence for causality.'''),
        code('''rng=np.random.default_rng(20261103)
rows=[]
for label,g in sample[sample.ew_ratio.gt(0)].groupby('experiment_group'):
    x=np.log10(g.ew_ratio.to_numpy())
    if len(x)<3:
        rows.append(dict(group=label,n=len(x),median_log_ratio=np.median(x),lo=np.nan,hi=np.nan));continue
    boot=np.median(rng.choice(x,size=(4000,len(x)),replace=True),axis=1)
    rows.append(dict(group=label,n=len(x),median_log_ratio=np.median(x),
        lo=np.quantile(boot,.025),hi=np.quantile(boot,.975)))
description=pd.DataFrame(rows)
display(description)
description.to_csv(A.OUT/'descriptive_group_bootstrap.csv',index=False)'''),
        md('''## Freeze the confirmatory experiment before fitting outcomes
* Define an object-level state transition using joint broad-line fits, meaningful
  upper limits, visual review and calibration checks. Keep ordinary variability,
  extreme variability and a changing-look transition distinct.
* Primary contrast: previously unconfirmed sources inside vs outside the frozen
  selected manifold regions. A/B/C/D labels are pre-run hypotheses, not inferred
  current states. Catalog matches need verified literature status.
* Match/adjust for redshift, archival brightness and host fraction, line
  detectability, archival-to-Palomar interval, optical/IR coverage and trigger
  selection. Check covariate overlap before fitting any model; a few comparison
  sources may not support a stable adjusted enrichment estimate.
* Recover the actual scheduling and backup decisions. Unknown inclusion
  probabilities prevent a defensible inverse-probability population rate.
  Restrict the claim to observed-sample enrichment if necessary.
* Inject broad lines over a width/flux grid into each epoch with real masks,
  noise and spectral resolution. Use recovery probabilities and upper limits,
  not a universal S/N cut, to model non-detections.
* Analyze state change as interval censored between observed epochs. Repeated
  transitions must be individually demonstrated. Do not turn the count of
  spectral pairs into the number of independent objects.
* Reserve Hβ as the primary predeclared diagnostic. Treat Hα, Mg II, He II,
  asymmetry, centroid shifts and continuum features as secondary searches;
  control false discoveries across the tested source×feature family after
  valid null simulations. Discovery and confirmation datasets should be distinct.

Additional opportunities: differential Hβ/Mg II response; broad-line breathing;
color–luminosity hysteresis; optical–IR state mismatch; rare asymmetric profiles;
and the location of these outcomes on the **frozen** manifold. Do not retrain UMAP
on the new labels and then present that same sample as independent validation.''')])

    save('05_proposal_mapping_and_next_measurements', [
        md('''# Proposal mapping and a concrete discovery work plan
The agreed scope is the Palomar science notes plus `jwst_proposal/proposal.tex`.
The original Palomar proposal was not used. The plan below separates measurements
already made here from tests requiring joint modeling or new observations.'''),code(SETUP),
        code('''display(T['jwst_coverage'][['id','internal_id','family','n_NGPS','n_archive','has_optical','has_IR']])
missing=T['jwst_coverage'][T['jwst_coverage'].n_NGPS.isna()]
display(Markdown('**Proposal targets without a delivered Palomar spectrum in this snapshot:** '+', '.join(missing.internal_id)))
'''),
        md('''## Science-to-measurement map
| Science question | What Palomar contributes now | What completes the test |
|---|---|---|
| Manifold predicts spectral state/transition | Frozen positions, observed denominator, multi-epoch line indices | Joint state fits; selection/completeness audit; controls |
| Dated photometric trigger predicts broad-line response | Optical/IR history and line changes with dates | Current total-flux photometry; calibration/host model |
| Individual recurrence and transition duration | All retained spectral epochs and five repeat targets | Broad-line state sequence; interval-censored transition dates |
| Differential line and continuum response | Hβ/Hα/Mg II wing indices and total continuum colors | Host + Fe II + narrow/broad decomposition; joint uncertainties |
| JWST warm-dust memory | New optical endpoints for 23/24 proposal targets | Actual MIRI spectra, optical driver extension, lag/host uncertainty |
| Dust replenishment vs fixed dust geometry | Optical/IR history constrains illumination alternatives | Forward models of continuum/silicates/PAHs; future MIRI data |
| Rare or blended objects | All-source atlas includes stellar/ancillary targets | Spatial association, stellar subtraction, ancillary archive work |

The JWST proposal's fading/rising labels describe selected histories. They are
not Palomar spectral state classifications. SPHEREx/WISE predictions and scenario
tables in the proposal are model products, not new observed MIRI data.'''),
        code('''# A machine-readable, per-source work queue with explicit missing evidence.
latest=pairs[pairs.baseline=='latest_archive']
flags=set(latest.loc[latest.review_large_change,'name'])
matches=T['photometry_matches']
contemporary=set(matches.loc[(matches.source=='NGPS')&(matches.n_points>0),'name'])
matched_visits=matches[(matches.source=='NGPS')&(matches.n_points>0)].groupby('name').epoch_id.nunique()
rows=[]
for obj in objects.itertuples(index=False):
    todo=[]
    if obj.population=='stellar/special':
        todo=['special-object classification / association; no AGN incidence inference']
    else:
        todo.append('joint host + FeII + broad/narrow fits; native resolution and covariance QA')
        if obj.n_archive==0:todo.append('retrieve / identify a valid archival spectrum')
        if matched_visits.get(obj.name,0)<obj.n_NGPS:todo.append('retrieve total-flux photometry for unmatched NGPS visits')
        if obj.n_NGPS>1:todo.append('paired month-scale NGPS fit with night-specific calibration')
        if obj.name in flags:todo.append('visually adjudicate robust wing-index trigger')
        if pd.notna(obj.jwst_id):todo.append('update JWST optical-state endpoint and dust-driver constraints')
    priority=0 if obj.name in flags or obj.n_NGPS>1 else (1 if pd.notna(obj.jwst_id) else 2)
    rows.append(dict(name=obj.name,priority=priority,jwst_id=obj.jwst_id,
                     review_trigger=obj.name in flags,contemporary_optical_visits=matched_visits.get(obj.name,0),
                     n_NGPS=obj.n_NGPS,
                     tasks='; '.join(todo)))
queue=pd.DataFrame(rows).sort_values(['priority','name'])
queue.to_csv(A.OUT/'work_queue.csv',index=False)
display(queue)'''),
        md('''## Ordered execution and stopping criteria
**1 — Data and calibration.** Recover missing native DESI uncertainties/resolution
and missing archival spectra; refresh endpoint photometry. Inspect individual
NGPS exposures, channel overlaps, sky/telluric residuals and extraction images for
every lead. Freeze the sample, identity ledger and known-object status. Stop
calling an effect intrinsic if it disappears under supported calibration choices.

**2 — Joint multi-epoch spectral fitting.** Galactic extinction; a consistent
rest frame; instrument-specific line-spread functions; host templates with
aperture-dependent amplitudes; Fe II; local/power-law continuum; linked narrow
doublets with physically justified ratios; flexible broad Balmer/Mg II profiles.
Fit grey and color calibration nuisance terms constrained by simultaneous
photometry and qualified narrow-line checks. Save model components, residuals,
posteriors, upper limits, covariance and fit failure flags. Validate with
injection/recovery and repeat-exposure tests. No forced broad component.

**3 — Discovery adjudication.** Inspect high-priority sources blind to manifold
labels where possible. Record unchanged/variable/extreme-variable/CL-candidate/
unusable states with evidence and reasons. For profile shifts, verify wavelength
calibration and narrow-line reference; for Balmer decrements, separate broad
components and account for cross-arm calibration. Literature audit comes before
claims of a new object or new transition.

**4 — Time-domain science.** Fit continuum and broad-line evolution jointly with
optical light curves. Measure response slopes, recurrence intervals and possible
breathing. Model sparse IR response only where simulations recover a lag; retain
nonidentifiable lags as such. Distinguish fading history from present-day state.

**5 — Population and proposal delivery.** Produce observed-sample enrichment,
detectability-adjusted outcomes if identifiable, and matched controls with
object-level uncertainty. Update the 24-target JWST table with the actual
Palomar endpoint and data gaps. Fixed-geometry vs evolving-dust conclusions await
appropriate IR spectroscopy and driver constraints.

**Deliverables to aim for:** a validated multi-epoch line/continuum catalog; a
reviewed variable-source catalog; an empirical manifold-selection test; and a
JWST optical/IR history package. A small set of especially unusual sources may
support separate case studies after the calibration and novelty checks.'''),
        md('''## Literature and project anchors
* [Hemmati et al. 2026, local manifold paper](../../../reference/papers/Hemmati_2026_ApJ_998_130.pdf): frozen manifold motivates the prospective spectral test.
* [Palomar science review](../../../observing/science_review/README.md): dated-change, blind-manifold and comparison groups; documented selection limitations.
* [Current JWST proposal](../../../jwst_proposal/proposal.tex): warm-dust memory and replenishment questions.
* [Zeltyn et al. 2024](https://arxiv.org/abs/2401.01933): multi-epoch spectroscopy and optical light curves distinguish CL behavior from other extreme variability.
* [Guo et al., DESI statistical study](https://arxiv.org/abs/2408.00402): narrow-line, pseudo-photometric and visual flux-calibration checks motivate multiple diagnostics.
* [Ren et al. 2024](https://arxiv.org/abs/2406.17598): prior-informed AGN/host decomposition addresses host-template degeneracy.

These sources guide the analysis design; none establishes the status or novelty
of an individual source in this new campaign. See `SCIENCE_PLAN.md` for the
consolidated scope and acceptance criteria.''')])
    print(f'Wrote {len(list(DEST.glob("*.ipynb")))} notebooks to {DEST}')


if __name__ == '__main__':
    main()
