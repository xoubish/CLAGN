# MIRI/MRS precedents for the proposal and Figure 2

Checked 2026-09-29. This is a focused reading set, not a complete MIRI/MRS bibliography. Five public programme records and six paper PDFs are saved here; URLs and checksums are in `downloads.json`. No proposal PDF was compiled. The proposal source and Figure 2 are unchanged by this literature review.

## Public programme examples

These STScI PDFs contain abstracts, observing descriptions and APT configurations. They are **not the full submitted science-justification documents**; no full narrative was verified in this search. Programme records can include later observing revisions.

| Programme | Relevance |
|---|---|
| [3696, Cycle 2](https://www.stsci.edu/jwst/phase2-public/3696.pdf), hidden tidal disruptions, PI De | Closest dust-echo precedent: four infrared transients; connects dust spectroscopy and ionization diagnostics to flare energetics. Resulting paper: Masterson below. |
| [8245, Cycle 4](https://www.stsci.edu/jwst/phase2-public/8245.pdf), SDSS1335+0728, PI Sánchez Sáez | Direct precedent for MRS following a 2019 nuclear brightening and testing a changing accretion state. Resulting paper: Ansky below. |
| [1670, Cycle 1](https://www.stsci.edu/jwst/phase2-public/1670.pdf), AGN outflow launching, PI Shimizu | Six AGN selected for a controlled comparison; abstracts explain why the sample tests the question. Science and offset-background configurations are public. |
| [2016, Cycle 1](https://www.stsci.edu/jwst/phase2-public/2016.pdf), ReveaLLAGN, PI Seth | Seven nearby low-luminosity AGN; links continuum and line diagnostics to accretion structure and host separation. |
| [1328, Cycle 1 ERS](https://www.stsci.edu/jwst/phase2-public/1328.pdf), GOALS, PI Armus | Nuclear and surrounding emission; useful extraction/decomposition precedent. Its merging-galaxy selection differs from ours. |

## Priority reading and figures

1. **[Masterson et al. (2025)](https://arxiv.org/abs/2503.08647)**, *JWST’s First View of Tidal Disruption Events*. Table 1 places the observations 5.5–9.3 years after disruption. Figure 5 compares real spectra and silicate strengths with AGN; Figure 6 varies physical parameters; Figure 7 connects assumed illumination to spectral components. The authors explicitly label the latter an illustrative, non-self-consistent model, not a direct fit. Appendix B compares nuclear and annular spectra. **Use:** empirical context and a physical history-to-spectrum illustration. **Limit:** these dormant-nucleus transients do not establish detectability in our AGN. Their spline continuum definition differs from our current local indices. Full relevant sections read; Figure 5 visually inspected in `masterson_page8.png`.

2. **[Ramos Almeida et al., arXiv:2512.02629](https://arxiv.org/abs/2512.02629)**, *Silicate emission in a type-2 quasar*. Figure 2 (PDF page 6) shows two dust-model families, posterior spectral bands and percentage residuals for J1010 at z≈0.1. Figure 3 contains model images, not resolved observations of the torus. **Use:** a clear template for displaying uncertainty and spectral leverage. The fitted 15% uncertainty allowance is specific to their analysis, not a universal MIRI calibration error. Geometry, composition and grain assumptions differ between libraries. Static libraries alone do not supply a time-dependent response. Methods and captions read; Figure 2 visually inspected in `ramos_almeida_page6.png`.

3. **[Sánchez-Sáez et al. (2026), Ansky](https://arxiv.org/abs/2607.00921)**, *Spatially resolved optical and mid-infrared spectroscopy of SDSS1335+0728*. The source brightened in December 2019. Figure 8 presents the nuclear MRS spectrum; Figure 17 explicitly compares continuum definitions for silicate measurements. Section 3.2.1 says silicate properties alone cannot distinguish pre-existing diffuse nuclear dust from recently heated material. **Use:** make the information supplied by our time histories explicit. The optical gas light-echo reconstruction concerns much longer timescales than the monitored flare. Abstract, relevant MRS section and captions read; do not treat the proposed physical interpretation as unique.

4. **[González-Martín et al. (2025)](https://arxiv.org/abs/2504.01103)**, *JWST reveals the diversity of nuclear obscuring dust in nearby AGN*. Twenty-one nearby AGN, point-source/extended decomposition and multiple dust-model libraries. Available models fit 12 objects acceptably and fail for the remainder. **Use:** empirical host constraints and realistic model flexibility. A poor static fit is not by itself evidence of time evolution: chemistry and missing components can matter. Abstract and relevant methods/results inspected. The [2026 correction](https://academic.oup.com/mnras/article/549/3/stag1045/8706970) adds a funding acknowledgement.

5. **[Donnan et al. (2024)](https://arxiv.org/abs/2402.17479)**, *Peeling Back the Layers of Extinction*. Differential extinction and dust-temperature distributions from joint NIRSpec/MIRI spectra. **Use:** check whether extinction can mimic a continuum-shape difference. Downloaded; abstract and selected methods inspected, not a full figure audit.

6. **[Armus et al. (2023)](https://arxiv.org/abs/2209.13125)**, *GOALS-JWST: Mid-infrared Spectroscopy of the Nucleus of NGC 7469*. **Use:** the additional information in nuclear ionic lines and the circumnuclear environment. Downloaded; supporting reference rather than the main forecasting template.

Additional pointers (not downloaded here): [GATOS V, PAHs (2024)](https://arxiv.org/abs/2409.05686); [ReveaLLAGN 0 (2024)](https://arxiv.org/abs/2307.01252); [ReveaLLAGN I (2026)](https://arxiv.org/abs/2601.16977).

## Consequences for our Figure 2 — proposed design, not a validated forecast

The figure should answer: **What does MIRI add to the histories and SPHEREx spectra, and which interpretations remain distinguishable after fitting nuisance parameters?**

- Use a real target's measured history and existing spectral constraints to generate allowed MIRI predictions. Separate an illustration of physical response from a statistical prediction.
- Show a simulated spectrum with realistic uncertainties, competing fitted model families, and residuals. Do not interpret an envelope of selected grid curves as a posterior credible interval.
- Test recovery across all 24 histories, including null outcomes and sensitivity to host emission, geometry, extinction, calibration and unobserved illumination. History dependence is the primary test; changing dust structure is a stronger, secondary inference.
- A real AGN/TDE comparison can explain the dust diagnostics, but cannot substitute for the history test. Recompute diagnostics consistently before combining published values with our predictions, and check rest-wavelength coverage target by target.
- Preserve the audit's counterexamples: separated hand-selected spectra do not demonstrate family-level identifiability. A paper's illustrative model is likewise not evidence that our model-selection forecast is validated.

## Toward a broader bibliography

[STScI's JWST bibliography](https://archive.stsci.edu/publishing/bibliography/) and [publication statistics](https://archive.stsci.edu/jwst/bibliography/pubstat.html) provide a starting index. MIRI includes imaging and LRS, so an instrument-wide count is not an MRS count. A systematic follow-up should deduplicate by DOI/arXiv ID and verify MRS usage in the methods, then separate AGN/dust/transients from unrelated targets.

The [eJWST AGN observation catalogue](https://arxiv.org/abs/2601.14905) provides a complementary programme-discovery route, not a complete MRS paper list. Its stated coverage has a cutoff. Use the [MIRI technical library](https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-technical-library) for instrument/calibration references and the [known-issues page](https://jwst-docs.stsci.edu/known-issues/miri-known-issues/miri-mrs-known-issues) when setting simulation uncertainties.
