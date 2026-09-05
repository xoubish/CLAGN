# Spectroscopic follow-up of CLAGN candidates: working plan

Started 2026-09-03. Working folder: `CLAGN/`. Legacy material from the AGNzoo repo is in
`from_AGNzoo/` (see its README for provenance). Published paper: `Hemmati_2026_ApJ_998_130.pdf`.

Strategy: **build wide, score, cut for observability last.** One master table with one row per
source; keep everything, filter at the end.

---

## 0. Where things stand

What exists and runs:

- `wise_manifold_zeltyn.ipynb` rebuilds the W1-only UMAP on Sample A, fetches WISE light curves
  for the Zeltyn+2024 CL-AGNs and EVQs at z<1, projects them, and defines a CLAGN-enriched region.
  Its last run gave:

  | Quantity | Value |
  |---|---|
  | Sample A objects after GP cut / NaN drop | 1960 |
  | Zeltyn CL-AGNs projected (z<1) | 111 |
  | Zeltyn EVQs projected (z<1) | 68 |
  | Enriched region (3x, 10x10 bins) | 8 of 100 bins |
  | Zeltyn CL-AGNs inside region | 41 / 111 (37%) |
  | Tier 1: EVQs inside region | 17 |
  | Tier 2: non-CLAGN Sample A inside region | 121 |
  | Tier 4: control pool in anti-CLAGN region | 361 |
  | CLAGN fraction at Zeltyn locations vs overall (Section 10) | **0.8x** |

- `data/sampleA_embedding.csv`, `data/zeltyn_embedding.csv`: umap_x, umap_y, label only.
- `data/zeltyn2024_table5.csv`: 204 rows, Name / Class / z (113 CL-AGN, 88 EVQ, 3 CL-AGN RM).

Added 2026-09-03 (scripts at top level, outputs in `data/`, large files gitignored):

- `01_rebuild_manifold.py` -> `sampleA_embedding_objectid.csv` (1960 rows with **objectid**, label, umap_x/y,
  fvar/mean/max W1, kNN `clagn_score`, `in_region_clagn`, `in_region_zeltyn`), `umap_w1_model.pkl` (fitted
  UMAP for `transform` of new pools), `zeltyn_embedding_full.csv` (all 204 Zeltyn rows with RA/Dec, 190 projected),
  `region_stats.txt`, `wise_cache/zeltyn/` (unWISE W1/W2 light curves, resumable chunks).
- `03_spectra_inventory.py <csv> <tag>` -> `spectra_epochs_<tag>.csv`, `spectra_summary_<tag>.csv`
  (SDSS DR19 allspec + DESI DR1). Run for `zeltyn`.
- `data/zeltyn_coords.csv` (J-names -> RA/Dec), `data/zeltyn_candidates_prelim.csv` (embedding + spectra join,
  sorted by in-region then clagn_score: the Tier 2/3 pool).
- `data/df_lc_020724.parquet.gzip` re-downloaded (243 MB, 8307 objects incl. 6209 SPIDER-only).
- 2026-09-03 evening: full pipeline run on the r<19 pool (`run_after_wise.sh`) + `rescore.sh` with NEOWISE:
  `targets_sep23.csv` 18 primaries / 20 backups (4.2 of 4.3 h), `targets_oct26.csv` 41 / 45 (9.4 of 9.5 h),
  `targets_oct27.csv` 41 / 45 (9.4 of 9.5 h); median exposure 8 min, median r ~ 18.0. Finder charts in `finders/`,
  cutouts in `data/cutouts/`, page at web/clagn_night_sheet.html (published artifact, About tab included).
  Pending: 19-19.5 mag WISE fetch (then re-run `run_after_wise.sh`), SDSS-V internal epochs, NGPS ETC check.

What is missing (gitignored, never committed anywhere):

- **Sample A coordinates.** Embedding row -> `objectid` (via `keepsw`) -> row of the sample table.
  The sample table was written by the Feb-2024 collector to `output/<name>_sample.ecsv`
  (`from_AGNzoo/light_curve_generator_feb2024.md`) and by `from_AGNzoo/build_sample.md` to
  `data/AGNsample_26Feb24.csv` / `data/sample.ecsv`. Look on Fornax or the Feb-2024 machine, or re-run.
- `data/df_lc_020724.parquet.gzip` (gdown id `1gb2vWn0V2unstElGTTrHIIWIftHbXJvz`), and the Zeltyn/EVQ
  WISE parquets. The notebook re-downloads / re-fetches them.
- Section 5 of the "P200 Target Selection" cell only prints counts; it never builds the
  objectid -> RA/Dec list.

Known caveat to fix before scoring: the region was defined from where the *Zeltyn* objects land,
yet their CLAGN enrichment came out at 0.8x. Define the region from Sample A Turn-on/Turn-off
(which the paper shows do separate), then test on Zeltyn as a held-out set.

---

## 1. Step zero: recover coordinates

- [ ] Locate the Feb-2024 sample table (ecsv/csv) on Fornax or wherever the collector ran.
- [ ] If not found: re-run `from_AGNzoo/build_sample.md` (run from inside `from_AGNzoo/`, it needs
      that folder's `code_src/`). Expect ~2100 rows with SkyCoord, redshift, label bitmask.
- [ ] Verify the mapping: `objectid` in the parquet == row index of the sample table. Spot-check a
      few Turn-on/Turn-off objects against the literature coordinates in
      `from_AGNzoo/ChangeLook_Sample_oct2023.ipynb`.
- [ ] Save `data/sampleA_coords.ecsv` (objectid, ra, dec, z, label_bits, labels_decoded).
- [ ] Finish Section 5 of the P200 cell: embedding row -> objectid -> coords, write
      `data/sampleA_embedding_coords.csv`.

---

## 2. Stage 1: master list (build wide)

Sources to union into one table, with a `source_catalog` column:

- [ ] Sample A (~1960 with embeddings).
- [ ] Zeltyn+2024 CL-AGNs and EVQs (z<1 subset already projected; keep the z>1 rows too, flagged).
- [x] Literature CLAGNs (`00_literature_clagn.py` -> `data/literature_clagn.csv`, 312 unique positions at 2":
      Graham19 145, Lyu22 54, Hon22 29, Yang18 20, MacLeod16 17, MacLeod19 16, Green22 16, Sheng20 7,
      Lopez-Navas22 4, Ruan16 3, LaMassa15 1). NED is flaky (timeouts / HTML errors): retry with 90 s pauses.
      Used in `04_score_tiers.py` to exclude published CLAGNs from Tier 1 (flag `is_literature_clagn`).
- [ ] Optional bigger parent pool if the in-region set thins out after cuts: all DR16Q quasars at
      z<0.8, r<19.5, Dec>-25, in the observable RA window -> `wise_get_lightcurves` -> same GP/normalise
      pipeline -> `mapp.transform`. Builds the candidate pool from observable objects in the first place.
      Note NEOWISE ended 2024-08, so WISE light curves stop there.

Output: `data/master_list.ecsv`, one row per unique position.

---

## 3. Stage 2: enrich every row

### 3a. Archival spectra inventory (the "have / have not" question)

Answer: **both, for different tiers.**
- Discovery targets = in CLAGN region + at least one archival spectrum as baseline + never reported
  as changing. Zero-spectrum objects drop to low priority: one new spectrum cannot show a transition.
- Already-changed objects (Zeltyn, literature CLAGNs, archival class-change flags) = revisit tier and
  validation set for whether the manifold region predicts change.

- [ ] SDSS DR19 (latest release; includes SDSS-V repeat spectra): per epoch MJD, plate/fiber or
      field/catalogid, pipeline `class`, `subclass`, z.
- [ ] DESI DR1 (public 2025; heavy overlap with SDSS at z<1): same per-epoch info. SDSS + DESI pair
      = free second epoch for many sources.
- [ ] LAMOST DR(latest): same.
- [ ] Derive: `n_spec`, `mjd_first`, `mjd_last`, `years_since_last_spec`, `class_per_epoch`,
      `class_change_flag` (pipeline class/subclass differs between epochs -> archival CLAGN candidate).

### 3b. Photometry and brightness

- [ ] PS1 or SDSS g r i (static).
- [ ] Current ZTF g, r (IRSA ZTF DR, or ZTF forced photometry) and ATLAS if needed for post-2024 epochs.
- [ ] W1, W2 mean and latest.
- [ ] Brightness at the epoch of the last spectrum (interpolate light curve at `mjd_last`).
- [ ] Galactic extinction E(B-V); ecliptic latitude (WISE cadence).

### 3c. Manifold position and score

- [ ] `umap_x`, `umap_y` for every row (Sample A from the fit; everything else via `transform`).
- [ ] Continuous `clagn_score`: fraction of Turn-on/Turn-off among the k nearest Sample A neighbours
      in embedding space (k ~ 30-50), instead of a binary bin flag.
- [ ] Distance to Turn-on centroid and Turn-off centroid separately.
- [ ] Held-out test: score distribution for Zeltyn CL-AGNs vs EVQs vs SDSS_QSO. Record the result.

### 3d. Variability and time ordering (the strongest prioritiser)

- [ ] W1 amplitude, ZTF amplitude over full baseline.
- [ ] **Change since last spectrum** in W1 and in ZTF g/r: flux(now) / flux(mjd_last).
- [ ] Current trend direction (rising / fading over last ~2 yr).
- [ ] Expected transition: type-1 baseline + fading -> turn-off candidate; galaxy/type-2 baseline +
      brightening -> turn-on candidate.

---

## 4. Stage 3: score, tier, exclude

Tiers (revise as the data come in):

| Tier | Definition | Purpose |
|---|---|---|
| 1 | EVQs / high-score objects with baseline spectrum and large change since last spectrum | Best discovery odds |
| 2 | Non-CLAGN Sample A in region with baseline spectrum | Discovery |
| 3 | Confirmed CLAGNs (Zeltyn, literature) — 5-10 brightest | Current state, validation |
| 4 | Control: anti-CLAGN region, matched in z and r to Tier 1+2 — 10-15 | Needed to claim the manifold predicts change |

Exclusions (`exclude_flag`, `exclude_reason`):
- [ ] Blazars / Fermi / radio-loud labels in Sample A.
- [ ] Spectrum within the last ~2 yr *unless* the light curve moved since.
- [ ] Published CLAGNs outside Tier 3.
- [ ] Stars / contaminants, bad z.

- [x] Priority score (implemented in `04_score_tiers.py`, 2026-09-03):
      M = clip(max(zeltyn_density_ratio/3, clagn_score/0.15), 0, 2)  (manifold CLAGN-likeness, 1 = threshold)
      P = clip(|log10(flux_now/flux_at_last_spec)| / log10(1.5), 0, 2)  (W1 from unWISE to 2020-12, ZTF r to 2025; larger of the two)
      S = 1 if years_since_last_spec >= 3 else 0.3;  B = 1.0 / 0.8 / 0.5 / 0.25 for r < 18.5 / 19 / 19.5 / fainter
      **priority = B * S * (M + P) + 0.5 * class_change_flag**;  per night **priority_night = priority * W_moon**
      with W_moon = 1.0 (>= 60 deg from the moon) / 0.7 (40-60) / 0.4 (30-40); < 30 deg excluded.
      Tiers: T3 Zeltyn CL-AGN; T2 Zeltyn EVQ with M >= 1 or in Zeltyn region; T1 pool object with M >= 1; T4 control
      (pool, clagn_score = 0, zeltyn_density_ratio <= 0.3, outside both regions). Allocation per night: hrs >= 1.5
      above airmass 2 inside the window, z <= 0.8 (T3 exempt); slots 14 / 34 / 34; caps T3 <= 3, T4 <= 4 per night;
      floors so every proposal tier is represented: Sep 23 T3 1 / T2 1 / T4 2, October nights T3 3 / T2 3 / T4 4.
      Controls (T4) use the opposite score, priority = B * S * (1 - min(P, 1)): bright, photometrically quiet objects
      off both CLAGN regions. Rank within a night = order of priority_night among the selected set (floors do not
      change ranks). P now takes the largest of: unWISE W1 (to 2020) vs W1 at the archival spectrum, NEOWISE-R W1
      (to 2024) vs at-spectrum, NEOWISE-R 2024 vs 2014, and ZTF r now vs r at the last spectrum.
      State-reversal candidates (2026-09-03): from the downloaded archival spectra, `spec_dir` = "turned off" if the
      rest-frame Hβ EW fell by > 50% between the first and last SDSS epoch (first EW > 8 Å), "turned on" if it more than
      doubled (last EW > 8 Å). If the current photometry (ZTF r since the last spectrum, NEOWISE W1 2014-24, W1 slope
      2022-24) points the other way, the target is a reversal candidate and gets +0.5 priority. Motivation: recurrent
      CLAGNs show dim-state plateaus of ~4-7 yr before re-brightening (Wang et al. 2025, ApJ 981, 129, eight recurrent
      CLAGNs; SDSS J1011+5442 back to type 1 in 2024 after turning off 2003-15), and turn-on flares can fade within
      months to years (1ES 1927+654). Objects that turned off in 2019-21 are therefore due, and 2018-22 turn-ons may
      already be fading.
      **No magnitude cut** (user decision 2026-09-03): brightness and moon distance only weight the ranking.
      The pool WISE fetch therefore covers r < 19.5 (r < 19 first, 19-19.5 added afterwards); 19.5-20 not fetched.
- [x] Write `data/master_list_scored.csv` and `data/targets_<night>.csv` (top list + backups).

---

## 5. Stage 4: observability (last)

- [x] Nights: **P200/NGPS, 2026-09-23 first half** (astro. twilight 20:05, window to 00:39 PDT, LST 19.5h-0.1h,
      moon 93% at RA 22.2h Dec -12, up all window) and **2026-10-26 + 10-27 full nights** (19:24-05:38 PST,
      LST 21h-7.2h; moon 98%/94% at RA 3-4h Dec +22..+26, rises ~3.5 h after twilight so the first third of
      each October night is moon-free). All bright time -> bright targets and/or > 40 deg from the moon.
      Computed by `05_observability.py` (astroplan); it adds hrs_<night>, minX_<night>, moonsep_<night> columns.
- [x] RA windows used for the parent pool (`02_parent_pool.py`): 15.5h-24h, 0h-4.5h, 6.5h-11h; Dec > -15.
- [x] Exposure model (2026-09-03, in `04_score_tiers.py`): t_exp = 9 min * (7/10)^2 * 10^(0.8 (r-18.5)) * 10^(0.4 dsky),
      dsky = 1 mag (Ha usable, z <= 0.55) or 2 mag (Hb needed) for the moon at > 60 deg, +0.5 at 40-60 deg, +1 inside
      40 deg; floor 8, cap 60 min; +5 min overhead. Anchor: proposal NGPS ETC, r=18.5 dark -> 9 min at S/N 10.
      Roughly: r=18 -> 4/h, r=18.5 -> 3/h (Ha) or 2/h (Hb), r=19 -> 2/h (Ha) or 1/h (Hb) under a bright moon.
      Allocation fills a time budget (4.3 h Sep 23; 9.5 h each October night) by priority per hour of telescope time;
      tier floors only take targets costing <= 30 min. Re-check with the NGPS ETC (moon phase set) before the run.
- [x] NGPS documentation read 2026-09-05 (https://caltechopticalobservatories.github.io/NGPS/). Facts that changed or
      confirmed our assumptions:
      * Coverage 3050-10400 Å in four simultaneous channels: U 3050-4430, G 4250-5960, R 5620-7950, I 7530-10400
        (was 3200-10400 in our code; corrected in 04 and the page ruler, which now shows the channels).
      * Slit: a 3-slice adjustable IFU, slices 50" long, 0.36-10" wide (0.37-10" tested); R ~ 4000-4500 at 0.4",
        R > 1500 at 1.5". The ETC page says SNR with all three slices is currently inaccurate because of "known issues
        with the current slicer hardware" and recommends the single central slice -> plan on ONE slit of the chosen
        width, not 3x capture. Target-list default slit is 1.3"; we use SET 1.3 with binning 2x3 (BINSPAT x BINSPEC,
        from the observing page's table: 2x2 for 1.0", 2x3 for 1.5" in 1-1.5" seeing; 4xN in > 2" seeing).
      * Acquisition: astrometric solve on the ACAM (4.4' x 4.15', 0.26"/px, offset ~8' from the slit, needs >= 4
        unsaturated stars; solves > 99% of fields), requirement < 120 s from slew end to exposure start, 90 s achieved
        in automatic mode. Guiding automatic; 30 s ACAM exposures OK in bright time. Our 5-min overhead per target
        stands (slew + acquisition + readout).
      * Readout (2x1 binning): U 64 s, G 31 s, R/I 55 s per exposure; cosmic-ray rates 1-3 %/hr -> "split long
        integrations into exposures of 900 s or less" (quick start). Our sub-exposure plan (>= 2 x <= 10 min) complies.
      * Calibrations: 3 ThAr + 3 FeAr arcs + 7 biases per channel and binning mode; >= 5 dome flats (7-10 for U, G) per
        slit-width + binning combination; taken in the afternoon; internal focus done by the Support Astronomer.
        Standards: at least one spectrophotometric standard near the start and one near the end of the night.
      * ETC: reports the average single-wavelength-bin SNR over a user window (matches the "per bin" S/N of our model);
        inputs are seeing at 6400 Å, V-band sky brightness (mag/arcsec^2), airmass, channel, binning, slit; no cloud
        term. Keep "No Slicer" checked. Our exposure model should be re-anchored with the ETC using the single-slice
        setting and a bright-moon V sky (~18.5-19.5 mag/arcsec^2) before the run.
      * Target list: CSV with a header (name, RA HH:MM:SS.S, DECL +DD:MM:SS, J2000) plus slitwidth ("SET 1.3" | "PSF X"
        | "SNR X"), exptime ("SET s" | "SNR X"), nexp, binspect, binspat, slitangle (deg | "PA"), airmass_max, and the
        ETC columns mag, magsystem, magfilter, channel, wrange ("a:b" Å) that let the sequencer solve exposure times;
        Note <= 24 chars, Comment <= 1024 chars. `06b_ngps_targetlist.py` writes finders/<night>/ngps_<night>_fixed.csv
        (our SET times) and ngps_<night>_snr.csv (exptime "SNR 7" with mag/channel/wrange on Hα or Hβ) for each night.
      * Data products: Quicklook DRP (MATLAB) reduces all four channels and three slices within tens of seconds of
        readout; spec2d/ multi-extension FITS and spec1d/ FITS saved from the GUI; flux calibration via CALSPEC
        sensitivity functions. The night sheet's spectrum slot currently expects a CSV; add a spec1d FITS reader once
        we have a sample file.
- [ ] Palomar: Dec > -25 comfortable.
- [ ] Which lines land in DBSP range: Hα to z~0.4; Hβ, [OIII] to z~0.8; MgII beyond.
- [ ] astroplan: hours above airmass 2 per night, moon separation, per target.
- [ ] Final list: `data/targets_<run>.ecsv` + finder charts. Keep a backup list per RA hour.

---

## 6. Master table schema

| Group | Columns |
|---|---|
| Identity | `objectid`, `labels_decoded`, `label_bits`, `source_catalog`, `sdss_objid`, `specobjid`, `desi_targetid` |
| Position | `ra`, `dec`, `ebv`, `ecl_lat` |
| Brightness | `g`, `r`, `i` (PS1/SDSS), `ztf_g_now`, `ztf_r_now`, `w1`, `w2`, `r_at_last_spec`, `w1_at_last_spec` |
| Redshift | `z`, `lines_in_range` |
| Manifold | `umap_x`, `umap_y`, `clagn_score`, `d_turnon`, `d_turnoff`, `in_region` |
| Spectra | `n_spec`, `mjd_list`, `survey_list`, `class_list`, `subclass_list`, `mjd_last`, `years_since_last_spec`, `class_change_flag` |
| Variability | `amp_w1`, `amp_ztf`, `dflux_w1_since_spec`, `dflux_ztf_since_spec`, `trend`, `expected_transition` |
| Bookkeeping | `tier`, `priority`, `exclude_flag`, `exclude_reason`, `notes` |

---

## 7. Proposed file layout

Actual layout as of 2026-09-03 (scripts, run with /opt/anaconda3/bin/python from this folder, in this order):

```
CLAGN/
  FOLLOWUP_PLAN.md                 <- this file
  wise_manifold_zeltyn.ipynb       <- original manifold + Zeltyn projection notebook
  00_literature_clagn.py           <- literature CLAGN coords (NED/VizieR/SIMBAD) -> data/literature_clagn.csv
  01_rebuild_manifold.py           <- Sample A W1 UMAP with objectids, Zeltyn projection, region stats, umap_w1_model.pkl
  02_parent_pool.py query|wise|project  <- DR16 z<0.8 r<20 quasars in the run RA windows; unWISE fetch (parallel
                                          workers: `wise 19.0 <worker> <nworkers>`); projection + scores
  03_spectra_inventory.py <csv> <tag>   <- SDSS DR19 allspec + DESI DR1 epochs per target (targets already in
                                          data/spectra_epochs_<tag>.csv are reused, only new ones are queried)
  03b_ztf_now.py <csv> <tag> ...        <- current ZTF g/r photometry and change since a reference MJD
  03d_fetch_spectra.py [targets csv]    <- downloads every SDSS spectrum epoch (sas_url from allspec: DR17 lite files for
                                          SDSS-I..IV, DR19 v6_1_3 lite incl. SDSS-V allepoch coadds) and DESI DR1 spectra via
                                          SPARCL (client >= 1.3.0 works; 1.2.5 did not); 6 Å rebinned JSON per target in
                                          data/spectra_dl/, model-free rest-frame EW(Hb, Ha) per epoch and the pipeline
                                          SPZLINE/ZLINE fits in data/spectra_lines.csv. Shown on the night sheet cards.
  03e_blending.py <csv> <tag>           <- WISE-beam blending check (user's point 2026-09-05: W1 PSF 6.1", 2.75" pixels):
                                          unWISE DR1 fracflux_w1/w2 (fraction of the W1 flux at the target that is the
                                          target) via Data Lab TAP, plus PS1 DR2 neighbours within 8" (count, separation,
                                          Δz) and a rough W1 contamination estimate (neighbours assumed 1.5 mag bluer in
                                          z-W1 than a quasar). blend_kind: "neighbour" (a PS1 source within 8" contributing
                                          > 15% of W1 or fracflux_w1 < 0.8 with a neighbour) -> priority scaled by the clean
                                          fraction (floor 0.5) unless the ZTF change confirms; "extended host" (fracflux < 0.8,
                                          no neighbour: the deblender split the host; dilution only) and "minor neighbour"
                                          (< 15%) are annotated only. 2026-09-05 result: of 1,061 candidates 75 neighbour blends,
                                          25 extended hosts, 48 minor; of the 100 primaries 8 / 6 / 4. Notable: P15022 has a
                                          star 2.5" away 2.2 mag brighter in z (W1 mostly the star's); P3478 owns 39% of its W1;
                                          P20067's neighbour 6.5" away is 0.5 mag brighter in z (its spectral change is real).
  03c_neowise_now.py <csv> <tag>        <- NEOWISE-R single-exposure W1/W2 to 2024-02 via IRSA Gator bulk upload
                                          (100 positions/request, ~35 s); per-visit medians; feeds P and the trend
  04_score_tiers.py                <- master list, priority, tiers, per-night allocation -> targets_<night>.csv
  08_time_axis.py                  <- test of the "region offset = time axis" idea (see Section 8)
  05_observability.py <in> <out>   <- astroplan hours/airmass/moon per night (run on Zeltyn and pool tables)
  06_finder_charts.py targets_<night>.csv  <- PS1 finder charts + text target list in finders/<night>/
  07b_cutouts.py                   <- 40" thumbnails per target (2026-09-05; was 64"): SDSS DR18 gri JPEG (SkyServer
                                      ImgCutout) and PS1 g, r stack cutouts (0.25"/pix, fitscut) -> data/cutouts/, embedded as
                                      data URIs. ZTF g/r reference cutouts (IRSA IBE; stored with CD1_1>0, CD2_2<0, flip both
                                      axes for N-up E-left) remain available with FETCH_ZTF=True but are too coarse at 40".
                                      The public copy of the page is docs/index.html -> https://xoubish.github.io/CLAGN/
  07_make_webpage.py               <- self-contained "CLAGN Night Sheet" (web/clagn_night_sheet.html): cards per target
                                      with light curves, manifold position, archival epochs, NGPS line coverage, and a
                                      slot for the observed spectrum (data/ngps_spectra/<name>.csv). Published as an artifact.
  run_after_wise.sh                <- runs 02 project -> subset -> 03 -> 03b -> 04 -> 06 -> 07
  data/                            <- all outputs (parquet/pkl/caches gitignored)
  from_AGNzoo/                     <- legacy notebooks, code, figures, paper-era data
```

Notebooks in `from_AGNzoo/` must be run from inside that folder (they import its `code_src/`).
The top-level `code_src/` is the newer Fornax code and is not import-compatible with them.

---

## 8. Open questions / decisions log

- 2026-09-03: Region definition. Zeltyn-defined region gave 0.8x enrichment; switch to
  Sample A Turn-on/Turn-off definition and hold Zeltyn out. **Decision pending.**
- 2026-09-03 (later): `01_rebuild_manifold.py` re-fit (unshuffled, objectids kept) and tested both definitions
  on the same 10x10 grid, 3x enrichment threshold (`data/region_stats.txt`):

  | Region | Sample A inside | known Turn-on/off inside | Zeltyn CL-AGN z<1 inside | Zeltyn EVQ z<1 inside |
  |---|---|---|---|---|
  | Zeltyn-defined (7 bins) | 55 (2.8%) | 2 / 167 | 33 / 111 (self-defined) | 12 / 68 = 18% -> **6.3x** |
  | Turn-on/off-defined (8 bins) | 174 (8.9%) | 60 / 167 (self-defined) | 12 / 111 = 11% -> **1.2x** | 17 / 68 = 25% -> 2.8x |

  kNN score (k=50, fraction of Turn-on/off neighbours): SDSS_QSO median 0.04, WISE_Variable 0.06, Turn-on 0.14,
  Turn-off 0.17, Zeltyn CL-AGN 0.10, Zeltyn EVQ 0.14. Reading: the SDSS-V CLAGNs do *not* land where the
  literature Turn-on/off objects sit (1.2x), and where the SDSS-V CLAGNs concentrate the literature CLAGNs are
  absent (2/167). The two CLAGN samples occupy different parts of the W1 manifold. EVQs follow the Zeltyn
  CL-AGNs (6.3x, a genuine held-out number since EVQs were not used to define the region). **Decision still
  pending**: which population the discovery tier should be modelled on. The proposal's 6x figure corresponds
  to the EVQ-vs-Sample-A number of the Zeltyn-defined region, not to a CL-AGN held-out test.
- 2026-09-03 (evening): DR16 parent pool projected (`02_parent_pool.py`): 25,913 z<0.8 r<20 quasars in the run RA
  windows; unWISE W1 fetched for the 10,660 with r<19, 10,575 projected. In Zeltyn-defined region 774 (7.3%, vs
  2.8% of Sample A), in Turn-on/off-defined region 398, kNN clagn_score>=0.15 1,397, Zeltyn density ratio>=3 744.
  The two region flags are disjoint (0 objects in both). M>=1: 1,941 (median z 0.42, median r 18.5);
  observable with moon sep >= 40 deg: 280 on Sep 23, ~910 on each October night -> Tier 1 supply is not the
  limiting factor; the photometric-change term P and brightness decide the ranking. Enrichment (spectra + ZTF)
  capped at the top 600 by M plus 150 controls.
- 2026-09-03 (late): Is the offset between the literature-CLAGN region and the Zeltyn region a *time* axis?
  (`08_time_axis.py`, `data/manifold_time_axis.{csv,png}`). Position along the axis between the two region centroids
  does not correlate with the epoch of the largest W1 step (Spearman rho ~ 0 for CLAGN+EVQ, n=344) nor with the share
  of the change after 2015. But the mean normalised W1 profiles differ in *shape*: Sample A Turn-off objects are
  bright in 2010-11 (1.25, 1.18 x median) and flat ~1.0 from 2013 on; Turn-on objects are faint in 2010-11 (0.83,
  0.80) and rise to 1.0 by 2014; Zeltyn CL-AGNs and EVQs decline gently and monotonically 1.06 -> 0.97 across
  2010-2020. So the two regions encode "large change at the start of the WISE window" vs "slow steady fade through
  2020", not a continuous clock. Consequence for 2026: do not extrapolate along the axis; the right signature for an
  object changing *now* is a change in the last epochs, which needs post-2020 photometry (ZTF to 2025 is in; NEOWISE
  single-exposure photometry to 2024 from IRSA would extend W1 by 3-4 yr and is the natural next enrichment).
- 2026-09-04: Multi-band manifold test (`09_multiband_test.py`, `data/multiband_test.csv`), held-out Zeltyn objects vs
  Sample A SDSS_QSO; ZTF binned to 3 d and cut at MJD 60067 to match Sample A:

  | bands | AUC CL-AGN | AUC EVQ | held-out Zeltyn-region enrichment | literature-region enrichment |
  |---|---|---|---|---|
  | W1 (DTW) | 0.656 | 0.707 | 10.6x (27% captured, 2.7% QSO inside) | 1.7x |
  | W1+W2 (DTW) | 0.627 | 0.702 | 4.2x | 1.8x |
  | ZTF g,r (manhattan) | 0.494 | 0.564 | 4.9x | 1.1x |
  | ZTF g,r + W1,W2 (manhattan) | 0.635 | 0.684 | 4.4x | 12x (11% captured, 0.9% QSO inside) |

  Reading: W1 alone is the best predictor of new CLAGNs; the optical bands add SED/redshift structure and short-timescale
  noise. The combined space's one virtue is that the literature and SDSS-V CLAGNs overlap there (12x in a tiny corner).
  Decision: keep the W1 manifold for selection; `10_combined_rescore.py` scores the enriched candidates in a ZTF g,r + W1
  space (manhattan): AUC 0.70 CL-AGN / 0.69 EVQ, held-out enrichment 7.1x; Spearman 0.18 with the W1-space M, so it is a
  near-independent second vote. Of 91 scored primaries only 4 (W1-position-only, P<1, M_combined<0.5) lacked support; 65
  non-selected T1 candidates score >= 1 in both spaces (28 with r<18.5). 2026-09-04: priority += 0.5 when M_combined >= 1,
  -= 0.3 when M_combined < 0.5 and P < 1 (T1 only); M_combined shown on the cards. Full-pool reselection in a combined
  space is not feasible before the runs (ZTF light curves for 16,500 more quasars).
- 2026-09-03: Step zero. `~/Dropbox/sample.ecsv` (2042 rows) is **not** the Sample A table: labels disagree
  for objectid >= 1000 and a light-curve fingerprint test (unWISE W1/W2 fetched at its coordinates vs the
  parquet light curve of the same objectid) shows different epochs/fluxes for objectids 0, 63, 153, 1000.
  Parquet objectid order is the build_sample.md order (SDSS_QSO 0-999, SPIDER 1000-7559, WISE_Variable
  7560-7939, Optical 7940-8061, Galex 8062-8106, Turn-on/off 8107-8273, TDE 8274-8306), so deterministic blocks
  could be rebuilt and fingerprint-verified if the original CSV never turns up.
- 2026-09-03: Archival spectra sources settled (`03_spectra_inventory.py`): SDSS DR19 `allspec` holds every
  SDSS-I..V spectrum incl. SDSS-V BOSS daily epochs to ~MJD 60000 (plate era 2020-21 and FPS era 2022);
  class/subclass/z come from `SpecObjAll` (SDSS-I..IV) and `mos_sdssv_boss_spall` (plate-era SDSS-V only; the
  FPS-era epochs have no class in CAS, would need spAll-lite v6_1_3 from SAS). DESI DR1 via NOIRLab Data Lab
  TAP `desi_dr1.zpix` (box query on mean_fiber_ra/dec; ADQL CONTAINS and q3c are not accepted; SPARCL client
  1.2.5 endpoint is dead). SkyServer REST rejects UNION ALL, so one cone query per target (~2 s each).
  Zeltyn test: 204/204 have SDSS spectra (median 8 epochs), 119 have DESI DR1, last spectrum 3.2-5.7 yr ago
  for all of them, 23 have a pipeline class change across epochs.
- W1 only vs W1+W2 vs ZTF+WISE manifold for scoring. Paper: W1 cleanest for CLAGN separation.
- Whether to build the bigger DR16Q parent pool (Section 2) or stay with Sample A.
- Priority score formula (Section 4).
