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
      **priority = B * S * (M + P) + 0.5 * class_change_flag**
      Tiers: T3 Zeltyn CL-AGN; T2 Zeltyn EVQ with M >= 1 or in Zeltyn region; T1 pool object with M >= 1; T4 control
      (pool, clagn_score = 0, zeltyn_density_ratio <= 0.3, outside both regions). Allocation per night: hrs >= 1.5,
      moon separation >= 40 deg, r <= 19.5, z <= 0.8 (T3 exempt); slots 14 / 34 / 34; caps T3 <= 3, T4 <= 4 per night.
- [x] Write `data/master_list_scored.csv` and `data/targets_<night>.csv` (top list + backups).

---

## 5. Stage 4: observability (last)

- [x] Nights: **P200/NGPS, 2026-09-23 first half** (astro. twilight 20:05, window to 00:39 PDT, LST 19.5h-0.1h,
      moon 93% at RA 22.2h Dec -12, up all window) and **2026-10-26 + 10-27 full nights** (19:24-05:38 PST,
      LST 21h-7.2h; moon 98%/94% at RA 3-4h Dec +22..+26, rises ~3.5 h after twilight so the first third of
      each October night is moon-free). All bright time -> bright targets and/or > 40 deg from the moon.
      Computed by `05_observability.py` (astroplan); it adds hrs_<night>, minX_<night>, moonsep_<night> columns.
- [x] RA windows used for the parent pool (`02_parent_pool.py`): 15.5h-24h, 0h-4.5h, 6.5h-11h; Dec > -15.
- [ ] Palomar: Dec > -25 comfortable. NGPS: r < 19.5-20 for ~15-30 min exposures (proposal ETC: r=18.5 -> 8-10 min).
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
  03_spectra_inventory.py <csv> <tag>   <- SDSS DR19 allspec + DESI DR1 epochs per target
  03b_ztf_now.py <csv> <tag> ...        <- current ZTF g/r photometry and change since a reference MJD
  04_score_tiers.py                <- master list, priority, tiers, per-night allocation -> targets_<night>.csv
  05_observability.py <in> <out>   <- astroplan hours/airmass/moon per night (run on Zeltyn and pool tables)
  06_finder_charts.py targets_<night>.csv  <- PS1 finder charts + text target list in finders/<night>/
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
