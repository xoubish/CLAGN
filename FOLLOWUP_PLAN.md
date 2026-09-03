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
- [ ] Literature CLAGNs from `ChangeLook_Sample_oct2023.ipynb` (LaMassa15, MacLeod16/19, Ruan16,
      Yang18, Sheng20, Green22) — these overlap with the Turn-on/Turn-off labels; de-duplicate at 1".
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

- [ ] Priority score = f(clagn_score, change_since_last_spec, years_since_last_spec, r mag, z).
      Write down the formula in this file when chosen.
- [ ] Write `data/master_list_scored.ecsv`.

---

## 5. Stage 4: observability (last)

- [ ] Nights / semester: ______  (RA window follows from this)
- [ ] Palomar: Dec > -25 comfortable. DBSP: r < 19.5-20 for ~30-60 min exposures.
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

```
CLAGN/
  FOLLOWUP_PLAN.md                 <- this file
  wise_manifold_zeltyn.ipynb       <- manifold + Zeltyn projection (exists)
  01_recover_coords.ipynb          <- step zero
  02_build_master_list.ipynb       <- stage 1
  03_enrich_spectra_phot.ipynb     <- stage 2 (archives, photometry)
  04_score_tiers.ipynb             <- stage 3
  05_observability.ipynb           <- stage 4
  data/
    sampleA_coords.ecsv
    master_list.ecsv
    master_list_scored.ecsv
    targets_<run>.ecsv
  from_AGNzoo/                     <- legacy notebooks, code, figures, paper-era data
```

Notebooks in `from_AGNzoo/` must be run from inside that folder (they import its `code_src/`).
The top-level `code_src/` is the newer Fornax code and is not import-compatible with them.

---

## 8. Open questions / decisions log

- 2026-09-03: Region definition. Zeltyn-defined region gave 0.8x enrichment; switch to
  Sample A Turn-on/Turn-off definition and hold Zeltyn out. **Decision pending.**
- W1 only vs W1+W2 vs ZTF+WISE manifold for scoring. Paper: W1 cleanest for CLAGN separation.
- Whether to build the bigger DR16Q parent pool (Section 2) or stay with Sample A.
- Priority score formula (Section 4).
