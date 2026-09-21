# Files copied from the AGNzoo repo (xoubish/AGNzoo, fork of fornax-navo demo notebooks)

Copied 2026-09-03 so that CLAGN can be the single working folder for the
spectroscopic follow-up target selection. Nothing here has been modified.

## Provenance

| File / dir | AGNzoo branch @ commit | Date | Notes |
|---|---|---|---|
| `../Hemmati_2026_ApJ_998_130.pdf` | working tree (untracked) | 2026-02-05 | Published paper. `*.pdf` is gitignored in CLAGN. |
| `build_sample.md` | main @ a6eb47f | 2024-02-26 edit | **Sample A recipe** (SDSS QSO, WISE/optical/Galex variables, Turn-on/off, SPIDERS, TDE, BOSS SF). Writes `data/AGNsample_26Feb24.csv` and `data/sample.ecsv` (both gitignored, never committed). |
| `build_sample_dualAGN_mar2025.md` | paper_sample_plots @ 5ddcac5 | 2025-03 | Same notebook repurposed for a dual-AGN sample (Charisi16, Chen20, Graham15, Liu19, Ward22, bigMAC, Rodriguez06). Not Sample A. |
| `ML_AGNzoo_main_apr2024.md` | main @ a6eb47f | 2024-04-14 | Manifold notebook, Sample A era. Loads `GP_ZTFWISE.npz`. |
| `ML_AGNzoo_paper_sample_plots_mar2025.md` | paper_sample_plots @ 5ddcac5 | 2025-03-10 | Latest version; made the `figures/*.png` plots. |
| `ML_AGNzoo-Copy1_mar2025.md` | paper_sample_plots @ 5ddcac5 | 2025-03 | Scratch copy. |
| `ML_AGNzoo-kauffmann.md`, `BPT.md` | main @ a6eb47f | 2024 | Sample B (Kauffmann narrow-line AGN, ~65k) and BPT / BOSS SF notebooks. |
| `light_curve_generator_feb2024.md` | main @ a6eb47f | 2024-02 | Version of the collector that produced `df_lc_020724.parquet.gzip`. Writes `output/<name>_sample.ecsv` with coordinates. |
| `GP_ZTFWISE.npz` | main @ a6eb47f | 2024 | `data`: (1544, 400) GP-interpolated ZTF+WISE features; `labc`: dict label -> row indices. No object ids or coordinates. |
| `requirements.txt` | main @ a6eb47f | | For the notebooks above (needs `ML_utils`, `sompy`, ...). |
| `code_src/` | main @ a6eb47f | 2024 | Feb-2024 Fornax code that these notebooks import (`ML_utils.py`, `sample_lc.py`, ...). **Do not** use the paper_sample_plots version of `ML_utils.py`: its `translate_bitwise_sum_to_labels` overrides the label list with dual-AGN labels. |
| `figures/` | paper_sample_plots @ 5ddcac5 (+ `Kauffmann.png` from main) | 2025-03 | PNG plots. |
| `ChangeLook_Sample_oct2023.ipynb` | shoobyFeb @ 0344fe8 | 2023-10-04 edit | Literature CLAGN compilation with RA/Dec, z, turn-on/off type (LaMassa15, MacLeod16/19, Ruan16, Yang18, Sheng20, Green22). Reads local `yangg.cat` / `table2.cat` that were never committed. |
| `spectroscopy/` | spectroscopy_faisst_jan1624 @ aeb47ab | 2024-02-06 | A. Faisst's `spectra_generator.md` + `code_src/` for pulling archival spectra (SDSS etc.). |

## What is NOT anywhere in either repo

- A Sample A table with coordinates. `build_sample.md` and `light_curve_generator_feb2024.md`
  write it to `data/` or `output/`, both gitignored. Look on Fornax / the machine used in Feb 2024,
  or re-run `build_sample.md`.
- The light-curve parquet `df_lc_020724.parquet.gzip` (Google Drive id `1gb2vWn0V2unstElGTTrHIIWIftHbXJvz`,
  downloaded by `../wise_manifold_zeltyn.ipynb` via `gdown`). It has `objectid` but no positions.

## Running the old notebooks

They do `sys.path.append('code_src/')`, so run them with this directory as the working directory.
The current CLAGN `../code_src/` is the newer Fornax code (renamed modules, `AGNzoo_functions.py`)
and is not import-compatible with these notebooks.
