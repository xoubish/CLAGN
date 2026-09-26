# JWST Cycle 6 proposal

The current source is [proposal.tex](proposal.tex); [proposal.pdf](proposal.pdf)
is its compiled draft. The working design is **24 AGN: 12 fading and 12 rising,
with one MIRI/MRS visit per target**. SPHEREx supplies the near-infrared spectra.
Recovery and reversal are retained in individual histories, not a third group.

## Current sample and status

[inputs/jwst_sample_cycle6.csv](inputs/jwst_sample_cycle6.csv) is the authoritative
working target list. Source identities, memberships and its existing estimates
have not been changed by the document cleanup. The user confirms SPHEREx data
and ZTF/WISE manifold coverage for all targets. The table records five completed
September NGPS observations and 19 scheduled for October 26–27, 2026.
Individual September exposure UTC dates can fall on September 24.

The intended category is Medium: more than 50 and at most 130 charged hours,
including overheads, under the [Cycle 6 rules](https://jwst-docs.stsci.edu/jwst-opportunities-and-policies/jwst-call-for-proposals-for-cycle-6/jwst-proposal-types-and-categories/jwst-general-observer-go-proposals).
No total time is yet validated. Target-specific ETC/APT calculations must determine it.

The next review must verify individual rising/fading assignments, event intervals,
host contributions and lag estimates. Three entries carry reversal flags. The
existing table assumes a warm-dust delay ten times the hot-dust delay throughout;
this is not a validated measurement. Figure 2, the target/epoch tables, model
separation, statistical power, visibility and duplication checks remain unfinished.
Red text marks unresolved quantities and the model figure. Procedural checks
are collected in [completion_notes.md](completion_notes.md).

The rebuilt draft has six pages: scientific justification and observations on
pages 1–4, supplemental information on page 5, and references on pages 5–6.
The required sections currently fit within the five-page Medium limit, leaving
space for the target table. Recheck the limit after completing figures and tables.

## Build

```sh
make
```

To refresh the sample summary and regenerate Figure 1, install `requirements.txt`
and run `make figures PYTHON=/path/to/python`, followed by `make`.
`make_selection.py` summarizes the authoritative CSV and checks its counts and
unique identifiers; it does not reselect targets or validate physical classifications.
`make clean` removes compilation files while retaining PDFs.
For Overleaf, upload this folder and select `proposal.tex` as the main document.

## Files

| File | Purpose |
| --- | --- |
| `proposal.tex` / `proposal.pdf` | Current text and compiled draft |
| `apt_title_abstract.txt` | Matching title and abstract for APT |
| `jwstproposaltemplate_v6.sty` | Official style, unchanged |
| `fig1_connected.pdf` / `.png` | P1823 measured light curves and spectra |
| `make_fig1_connected.py` | Figure generator |
| `make_selection.py` | Summarize the agreed 24-source CSV |
| `spherex_data.py` | Figure SPHEREx input loading and epoch grouping |
| `inputs/jwst_sample_cycle6.csv` | Authoritative working sample |
| `inputs/selection_summary.json` | Counts, observing status and pending audits |
| `inputs/provenance.json` | Input origins and checksums |
| `validation.json` | Current document build checks and remaining work |
| `completion_notes.md` | Outstanding checks kept outside the proposal narrative |

The superseded selection/timing scripts, sample tables and summaries are preserved
in [the project archive](../archive/jwst_before_24_target_miri/README.md).
Older proposal drafts under `observing/` are historical and are not build inputs.

## Figure 1

Figure 1 contains P1823 only. It preserves the P330E calibration and the previously
requested exclusion of one SPHEREx point; the input CSV is intact. Its optical
spectrum spans 0–0.6 mJy. P1823's 18-micron peak lies beyond MRS; its test uses
the warm continuum and accessible 9.7-micron feature. The figure does not itself
establish significant SPHEREx variability or separation of competing dust models.
Only P1823's figure data are bundled here; availability across the full sample is
confirmed by the user and does not imply that all spectra are in this folder.

Gold dashed lines mark expected redshifted emission-line positions, not fitted
detections. Vacuum wavelengths follow the [SDSS optical line table](https://classic.sdss.org/dr6/algorithms/linestable.php)
and [STScI infrared line lists](https://www.stsci.edu/instruments/nicmos/documents/handbooks/instrument/v5/Appendix_26.html).
The figure provenance records every plotted line position; a few labels are
offset horizontally to keep neighboring lines readable.

Date keys above the optical and SPHEREx panels use the same colors as their
spectra and the time markers on both light curves. Archival optical spectra
are grouped by year in the key; the latest spectrum is labelled 2026.
NGPS observing status, months and program identifiers are omitted from the
proposal and its target table, following the user's presentation preference.
Each of the 20 optical epochs since 2015 has a dashed marker at its saved MJD.
The labelled 2002 spectrum is outside the displayed window. SPHEREx keys show
each visit's month range; dotted lines mark median observation dates, with
shading covering the full visit. Exact dates and colors are recorded in
`inputs/fig1_connected_provenance.json`.

Both ZTF and WISE panels start at 2015. W1 is shown as purple circles and W2 as
brown squares, on the same mJy scale without offsets. The original W1 data are
retained (filled circles: saved unWISE series; open circles: NEOWISE visits).
W2 uses the bundled target-only NEOWISE exposure snapshot in
`inputs/figure_neowise_exposures.csv`. The figure generator applies the existing
W1 frame-quality and six-month grouping rules, additionally requiring finite,
positive W2 uncertainties, and keeps visits with at least three measurements.
Visit magnitudes are medians; their approximate statistical errors are
`1.2533 * sample_std / sqrt(n)`. Conversion to mJy uses a W2 Vega zero point of
171.787 Jy from the [WISE calibration table](https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/sec4_4h.html).
These plotted errors exclude absolute calibration uncertainty; no color
correction is applied. Figure provenance records all W2 visits and quality cuts.
Earlier measurements remain in the inputs, outside the displayed time window.
