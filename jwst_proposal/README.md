# JWST Cycle 6 proposal

The current draft is **[proposal.tex](proposal.tex)** for editing in Overleaf.
Rebuild `proposal.pdf` after every proposal edit, as requested by the user on
September 30. Run `make` from this folder. Figure PDFs remain required Overleaf assets.
It requests 67.0 hours for 24 AGN with rises, declines, flares and recoveries,
with one MIRI/MRS science visit and a linked sky visit per target.
A September 30 no-PDF LaTeX check retained five core pages after the co-I
clarity edits and two-row Figure 2; subsequent edits now receive a full PDF build.

## Working files

| File | Purpose |
| --- | --- |
| `proposal.tex` | Current source for Overleaf |
| `fig1_targets.pdf` | Three measured histories: F06, F04 and R01 |
| `fig2_miri.pdf` | Two-row literature figure: spectrum and MRS channels above spatial separation and delayed response |
| [review/fig2_three_panel/](review/fig2_three_panel/README.md) | Adopted Figure 2 preview, caption, sources and reproduction instructions |
| [inputs/jwst_sample_cycle6.csv](inputs/jwst_sample_cycle6.csv) | Authoritative sample |
| [apt/clagn24_miri_draft.aptx](apt/clagn24_miri_draft.aptx) | APT project; 48 visits |
| [apt_title_abstract.txt](apt_title_abstract.txt) | Revised title/abstract draft for transfer to APT; saved APT text not yet synchronized |
| [review/README.md](review/README.md) | Technical checks and supporting analyses |
| [validation.json](validation.json) | Current build status and scientific validation scope |
| [older/README.md](older/README.md) | Superseded drafts, figures and revision snapshots |

## Figure regeneration and Overleaf

From this folder:

```sh
python make_fig2_miri.py
```

This regenerates the adopted Figure 2 without compiling the proposal.
The official style file is unchanged. For Overleaf, upload `proposal.tex`,
`jwstproposaltemplate_v6.sty`, and the two current figure PDFs above.

To regenerate both figures from saved inputs, install `requirements.txt` and run
`make figures PYTHON=/path/to/python`. The downloaded source-paper PDFs must
already be present for Figure 2. Check page layout in Overleaf; the local
`review/audit_consistency.py` PDF audit requires a compiled PDF and build log
and has not been rerun for this source-only revision.
The Figure 1 generator retains its historical filename `make_fig1_five.py` but
now writes the three-target `fig1_targets.pdf` directly.

## Scientific inputs

Figure 1 combines SDSS images, ZTF/alert and WISE histories, archival SDSS spectra,
NGPS P330E spectra and SPHEREx for F06, F04 and R01. Its display inventory is in
`inputs/fig1_targets_provenance.json`. The SPHEREx inputs and fits are in
`inputs/spherex_extracted/`, `inputs/spherex_fits/` and the supplied R01 CSV.
All previously selected optical epochs and SPHEREx points remain plotted.

Figure 2 presents three published examples demonstrating the methods and physical
basis of the experiment: a digitized spectrum and decomposition redrawn at full
width (NGC 7469), central/surrounding
separation (NGC 6552), and a delayed infrared response (NGC 4151). The accepted
PDF is installed as `fig2_miri.pdf`, and its caption is in `proposal.tex`.
The examples do not predict component fractions, resolved scales or lags for
our targets. See [the figure documentation](review/fig2_three_panel/README.md).

`make_fig2_miri.py` regenerates and installs the adopted figure; `make figures`
uses this entry point and no longer restores the superseded mock figure.
The source paper PDFs must be available under `review/miri_literature/`.
The previous mock figure, generator, provenance and proposal source are archived
in `older/before_literature_figure2_20260929.zip`. Its numerical checks remain in
[review/figure2_recovery/](review/figure2_recovery/README.md) as support for the
conditional precision estimates retained in the text.

The old conditional ensemble forecast remains in the saved inputs for audit,
but its 4.5σ/3.7σ estimates were removed from the proposal because they do not
establish model discrimination. Full-sample physical-model recovery remains
unvalidated. The source now distinguishes SPHEREx measurements for 23 targets
from the remaining crowded target's WISE constraints. Figure 1 illustrates the
sample's data types without claiming identical complete data for every target.
Public SPHEREx retrieval is complete and preserved in the Git-ignored
`local_data/spherex/sample/` cache. See [retrieval and restart instructions](review/SPHEREX.md).
Downloaded diagnostic FITS are in `local_data/spherex/rejection_checks/`.

The forecast calculations remain reproducible with `forecast_warm_response.py`
and `forecast_scenarios.py`. Their saved numerical results were not refitted in
this editorial revision. Shared shell and fitting scripts remain here because
the forecast and technical audits import them.

## Earlier work

Recent rejected Figure 2 layouts and their generators are archived in the
local-only `older/figure2_trials_20260929.zip`, with checksums and original paths
in [the cleanup manifest](review/cleanup_20260929.json). The accepted three-panel
figure is in [review/fig2_three_panel/](review/fig2_three_panel/README.md).

The `older` folder preserves both incoming drafts and their figure dependencies,
plus superseded single-target and boundary-evolution figures. Earlier project
archives remain under [../archive/](../archive/README.md). The September 26
technical review and September 27 history audit retain their dated scope.
