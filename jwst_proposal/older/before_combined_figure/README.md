# JWST Cycle 6 proposal

The current draft is **[proposal_reframed.pdf](proposal_reframed.pdf)**;
edit **[proposal_reframed.tex](proposal_reframed.tex)**.
It requests 67.0 hours for 24 AGN (12 rising and 12 fading selection histories),
with one MIRI/MRS science visit and a linked sky visit per target.
The attachment has five core pages and seven pages overall.

## Working files

| File | Purpose |
| --- | --- |
| `proposal_reframed.tex` / `proposal_reframed.pdf` | Current source and compiled proposal |
| `fig1_targets.pdf` | Three measured histories: F06, F04 and R01 |
| `fig2_diagnostics.pdf` | Delayed-heating physics and spatial decomposition |
| `fig2_forecast.pdf` | Conditional forecasts tied to the recorded histories |
| [inputs/jwst_sample_cycle6.csv](inputs/jwst_sample_cycle6.csv) | Authoritative sample |
| [apt/clagn24_miri_draft.aptx](apt/clagn24_miri_draft.aptx) | APT project; 48 visits |
| [apt_title_abstract.txt](apt_title_abstract.txt) | Title and abstract matching saved APT |
| [review/README.md](review/README.md) | Technical checks and supporting analyses |
| [validation.json](validation.json) | Current build status and scientific validation scope |
| [older/README.md](older/README.md) | Superseded drafts, figures and revision snapshots |

## Build

From this folder:

```sh
make
python review/audit_consistency.py
make clean
```

The audit checks the current reframed PDF and its build log. The official style
file is unchanged. For Overleaf, upload `proposal_reframed.tex`,
`jwstproposaltemplate_v6.sty`, and the three current figure PDFs above.

To regenerate figures from the saved inputs, install `requirements.txt` and run
`make figures PYTHON=/path/to/python`, then `make`. No downloads are needed.
The Figure 1 generator retains its historical filename `make_fig1_five.py` but
now writes the three-target `fig1_targets.pdf` directly.

## Scientific inputs

Figure 1 combines SDSS images, ZTF/alert and WISE histories, archival SDSS spectra,
NGPS P330E spectra and SPHEREx for F06, F04 and R01. Its display inventory is in
`inputs/fig1_targets_provenance.json`. The SPHEREx inputs and fits are in
`inputs/spherex_extracted/`, `inputs/spherex_fits/` and the supplied R01 CSV.
All previously selected optical epochs and SPHEREx points remain plotted.

Figure 2 uses fixed optically thin shells under step illumination. Figure 3 uses
the same grain physics with the sample W1 histories and assumed host fractions,
size scales and extrapolations. These are conditional planning calculations,
not validated population-level detections. The forecast uses monochromatic
observed-frame 12-micron flux as a W3 proxy; the science analysis will compare
bandpass-integrated models and use the rest-frame 8–13-micron spectrum to test
structural evolution. `review/reframed_revision.md` documents the revisions,
error assumptions and figure consistency checks.

The forecast calculations remain reproducible with `forecast_warm_response.py`
and `forecast_scenarios.py`. Their saved numerical results were not refitted in
this editorial revision. Shared shell and fitting scripts remain here because
the forecast and technical audits import them.

## Earlier work

The `older` folder preserves both incoming drafts and their figure dependencies,
plus superseded single-target and boundary-evolution figures. Earlier project
archives remain under [../archive/](../archive/README.md). The September 26
technical review and September 27 history audit retain their dated scope.
