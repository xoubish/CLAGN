# JWST Cycle 6 proposal

The current draft is **[proposal.pdf](proposal.pdf)**;
edit **[proposal.tex](proposal.tex)**.
It requests 67.0 hours for 24 AGN (12 rising and 12 fading selection histories),
with one MIRI/MRS science visit and a linked sky visit per target.
The attachment has five core pages and seven pages overall.

## Working files

| File | Purpose |
| --- | --- |
| `proposal.tex` / `proposal.pdf` | Current source and compiled proposal |
| `fig1_targets.pdf` | Three measured histories: F06, F04 and R01 |
| `fig2_miri.pdf` | What MIRI adds: spectral memory and dust-boundary evolution |
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

The audit checks the current PDF and its build log. The official style
file is unchanged. For Overleaf, upload `proposal.tex`,
`jwstproposaltemplate_v6.sty`, and the two current figure PDFs above.

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

Figure 2 follows the analysis of the proposed MRS data: separate nucleus and
host; measure the continuum, both silicate profiles, PAHs and gas lines; connect
these observables to each target's history; then test spectral memory and the
need for dust evolution. Supporting branches identify the obscuration,
excitation and host science. The extraction profile is a labelled schematic,
and the nuclear spectrum is an illustration from the existing shell model.
No forecast detections or model-separation amplitudes are plotted.

The generator is `make_fig2_miri.py`. Its current numerical illustration and
provenance are `inputs/fig2_miri_workflow_spectrum.csv` and
`inputs/fig2_miri_provenance.json`. Earlier response-curve and spectrum CSVs
remain as numerical records of the superseded figure, with copies in the archive.
[review/mrs_measurement_plan.md](review/mrs_measurement_plan.md) maps each
measurement to its comparison, scientific question and interpretation limits.

The conditional, history-driven forecast remains in the proposal text and saved
inputs. It uses monochromatic observed-frame 12-micron flux as a W3 proxy;
the science analysis will compare bandpass-integrated models and the rest-frame
8–13-micron spectrum. `review/reframed_revision.md` documents its assumptions;
`review/combined_figure_revision.md` records the current presentation changes.

The forecast calculations remain reproducible with `forecast_warm_response.py`
and `forecast_scenarios.py`. Their saved numerical results were not refitted in
this editorial revision. Shared shell and fitting scripts remain here because
the forecast and technical audits import them.

## Earlier work

The `older` folder preserves both incoming drafts and their figure dependencies,
plus superseded single-target and boundary-evolution figures. Earlier project
archives remain under [../archive/](../archive/README.md). The September 26
technical review and September 27 history audit retain their dated scope.
