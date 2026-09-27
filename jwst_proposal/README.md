# JWST Cycle 6 proposal

Open **[proposal.pdf](proposal.pdf)** to read the current draft or edit
**[proposal.tex](proposal.tex)**. The programme requests **67.0 hours for 24 AGN**
(12 rising and 12 fading), with one MIRI/MRS science visit and a linked sky visit
per target. The PDF has five core pages and seven pages overall.

## Working files

| Location | Contents |
| --- | --- |
| [proposal.tex](proposal.tex) / [proposal.pdf](proposal.pdf) | Current source and compiled proposal |
| [apt/clagn24_miri_draft.aptx](apt/clagn24_miri_draft.aptx) | Current APT project; 48 visits |
| [apt_title_abstract.txt](apt_title_abstract.txt) | Title and abstract matching APT |
| [inputs/jwst_sample_cycle6.csv](inputs/jwst_sample_cycle6.csv) | Authoritative sample |
| `fig1_connected.pdf`, `fig2_diagnostics.pdf`, `fig3_dust_evolution.pdf` | Figures used by LaTeX |
| [inputs/](inputs/) | Saved figure data, models and provenance |
| [review/README.md](review/README.md) | Technical checks and supporting analyses |
| [apt/README.md](apt/README.md) | APT setup, accepted sky fields and timing reports |
| [etc/LOCAL_SETUP.md](etc/LOCAL_SETUP.md) | ETC calculations and local runtime |
| [duplication/README.md](duplication/README.md) | Duplication search |
| [validation.json](validation.json) | Build status and outstanding scientific/technical checks |

## Build

From this folder:

```sh
make
python review/audit_consistency.py
make clean
```

`make clean` removes temporary LaTeX files and retains the PDFs. Run the audit
before cleaning, since it checks the build log. The official style file remains
unchanged.

To regenerate the sample summary and all three figures, install
`requirements.txt`, then run `make figures PYTHON=/path/to/python` and `make`.
The generators save only the figure PDFs needed by the proposal.
`make_selection.py` summarizes the agreed CSV; it does not reselect targets or
validate their physical classifications.

For a compact Overleaf upload, use `proposal.tex`, `jwstproposaltemplate_v6.sty`
and the three figure PDFs. Select `proposal.tex` as the main document.

## Scientific inputs

Figure 1 uses the P9694 light curves, optical spectra and all 289 supplied
SPHEREx measurements. The dashed curve is the pooled host + disc + blackbody
fit, with a model-dependent colour temperature near 1200 K. Its inputs and
masks are recorded in `inputs/fig1_connected_provenance.json` and
[the fit report](review/p9694_hot_dust/README.md).

Figure 2 illustrates history dependence relative to equilibrium at current
nuclear power. Figure 3 compares delayed fixed dust with an evolving boundary
under the same illumination history. Both are illustrations, not validated
sample-level detection forecasts. The exploratory
[P9694 predictivity analysis](review/p9694_predictivity/README.md) and
[validation status](validation.json) retain the limits of the current checks.

## Archived material

Old proposal and figure versions, APT snapshots, inactive logs and editing notes
are consolidated in [one ZIP archive](../archive/jwst_cleanup_20260926.zip),
along with the duplicate PNG previews. Earlier project
archives remain under [archive/](../archive/README.md).
