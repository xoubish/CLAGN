# Palomar CLAGN observing project

Start with **[Observing files](observing/README.md)** for the September 23 sequence, NGPS CSVs, candidate pages, and visibility plots for all three nights.

| Folder | Contents |
| --- | --- |
| [jwst_proposal/](jwst_proposal/README.md) | Editable JWST Cycle 6 proposal, compiled PDF, figures, and portable figure inputs. |
| [observing/](observing/README.md) | Current September primary/backup packet and links to the three visibility plots. Local observing products are git-ignored. |
| [scripts/](scripts/README.md) | Selection, acquisition, geometry, sensitivity, and page-building code. |
| [data/](data/README.md) | Catalogues, manifold model, spectra, light curves, imaging caches, and private working tables. |
| [docs/](docs/README.md) | Public website served by GitHub Pages. |
| [web/](web/README.md) | Page templates, shared spectral plotting code, and generated public page copies. |
| [reference/](reference/README.md) | Manifold paper, data-access document, and links to cached NGPS documentation. |
| [research/](research/README.md) | Original AGNzoo work, science notebooks, and background notes. |
| code_src/ | Shared light-curve/manifold library; retained at this path for the saved model and existing imports. |
| requirements/ | Original dependency lists for the research workflows. |
| tests/ | Observability regression checks. |
| [archive/](archive/README.md) | Superseded lists, pages, reports, drivers, and inactive logs; retained locally. |

Python entry points now live in `scripts/`. From the project directory, for example:

```sh
/opt/anaconda3/bin/python scripts/16_candidate_webpage.py
/opt/anaconda3/bin/python -m unittest discover -s tests -v
```

The page command rebuilds the current candidate explorer from cached tables. It does not commit, push, or publish. Rebuilding the September packet is a separate action; see the script guide.

The full bright-parent W1 acquisition and continuation pipeline were still running during the September 20 cleanup. Check the [live pipeline status](data/reselection_2026-09-20/completion_pipeline_status.json); the curated September sequence remains separate from that automatic rebuild.

Folder cleanup did not change the source selection, exposure settings, or observing order. Collaboration data stay in ignored local directories. See the archive guide for the move manifest and pre-cleanup copies.
