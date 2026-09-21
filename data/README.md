# Data and caches

These files are retained in place because the selection pipeline and saved manifold depend on them. The cleanup moved inactive top-level logs into the archive; it did not discard spectra, photometry, catalogue inputs, or fitted models.

| Location | Use |
| --- | --- |
| `reselection_2026-09-20/` | Current private parent search, prepared candidates, geometry, sensitivity plans, acquisition audits, and local explorer. |
| `reselection_2026-09-20/completion_pipeline_status.json` | Live continuation-job status. |
| `reselection_2026-09-20/completion_pipeline.log` | Continuation log, including the restart after folder cleanup. |
| `reselection_2026-09-20/three_night_wise_focused.log` | Active bright-parent W1 acquisition log. |
| `reselection_2026-09-20/three_night_review/` | Science/sensitivity tables and current visibility plots. |
| `reselection_2026-09-20/ngps_etc/` | Cached NGPS exposure-time calculator used by the sensitivity scripts. |
| `spectra_cache/`, `spectra_dl/` | Downloaded public spectral files and parsed spectra. |
| `sdssv_internal/`, `reselection_2026-09-20/sdssv_spectra/` | Collaboration spectral data. |
| `wise_cache/`, `neowise_cache/`, `ztf_cache/`, `ztf_objects_cache/` | Reusable light-curve and query caches. |
| `cutouts/` | SDSS imaging used by the current pages. |
| `external/` | External catalogue inputs. |
| `ngps_spectra/` | Local NGPS follow-up spectra. |
| `umap_w1_model.pkl`, `sampleA_embedding*.csv`, `zeltyn_embedding*.csv` | Saved manifold model and reference embeddings. |

The historical root-level `targets_*.csv` and `schedule_*.csv` remain as inputs to earlier analyses and regression checks. **The current September telescope CSVs are in [observing/sep23/](../observing/sep23/), not these historical tables.**

Private derived tables and large caches are covered by `.gitignore`. Keep those rules when adding or moving data. New code should use the project root explicitly rather than infer a data directory from the shell's working directory.
