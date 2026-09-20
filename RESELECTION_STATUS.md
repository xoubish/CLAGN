# Manifold-first reselection

Started 2026-09-20 following the Moon-coordinate audit. The intended order is:

1. Reopen the broad parent, preserving the paper's W1 manifold as the reference.
2. Calculate Palomar observing windows using apparent Moon geometry and retain brightness/airmass alternatives.
3. Complete W1 coverage for newly feasible sources.
4. Review dated spectra, including available SDSS-V, to decide scientific weights.
5. Produce a new observing sequence and matching webpage.

The old scheduled lists and currently published webpage contain the previously identified Moon-distance errors. They are not the result of this rebuild. A replacement candidate-review page has now been built and checked locally; publication is awaiting explicit public-export approval after automatic review rejected the push.

## Candidate webpage update

Local commit `23fe570` replaces `docs/index.html` and `web/clagn_night_sheet.html` with the provisional 301-object review snapshot. It includes night/region/history filters, corrected continuous observing windows, a clickable manifold, cached light curves, public spectral overlays, and CSV export. The public payload has 104 objects with repeat spectral dates; the ignored [local collaboration copy](data/reselection_2026-09-20/candidate_review_local.html) includes the fuller 175-object repeat-coverage count and SDSS-V metadata. Internal metadata is excluded from the public payload. The publication destination is `https://xoubish.github.io/CLAGN/`, repository `xoubish/CLAGN`, branch `main`.

All 301 cards, night counts, filters, CSV export, and desktop/mobile layouts passed browser checks. All three Moon/observability regression tests passed. Build with `16_candidate_webpage.py` after refreshing a consistent compact pool and spectral audit. The builder checks publication provenance and stops if future candidates come from additional parent origins that have not yet been reviewed for public export. The background projection pipeline does not automatically publish or change this webpage snapshot.

The parent and first complete geometry pass are available in [the pool report](data/reselection_2026-09-20/POOL_REBUILD.md). Scientific weights remain unset. [The current manifold review queue](data/reselection_2026-09-20/manifold_observable_review_queue.csv) contains all currently projected, observable objects with historical r≤19.5, without an optical-trigger requirement.

Following the request to reduce the working pool, use the [compact review](data/reselection_2026-09-20/COMPACT_REVIEW.md) and its [one-row-per-source table](data/reselection_2026-09-20/compact_review_objects.csv) for the provisional first spectral pass: either existing candidate manifold region, historical r≤18.5, and ≥60 continuous minutes simultaneously at X≤1.5 and Moon separation ≥60°. This is an operational subset, with both regions separately labeled and no finalized science ranking. Its initial 301 objects come from currently available projections; the report automatically reapplies the cuts after expanded projections finish. The broad parent remains a search catalog, not the working observing list.

The unWISE expansion is a substantial batch: 92,757 newly feasible objects with estimated r≤20.5 across 3,596 sky partitions. Acquisition is resumable. Fainter and unknown-brightness objects remain explicitly in the parent. The batch is not declared complete until acquisition and projection finish.

- [Automatic pipeline status](data/reselection_2026-09-20/pipeline_status.json)
- Download log: `data/reselection_wise_fetch.log`
- Continuation log: `data/reselection_pipeline.log`
- Scripts: `15a_expand_parent.py`, `15b_rebuild_observing_pool.py`, `15c_project_expanded_pool.py`, `15d_finish_reselection.py`, `15e_pool_report.py`.

The continuation process waits for acquisition to finish, projects the new light curves onto the saved manifold, refreshes the complete queues, and rewrites the pool report. It does not assign scientific weights or publish anything. Failures are recorded in the status file; a failed acquisition can be resumed with `15c_project_expanded_pool.py fetch --rmax 20.5 --workers 32`, followed by rerunning the continuation.

Two coverage issues remain explicit: the paper's original object-ID-to-coordinate lookup is missing, and the SDSS-V rolling master is newer than the cached tagged catalog used for this first pass. Resolve these before claiming a complete parent or assigning final spectral-recency weights. All joined research products stay local and git-ignored because they include collaboration metadata.
