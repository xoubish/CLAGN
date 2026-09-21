# Local historical archive

`2026-09-20/` contains material removed from the active project layout during cleanup:

- `finders/`: earlier finder charts, target lists, and observing CSVs.
- `plots/`: earlier September-only visibility plots.
- `legacy_pages/`: superseded v2 pages.
- `legacy_drivers/`: six original shell drivers.
- `planning/`: previous follow-up plan and status narrative.
- `review_2026-09-20/`: the initial detailed review, audit figures, and cached reference documents.
- `logs/`: inactive top-level data logs.
- `before_cleanup/`: original Python sources and affected page/template copies before path edits.
- `move_manifest.json`: old/new paths and SHA-256 hashes of the moved files.

The archive is local and git-ignored because it includes collaboration-derived products. Files were moved, not discarded. Old scripts/reports retain their historical paths and should be treated as snapshots; use `scripts/` for the maintained entry points.

To recover an individual historical file, use the manifest to locate its destination and copy it where needed. Its hash records the content before any subsequent path edits. Avoid restoring the whole tree over newer observing products.
