# Local downloaded data

Everything below this directory except this README is ignored by Git.

- `spherex/sample/`: QR2 query responses, image indexes, per-image extraction
  records, extracted tables, progress, and resumable-worker state.
- `spherex/rejection_checks/`: downloaded FITS used to diagnose crowded R12.

These are preserved working data, not disposable failed downloads. In particular,
the per-image JSON records allow retrieval to resume without repeating completed
work. See [the retrieval instructions](../review/SPHEREX.md).

Curated Figure 1 spectra remain in `../inputs/`. Scientific review tables and
plots remain in `../review/spherex_sample_review/`; the compact R12 diagnostics
remain there under `rejection_checks/diagnostics.json`.
