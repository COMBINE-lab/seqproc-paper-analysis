# Journal rerun campaign (2026-08-17)

This directory contains curated, commit-sized artifacts from the staged journal
rerun. The full materialized outputs live under
`/scratch1/seqproc-ecosystem/campaigns/journal-rerun-2026-08-17/`.

The accuracy manifest contains all current technology blocks. Blocks are run
incrementally with `run_frozen_schedule.py --dataset ...`, and the reusable
`scripts/summarize_publication_accuracy.py` evaluator regenerates accession
bitmaps, validates splitcode component ID equivalence, unions LR orientations,
and computes accuracy against the revised conservative structural references.

See `PUBLICATION_RERUN_STATUS.md` at the repository root for the active matrix
of completed and pending accuracy and performance results.
