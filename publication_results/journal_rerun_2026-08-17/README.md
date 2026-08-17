# Journal rerun campaign (2026-08-17)

This directory contains curated, commit-sized artifacts from the staged journal
rerun. The full materialized outputs live under
`/scratch1/seqproc-ecosystem/campaigns/journal-rerun-2026-08-17/` and its
configuration-correction revision
`/scratch1/seqproc-ecosystem/campaigns/journal-rerun-2026-08-17-r2/`.

The accuracy manifest contains all current technology blocks. Blocks are run
incrementally with `run_frozen_schedule.py --dataset ...`, and the reusable
`scripts/summarize_publication_accuracy.py` evaluator regenerates accession
bitmaps, validates splitcode component ID equivalence, unions LR orientations,
and computes accuracy against the revised conservative structural references.

See `PUBLICATION_RERUN_STATUS.md` at the repository root for the active matrix
of completed and pending accuracy and performance results.

`lr_matchbox_hamming_expansion_sensitivity.json` evaluates the deterministic
August 11 Matchbox Hamming-expanded/exact-linker output against the revised
current LR-SPLiT-seq reference and current seqproc/splitcode outputs. Complete
binary, input, configuration, expansion-list, timing, and ambiguity provenance
for that output is retained in
`../lr_splitseq_hamming_expansion_grid_2026-08-11/results.json`.

The `scirnaseq3_*` artifacts and the `*-r2` frozen manifest/schedule record the
official sci-RNA-seq3 accuracy block. Revision 2 replaces an ordered two-branch
Matchbox expression that accidentally shadowed the 10-nt BC1 geometry with a
single variable-prefix capture guarded to the documented 9- or 10-nt lengths.
