# SPLiT-seq PE guarded-fuzzy Matchbox update (2026-08-29)

This directory freezes the compact outputs used to update the manuscript after
replacing Matchbox's exact-linker SPLiT-seq PE configuration with a
boundary-safe configuration that combines:

- exact membership in the canonical barcode lists;
- edit distance at most three for each 30-nt linker; and
- explicit length guards for the UMI and three captured barcode fields.

The guards prevent approximate linker matches from silently shifting captured
field boundaries.  The primary configuration therefore allows the same linker
error budget as the protocol's conservative structural reference without
requiring an externally expanded barcode list.  The Hamming-1-expanded variant
changes only the barcode lists and remains a sensitivity analysis.

## Primary preprocessing results

Against the 57,437,503-read conservative structural reference, the canonical
configuration emits 50,266,141 pairs (64.7583% of input), with 99.7777%
precision, 87.3200% recall, and F1 0.931341.  The expanded sensitivity emits
57,252,325 pairs (73.7586%), with 99.5760% precision, 99.2550% recall, and F1
0.994153.  The canonical set is a strict subset of the expanded set.

The three-replicate `/dev/null` primary timing means are 4,819.408, 1,273.112,
470.386, and 358.173 seconds at 1, 4, 16, and 32 threads, respectively.  The
largest peak RSS at the fastest tested point (32 threads) is 11,085.5 MiB.  A
separate 32-thread materialized-output expanded sensitivity run took 6,312.330
seconds and 12,700.9 MiB; it is not a replicate-matched `/dev/null` benchmark.

## Downstream configurations

The four subdirectories hold STARsolo and downstream summaries for canonical
or Hamming-1-expanded Matchbox input under `EditDist_2` or `1MM`.  Only the
Matchbox STARsolo matrix was recomputed.  The unchanged seqproc and splitcode
matrices were reused from the pinned runs with the corresponding correction
mode.  Each subdirectory includes input provenance, structural-reference and
pairwise read-set metrics, STARsolo resources, count concordance, marker-based
cell typing, and clustering concordance.

`canonical-editdist2` is the primary downstream analysis.  The other three
subdirectories are controlled sensitivity analyses.

## Source campaign

Large FASTQs, bitmaps, and full run records remain outside Git under:

`/scratch1/seqproc-ecosystem/campaigns/journal-matchbox-splitseq-pe-fuzzy-2026-08-28`

The source configuration is `configs/matchbox/publication_splitseq_pe.mb`; the
expanded sensitivity is
`configs/matchbox/sensitivity_splitseq_pe_ham1_expanded.mb`.  The reusable
Matchbox-only downstream driver is
`biological_analysis/run_matchbox_variant_downstream.sh`.
