# Journal performance campaign (2026-08-17)

This directory contains commit-sized artifacts from the frozen full-data
performance campaign. The master grid uses three randomized replicates for
each tool at 1, 4, 16, and 32 threads, pins requested threads to physical CPUs
1 through N, consumes uncompressed inputs, and directs requested biological
sequence output to `/dev/null`.

`frozen_manifest.yaml` records exact inputs, configurations, binaries,
repositories, commands, hashes, hardware, and measurement policy.
`frozen_schedule.json` and its digest record the immutable randomized order.
`execution-log.jsonl` contains coordinator events for the completed blocks.
Each completed block directory contains replicate-level `runs.csv`, aggregated
`summary.csv`, machine-readable `aggregate.json`, and `correctness.json`.

Completed in this checkpoint:

- `scirnaseq3`: 36/36 successful conditions.
- `lr_splitseq_forward`: 36/36 successful conditions.

Timing runs intentionally do not materialize FASTQ output. Exact read-set and
output-shape validation is reported separately in the journal accuracy
campaign. In particular, the LR Matchbox timings use the accuracy-audited,
boundary-safe exact canonical-list configuration; they must not be compared to
the much faster superseded fuzzy-capture configuration.

The full resumable run tree remains at
`/scratch1/seqproc-ecosystem/campaigns/journal-performance-2026-08-17/`.
