# Aligned performance campaign (2026-08-10)

This directory archives the frozen inputs and machine-generated aggregates for
the revised manuscript performance table. The primary grid contains three
randomized full-data replicates for each tool at 1, 4, 16, and 32 threads.
Inputs were uncompressed and requested sequence output was directed to
`/dev/null`.

- `frozen_manifest.yaml`: commands, qualitatively aligned configurations,
  executable and input hashes, hardware, affinity, and measurement policy.
- `frozen_schedule.json`: randomized execution order.
- `execution-log.jsonl`: condition-level execution records.
- `<block>/runs.csv`: replicate-level wall time and peak RSS.
- `<block>/summary.csv`: per-condition aggregates used in the manuscript.
- `<block>/aggregate.json`: machine-readable aggregate plus provenance.

`lr_splitseq_dual` is the primary LR comparison. seqproc and matchbox search
both orientations natively; splitcode is the sum of two sequential passes over
forward and precomputed reverse-complement inputs. Reverse-complement creation
and duplicate reconciliation are not timed. `lr_splitseq_forward` is the
controlled forward-only supplementary block.

The pinned executables are seqproc commit
`edbf9efcfa69d588a42d3f1fb714742b675502eb`, matchbox 0.3.2, and splitcode
0.31.6. Exact executable SHA-256 values are in the frozen manifest.
