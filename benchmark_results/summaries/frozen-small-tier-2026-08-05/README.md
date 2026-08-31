# Frozen scheduler integration tier — 2026-08-05

This is a development-mode acceptance test of the publication benchmark
harness, not manuscript performance evidence. Each tool invocation processes
only 100,000 records and many conditions finish well below one second; their
timings are retained for diagnostics but are too short for comparative claims.

## Frozen design

- seqproc binary SHA-256:
  `8a57dbc9c502cdffe5a7d5f8e6cf50d326a4999090b7428efac76564b68978e8`
- seqproc commit: `c9d71476527eea050cedd0431d29d39043e4e517`
- antisequence commit: `19eef2c00bd19ec1fde8b8bc041db3371e1a9875`
- manifest SHA-256:
  `517aa2f9fa1ec93f388e0ea52dbd846e60f8e87e850bedfd8c4269b5ec49f239`
- schedule SHA-256:
  `572ab9601c56d5fc0e4c5f60b8bafc93a6597ba91ca31d63fca0227d0b50a359`
- protocols: SPLiT-seq paired-end, LR-SPLiT-seq, 10x v2, sci-RNA-seq3
- modes: default worker path and bounded staged pipeline
- threads: 1, 8, 32
- replicates: 3
- scheduled conditions: 72

The checked-in manifest hashes the binary, inputs, geometries, mapping files,
and the generator, runner, harness, and aggregator scripts. The schedule was
created once from seed `741211` and protected by a SHA-256 sidecar.

## Acceptance result

- 72/72 scheduled conditions produced one valid successful attempt.
- 0 conditions were excluded by aggregation.
- All output FASTQs parsed fully and contained 100,000 records per mate.
- For every protocol, the normalized multiset digest was identical across all
  modes, thread counts, and replicates.
- An immediate resume completed 0 runs and skipped all 72 successful runs.
- The scheduler now holds an OS-backed output-root lock; its test suite refuses
  a concurrent coordinator while allowing recovery after process exit.

`aggregate.json` is the canonical machine-readable result. `runs.csv` contains
one row per valid attempt, `summary.csv` contains diagnostic timing/RSS groups,
and `correctness.json` isolates the output-equivalence gate. The 1.7 GB raw run
tree remains outside Git under
`benchmark_results/development/current-small-tier-locked-2026-08-05/`.

## Interpretation

This tier establishes that the current default and staged implementations are
semantically identical on the four paper protocols under the tested thread
counts. It also validates randomized scheduling, frozen input/config checks,
FASTQ structural checks, output counting, order-independent complete-record
digests, failure preservation, and idempotent resume.

It does **not** replace the full five-replicate, full-dataset, exclusive-node,
cross-tool rerun. In particular, sub-second timings here should not be cited in
the manuscript or used to choose between default and staged execution.
