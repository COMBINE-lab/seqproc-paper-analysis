# Statistics collection overhead (development checkpoint)

This experiment measures the runtime cost of optional ANTISEQUENCE statistics
after moving hot counters into worker-local accumulators. It is a development
checkpoint, not a substitute for the randomized full-paper benchmark rerun.

## Workload

- 5,000,000 synthetic FASTQ records per repetition
- short-pattern Hamming lookup, which is deliberately cheap and therefore a
  conservative stress test for instrumentation overhead
- 8 transform workers, bounded pipeline, batch size 256, null output
- 30 internal repetitions per condition
- current implementation: ANTISEQUENCE commit `96c9c9ffe1b516d5a76b469adecbc05058b91ae0`
- pre-refactor implementation: commit `2f930767`; archived binary SHA-256
  `463244047de2aecd9db51fb6a42d4ce330e8d6f66759bdd26012677d641af831`

The benchmark's `seconds` interval includes FASTQ parsing, transformation,
statistics aggregation, and graph teardown, but excludes synthetic FASTQ
construction and graph construction. GNU time measures the complete 30-repeat
process and supplies peak RSS.

## Results

| Implementation | Statistics | Mean seconds | Throughput (M reads/s) | Overhead vs same-build off | Peak RSS (KiB) |
|---|---:|---:|---:|---:|---:|
| pre-refactor | off | 0.447020 | 11.185 | reference | 330,992 |
| pre-refactor | detailed | 0.457912 | 10.919 | 2.44% | 331,272 |
| worker-local | off | 0.452167 | 11.058 | reference | 331,068 |
| worker-local | basic | 0.452187 | 11.057 | 0.004% | 331,000 |
| worker-local | detailed | 0.459046 | 10.892 | 1.52% | 331,000 |

The paired overhead fell from 2.44% to 1.52% for detailed collection, a 38%
reduction in the instrumentation penalty. Basic collection was indistinguishable
from off in this run. End-of-run aggregation itself averaged 2.47 microseconds
for basic and 4.96 microseconds for detailed statistics. Peak-RSS differences
were below 0.1% and should be treated as measurement noise.

Absolute current-versus-pre-refactor times should not be interpreted as an A/B
speed claim: the two paired blocks were not interleaved or randomized. The
within-build on/off contrasts are the intended estimands. The manuscript's
headline performance results should remain statistics-off, accompanied by a
separate paired instrumentation-overhead experiment.

## Reproduction and provenance

The exact current run specifications in this repository are:

- `benchmark_specs/development/hamming-statistics-off.json`
- `benchmark_specs/development/hamming-statistics-basic.json`
- `benchmark_specs/development/hamming-statistics-on.json`

The pre-refactor specifications are:

- `benchmark_specs/development/hamming-statistics-pre-refactor-off.json`
- `benchmark_specs/development/hamming-statistics-pre-refactor-detailed.json`

Run a specification with `scripts/benchmark_harness.py`, then aggregate the
result directories with `scripts/summarize_hot_path.py`. The complete raw run
directories are retained in the development workspace at
`journal-readiness-review/results/statistics-overhead-2026-08-05`; each contains
the exact argv, environment, binary digest, host snapshot, raw benchmark JSON,
GNU-time output, and run ID. `current-summary.json` and
`pre-refactor-summary.json` are the committed generated artifacts and numerical
sources for the table above. The raw run directories should be included in the
eventual DOI reproduction archive.
