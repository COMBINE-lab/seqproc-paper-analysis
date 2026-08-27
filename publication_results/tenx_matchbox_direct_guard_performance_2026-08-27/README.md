# 10x v2 Matchbox direct-guard performance rerun (2026-08-27)

This targeted frozen rerun corrects an avoidable Matchbox configuration cost in
the 10x Chromium v2 performance block.  The previous configuration decomposed
R1 into a nested `bcumi` capture containing 16-base `bc` and 10-base `umi`
captures, even though none of those captured values were used.  The benchmark's
declared workload is only to reject pairs whose R1 is shorter than 26 bases and
otherwise pass both reads through unchanged.  The revised configuration states
that operation directly:

```text
if read.r1.seq.len() >= 26 {
    read.r1.out!('mb_r1.fq')
    read.r2.out!('mb_r2.fq')
}
```

The archived full-data correctness audit in
`../journal_rerun_2026-08-17/tenx_v2_accuracy_metrics.json` establishes that
the previous Matchbox output contained all 234,382,218 input pairs, that every
emitted R1 was exactly 26 bases, and that both output FASTQs were byte-identical
to the inputs.  A separate 10,000-pair A/B check also produced byte-identical
old/new outputs.  The change therefore removes unused capture construction
without changing the accepted set or output bytes for this dataset.

The rerun used the same Matchbox v0.3.2 release binary, full uncompressed
inputs, CPU-affinity policy, `/dev/null` biological-output policy, thread grid
(1, 4, 16, and 32), and three-replicate randomized schedule as the publication
campaign.  The frozen full schedule contains all publication conditions; the
new `--dataset tenx_v2 --tool matchbox` selectors preserve the selected
conditions' relative frozen order while executing and summarizing only these
12 conditions.

| Threads | Old mean (s) | New mean ± SD (s) | Runtime reduction | New peak RSS cap (MiB) |
|---:|---:|---:|---:|---:|
| 1 | 1468.733 | 1075.199 ± 3.455 | 26.8% | 37.2 |
| 4 | 582.965 | 476.547 ± 3.052 | 18.3% | 38.4 |
| 16 | 534.499 | 473.005 ± 2.159 | 11.5% | 40.5 |
| 32 | 575.029 | 522.985 ± 5.658 | 9.1% | 42.0 |

All 12 scheduled conditions completed successfully and the aggregate has zero
exclusions.  Runtime improves substantially, especially at one thread, while
the plateau above four threads remains.  The latter is attributable to the
effective parallelism of Matchbox's paired FASTQ streaming/output path rather
than the now-trivial 10x predicate.

Files:

- `frozen_manifest.yaml`: exact inputs, configurations, binaries, repositories,
  commands, hashes, hardware, and measurement policy.
- `frozen_schedule.json` and `.sha256`: immutable randomized full schedule.
- `execution-log.jsonl`: start/finish records for the 12 selected conditions.
- `tenx_v2/runs.csv`: replicate-level measurements.
- `tenx_v2/summary.csv`: machine-generated thread-level aggregates.
- `tenx_v2/aggregate.json`: selection provenance, validation counts, summaries,
  and exclusions.

The complete resumable run tree remains at
`/scratch1/seqproc-ecosystem/campaigns/tenx-matchbox-direct-2026-08-27/`.
