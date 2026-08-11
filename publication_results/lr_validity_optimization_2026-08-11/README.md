# LR-SPLiT-seq validity optimization (exploratory)

These are one-run diagnostic measurements, not final benchmark estimates.
Every seqproc row used 32 threads pinned to CPUs 1-32, the optimized seqproc
binary at commit `edbf9ef`, the complete 5,764,421-read
`SRR13948564_full.fastq`, staged execution, and materialized FASTQ output.
Wall time and maximum RSS come from GNU `time -v`.

The repaired structural reference contains 588,233 reads (10.20% of input).
Its full Python/edlib scan took 108.91 s and 97,112 KiB maximum RSS. Relative
to the prior 601,603-read reference, 564,063 IDs overlap, 37,540 old IDs are
removed, and 24,170 IDs are newly accepted.

| Configuration | Wall (s) | RSS (MiB) | Emitted | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| Previous dual geometry, rescored | 3.74 | 129.0 | 2,865,883 | 20.07% | 97.80% | 0.3331 |
| Complete leading components | 3.73 | 128.0 | 2,507,488 | 22.94% | 97.78% | 0.3716 |
| Complete components, searched L2 (6,4) | 3.71 | 124.0 | 2,297,584 | 24.41% | 95.36% | 0.3888 |
| Complete components, searched L2 (3,3) | 4.14 | 128.0 | 2,258,682 | 24.40% | 93.67% | 0.3871 |
| Complete components, searched L2 (4,3) | 3.93 | 125.0 | 2,269,777 | 24.34% | 93.93% | 0.3867 |
| Complete components, searched L2 (5,4) | 3.81 | 128.0 | 2,293,874 | 24.44% | 95.31% | 0.3891 |
| Searched L2 (5,4) + whitelist | 3.93 | 110.0 | 546,324 | 98.86% | 91.81% | 0.9520 |
| Adjacent L2 (5,4), no whitelist | 3.73 | 126.0 | 2,474,202 | 23.12% | 97.23% | 0.3735 |
| Adjacent L2 (5,4) + whitelist | 4.03 | 110.0 | 553,840 | 100.00% | 94.15% | 0.9699 |

All outputs in the complete-component rows have exactly 32 sequence bases.
The best permissive searched-linker threshold in this small sweep is (5,4).
The recommended optional high-specificity mode is adjacent L2 (5,4) with
whitelist filtering.

The high-specificity metric is agreement with a conservative structural
reference, not biological ground truth. The seqproc mode and reference share
the same linker sequences, barcode whitelists, and tolerance family, although
the reference is implemented independently and searches multiple candidate
placements. Perfect observed precision must therefore not be presented as an
independent accuracy estimate.

Discordance adjudication scanned the full FASTQ in 8.80 s with 58,984 KiB RSS.
The permissive (5,4) mode had 40,821 calls absent from both matchbox and the
splitcode dual-pass union; 7,026 passed the repaired reference and 33,795 did
not. The adjacent-whitelist mode had 6,869 competitor-unique calls, all of
which passed the reference; 4,411 (64.2%) were reverse-orientation reads.

Exact machine-readable values are in `results.json`. Full materialized outputs,
bitmaps, GNU-time records, validator output, and the adjudication summary are
under `/scratch1/seqproc-ecosystem/lr-validity-optimization-2026-08-11/`.
