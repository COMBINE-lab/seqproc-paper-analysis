# Matchbox SPLiT-seq PE anchor exploration (2026-08-29)

This directory records the output-equivalence and performance gates used to
evaluate fixed-window/post-placement whitelist checks for the Matchbox
SPLiT-seq PE configurations.

## Decision

- The primary canonical-list configuration remains unchanged. A safe
  terminal-BC1-anchor transcription was byte-identical on all 77,621,181
  input pairs, but its full 32-thread runtime and memory were effectively
  neutral relative to the current primary configuration.
- The Hamming-one-expanded sensitivity configuration now retains terminal BC1
  as an exact pattern anchor, captures the fixed-position BC3 and BC2 windows,
  and checks those two windows for list membership after linker placement.
  Canonical entries precede noncanonical radius-one variants because Matchbox
  0.3.2 implements `contains()` as an ordered linear scan.
- The adopted expanded configuration is byte-identical to the former
  three-anchor configuration across the complete transformed R2 FASTQ:
  57,252,325 records and SHA-256
  `9d8c8f68e2c395454ee10f7d2c7fd9f1eab5aedd3c7b4ce10924214ace504aab`.
- The frozen, uncontended, 32-thread `/dev/null` timing is 4,272.491561808 s
  with peak RSS 12,969,028 KiB (12,665.1 MiB). The former value was
  6,312.330 s and 12,700.9 MiB, so the safe configuration is 32.3% faster
  while retaining the same output.

## Rejected transcriptions

- Moving all three barcodes to body membership checks emitted only 3,420,289
  records on the full canonical-list input rather than 50,266,141. The fuzzy
  linker placement is selected before body guards and is not reconsidered
  when a captured window fails membership.
- On the first 100,000 pairs, retaining BC2 rather than terminal BC1 emitted
  only 4,195 records rather than 63,012.
- A nested exact-matcher formulation either reached Matchbox's unimplemented
  both-ends-pinned matcher case or, with an unbounded terminal hole, emitted
  63,134 records rather than 63,012 on the first 100,000 pairs.

These failures show why a full output-equivalence gate is required even for a
fixed-orientation, fixed-position PE layout.

## Frozen evidence

- `frozen_timing_manifest.yaml` and `frozen_timing_schedule.json`: immutable
  single-condition timing specification.
- `timing_run.json` and `timing_time.txt`: harness result and raw GNU time
  record for run ID `74c950bdb14b944032b7d6117bd3369d7fb3d168b95f1a1c8353bf17ffdc599d`.
- `equivalence.json`: complete and subset equivalence outcomes, including
  diagnostic configurations that were rejected.

The timing run records analysis commit
`81ae0593739eaea41e0d189506cfeff231fc394e`.
