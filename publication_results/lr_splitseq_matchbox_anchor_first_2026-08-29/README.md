# Matchbox LR-SPLiT-seq anchor-first update (2026-08-29)

This directory archives the compact, machine-readable results and frozen
specifications used to update the primary dual-orientation and supplementary
forward-only Matchbox LR-SPLiT-seq results. Raw attempts remain under:

- `/scratch1/seqproc-ecosystem/campaigns/journal-matchbox-lr-anchor-first-2026-08-29-v2`
- `/scratch1/seqproc-ecosystem/campaigns/journal-matchbox-lr-forward-anchor-first-2026-08-29`

The final configuration searches an exact linker over the long read before
testing the captured terminal barcode against its canonical list. It keeps an
orientation-discriminating barcode in each pattern. This preserves the exact
barcode/linker constraints while avoiding a full long-read scan for every
terminal-barcode candidate.

## Primary dual-orientation results

The `/dev/null` timing means +/- sample SD over three full-input replicates
are 559.617 +/- 2.268, 196.187 +/- 11.115, 44.183 +/- 2.453, and
25.442 +/- 0.302 seconds at 1, 4, 16, and 32 threads. Peak RSS at the fastest
tested point is 68,556 KiB (67.0 MiB). The 1-thread values come from a complete
replacement block because two original attempts consumed the expected CPU work
but were externally stalled for substantial wall time; the original attempts
are preserved and were not selectively mixed into the replacement summary.

The materialized 32-thread run emitted 8,273 reads. All are in the 560,699-read
conservative structural reference: precision 1.0, recall
0.014754797137144886, and F1 0.029080517143198613.

## Four-read adjudication

Relative to the former direct barcode-first transcription, the anchor-first
`one-best` configuration omits exactly four reads, listed in
`adjudication.json`. Each contains one exact, whitelist-valid cassette and at
least one additional exact-linker cassette with an invalid terminal barcode.
Matchbox chooses among pattern matches before a body whitelist guard is
evaluated, so these molecules become ambiguous under the anchor-first pattern.
`--match-mode all` recovers them but also emits duplicate FASTQ records for
multi-cassette molecules; the deterministic `one-best` policy is retained.

## Forward-only results

The forward anchor-first configuration emits a byte-identical 4,162-record
FASTQ relative to the former barcode-first configuration on the complete
input (SHA-256
`89980c6083754d654a7548df1289e7e2d759487b065d4145dea2f4c0d5eb5164`).
The `/dev/null` timing means +/- sample SD are 420.893 +/- 0.703,
110.146 +/- 0.281, 34.119 +/- 0.151, and 21.671 +/- 0.307 seconds at
1, 4, 16, and 32 threads. Peak RSS at the fastest point is 36.0 MiB.

