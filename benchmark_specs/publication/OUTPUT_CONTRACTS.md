# Publication benchmark workload and output contracts

## Benchmark policy

The benchmark uses two separate full-input tracks.

1. **Primary timing/RSS track.** Every biological sequence output is directed
   to `/dev/null`. The tools still parse, match, filter, correct, project, and
   serialize the configured products, but storage throughput and retained-file
   size do not determine the headline runtime. Small harness metadata (GNU time,
   logs, and the splitcode dual-pass report) remains materialized.
2. **Correctness track.** One run per dataset and tool materializes the FASTQ
   products on the same reported input and validates FASTQ structure, record
   counts, accepted numeric-accession ID sets, and product-specific semantics.
   These runs are retained for the downstream validity/concordance checks and
   are not used as headline timing observations. Run them one technology block
   at a time and remove a block only after its derived checksummed reports have
   been generated, because the complete materialized set is large.

Configurations are **qualitatively aligned**: each tool must perform the same
kind of biological operation. Byte-identical cross-tool output is not required
where the tools expose different representations or matching semantics. Within
each tool and mode, deterministic read sets remain required.

## Per-technology contracts

| Technology | Semantic workload | seqproc / matchbox product | splitcode product | Known qualification |
|---|---|---|---|---|
| SPLiT-seq paired-end | Filter on both linkers and three whitelist barcodes; project UMI+BC3+BC2+BC1 while retaining R1 | R1 plus a compact read: seqproc emits canonical whitelist values; matchbox emits the matched segments | Selected R1 with assignment-prefix output disabled, plus `umi_bc3_bc2_bc1.fastq`, one 30-nt UMI/canonical-barcode projection | Matching and ambiguity policies are best-practical rather than identical; indel-tolerant matchbox captures can have non-nominal lengths, which are counted |
| LR-SPLiT-seq, dual orientation | Filter both orientations on two approximate linkers and project UMI+BC3+BC2+BC1 (nominally 10+8+8+6 nt) | One compact read; native dual-orientation processing | Two sequential passes over forward and precomputed reverse-complement inputs; `prefix`, `bc2`, and `bc1` component FASTQs jointly form the compact product | Approximate-anchor edge cases can yield non-nominal component lengths and are counted as validity outcomes; no duplicate reconciliation is performed; splitcode matching is substitution-only |
| LR-SPLiT-seq, forward-only | Same projection, forward orientation only | One compact read | `prefix`, `bc2`, and `bc1` component FASTQs | Supplementary controlled comparison |
| 10x Chromium v2 | Reject a pair when R1 is shorter than 26 nt; otherwise pass R1/R2 through | R1/R2 | R1/R2 under `--trim-only` with a 26-nt R1 minimum; no `@extract` | On the reported data, R1 is 26 nt and R2 is 98 nt |
| sci-RNA-seq3 | Match the approximate anchor and project BC1+BC2+UMI while retaining R2 | Normalized BC1+BC2+UMI and R2 | Selected R2 with assignment-prefix output disabled, plus raw extracted BC1+BC2+UMI | The current splitcode configuration does not reproduce seqproc's variable-length `norm` encoding; validate logical components separately |

## Reporting rules

- Runtime and peak RSS tables use only records whose metadata says
  `measurement_track: timing` and `sequence_output_policy: dev-null`.
- Accuracy, validity, concordance, and transformed-sequence checks use only the
  materialized correctness track.
- Product-length histograms are recorded for every correctness FASTQ. Products
  whose semantics require fixed lengths fail the run if a non-nominal length
  appears. Approximate-match-derived projections (LR outputs and the matchbox
  SPLiT-seq PE compact output) retain such records but report their non-nominal
  count as a validity result rather than silently treating them as well-formed
  barcodes.
- Output-to-disk and gzip measurements, if reported, are a separate I/O
  supplement and must not be mixed into the primary `/dev/null` table.
- The LR dual splitcode row includes both sequential passes. Reverse-complement
  staging occurs outside measurement; duplicate reconciliation is described as
  an unsupported post-processing requirement and is neither run nor timed.
