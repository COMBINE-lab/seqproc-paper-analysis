# Publication rerun status

Last updated: 2026-08-17

This is the active checklist for results computed with the revised structural
references and semantically aligned publication configurations. Historical
results produced with preprint-era or August diagnostic configurations do not
count as complete here.

## Completion matrix

| Technology block | Structural reference | Accuracy: 1 run/tool, 32 threads | Final speed/RSS: 3 replicates at 1/4/16/32 threads |
|---|---|---|---|
| SPLiT-seq PE, primary | Complete | **Complete** | Pending |
| LR-SPLiT-seq dual orientation, primary | Complete | **Complete** | Pending |
| LR-SPLiT-seq forward-only, supplementary | Definition/configs complete; forward-only reference subset pending | Pending | **Complete** |
| 10x Chromium v2, primary | Complete | **Complete** | Pending |
| sci-RNA-seq3, primary | Complete | **Complete** | **Complete** |

Final performance runs direct biological outputs to `/dev/null`; correctness
runs materialize and validate outputs. Therefore, diagnostic time/RSS recorded
during an accuracy run is not substituted for the final performance cells.

## Newly completed: LR-SPLiT-seq dual-orientation accuracy

Input: 5,764,421 reads. Conservative structural reference: 560,699 reads.
All tools used one frozen 32-thread run pinned to physical CPUs 1--32.

| Tool | Emitted | Intersection | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|
| seqproc | 545,702 | 545,702 | 1.000000 | 0.973253 | 0.986445 |
| matchbox | 8,277 | 8,277 | 1.000000 | 0.014762 | 0.029094 |
| splitcode, union of two orientations | 428,941 | 427,895 | 0.997561 | 0.763146 | 0.864749 |

The low matchbox recall is not an accidental use of the old script. It is the
known limitation of the boundary-safe, canonical-list, no-preprocessing mode:
native approximate captures can move adjacent cassette boundaries, while exact
matching requires every barcode and linker to be error-free. Explicit Hamming
expansion improves recall but is retained only as a disclosed sensitivity
analysis because it requires user-visible preprocessing and ambiguity policy.

Splitcode emitted 218,167 forward and 211,490 reverse-pass reads; their union is
428,941, so 716 read IDs occur in both orientations. Its prefix, BC2, and BC1
component outputs had identical read-ID sets within each pass.

### Accuracy-run resource diagnostics (not final performance numbers)

| Tool | Wall time | Peak RSS |
|---|---:|---:|
| seqproc | 4.219 s | 111.0 MiB |
| matchbox | 203.005 s | 56.1 MiB |
| splitcode two-pass wrapper | 17.285 s | 531.4 MiB |

### Matchbox Hamming-expansion sensitivity analysis

Expanding each canonical barcode to all unique sequences within Hamming
distance one, while retaining **exact linker matching**, substantially rescues
Matchbox without moving cassette boundaries:

| Matchbox configuration | Emitted | Intersection | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|
| Canonical lists, exact matching (primary) | 8,277 | 8,277 | 1.000000 | 0.014762 | 0.029094 |
| Hamming-expanded barcodes, exact linkers (sensitivity) | 406,257 | 406,257 | 1.000000 | 0.724555 | 0.840280 |

This is a 49.1-fold increase in retained reads and a 70.98-percentage-point
increase in recall. Against the same reference, it approaches but remains below
splitcode (recall 0.763146, F1 0.864749), and remains below seqproc (recall
0.973253, F1 0.986445). The expanded Matchbox set overlaps 405,829 seqproc
reads and 386,396 splitcode reads; the corresponding Jaccard indices are
0.743100 and 0.860950.

The expansion contains 2,400 unique BC2/BC3 variants from 96 canonical 8-mers,
with no collisions, and 1,614 unique BC1 variants from 96 canonical 6-mers.
For BC1, 206 expanded sequences have two or three possible canonical owners.
Matchbox uses these lists as an acceptance filter and emits the observed
sequence, so this configuration means "accept once" rather than performing an
unambiguous correction to a canonical barcode. That ambiguity policy and the
need for an external expansion step are material usability differences from
seqproc's native approximate-whitelist matching.

Approximate-linker matching must not be combined with this expansion: the
tested configuration shifted capture endpoints and emitted 237,655 records of
33--47 nt instead of the intended 32 nt. Exact linkers produced 406,257/406,257
correct-length outputs. The expanded exact-linker result was recomputed against
the current 560,699-read reference from a deterministic August 11 output whose
input, Matchbox binary, configuration, and expanded-list hashes are unchanged.
Its original one-run diagnostic cost was 4,308 s (71.8 min) and 434.1 MiB peak
RSS; this is not a final `/dev/null` performance result.

## Newly completed: sci-RNA-seq3 accuracy

Input: 22,088,821 read pairs. Conservative structural reference: 19,864,110
pairs. All tools used one frozen 32-thread run pinned to physical CPUs 1--32
and performed the aligned anchor-filter plus BC1+BC2+UMI projection workload.

| Tool | Emitted | Intersection | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|
| seqproc | 19,879,895 | 19,853,935 | 0.998694 | 0.999488 | 0.999091 |
| matchbox | 19,559,040 | 19,533,875 | 0.998713 | 0.983375 | 0.990985 |
| splitcode | 19,695,955 | 19,669,995 | 0.998682 | 0.990228 | 0.994437 |

Matchbox required a configuration correction before this result was frozen.
The first attempt expressed the documented 9- and 10-nt BC1 geometries as two
ordered approximate-match branches. A match of the first branch shadows the
second branch before its length guard is evaluated, leaving only 4,571 28-nt
products and reducing recall to 0.600593. The corrected script captures one
variable-length prefix and then accepts only lengths 9 or 10. It emits
7,901,629 28-nt and 11,657,411 27-nt products, raising recall to 0.983375
without reducing precision. The superseded result is diagnostic and is not a
publication result.

All Matchbox and splitcode accepted IDs are subsets of seqproc's accepted IDs.
seqproc and splitcode overlap on 19,695,955 IDs (Jaccard 0.990747); seqproc and
Matchbox overlap on 19,559,040 (Jaccard 0.983860); Matchbox and splitcode
overlap on 19,489,430 (Jaccard 0.986029).

The reference rejects 25,960 reads whose equally good edit-one anchor matches
occur at both allowed BC1 offsets. seqproc and splitcode each retain all 25,960
of these, while Matchbox retains 25,165. This accounts for essentially all
nominal false positives and should be described as an ambiguity-policy
difference rather than evidence of biologically invalid reads.

### Accuracy-run resource diagnostics (not final performance numbers)

| Tool | Wall time | Peak RSS |
|---|---:|---:|
| seqproc | 6.772 s | 74.0 MiB |
| matchbox | 78.061 s | 3,102.2 MiB |
| splitcode | 11.179 s | 1,834.9 MiB |

These runs materialized and validated outputs. In particular, Matchbox's
corrected variable-prefix search changes its performance profile, so the old
Matchbox sci-RNA-seq3 `/dev/null` timing results must be replaced in the final
performance campaign.

## Newly completed: paired-end SPLiT-seq accuracy

Input: 77,621,181 read pairs. Conservative structural reference: 57,437,503
pairs. All tools used one frozen 32-thread run pinned to physical CPUs 1--32
and performed the aligned barcode-matching, replacement/projection, and
retention workload.

| Tool | Emitted | Intersection | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|
| seqproc | 56,991,381 | 56,796,243 | 0.996576 | 0.988836 | 0.992691 |
| matchbox | 35,160,366 | 35,160,366 | 1.000000 | 0.612150 | 0.759421 |
| splitcode | 53,495,595 | 53,472,550 | 0.999569 | 0.930969 | 0.964050 |

Matchbox's accepted IDs are a strict subset of both other tools' accepted IDs.
It uses the same boundary-safe, canonical-list, no-preprocessing policy as the
primary long-read comparison: exact whitelist matching avoids shifting
adjacent barcode boundaries, but rejects reads containing barcode errors. The
result is therefore high precision but substantially lower recall. seqproc and
splitcode overlap on 53,232,995 IDs (Jaccard 0.929769); seqproc contains
3,758,386 additional IDs and splitcode contains 262,600 additional IDs.

### Accuracy-run resource diagnostics (not final performance numbers)

| Tool | Wall time | Peak RSS |
|---|---:|---:|
| seqproc | 27.897 s | 92.0 MiB |
| matchbox | 281.690 s | 7,751.1 MiB |
| splitcode | 50.249 s | 1,282.5 MiB |

## Newly completed: 10x Chromium v2 accuracy

Input: 234,382,218 read pairs. Every R1 in this accession satisfied the
26-nucleotide structural length requirement, so the conservative reference
also contains 234,382,218 pairs. All tools retained exactly that complete set:

| Tool | Emitted | Intersection | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|
| seqproc | 234,382,218 | 234,382,218 | 1.000000 | 1.000000 | 1.000000 |
| matchbox | 234,382,218 | 234,382,218 | 1.000000 | 1.000000 | 1.000000 |
| splitcode | 234,382,218 | 234,382,218 | 1.000000 | 1.000000 | 1.000000 |

Every pairwise intersection is 234,382,218 and every pairwise Jaccard index is
1.0. Output-contract audits also agree: all R1 products are 26 nt and all R2
products are 98 nt. This block is a useful exact-equivalence and output-shape
control, but this particular input is not a discriminative accuracy challenge
because no pair fails the sole retention criterion.

### Accuracy-run resource diagnostics (not final performance numbers)

| Tool | Wall time | Peak RSS |
|---|---:|---:|
| seqproc | 43.562 s | 76.0 MiB |
| matchbox | 601.264 s | 40.1 MiB |
| splitcode | 88.674 s | 1,234.8 MiB |

## Newly completed: sci-RNA-seq3 final performance

The final timing block used uncompressed inputs and directed requested
biological sequence output to `/dev/null`. Values are mean wall time +/- sample
standard deviation over three randomized full-data replicates.

| Threads | seqproc | splitcode | matchbox |
|---:|---:|---:|---:|
| 1 | 16.118 +/- 0.029 s | 26.966 +/- 0.226 s | 233.997 +/- 2.616 s |
| 4 | 6.423 +/- 0.150 s | 12.881 +/- 0.000 s | 88.241 +/- 0.958 s |
| 16 | 6.456 +/- 0.029 s | 11.664 +/- 0.031 s | 72.204 +/- 0.655 s |
| 32 | 7.256 +/- 0.838 s | 10.379 +/- 0.050 s | 75.156 +/- 0.834 s |

seqproc is fastest at every thread count. At one thread it is 1.67-fold faster
than splitcode and 14.52-fold faster than Matchbox; at 32 threads those ratios
are 1.43-fold and 10.36-fold. seqproc reaches its best mean at four threads,
consistent with an I/O or memory-bandwidth ceiling on this uncompressed input.
Maximum observed 32-thread RSS was 73.0 MiB for seqproc, 1,833.7 MiB for
splitcode, and 3,104.4 MiB for Matchbox.

These Matchbox timings use the corrected variable-prefix expression that
accepts both documented 9- and 10-nt BC1 geometries. Timings from the
superseded shadowing two-branch script are not publication results.

## Newly completed: LR-SPLiT-seq forward-only final performance

This supplementary controlled comparison uses the same complete-component,
canonical-list constraints as the primary dual-orientation block. Matchbox is
boundary-safe exact, splitcode linker matching is substitution-only, and all
requested biological outputs are directed to `/dev/null`.

| Threads | seqproc | splitcode | matchbox |
|---:|---:|---:|---:|
| 1 | 14.650 +/- 0.058 s | 221.907 +/- 1.950 s | 3,075.491 +/- 4.018 s |
| 4 | 3.753 +/- 0.029 s | 57.324 +/- 0.603 s | 780.362 +/- 2.330 s |
| 16 | 2.801 +/- 0.029 s | 16.554 +/- 0.307 s | 206.140 +/- 1.468 s |
| 32 | 3.118 +/- 0.050 s | 8.610 +/- 0.175 s | 106.246 +/- 0.347 s |

seqproc is fastest at every thread count. At one thread it is 15.15-fold faster
than splitcode and 209.93-fold faster than Matchbox; at 32 threads those ratios
are 2.76-fold and 34.07-fold. Maximum observed 32-thread RSS was 71.0 MiB for
seqproc, 532.9 MiB for splitcode, and 30.8 MiB for Matchbox.

The Matchbox result is intentionally very different from the old performance
table. The conservative-reference audit showed that native fuzzy adjacent
captures can move biological field boundaries, so the semantically aligned
configuration uses exact canonical barcode and linker lists. That corrected
configuration is highly reproducible and scales almost linearly, but takes
51.26 minutes at one thread and 106.25 seconds at 32 threads. The old fast
fuzzy configuration must not be used for either timing or accuracy claims.

## Frozen artifacts

- `publication_results/journal_rerun_2026-08-17/publication-core-correctness-t32.yaml`
- `publication_results/journal_rerun_2026-08-17/publication-core-correctness-t32.schedule.json`
- `publication_results/journal_rerun_2026-08-17/lr_splitseq_dual_accuracy_metrics.json`
- `publication_results/journal_rerun_2026-08-17/lr_splitseq_dual_accuracy_metrics.csv`
- `publication_results/journal_rerun_2026-08-17/lr_splitseq_dual_pairwise_metrics.csv`
- `publication_results/journal_rerun_2026-08-17/lr_matchbox_hamming_expansion_sensitivity.json`
- `publication_results/journal_rerun_2026-08-17/publication-core-correctness-t32-r2.yaml`
- `publication_results/journal_rerun_2026-08-17/publication-core-correctness-t32-r2.schedule.json`
- `publication_results/journal_rerun_2026-08-17/scirnaseq3_accuracy_metrics.json`
- `publication_results/journal_rerun_2026-08-17/scirnaseq3_accuracy_metrics.csv`
- `publication_results/journal_rerun_2026-08-17/scirnaseq3_pairwise_metrics.csv`
- `publication_results/journal_rerun_2026-08-17/splitseq_pe_accuracy_metrics.json`
- `publication_results/journal_rerun_2026-08-17/splitseq_pe_accuracy_metrics.csv`
- `publication_results/journal_rerun_2026-08-17/splitseq_pe_pairwise_metrics.csv`
- `publication_results/journal_rerun_2026-08-17/tenx_v2_accuracy_metrics.json`
- `publication_results/journal_rerun_2026-08-17/tenx_v2_accuracy_metrics.csv`
- `publication_results/journal_rerun_2026-08-17/tenx_v2_pairwise_metrics.csv`
- `publication_results/journal_performance_2026-08-17/frozen_manifest.yaml`
- `publication_results/journal_performance_2026-08-17/frozen_schedule.json`
- `publication_results/journal_performance_2026-08-17/execution-log.jsonl`
- `publication_results/journal_performance_2026-08-17/scirnaseq3/`
- `publication_results/journal_performance_2026-08-17/lr_splitseq_forward/`

Large materialized FASTQs and compact accession bitmaps remain in the external
campaign directory recorded by the metrics artifact.
