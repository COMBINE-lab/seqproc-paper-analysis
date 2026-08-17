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
| LR-SPLiT-seq forward-only, supplementary | Definition/configs complete; forward-only reference subset pending | Pending | Pending |
| 10x Chromium v2, primary | Complete | **Complete** | Pending |
| sci-RNA-seq3, primary | Complete | **Complete** | Pending |

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

Large materialized FASTQs and compact accession bitmaps remain in the external
campaign directory recorded by the metrics artifact.
