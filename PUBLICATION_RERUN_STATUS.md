# Publication rerun status

Last updated: 2026-08-17

This is the active checklist for results computed with the revised structural
references and semantically aligned publication configurations. Historical
results produced with preprint-era or August diagnostic configurations do not
count as complete here.

## Completion matrix

| Technology block | Structural reference | Accuracy: 1 run/tool, 32 threads | Final speed/RSS: 3 replicates at 1/4/16/32 threads |
|---|---|---|---|
| SPLiT-seq PE, primary | Complete | Pending | Pending |
| LR-SPLiT-seq dual orientation, primary | Complete | **Complete** | Pending |
| LR-SPLiT-seq forward-only, supplementary | Definition/configs complete; forward-only reference subset pending | Pending | Pending |
| 10x Chromium v2, primary | Complete | Pending | Pending |
| sci-RNA-seq3, primary | Complete | Pending | Pending |

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

## Frozen artifacts

- `publication_results/journal_rerun_2026-08-17/publication-core-correctness-t32.yaml`
- `publication_results/journal_rerun_2026-08-17/publication-core-correctness-t32.schedule.json`
- `publication_results/journal_rerun_2026-08-17/lr_splitseq_dual_accuracy_metrics.json`
- `publication_results/journal_rerun_2026-08-17/lr_splitseq_dual_accuracy_metrics.csv`
- `publication_results/journal_rerun_2026-08-17/lr_splitseq_dual_pairwise_metrics.csv`
- `publication_results/journal_rerun_2026-08-17/lr_matchbox_hamming_expansion_sensitivity.json`

Large materialized FASTQs and compact accession bitmaps remain in the external
campaign directory recorded by the metrics artifact.
