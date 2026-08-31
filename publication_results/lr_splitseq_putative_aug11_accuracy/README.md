# Putative new LR-SPLiT-seq configuration (Aug 11)

Status: putative configuration under accuracy evaluation; not a frozen
publication configuration.

This evaluation aligns the three tools as closely as their native interfaces
permit, without pre-expanding or deduplicating user whitelist files:

- seqproc: complete adjacent cassette, canonical Hamming-one whitelist
  filtering at BC3/BC2/BC1, 5/4-edit linker limits, native dual orientation,
  and a 32-nt UMI+BC3+BC2+BC1 projection.
- matchbox: canonical, unexpanded whitelists and native dual orientation.
  Barcode and linker terms are exact in the defensible no-expansion mode.
  Approximate terms can slide within their search window and select alternate
  cassette boundaries, so they cannot preserve the intended adjacency.
- splitcode: canonical, unexpanded whitelists and exact adjacency in two passes
  over forward and precomputed reverse-complement input, with three component
  projections. BC3 is an exact chain initiator; BC2 and BC1 allow Hamming one.
  Its existing three-edit linker limit is retained because splitcode
  materializes every allowed edit variant during initialization; the seqproc
  5/4 budgets are computationally impractical in this backend.

Configurations:

- `configs/seqproc/publication_lr_splitseq_dual_complete_adjacent_whitelist_e5_e4.geom`
- `configs/matchbox/putative_aug11_lr_splitseq_dual.mb`
- `configs/splitcode/putative_aug11_lr_splitseq_forward.txt`

All full-data accuracy runs use one execution at 32 threads. Timing and peak
RSS are recorded as useful diagnostics, not final timing-benchmark estimates.
The repaired structural reference remains a conservative reference rather
than biological ground truth.

The final result artifact is `results.json`. Diagnostic runs retained outside
this repository show why the stricter competitor modes are necessary:

- matchbox approximate barcodes plus 5/4-edit linkers: positional drift,
  precision 0.2671, recall 0.0972, F1 0.1425.
- matchbox exact barcodes plus 5/4-edit linkers: alternate linker endpoints,
  precision 0.1434, recall 0.0228, F1 0.0393.
- splitcode with unbounded `0` chain gaps: not exact adjacency, precision
  0.2521, recall 0.7975, F1 0.3831.
- splitcode with exact `0-0` gaps and Hamming-one BC3 initiators: random
  neighbor hits consume its single forward chain, precision 0.9962, recall
  0.5243, F1 0.6870.
