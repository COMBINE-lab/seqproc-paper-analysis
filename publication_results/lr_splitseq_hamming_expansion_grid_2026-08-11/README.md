# LR-SPLiT-seq Hamming-expansion diagnostic grid (Aug 11)

Status: diagnostic sensitivity analysis, not a proposed user-facing workflow
and not a replacement for the no-preprocessing putative Aug. 11 configuration.

This grid asks what explicit Hamming-distance-one whitelist expansion buys
matchbox and splitcode. It deliberately performs the pre-preprocessing that we
do not recommend requiring from users, so that limitations of tool semantics
can be separated from limitations of the available configurations.

## Expansion semantics

The generator `scripts/build_aug11_hamming_expansion_grid.py` reads the
canonical seqproc whitelists, generates every sequence within Hamming distance
one, sorts and deduplicates the result, and writes tool inputs deterministically.

- BC2/BC3: 96 unique canonical 8-mers become 2,400 unique variants; no
  expanded sequence has multiple canonical owners.
- BC1: 96 unique canonical 6-mers become 1,614 unique variants; 206 expanded
  sequences have multiple canonical owners, with maximum multiplicity three.
- A collision is accepted once. This reproduces permissive whitelist filtering,
  not unique barcode correction. The complete collision map is in
  `expanded_whitelist_metadata.json`.

The generated files are diagnostic assets, not required seqproc inputs:

- `configs/diagnostics/aug11_lr_splitseq_bc23_ham1.csv`
- `configs/diagnostics/aug11_lr_splitseq_bc1_ham1.csv`
- `configs/splitcode/diagnostic_aug11_lr_splitseq_expanded_exact.txt`

## Results

All rows use one execution at 32 threads on the same AMD EPYC 9555 node and
the same 5,764,421-read input. Accuracy is evaluated against the repaired
588,233-read conservative structural reference, not biological ground truth.
Times and RSS are single-run diagnostics. The seqproc row is the existing
Aug. 11 comparator and was not rerun for this grid.

| Tool/configuration | Calls | Precision | Recall | F1 | Wall time | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|
| seqproc native Hamming 1, linkers e5/e4 | 553,840 | 1.0000 | 0.9415 | 0.9699 | 3.90 s | 108.0 MiB |
| matchbox expanded barcodes, linkers e5/e4 | 695,373 | 0.4308 | 0.5093 | 0.4668 | 1:15:26 | 520.6 MiB |
| matchbox expanded barcodes, exact linkers | 406,257 | 1.0000 | 0.6906 | 0.8170 | 1:11:48 | 434.1 MiB |
| splitcode expanded exact barcodes, exact gaps | 309,601 | 0.9962 | 0.5243 | 0.6870 | 6.79 s | 533.5 MiB |

### Annotations

1. **The appropriate matchbox sensitivity control is expansion plus exact
   linkers.** It retains perfect nominal precision, raises recall from 0.0141
   in the no-expansion exact configuration to 0.6906, and emits only the
   intended 32-nt projection. It overlaps 405,850 seqproc calls, has 407
   matchbox-only calls, and misses 147,990 seqproc calls.

2. **Approximate matchbox linkers remain semantically unsafe.** With the same
   expanded barcode sets and 5/4-edit linker budgets, 237,655 of 695,373
   outputs are 33--47 nt rather than 32 nt. Precision falls to 0.4308. The
   matcher is selecting alternate linker endpoints, not merely allowing
   substitutions/indels while preserving cassette boundaries.

3. **Explicit expansion is operationally poor for matchbox.** Both matchbox
   runs require over 71 minutes and 434--521 MiB, versus 3.90 seconds and
   108 MiB for seqproc's native Hamming filtering. Expansion also requires an
   extra user-visible preprocessing/normalization step with a collision policy.

4. **Expansion buys splitcode nothing.** Its forward and reverse bitmap files
   are byte-identical to the earlier native Hamming-one, exact-`0-0` diagnostic.
   The 0.5243 recall limitation is therefore not caused by failure to enumerate
   Hamming neighbors. The result is consistent with its single forward-chain
   selection consuming candidate placements differently from seqproc's complete
   adjacent-cassette search.

5. **This does not change the proposed primary comparison.** The expansion
   grid is useful as a sensitivity analysis showing what competitors can do
   with substantial manual assistance. It should not be substituted for the
   no-preprocessing, qualitatively aligned configurations without explicitly
   disclosing the additional whitelist construction and ambiguity handling.

Machine-readable values, hashes, annotations, timing details, and raw campaign
paths are in `results.json`; the evaluator's complete pairwise matrix is in
`accuracy-grid.json`. Raw output bitmaps and time reports are under
`/scratch1/seqproc-ecosystem/lr-splitseq-hamming-expansion-grid-2026-08-11`.
