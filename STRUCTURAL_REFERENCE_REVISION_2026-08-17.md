# Revised configurations and conservative structural references (2026-08-17)

## Scope and interpretation

These retained sets are protocol-informed, conservative structural references;
they are not experimental ground truth. They answer whether a complete read has
a defensible protocol structure under a declared error model. Accuracy against
these sets should therefore be described as *agreement with the structural
reference*, accompanied by pairwise concordance and sensitivity analyses.

All tools are asked to do the same qualitative job: filter on the relevant
structure and, when the seqproc geometry projects biological components, write
those same components. Backend-specific error semantics that cannot be aligned
without preprocessing are reported as capability differences. No Hamming-
expanded or deduplicated whitelist is required from the user.

## Technology-by-technology changes from the original preprint

### Paired-end SPLiT-seq

The original configurations modeled R2 as a 2-nt skip plus an 8-nt UMI, used
`GTGGCCGCTGTTTCGCATCGGCGTACGACT` (C at the eighth position) for Linker1, and
used a 6-nt truncation of the round-one barcode. The documented physical layout
is instead:

`UMI(10) + BC3(8) + Linker1(30) + BC2(8) + Linker2(30) + BC1(8)`.

The revised seqproc, matchbox, and splitcode configurations use the 10-nt UMI,
the A-containing Linker1, and the full canonical 8-nt round-one list, and emit a
34-nt `UMI+BC3+BC2+BC1` projection. Splitcode now encodes the complete adjacent
component chain instead of independent windowed hits. Seqproc uses edit distance
3 for linkers and Hamming distance 1 for barcodes. Matchbox uses exact canonical
barcodes and linkers because its approximate regions can shift an adjacent
capture boundary even for an exact cassette; this is a material capability
difference, not silently compensated for by expanding inputs.

Verification now requires complete paired FASTQs with matching IDs, full
8-nt barcode membership within Hamming distance 1, the corrected linker
sequences, and adjacent components. The primary linker budget is 3/3. The
reference changes from 60,675,548 reads (78.17%) in the original paper to
57,437,503 of 77,621,181 pairs (73.9972%). These numbers are not a simple
recalculation: the UMI, linker, barcode length, whitelist, error budget, and
paired-integrity criteria all changed.

### Long-read SPLiT-seq

The original reference used a broad 6/6 linker budget, stopped after the first
successful orientation, and capped Linker1 candidate enumeration at eight. The
new primary core requires a complete adjacent cassette in either orientation,
linker edit distance at most 3/3, Hamming-distance-one membership in canonical
8-nt BC2/BC3 and 6-nt BC1 lists, and no semantic candidate cap. Both orientations
are always evaluated and reads accepted in both are reported explicitly.

Seqproc implements the full dual-orientation, whitelist-aware 3/3 geometry
natively. Matchbox processes both orientations natively but uses exact anchors
and canonical barcodes because fuzzy captures can move field boundaries.
Splitcode uses its canonical unexpanded lists, 3-mismatch linkers, an exact BC3
chain initiator, and Hamming-one BC2/BC1; the primary dual-orientation benchmark
is two sequential passes over forward and precomputed reverse-complement reads.
As agreed, duplicate reconciliation is discussed but is not added to the timed
splitcode workflow. Forward-only remains the controlled supplementary block.

The new core is 560,699 of 5,764,421 reads (9.7269%). Wider sensitivity sets are
573,518 at 5/4 (9.9493%) and 590,579 at 6/6 (10.2452%). The original paper's
6/6-style set was 601,603 (10.44%) and should no longer be used. In the new core,
1,793 reads pass in both orientations; 129,024 accepted reads have an ambiguous
six-base barcode correction and are retained as structurally compatible but
reported separately. The unique-correction subset contains 431,675 reads.

### 10x Chromium v2

The benchmarked operation remains intentionally narrow: reject a pair when R1
is shorter than 26 nt, otherwise pass both reads through. The original validator
required R1 to be *exactly* 26 nt and did not verify the mate. The revised
reference uses `length >= 26`, validates complete paired FASTQs and matching IDs,
and does not introduce a barcode whitelist into a length-only benchmark.

The seqproc, matchbox, and splitcode publication configurations remain a
length-filter/passthrough task; splitcode is not asked to perform an unused
extraction. All 234,382,218 pairs pass, as before, because this input's R1 reads
already satisfy the minimum. The definition is now correct for future inputs
with longer technical reads and detects paired-file corruption.

### sci-RNA-seq3

The original validator found one global best CAGAGC match in the first 18 bases.
A better decoy at an illegal offset could therefore hide a valid edit-one anchor
at the documented BC1 boundary. The revised validator evaluates offsets 9 and 10
independently, requires room for UMI(8)+BC2(10), selects a unique minimum-edit
placement, and rejects/reports equal-best offset ambiguity. Boundary indels are
inherently observationally equivalent to an exact anchor one base outside the
nominal boundary; that limitation is explicit rather than silently classified.

Seqproc no longer pads the variable 9/10-nt BC1 and now emits the natural
27/28-nt `BC1+BC2+UMI` projection. Matchbox retains edit-one anchors but guards
all component lengths before output, eliminating the previous 29/30-nt boundary-
shifted products. Splitcode already emitted 27/28 nt and is unchanged.

The revised reference contains 19,864,110 of 22,088,821 pairs (89.9283%), versus
19,767,975 (89.49%) in the original paper. It rejects 25,960 equal-best offset
ambiguities; 13,840,651 accepted reads have an exact anchor and 6,023,459 have an
edit-one anchor.

## Validator implementation and performance

`scripts/edit_tolerant_validity.py` is now the single reason-coded validator for
all four technologies. It uses strict binary FASTQ parsing, optional mate-ID and
record-count checks, owner-aware Hamming lookup tables, exact-structure fast
paths, ordered multiprocessing, and streamed ID output. Streaming removes the
previous need to retain tens or hundreds of millions of Python strings. Accepted
IDs are deterministic and written in input order. A separate finalizer verifies
their line counts, hashes outputs, and attaches input checksums and exact commands.

On 100,000 reads, LR throughput increased from 55,103 reads/s at one process to
286,222 at eight; PE increased from 171,126 to 663,290 reads/s. Full 32-process
runs on the dedicated benchmark node were:

| Reference | Reads scanned | Retained | Wall time | Throughput | Controller / max-worker RSS |
|---|---:|---:|---:|---:|---:|
| PE SPLiT-seq core | 77,621,181 | 57,437,503 | 106.77 s | 726,964/s | 43.8 / 26.7 MiB |
| LR SPLiT-seq core 3/3 | 5,764,421 | 560,699 | 12.59 s | 457,791/s | 47.9 / 52.4 MiB |
| LR sensitivity 5/4 | 5,764,421 | 573,518 | 13.07 s | 441,150/s | 46.2 / 50.0 MiB |
| LR sensitivity 6/6 | 5,764,421 | 590,579 | 12.55 s | 459,342/s | 46.2 / 50.0 MiB |
| sci-RNA-seq3 core | 22,088,821 | 19,864,110 | 29.19 s | 756,802/s | 35.6 / 19.0 MiB |
| 10x v2 core | 234,382,218 | 234,382,218 | 307.27 s | 762,784/s | 31.2 / 19.0 MiB |

The 10x run is dominated by reading and cross-checking both 20-GiB/51-GiB mate
files and writing 4.5 GiB of accepted IDs, rather than by the length predicate.

Machine-readable summaries and provenance are under
`publication_results/structural_reference_2026-08-17/`. The multi-gigabyte ID
lists remain local/ignored and should be archived separately with their recorded
SHA-256 digests.
