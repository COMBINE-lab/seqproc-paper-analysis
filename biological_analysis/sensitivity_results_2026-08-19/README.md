# SPLiT-seq PE downstream sensitivities (2026-08-19)

These compact results use the final publication-campaign FASTQ sets, STAR
2.7.11b, GRCm38/Ensembl release 102, 32 threads, no STAR cell filter, and a
fixed 200-UMI analysis threshold. They separate two questions:

1. Does STARsolo's barcode-correction mode place the same burden on each
   upstream product?
2. Does an externally generated radius-one barcode expansion materially rescue
   Matchbox's exact, boundary-safe PE configuration?

The primary exact/`1MM` result is archived separately in
`../final_run_results_2026-08-19/`. The three subdirectories here are
sensitivities and must not be silently substituted into the upstream accuracy
or performance tables.

## Results

| Matchbox input | STAR correction | Matchbox reads | STAR-valid Matchbox barcodes | Called cells (seqproc/splitcode/Matchbox) | All-tool cell-type agreement | Per-gene Pearson, seqproc/Matchbox | Per-gene Pearson, splitcode/Matchbox |
|---|---|---:|---:|---:|---:|---:|---:|
| Exact (primary) | `1MM` | 35,160,366 | 99.9747% | 225 / 220 / 211 | 0.8673 | 0.9860 | 0.9885 |
| Exact (primary) | `EditDist_2` | 35,160,366 | 99.9747% | 225 / 220 / 211 | 0.8673 | 0.9856 | 0.9885 |
| Radius-1 expanded | `1MM` | 39,758,892 | 99.2976% | 225 / 220 / 215 | 0.9023 | 0.9889 | 0.9919 |
| Radius-1 expanded | `EditDist_2` | 39,758,892 | 99.9744% | 225 / 220 / 217 | 0.8848 | 0.9888 | 0.9922 |

For the exact final read sets, independent per-piece correction raises
seqproc's STAR-valid barcode fraction from 99.1886% to 99.9731% and improves
seqproc--splitcode per-barcode Pearson correlation from 0.9722 to 0.9809. Cell
counts and the all-tool cell-type agreement are unchanged. The non-monotonic
cell-type agreement in the expanded rows reflects threshold/marker sensitivity
among low-signal cell types; count correlations are the more stable comparison.

## Expanded Matchbox structural-reference result

The 96 eight-base barcodes have minimum Hamming distance four. Their radius-one
expansion therefore contains 2,400 distinct observations with no multiply
owned variants. With exact linkers, expanded Matchbox is a strict superset of
exact Matchbox and a strict subset of both other tools:

- emitted: 39,758,892 (51.22% of input), +4,598,526 over exact Matchbox;
- conservative-reference intersection: 39,758,892;
- precision: 1.0000; recall: 0.6922; F1: 0.8181;
- exact Matchbox recall/F1: 0.6121/0.7594;
- remaining gap to seqproc: 17,232,489 reads.

The materialized-output diagnostic took 5,454 seconds (1:30:54) and peaked at
9,067,304 KiB RSS. The corresponding exact Matchbox materialization took
281.69 seconds and 7,937,164 KiB in the final campaign. Thus the external
expansion was about 19.4-fold slower for 13.1% more retained reads. These are
diagnostic materialization measurements, not replacements for the paper's
randomized `/dev/null` performance grid.

## Interpretation

- Keep exact, boundary-safe Matchbox as the primary configuration unless the
  expansion is also promoted and rerun throughout the accuracy/performance
  campaign. Report expansion as a capability/workaround sensitivity.
- Prefer `EditDist_2` for the controlled downstream comparison. In STARsolo
  `CB_UMI_Complex`, `1MM` rejects a read when more than one barcode piece needs
  correction. splitcode has already canonicalized its extracted tag sequences,
  Matchbox exact is canonical by construction, and seqproc retains observed
  bases; independent correction removes that unequal downstream burden without
  changing upstream set membership.
- Do not compare STAR's CB+UMI Q30 values as if they were equivalent upstream
  measurements. splitcode retains the observed UMI qualities but assigns `K`
  qualities to its extracted/canonicalized barcode fields; seqproc and Matchbox
  retain observed qualities. Neither `1MM` nor `EditDist_2` uses those qualities
  to resolve these minimum-distance-four barcode matches.
