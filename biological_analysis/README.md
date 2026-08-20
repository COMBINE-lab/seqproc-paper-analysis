# Downstream biological validation

This directory reproduces the SPLiT-seq PE downstream comparison from the
final publication-campaign read sets. Each tool supplies a paired 66-nt cDNA
FASTQ and a 34-nt `UMI10+BC3_8+BC2_8+BC1_8` FASTQ. One pinned STARsolo
configuration then generates count matrices for cell filtering, clustering,
marker-based cell typing, and count-concordance analysis.

The archived outputs in `full_run_results/` were generated from the preprint
configurations and are retained only for history. They must not be combined
with the final campaign accuracy or performance results.

The curated, compact outputs from the final-read-set rerun are in
`final_run_results_2026-08-19/`. Large FASTQs, STAR matrices, and indices are
gitignored but are fully described by the committed provenance artifacts.
The four controlled correction/Matchbox-expansion comparisons are summarized
in `sensitivity_results_2026-08-19/README.md` with compact machine-readable
artifacts for each configuration.

## Reproduction

```bash
# Locked Python 3.12 environment.
bash biological_analysis/setup_env.sh

# STAR 2.7.11b plus checksummed GRCm38/Ensembl-release-102 reference.
bash biological_analysis/reference/prepare_reference.sh

# Recreate splitcode's accepted cDNA mate. The final accuracy run retained its
# extracted barcode FASTQ but intentionally sent ordinary output to null.
biological_analysis/materialize_splitcode_downstream.sh \
  --r1 <RAW_R1.fastq> --r2 <RAW_R2.fastq> --outdir <SPLITCODE_PAIR> --threads 32

# Quantify the three final pairs and regenerate all downstream outputs.
biological_analysis/run_current_downstream.sh \
  --genome biological_analysis/reference/star_GRCm38_ensembl102 \
  --outdir <OUT> \
  --seqproc-r1 <SEQPROC_R1> --seqproc-r2 <SEQPROC_R2> \
  --splitcode-r1 <SPLITCODE_R1> --splitcode-r2 <SPLITCODE_R2> \
  --matchbox-r1 <MATCHBOX_R1> --matchbox-r2 <MATCHBOX_R2> \
  --threads 32 --min-umi 200 --cb-match EditDist_2
```

Optional `--seqproc-bitmap`, `--splitcode-bitmap`, and `--matchbox-bitmap`
arguments regenerate read-set Jaccard values from the compact final-campaign
bitmaps without loading tens of millions of identifiers into Python sets.

## Fixed analysis choices

- STAR 2.7.11b.
- Mus musculus GRCm38 primary assembly and Ensembl release-102 GTF.
- `sjdbOverhang=65` for the 66-nt cDNA read.
- STARsolo `CB_UMI_Complex`, UMI position `0_0_0_9`, and barcode positions
  `0_10_0_17`, `0_18_0_25`, and `0_26_0_33`.
- The canonical 96-member whitelist for all three barcode rounds and
  `--soloCBmatchWLtype EditDist_2` for the primary analysis. The historical
  exact-whitelist `1MM` result is retained as a sensitivity analysis.
- No STARsolo cell filter; analysis calls cells at at least 200 total UMIs.
- Random seed 0 for PCA, neighbors, Leiden, and marker scoring.

The tools retain their actual quality behavior. seqproc and Matchbox retain
observed qualities, while splitcode assigns synthetic qualities to extracted
barcode fields. This is a tool capability difference, not normalized away by
the workflow. Under STARsolo's `CB_UMI_Complex`/`1MM` implementation, a unique
one-mismatch whitelist hit is accepted without quality-based scoring. The
synthetic values therefore affect the reported barcode Q30 statistic but not
correction for this minimum-distance-four whitelist.

The exact-whitelist `1MM` result is a STARsolo sensitivity worth distinguishing
from the upstream-tool comparison. In `CB_UMI_Complex` mode, `1MM` rejects a read when
more than one of the three barcode pieces requires correction. splitcode's
extracted tag sequences are already canonical, Matchbox's primary PE result is
exact, and seqproc deliberately retains the observed sequence after filtering;
therefore this rule places different correction burdens on the final products.
The driver accepts `--cb-match EditDist_2`, which corrects each barcode piece
independently. Because all upstream configurations already constrain accepted
barcode observations to Hamming distance at most one (or exact), that setting
does not broaden upstream read acceptance and is the appropriate controlled
primary downstream configuration. It also does not use quality-based
tie-breaking.

## Matchbox expanded-whitelist sensitivity

`configs/matchbox/sensitivity_splitseq_pe_ham1_expanded.mb` is deliberately
not the primary benchmark configuration. It matches an automatically
generated, exact radius-one expansion of the 96-member barcode whitelist.
The barcode code has minimum Hamming distance four, so all 2,400 expanded
sequences have one owner. The workaround avoids Matchbox's approximate-field
boundary shifts, but requires external configuration generation and retains
exact linker matching. If promoted to a primary configuration, Matchbox's
accuracy and performance benchmarks must be rerun as well.
`run_matchbox_ham1_sensitivity.sh` regenerates its accession accuracy,
STARsolo matrix, downstream concordance, and provenance while reusing the
primary seqproc and splitcode matrices for the same STAR matching mode.

## Outputs

The current driver writes:

- validated FASTQ counts, lengths, paired identifiers, and SHA-256 digests;
- STAR logs, summaries, timing, and peak RSS for every tool;
- `analysis/biological_metrics.json`;
- `analysis/count_concordance.json` and its PNG/PDF figure;
- `analysis/jaccard_confusion.json`;
- optional `analysis/read_set_jaccard.json` from campaign bitmaps;
- `downstream_provenance.json`, including the reference-index manifest and
  exact Python environment; and
- `downstream_bundle.tar.gz` with the compact reviewable artifacts.

## Regenerate manuscript tables and sensitivity figure

The manuscript tables are generated directly from the four frozen summary
JSON files rather than transcribed by hand:

```bash
PY=biological_analysis/.venv_downstream/bin/python

$PY biological_analysis/scripts/plot_downstream_sensitivities.py \
  --exact-editdist2 biological_analysis/sensitivity_results_2026-08-19/primary_exact_editdist2 \
  --exact-1mm biological_analysis/final_run_results_2026-08-19 \
  --expanded-editdist2 biological_analysis/sensitivity_results_2026-08-19/matchbox_ham1_editdist2 \
  --expanded-1mm biological_analysis/sensitivity_results_2026-08-19/matchbox_ham1_1mm \
  --output biological_analysis/sensitivity_results_2026-08-19/downstream_sensitivity_comparison

$PY biological_analysis/scripts/write_downstream_latex_tables.py \
  --exact-editdist2 biological_analysis/sensitivity_results_2026-08-19/primary_exact_editdist2 \
  --exact-1mm biological_analysis/final_run_results_2026-08-19 \
  --expanded-editdist2 biological_analysis/sensitivity_results_2026-08-19/matchbox_ham1_editdist2 \
  --expanded-1mm biological_analysis/sensitivity_results_2026-08-19/matchbox_ham1_1mm \
  --accuracy biological_analysis/sensitivity_results_2026-08-19/matchbox_ham1_editdist2/accuracy.json \
  --output-dir ../seqproc-paper/sections
```
