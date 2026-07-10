#!/usr/bin/env bash
# Build a chr19-only STAR index for box-side Phase 2A pipeline development.
# The real cluster run uses a full GRCm38 index instead. chr19 is too small
# for trustworthy cell calling and is only for proving the pipeline runs.
set -euo pipefail
FASTA=${1:?usage: build_chr19_index.sh <chr19.fa> <chr19.gtf> [outdir]}
GTF=${2:?usage: build_chr19_index.sh <chr19.fa> <chr19.gtf> [outdir]}
OUT=${3:-biological_analysis/refs/star_chr19}
mkdir -p "$OUT"
# chr19 ~61 Mb, so genomeSAindexNbases is reduced per STAR's small-genome rule.
STAR --runMode genomeGenerate \
     --genomeDir "$OUT" \
     --genomeFastaFiles "$FASTA" \
     --sjdbGTFfile "$GTF" \
     --sjdbOverhang 99 \
     --genomeSAindexNbases 11 \
     --runThreadN 8
echo "chr19 STAR index built at $OUT"
