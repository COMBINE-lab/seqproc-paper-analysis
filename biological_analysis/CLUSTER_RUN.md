# Full downstream rerun

Run on a dedicated node with at least 48 GB RAM. STARsolo uses approximately
30 GB with the full GRCm38 index. The final campaign uses 32 physical cores.

## One-time setup

```bash
cd /path/to/seqproc-paper-analysis
bash biological_analysis/setup_env.sh
bash biological_analysis/reference/prepare_reference.sh
```

`prepare_reference.sh` downloads and verifies the GRCm38 primary assembly and
Ensembl release-102 GTF, builds STAR 2.7.11b from its tagged source archive,
generates the index with `sjdbOverhang=65`, and records SHA-256 provenance for
all inputs and generated index files.

## Materialize the final splitcode pair

The accuracy campaign retained splitcode's extracted barcode product but sent
its ordinary read output to `/dev/null`. Re-run the identical frozen
configuration once while retaining R1:

```bash
biological_analysis/materialize_splitcode_downstream.sh \
  --r1 /path/to/SRR6750041_R1.fastq \
  --r2 /path/to/SRR6750041_R2.fastq \
  --outdir "$WORK/splitcode-final-downstream" --threads 32
```

The materializer validates all paired IDs and writes checksums. Its generated
34-nt barcode FASTQ must be byte-identical to the final campaign artifact.

## Quantify and analyze

```bash
BMAP=/path/to/final-campaign/aggregates/accuracy/bitmaps
biological_analysis/run_current_downstream.sh \
  --genome biological_analysis/reference/star_GRCm38_ensembl102 \
  --outdir "$WORK/downstream-final" \
  --seqproc-r1 <final-seqproc-R1> --seqproc-r2 <final-seqproc-R2> \
  --splitcode-r1 "$WORK/splitcode-final-downstream/sc_cdna.fastq" \
  --splitcode-r2 "$WORK/splitcode-final-downstream/umi_bc3_bc2_bc1.fastq" \
  --matchbox-r1 <final-matchbox-R1> --matchbox-r2 <final-matchbox-R2> \
  --seqproc-bitmap "$BMAP/splitseq_pe.seqproc.accepted.1.bitmap" \
  --splitcode-bitmap "$BMAP/splitseq_pe.splitcode.accepted.1.bitmap" \
  --matchbox-bitmap "$BMAP/splitseq_pe.matchbox.accepted.1.bitmap" \
  --threads 32 --min-umi 200 --cb-match 1MM
```

The driver verifies every FASTQ before STARsolo and produces a compact
`downstream_bundle.tar.gz`. The large STAR matrices remain in the output
directory and are not committed.

The command above exactly reproduces the historical STARsolo `1MM` matching
choice on the final upstream outputs. For a controlled correction sensitivity,
rerun to a separate output directory with `--cb-match EditDist_2`. STARsolo's
complex `1MM` mode permits a mismatch in only one barcode piece, whereas
`EditDist_2` corrects the three pieces independently. This matters because
seqproc retains observed barcode bases, while splitcode's extracted tag output
is already canonicalized and Matchbox's primary PE output is exact.

After materializing the expanded-whitelist Matchbox FASTQs, reproduce the
separate sensitivity analysis with:

```bash
biological_analysis/run_matchbox_ham1_sensitivity.sh \
  --genome biological_analysis/reference/star_GRCm38_ensembl102 \
  --primary-run "$WORK/downstream-final" \
  --outdir "$WORK/downstream-matchbox-ham1" \
  --matchbox-r1 <expanded-matchbox-R1> --matchbox-r2 <expanded-matchbox-R2> \
  --reference-bitmap "$BMAP/splitseq_pe.reference.raw" \
  --seqproc-bitmap "$BMAP/splitseq_pe.seqproc.accepted.1.bitmap" \
  --splitcode-bitmap "$BMAP/splitseq_pe.splitcode.accepted.1.bitmap" \
  --matchbox-exact-bitmap "$BMAP/splitseq_pe.matchbox.accepted.1.bitmap" \
  --threads 32 --min-umi 200 --cb-match 1MM
```

The script refuses to reuse primary matrices generated with a different
barcode-matching mode.

## Refresh analysis without realignment

```bash
PY=biological_analysis/.venv_downstream/bin/python
OUT="$WORK/downstream-final"
ARGS=(seqproc:$OUT/sp_Solo.out/Gene splitcode:$OUT/sc_Solo.out/Gene matchbox:$OUT/mb_Solo.out/Gene)
$PY biological_analysis/scripts/biological_analysis.py "$OUT/analysis" 200 "${ARGS[@]}"
$PY biological_analysis/scripts/count_concordance.py "$OUT/analysis" "${ARGS[@]}"
$PY biological_analysis/scripts/jaccard_confusion.py "$OUT/analysis" 200 "${ARGS[@]}"
```

The barcode-rank inflection comes from the committed Python port of
`DropletUtils::barcodeRanks`; no `kneed` dependency or ad hoc manual step is
used.
