# Phase 2A downstream biological validation (seqproc, splitcode, matchbox)

End-to-end, all three tools emit **symmetric** SPLiT-seq barcode reads (UMI 10 + bc3 8 + bc2 8 +
bc1 8 = 34bp, observed barcodes), then the **same** STARsolo `CB_UMI_Complex` config quantifies all
three so barcode error-correction is identical. Downstream we compare count matrices and biology
(cell-calling, clustering, mouse-brain cell typing), plus runtime and peak memory across the tools.

## One-time setup
```bash
# 1. python env (scanpy + clustering + notebook tooling)
bash biological_analysis/setup_env.sh

# 2. STAR index for the genome (cluster uses full GRCm38). Skip if you already have one.
#    --genomeSAsparseD 2 keeps it ~14GB / ~16GB RAM; drop it if you have >32GB RAM.
STAR --runMode genomeGenerate --genomeDir <INDEX_DIR> \
     --genomeFastaFiles <GRCm38.primary_assembly.fa> \
     --sjdbGTFfile <GRCm38.102.gtf> --sjdbOverhang <readlen-1> \
     --genomeSAsparseD 2 --runThreadN 8
```

## Run (one command)
```bash
biological_analysis/run_phase2a.sh \
  --r1 <cDNA_R1.fastq> --r2 <barcode_R2.fastq> \
  --genome <INDEX_DIR> --outdir <OUT> --threads 8
```
`R1` is the cDNA read, `R2` is the 94bp barcode read. Binaries are found via the repo layout or the
`SEQPROC_BIN`, `SPLITCODE_BIN`, and `MATCHBOX_BIN` env vars.

## Outputs (and what to copy back)
- `<OUT>/analysis/count_concordance.png` shows the barcode-rank knee and the per-barcode and per-gene
  correlations (Pearson on log1p, plus Spearman).
- `<OUT>/analysis/biological_analysis.png` shows the joint-embedding UMAPs, cell-type fractions, and
  the robust concordance metrics.
- `<OUT>/analysis/resource_usage.png` shows read-processing runtime and peak memory per tool, with
  `resource_table.md` and `resources.csv` carrying the numbers.
- `<OUT>/analysis/biological_metrics.json` and `count_concordance.json` hold the metrics.
- `<OUT>/phase2a_bundle.tar.gz` bundles all of the above. Copy it back to the box and open the PNGs
  or the notebook.

## Notebook
`notebooks/phase2a_tool_concordance.ipynb` reproduces the count-level concordance interactively
(point `RES` at `<OUT>`); run it with the venv kernel from `setup_env.sh`.

## Validated on
Sandbox box, full GRCm38, SRR6750041 1M-read subsample. The cluster run is the same command on the
full dataset; it additionally enables the split-pipe vendor comparison (vendor matrix present on the
cluster).
