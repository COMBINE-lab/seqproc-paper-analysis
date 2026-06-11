# Full Dataset Benchmark Walkthrough

Step-by-step guide for running the complete paper benchmarks on full SRA
datasets instead of the 1M-read subsets.

## Prerequisites

1. **Tool binaries** must be compiled and accessible:
   ```bash
   export SEQPROC_BIN=/path/to/seqproc/target/release/seqproc
   export MATCHBOX_BIN=/path/to/matchbox/target/release/matchbox
   export SPLITCODE_BIN=/path/to/splitcode/build/src/splitcode
   ```

2. **Python environment**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Server requirements** (for full datasets):
   - RAM: 32 GB minimum (splitcode on SPLiT-seq PE uses ~250 MB,
     but matchbox can use 300+ MB)
   - Disk: ~50 GB free for FASTQ data + intermediate outputs
   - CPU: 4+ cores recommended

## Step 1: Check data availability

```bash
python scripts/data_config.py --reads full
```

This prints which datasets are present and which need downloading.

## Step 2: Download missing datasets

For any datasets marked `[MISSING]`, download from SRA:

```bash
cd data/

# SPLiT-seq PE (paired-end, ~20 GB)
fasterq-dump --split-files SRR6750041
mv SRR6750041_1.fastq SRR6750041_R1.fastq
mv SRR6750041_2.fastq SRR6750041_R2.fastq

# 10x Chromium v2 (paired-end, ~10 GB)
fasterq-dump --split-files SRR8315379
mkdir -p 10x_short
mv SRR8315379_1.fastq 10x_short/SRR8315379_R1.fastq
mv SRR8315379_2.fastq 10x_short/SRR8315379_R2.fastq

# LR-SPLiT-seq -- already present as SRR13948564_full.fastq
# sci-RNA-seq3 -- already present as SRR7827254_{1,2}.fastq
```

Verify all data is available:
```bash
python scripts/data_config.py --reads full
# All four datasets should show [OK]
```

## Step 3: Run the full pipeline (single command)

```bash
./scripts/run_all.sh --reads full --threads 4 --replicates 3
```

This runs all five steps in sequence:
1. Performance benchmarks (Table 2) -- 3 replicates per tool per dataset
2. Concordance analysis (Figures 3-5) -- pairwise Jaccard, hamming vs edit
3. Discordant read validation (Supp. Figure S1) -- structural validation
4. Figure generation -- all publication PDFs and PNGs

**Estimated runtime** (4 threads, full data):
- SPLiT-seq PE (~87M reads): ~30-60 min per tool per replicate
- LR-SPLiT-seq (~4.2M reads): ~2-5 min per tool per replicate
- 10x Short (~56M reads): ~10-30 min per tool per replicate
- sci-RNA-seq3 (~10M reads): ~5-15 min per tool per replicate
- Total: approximately 3-6 hours

## Step 4: Verify results

After the pipeline completes, check the outputs:

```bash
# Results JSON
cat results/paper_figures/benchmark_results.json | python -m json.tool | head -50

# Concordance JSON
cat results/concordance/concordance_results.json | python -m json.tool | head -50

# Figures
ls -la results/paper_figures/*.pdf
```

## Step 5: Copy figures to paper

```bash
cp results/paper_figures/*.pdf ../paper/Figures/
```

## Running with 1M subsets (default)

For quick iteration or testing, use the 1M-read subsets:

```bash
./scripts/run_all.sh                          # defaults to --reads 1m
./scripts/run_all.sh --reads 1m --threads 4   # explicit
```

## Running individual scripts

Each script also accepts the `--reads` flag independently:

```bash
# Performance benchmarks only
python scripts/run_paper_benchmarks.py --threads 4 --replicates 3 --reads full

# Concordance analysis only
python scripts/concordance_analysis.py --threads 4 --reads full

# Discordant analysis (uses cached concordance results, no --reads flag needed)
python scripts/discordant_analysis.py

# Figure generation (reads from JSON, no --reads flag needed)
python scripts/generate_figures.py
```

## Datasets

| Dataset | SRA Accession | Full Reads | 1M Subset |
|---------|---------------|-----------|-----------|
| SPLiT-seq PE | SRR6750041 | 86,820,578 | 1,000,000 |
| LR-SPLiT-seq | SRR13948564 | 4,229,250 | 1,000,000 |
| 10x Chromium v2 | SRR8315379 | 56,514,800 | 1,000,000 |
| sci-RNA-seq3 | SRR7827254 | 10,177,866 | 1,000,000 |

## Troubleshooting

- **splitcode OOM**: If splitcode runs out of memory on SPLiT-seq PE full data,
  the config already uses `dist 3` (not `dist 6`) to keep memory under control.
- **matchbox slow**: matchbox on SPLiT-seq PE is CPU-intensive (~95s per 1M reads).
  On full data (~87M reads) this can take several hours.
- **Missing tool binaries**: Set `SEQPROC_BIN`, `MATCHBOX_BIN`, `SPLITCODE_BIN`
  environment variables to point to the correct binaries.
- **Data files not found**: Run `python scripts/data_config.py --reads full` to
  see exactly which files are missing and their expected paths.
