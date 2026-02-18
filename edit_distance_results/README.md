# Edit Distance Benchmark Results

This directory contains all scripts, configurations, and results for running the seqproc paper benchmarks with edit distance (Levenshtein) instead of Hamming distance.

## Overview

The benchmark compares three tools for barcode extraction from single-cell sequencing data:
- **seqproc**: Using the new edit distance feature
- **matchbox**: Using edit distance (configured with -e 0.2)
- **splitcode**: Using edit distance (configured with distance ratios)

## Prerequisites

### 1. Data Requirements

The benchmark expects the following data files in the `data/` directory (relative to the repository root):

```
data/
├── SRR6750041_1M_R1.fastq          # SPLiT-seq PE R1 (1M reads)
├── SRR6750041_1M_R2.fastq          # SPLiT-seq PE R2 (1M reads)
├── SRR13948564_1M.fastq           # SPLiT-seq Long Read (1M reads)
├── 10x_short/
│   └── SRR8315379_1M_R1.fastq     # 10x v2 Short Read R1 (1M reads)
│   └── SRR8315379_1M_R2.fastq     # 10x v2 Short Read R2 (1M reads)
├── SRR7827254_1M_1.fastq          # Sci-Seq 3 R1 (1M reads)
├── SRR7827254_1M_2.fastq          # Sci-Seq 3 R2 (1M reads)
├── 10x/
│   └── ERR9958134_1M.fastq        # 10x GridION (1M reads)
│   └── ERR9958135_1M.fastq        # 10x PromethION (1M reads)
└── 3M-february-2018.txt.gz        # 10x barcode whitelist
```

### 2. Tool Binaries

The benchmark expects the following binaries (or set environment variables):

```bash
export SEQPROC_BIN="/path/to/seqproc/target/release/seqproc"
export MATCHBOX_BIN="/path/to/matchbox/target/release/matchbox"  
export SPLITCODE_BIN="/path/to/splitcode/build/src/splitcode"
```

Default paths (relative to repository parent):
- `../combine-lab/seqproc/target/release/seqproc`
- `../matchbox/target/release/matchbox`
- `../splitcode/build/src/splitcode`

### 3. Python Dependencies

```bash
# On Ubuntu/Debian
sudo apt install python3-numpy python3-matplotlib

# Or with pip (may require virtual environment)
pip install numpy matplotlib
```

## Running the Benchmark

### Single Replicate (Quick Test)

```bash
cd /path/to/seqproc-paper-analysis-clean
python3 edit_distance_results/run_paper_benchmarks.py --threads 4 --replicates 1
```

### Full Benchmark (3 Replicates)

```bash
python3 edit_distance_results/run_paper_benchmarks.py --threads 4 --replicates 3
```

### Output

Results will be saved to `results/paper_figures/` (relative to repository root):
- `benchmark_results.json`: Raw results in JSON format
- `fig1_performance_distribution.png`: Performance distribution plot
- `fig2_recovery_table.png`: Recovery rate comparison table
- `fig3_summary_table.png`: Summary statistics table

## Key Changes from Original Benchmark

### 1. Seqproc Geometry Files

All seqproc geometry files have been updated to use edit distance:

- `splitseq_filter_edit.geom`: Uses `edit()` instead of `hamming()` for linkers
- `splitseq_replacement_edit.geom`: Uses `map_with_edit()` for barcode mapping
- `splitseq_singleend_primer_edit.geom`: Uses `edit()` for linker matching
- `sciseq3_edit.geom`: Uses `edit()` for anchor matching
- `10x_longread_fwd_edit.geom`: Forward orientation for 10x long reads
- `10x_longread_rev_edit.geom`: Reverse orientation for 10x long reads

### 2. Datasets Added

The following datasets were added to the benchmark:
- Sci-Seq 3: Single-cell RNA-seq with anchor-based barcoding
- 10x GridION Long-Read: Oxford Nanopore long reads with 10x chemistry
- 10x PromethION Long-Read: Oxford Nanopore long reads with 10x chemistry

### 3. Whitelist Files

Whitelist files are included in `configs/`:
- `splitseq_bc23_whitelist.txt`: Combined BC2 and BC3 whitelists (8bp barcodes)
- `splitseq_bc1_whitelist.txt`: BC1 whitelist (8bp barcodes)

## Expected Results

Based on the benchmark run with edit distance:

| Dataset | Tool | Runtime (s) | Memory (MB) | Recovery (%) |
|---------|------|-------------|-------------|--------------|
| SPLiT-seq PE | seqproc | 1.06 | 28.7 | 81.32 |
| SPLiT-seq PE | matchbox | 97.28 | 331.4 | 81.75 |
| 10x GridION | seqproc | 4.29 | 33.2 | 92.71 |
| 10x GridION | matchbox | 7.36 | 355.1 | 92.71 |

Key observations:
- Seqproc is 91.8x faster than Matchbox on SPLiT-seq PE
- Seqproc uses 11.6x less memory than Matchbox
- Recovery rates are comparable between tools when using edit distance

## Troubleshooting

### 1. Seqproc produces 0 reads

Check that:
- The SEQPROC_BIN environment variable is set correctly
- The binary is executable and at the correct path
- Whitelist files exist at `configs/seqproc/splitseq_bc23_whitelist.txt` and `configs/seqproc/splitseq_bc1_whitelist.txt`
- Run `./edit_distance_results/install_whitelists.sh` from the repo root to generate them

### 2. Permission errors

Make sure:
- The data files are readable
- The output directory is writable
- Tool binaries have execute permissions

### 3. Missing dependencies

Install required Python packages:
```bash
sudo apt install python3-numpy python3-matplotlib
```

## Reproducing the Exact Results

To reproduce the exact results from the benchmark:

1. Ensure you have the exact same 1M read subsets of the data
2. Use the same tool versions (seqproc with edit distance support)
3. Run with the same number of threads (recommended: 4)
4. The benchmark uses deterministic algorithms, so results should be identical

## Citation

If you use these results, please cite the seqproc paper and mention that the benchmarks were performed using edit distance (Levenshtein) for barcode matching across all tools.
