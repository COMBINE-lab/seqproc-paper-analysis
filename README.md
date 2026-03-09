# seqproc Paper Analysis

Benchmarking and accuracy analysis for the seqproc paper, comparing seqproc
against matchbox and splitcode across four single-cell RNA-seq chemistries.

## Paper Artifacts

Each paper figure and table has a dedicated directory containing a README with
reproduction instructions, links to configs, scripts, and output files.

| Directory | Paper artifact | Description |
|-----------|---------------|-------------|
| [`table2/`](table2/) | Table 1 | Benchmark summary: recovery, runtime, memory |
| [`fig_recovery/`](fig_recovery/) | Figure 3 | Read recovery comparison bar chart |
| [`fig_concordance/`](fig_concordance/) | Figure 4 | Pairwise concordance heatmaps (Jaccard) |
| [`fig_hamming_vs_edit/`](fig_hamming_vs_edit/) | Figure 5 | Hamming vs edit distance comparison |
| [`fig_discordant/`](fig_discordant/) | Supp. Figure S1 | Discordant read recovery breakdown |

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Set tool paths (if not in PATH)
export SEQPROC_BIN=/path/to/seqproc
export MATCHBOX_BIN=/path/to/matchbox
export SPLITCODE_BIN=/path/to/splitcode
```

## Full Reproduction

To reproduce all paper results from scratch:

```bash
# Step 1: Run performance benchmarks (Table 1)
python scripts/run_paper_benchmarks.py --threads 4 --replicates 3

# Step 2: Run concordance analysis (Figs 3-5 data)
python scripts/concordance_analysis.py --threads 4

# Step 3: Run discordant read structural validation (Supp. Fig S1 data)
python scripts/discordant_analysis.py

# Step 4: Generate all figures and combined JSON
python scripts/generate_figures.py

# Step 5: (Optional) Re-run LR-SPLiT-seq performance with current config
python scripts/lr_perf_rerun.py --threads 4 --reps 3
```

### Datasets

| Dataset | SRA Accession | Reads | Mode |
|---------|---------------|-------|------|
| SPLiT-seq PE | SRR6750041 | 1,000,000 | paired-end |
| LR-SPLiT-seq | SRR13948564 | 1,000,000 | single-end |
| 10x Chromium v2 | SRR8315379 | 1,000,000 | paired-end |
| sci-RNA-seq3 | SRR7827254 | 1,000,000 | paired-end |

FASTQ files are stored in `data/` (gitignored).

## Testing

```bash
# Run the full regression test suite (73 tests)
python -m pytest tests/ -v
```

Tests cover script compilation, module imports, path resolution, banned
terminology guards, result integrity, figure determinism, and bug-fix
regression checks. A pre-commit hook runs the suite automatically before
every commit.

## Directory Structure

```
table2/                     # --> Paper Table 1 reproduction README
fig_recovery/               # --> Paper Figure 3 reproduction README
fig_concordance/            # --> Paper Figure 4 reproduction README
fig_hamming_vs_edit/        # --> Paper Figure 5 reproduction README
fig_discordant/             # --> Paper Supp. Figure S1 reproduction README
configs/
  seqproc/                  # EFGDL geometry files (.geom) + whitelists
  matchbox/                 # Matchbox pattern files (.mb) + support CSVs
  splitcode/                # Splitcode config files (.config)
scripts/
  run_paper_benchmarks.py   # Main benchmark runner (Table 1)
  concordance_analysis.py   # Pairwise concordance + hamming vs edit
  discordant_analysis.py    # Structural validation of tool-unique reads
  generate_figures.py       # Publication figure generation
  lr_perf_rerun.py          # LR-SPLiT-seq performance re-run
  splitcode_dual.py         # Splitcode dual-orientation helper
data/                       # Input FASTQ datasets (gitignored)
results/
  paper_figures/            # Final figures (PDF + PNG), benchmark JSON
  concordance/              # Pairwise concordance analysis outputs
  orientation/              # Orientation benchmark outputs
  lr_perf/                  # LR-SPLiT-seq performance re-run results
  splitseq_pe_perf/         # SPLiT-seq PE performance re-run results
tests/
  test_regression.py        # Comprehensive regression test suite
  test_benchmark_pipeline.py # Bug-fix unit tests
  conftest.py               # Shared pytest fixtures
```

## Configurations

### seqproc (`configs/seqproc/`)

| Config | Dataset | Notes |
|--------|---------|-------|
| `splitseq_filter_edit.geom` | SPLiT-seq PE | Edit distance, whitelist filtering |
| `splitseq_filter_hamming6.geom` | SPLiT-seq PE | Hamming baseline for comparison |
| `splitseq_replacement_edit.geom` | SPLiT-seq PE | Replacement mode variant |
| `splitseq_singleend_edit_ann.geom` | LR-SPLiT-seq | Annotation + edit distance |
| `splitseq_singleend_ann.geom` | LR-SPLiT-seq | Annotation + hamming |
| `splitseq_singleend_edit.geom` | LR-SPLiT-seq | Forward-only + edit |
| `splitseq_singleend.geom` | LR-SPLiT-seq | Forward-only + hamming |
| `splitseq_singleend_primer_edit.geom` | LR-SPLiT-seq | Primer-based variant |
| `10x_v2.geom` | 10x Chromium v2 | Fixed-position extraction |
| `10x_longread_fwd_edit.geom` | 10x Long Read | Forward orientation |
| `10x_longread_rev_edit.geom` | 10x Long Read | Reverse orientation |
| `sciseq3.geom` | sci-RNA-seq3 | Hamming baseline |
| `sciseq3_edit.geom` | sci-RNA-seq3 | Edit distance |

Support files: `splitseq_bc1_whitelist_6bp.txt`, `splitseq_bc23_whitelist.txt`

### matchbox (`configs/matchbox/`)

| Config | Dataset |
|--------|---------|
| `splitseq_replacement.mb` | SPLiT-seq PE |
| `splitseq_singleend.mb` | LR-SPLiT-seq (forward) |
| `splitseq_singleend_dual.mb` | LR-SPLiT-seq (dual orientation) |
| `10x_v2.mb` | 10x Chromium v2 |
| `10x_longread.mb` | 10x Long Read |
| `sciseq3.mb` | sci-RNA-seq3 |

Support files: `rt.csv`, `r2_r3.txt`

### splitcode (`configs/splitcode/`)

| Config | Dataset |
|--------|---------|
| `splitseq_paper.config` | SPLiT-seq PE |
| `splitseq_singleend.config` | LR-SPLiT-seq |
| `10x_v2.config` | 10x Chromium v2 |
| `10x_longread.config` | 10x Long Read |
| `sciseq3.config` | sci-RNA-seq3 |

## Citation

If you use this analysis, please cite the seqproc paper.
