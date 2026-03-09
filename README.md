# seqproc Paper Analysis

Benchmarking and accuracy analysis scripts for the seqproc paper, comparing
seqproc against matchbox and splitcode across four single-cell RNA-seq
chemistries.

## Results

Results are generated into `results/paper_figures/` and include:
- `benchmark_results.json` -- all benchmark numbers (recovery, runtime, memory)
- Concordance heatmaps, recovery comparison, hamming vs edit, discordant summary figures

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Set tool paths (if not in PATH)
export SEQPROC_BIN=/path/to/seqproc
export MATCHBOX_BIN=/path/to/matchbox
export SPLITCODE_BIN=/path/to/splitcode
```

## Reproduction Steps

### 1. Main Paper Benchmarks (Table 2)

Runs all three tools on all four datasets with 3 replicates each.

**Script:** `scripts/run_paper_benchmarks.py`

```bash
python scripts/run_paper_benchmarks.py --threads 4 --replicates 3
```

**Datasets:**
- SPLiT-seq PE (`SRR6750041`, 1M reads)
- LR-SPLiT-seq (`SRR13948564`, 1M reads)
- 10x Chromium v2 (`SRR8315379`, 1M reads)
- sci-RNA-seq3 (`SRR7827254`, 1M reads)

### 2. Concordance Analysis

Pairwise Jaccard concordance, discordant read characterization, and
hamming vs edit distance comparison across all datasets.

**Script:** `scripts/phase4_concordance.py`

```bash
python scripts/phase4_concordance.py --threads 4
```

### 3. Discordant Read Analysis

Structural validation of tool-unique reads (e.g., splitcode false positive
characterization on SPLiT-seq PE).

**Script:** `scripts/phase4_discordant_analysis.py`

```bash
python scripts/phase4_discordant_analysis.py
```

### 4. Figure Generation

Generates publication figures from the concordance and benchmark results.

**Script:** `scripts/phase4_figures.py`

```bash
python scripts/phase4_figures.py
```

### 5. LR-SPLiT-seq Performance Re-run

Fresh 3-replicate performance measurements for the LR-SPLiT-seq dataset.

**Script:** `scripts/phase5_lr_perf_rerun.py`

```bash
python scripts/phase5_lr_perf_rerun.py --threads 4 --reps 3
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

## Directory Structure

```
configs/                    # Tool configurations (.geom, .mb, .config)
data/                       # Input datasets (gitignored)
results/
  paper_figures/            # Final benchmark results, figures, and JSON
  phase4_concordance/       # Pairwise concordance analysis outputs
  phase3_orientation/       # Orientation benchmark outputs
  phase5_lr_perf/           # LR-SPLiT-seq performance re-run
  phase5_splitseq_pe_perf/  # SPLiT-seq PE performance re-run
scripts/                    # Analysis and figure generation scripts
tests/                      # Pipeline regression tests
```

## Citation

If you use this analysis, please cite the seqproc paper.
