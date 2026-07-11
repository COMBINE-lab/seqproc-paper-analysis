# Table 2: Benchmark Summary (Recovery, Runtime, Memory)

**Paper reference:** Table 1 in Section 4.2 ("Efficient pre-processing of scRNA
sequencing data")

## What this table shows

Average recovery rate, runtime, and memory usage for seqproc, matchbox, and
splitcode across four single-cell RNA-seq chemistries (3 replicates each,
1,000,000 reads per dataset).

## Reproduction

```bash
# 1. Run all benchmarks (requires tool binaries and input FASTQ data)
python scripts/run_paper_benchmarks.py --threads 4 --replicates 3

# 2. Generate the combined JSON with merged concordance + performance data
python scripts/generate_figures.py
```

The table in the paper is typeset directly in LaTeX from the numbers in
`benchmark_results.json`. The script also produces a rendered PNG table
(`fig_performance_table.png`) for quick reference.

## Input data

| Dataset | SRA Accession | Reads | Mode |
|---------|---------------|-------|------|
| SPLiT-seq PE | SRR6750041 | 1,000,000 | paired-end |
| LR-SPLiT-seq | SRR13948564 | 1,000,000 | single-end |
| 10x Chromium v2 | SRR8315379 | 1,000,000 | paired-end |
| sci-RNA-seq3 | SRR7827254 | 1,000,000 | paired-end |

FASTQ files are stored in `../data/` (gitignored).

## Tool configurations

### seqproc
| Dataset | Config |
|---------|--------|
| SPLiT-seq PE | [`configs/seqproc/splitseq_filter_edit.geom`](../configs/seqproc/splitseq_filter_edit.geom) |
| LR-SPLiT-seq | [`configs/seqproc/splitseq_singleend_edit_ann.geom`](../configs/seqproc/splitseq_singleend_edit_ann.geom) |
| 10x Chromium v2 | [`configs/seqproc/10x_v2.geom`](../configs/seqproc/10x_v2.geom) |
| sci-RNA-seq3 | [`configs/seqproc/sciseq3_edit.geom`](../configs/seqproc/sciseq3_edit.geom) |

Support files: [`splitseq_bc1_whitelist_6bp.txt`](../configs/seqproc/splitseq_bc1_whitelist_6bp.txt), [`splitseq_bc23_whitelist.txt`](../configs/seqproc/splitseq_bc23_whitelist.txt)

### matchbox
| Dataset | Config |
|---------|--------|
| SPLiT-seq PE | [`configs/matchbox/splitseq_replacement.mb`](../configs/matchbox/splitseq_replacement.mb) |
| LR-SPLiT-seq | [`configs/matchbox/splitseq_singleend.mb`](../configs/matchbox/splitseq_singleend.mb) |
| 10x Chromium v2 | [`configs/matchbox/10x_v2.mb`](../configs/matchbox/10x_v2.mb) |
| sci-RNA-seq3 | [`configs/matchbox/sciseq3.mb`](../configs/matchbox/sciseq3.mb) |

Support files: [`rt.csv`](../configs/matchbox/rt.csv), [`r2_r3.txt`](../configs/matchbox/r2_r3.txt)

### splitcode
| Dataset | Config |
|---------|--------|
| SPLiT-seq PE | [`configs/splitcode/splitseq_paper.config`](../configs/splitcode/splitseq_paper.config) |
| LR-SPLiT-seq | [`configs/splitcode/splitseq_singleend.config`](../configs/splitcode/splitseq_singleend.config) |
| 10x Chromium v2 | [`configs/splitcode/10x_v2.config`](../configs/splitcode/10x_v2.config) |
| sci-RNA-seq3 | [`configs/splitcode/sciseq3.config`](../configs/splitcode/sciseq3.config) |

## Output files

- [`results/paper_figures/benchmark_results.json`](../results/paper_figures/benchmark_results.json) -- combined results (recovery, runtime, memory, concordance)
- [`results/paper_figures/benchmark_results_perf.json`](../results/paper_figures/benchmark_results_perf.json) -- performance-only backup
- [`results/paper_figures/fig_performance_table.png`](../results/paper_figures/fig_performance_table.png) -- rendered table image

## Scripts

- [`scripts/run_paper_benchmarks.py`](../scripts/run_paper_benchmarks.py) -- main benchmark runner
- [`scripts/generate_figures.py`](../scripts/generate_figures.py) -- figure and JSON generation
