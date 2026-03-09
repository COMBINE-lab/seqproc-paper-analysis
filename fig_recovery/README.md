# Figure 3: Read Recovery Comparison

**Paper reference:** Figure 3 in Section 4.2 ("Efficient pre-processing of
scRNA sequencing data")

## What this figure shows

Grouped bar chart of read recovery rates (percentage of input reads successfully
parsed) for seqproc, matchbox, and splitcode across all four benchmark datasets.

## Key results

| Dataset | seqproc | matchbox | splitcode |
|---------|---------|----------|-----------|
| SPLiT-seq PE | 84.1% | 77.8% | 91.6%* |
| LR-SPLiT-seq | 49.9% | 39.7% | 28.0% |
| 10x Short | 100.0% | 100.0% | 100.0% |
| sci-RNA-seq3 | 89.3% | 89.8% | 88.4% |

*99.2% of splitcode-unique SPLiT-seq PE reads lack valid linker structure
(false positives).

## Reproduction

```bash
# 1. Run concordance analysis (computes recovery rates)
python scripts/concordance_analysis.py --threads 4

# 2. Generate figures
python scripts/generate_figures.py
```

## Output files

- [`results/paper_figures/fig_recovery_comparison.pdf`](../results/paper_figures/fig_recovery_comparison.pdf) -- publication PDF
- [`results/paper_figures/fig_recovery_comparison.png`](../results/paper_figures/fig_recovery_comparison.png) -- PNG preview

## Scripts

- [`scripts/concordance_analysis.py`](../scripts/concordance_analysis.py) -- recovery rate computation
- [`scripts/generate_figures.py`](../scripts/generate_figures.py) -- figure rendering
