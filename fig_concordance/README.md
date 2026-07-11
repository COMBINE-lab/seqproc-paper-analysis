# Figure 4: Pairwise Concordance Heatmaps

**Paper reference:** Figure 4 in Section 4.2 ("Efficient pre-processing of
scRNA sequencing data")

## What this figure shows

Pairwise Jaccard index between seqproc, matchbox, and splitcode for each of the
four benchmark datasets. High values (dark red) indicate strong agreement
between tools on which reads are successfully parsed.

## Key results

| Dataset | seqproc vs matchbox | seqproc vs splitcode | matchbox vs splitcode |
|---------|--------------------|--------------------|---------------------|
| SPLiT-seq PE | 0.925 | 0.913 | 0.846 |
| LR-SPLiT-seq | 0.797 | 0.501 | 0.451 |
| 10x Short | 1.000 | 1.000 | 1.000 |
| sci-RNA-seq3 | 0.994 | 0.991 | 0.985 |

## Reproduction

```bash
# 1. Run concordance analysis (requires tool binaries and FASTQ data)
python scripts/concordance_analysis.py --threads 4

# 2. Generate figures from concordance results
python scripts/generate_figures.py
```

## Output files

- [`results/paper_figures/fig_concordance_heatmaps.pdf`](../results/paper_figures/fig_concordance_heatmaps.pdf) -- publication PDF
- [`results/paper_figures/fig_concordance_heatmaps.png`](../results/paper_figures/fig_concordance_heatmaps.png) -- PNG preview
- [`results/concordance/concordance_results.json`](../results/concordance/concordance_results.json) -- raw concordance data

## Scripts

- [`scripts/concordance_analysis.py`](../scripts/concordance_analysis.py) -- pairwise concordance computation
- [`scripts/generate_figures.py`](../scripts/generate_figures.py) -- figure rendering

## Configs

Same tool configurations as Table 2 (see [`table2/README.md`](../table2/README.md)).
