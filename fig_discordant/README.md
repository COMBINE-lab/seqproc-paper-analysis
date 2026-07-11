# Supplementary Figure S1: Discordant Read Recovery Breakdown

**Paper reference:** Supplementary Figure S1 ("Supplementary Information")

## What this figure shows

Stacked bar chart showing read recovery broken down by tool concordance: reads
agreed upon by all three tools, reads shared by exactly two tools, and reads
unique to each tool. Dashed lines indicate total input reads per dataset.

## Key finding

On SPLiT-seq PE, the large splitcode-only segment (77,395 reads) consists
almost entirely (99.2%) of structurally invalid reads that lack valid linker
sequences at expected positions, indicating false-positive recovery.

## Reproduction

```bash
# 1. Run concordance analysis
python scripts/concordance_analysis.py --threads 4

# 2. Run discordant read structural validation
python scripts/discordant_analysis.py

# 3. Generate figures
python scripts/generate_figures.py
```

## Output files

- [`results/paper_figures/fig_discordant_summary.pdf`](../results/paper_figures/fig_discordant_summary.pdf) -- publication PDF
- [`results/paper_figures/fig_discordant_summary.png`](../results/paper_figures/fig_discordant_summary.png) -- PNG preview

## Scripts

- [`scripts/concordance_analysis.py`](../scripts/concordance_analysis.py) -- identifies discordant reads
- [`scripts/discordant_analysis.py`](../scripts/discordant_analysis.py) -- structural validation of tool-unique reads
- [`scripts/generate_figures.py`](../scripts/generate_figures.py) -- figure rendering
