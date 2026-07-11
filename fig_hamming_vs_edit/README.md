# Figure 5: Hamming vs Edit Distance Comparison

**Paper reference:** Figure 5 in Section 4.3 ("Validation of read recovery
accuracy")

## What this figure shows

Two-panel comparison of seqproc read recovery using Hamming distance vs edit
distance matching. Left panel: absolute read counts. Right panel: percentage
gain from edit distance. Edit distance is not applicable to 10x Chromium v2
(no anchor sequences).

## Key results

| Dataset | Hamming reads | Edit reads | Gain |
|---------|--------------|------------|------|
| SPLiT-seq PE | 802,750 | 841,089 | +4.8% |
| LR-SPLiT-seq | 430,879 | 498,723 | +15.8% |
| sci-RNA-seq3 | 884,389 | 892,626 | +0.9% |

## Reproduction

```bash
# 1. Run concordance analysis (includes hamming vs edit comparison)
python scripts/concordance_analysis.py --threads 4

# 2. Generate figures
python scripts/generate_figures.py
```

## Output files

- [`results/paper_figures/fig_hamming_vs_edit.pdf`](../results/paper_figures/fig_hamming_vs_edit.pdf) -- publication PDF
- [`results/paper_figures/fig_hamming_vs_edit.png`](../results/paper_figures/fig_hamming_vs_edit.png) -- PNG preview

## Configs

Hamming vs edit comparison uses two seqproc configs per dataset:

| Dataset | Hamming config | Edit config |
|---------|---------------|-------------|
| SPLiT-seq PE | [`splitseq_filter_hamming6.geom`](../configs/seqproc/splitseq_filter_hamming6.geom) | [`splitseq_filter_edit.geom`](../configs/seqproc/splitseq_filter_edit.geom) |
| LR-SPLiT-seq | [`splitseq_singleend_ann.geom`](../configs/seqproc/splitseq_singleend_ann.geom) | [`splitseq_singleend_edit_ann.geom`](../configs/seqproc/splitseq_singleend_edit_ann.geom) |
| sci-RNA-seq3 | [`sciseq3.geom`](../configs/seqproc/sciseq3.geom) | [`sciseq3_edit.geom`](../configs/seqproc/sciseq3_edit.geom) |

## Scripts

- [`scripts/concordance_analysis.py`](../scripts/concordance_analysis.py) -- runs both hamming and edit configs
- [`scripts/generate_figures.py`](../scripts/generate_figures.py) -- figure rendering
