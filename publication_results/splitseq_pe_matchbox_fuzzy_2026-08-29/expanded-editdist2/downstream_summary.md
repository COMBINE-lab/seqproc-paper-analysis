# Final SPLiT-seq PE downstream summary

| Tool | Input pairs | Valid barcodes | Called cells | Inflection rank | STARsolo time (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| seqproc | 56,991,381 | 99.9731% | 225 | 258 | 46.54 | 29311.4 |
| splitcode | 53,495,595 | 99.9723% | 220 | 255 | 41.78 | 29294.2 |
| matchbox | 57,252,325 | 99.9734% | 225 | 259 | 47.05 | 29312.8 |

| Tool pair | Read-set Jaccard | Per-gene Pearson | Per-barcode Pearson | Cell-type agreement | Mean type Jaccard | Cluster ARI |
|---|---:|---:|---:|---:|---:|---:|
| seqproc / splitcode | 0.930 | 0.995 | 0.981 | 0.936 | 0.744 | 0.634 |
| seqproc / matchbox | 0.982 | 0.999 | 0.985 | 0.936 | 0.766 | 0.724 |
| splitcode / matchbox | 0.922 | 0.995 | 0.975 | 0.923 | 0.731 | 0.545 |

Shared called cells: **220**.  All-tool cell-type agreement: **0.900**.  Mean per-type Jaccard: **0.747**.
