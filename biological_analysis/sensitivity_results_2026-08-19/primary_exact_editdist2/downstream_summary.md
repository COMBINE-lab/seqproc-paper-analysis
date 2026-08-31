# Final SPLiT-seq PE downstream summary

| Tool | Input pairs | Valid barcodes | Called cells | Inflection rank | STARsolo time (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| seqproc | 56,991,381 | 99.9731% | 225 | 258 | 46.54 | 29311.4 |
| splitcode | 53,495,595 | 99.9723% | 220 | 255 | 41.78 | 29294.2 |
| matchbox | 35,160,366 | 99.9747% | 211 | 252 | 30.36 | 29281.3 |

| Tool pair | Read-set Jaccard | Per-gene Pearson | Per-barcode Pearson | Cell-type agreement | Mean type Jaccard | Cluster ARI |
|---|---:|---:|---:|---:|---:|---:|
| seqproc / splitcode | 0.930 | 0.995 | 0.981 | 0.943 | 0.759 | 0.648 |
| seqproc / matchbox | 0.617 | 0.986 | 0.948 | 0.872 | 0.588 | 0.576 |
| splitcode / matchbox | 0.657 | 0.989 | 0.956 | 0.900 | 0.669 | 0.598 |

Shared called cells: **211**.  All-tool cell-type agreement: **0.867**.  Mean per-type Jaccard: **0.672**.
