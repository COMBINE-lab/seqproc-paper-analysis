# Final SPLiT-seq PE downstream summary

| Tool | Input pairs | Valid barcodes | Called cells | Inflection rank | STARsolo time (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| seqproc | 56,991,381 | 99.1886% | 225 | 257 | 47.57 | 29309.3 |
| splitcode | 53,495,595 | 99.9723% | 220 | 255 | 43.45 | 29292.8 |
| matchbox | 35,160,366 | 99.9747% | 211 | 252 | 30.21 | 29282.9 |

| Tool pair | Read-set Jaccard | Per-gene Pearson | Per-barcode Pearson | Cell-type agreement | Mean type Jaccard | Cluster ARI |
|---|---:|---:|---:|---:|---:|---:|
| seqproc / splitcode | 0.930 | 0.995 | 0.972 | 0.919 | 0.699 | 0.620 |
| seqproc / matchbox | 0.617 | 0.986 | 0.954 | 0.891 | 0.631 | 0.614 |
| splitcode / matchbox | 0.657 | 0.989 | 0.956 | 0.900 | 0.669 | 0.598 |

Shared called cells: **211**.  All-tool cell-type agreement: **0.867**.  Mean per-type Jaccard: **0.666**.
