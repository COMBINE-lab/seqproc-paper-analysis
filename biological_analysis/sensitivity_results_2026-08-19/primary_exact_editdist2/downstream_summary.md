# Final SPLiT-seq PE downstream summary

| Tool | Input pairs | Valid barcodes | Called cells | Inflection rank | STARsolo time (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| seqproc | 56,991,381 | 99.9731% | 225 | 258 | 46.54 | 29311.4 |
| splitcode | 53,495,595 | 99.9723% | 220 | 255 | 41.78 | 29294.2 |
| matchbox | 35,160,366 | 99.9747% | 211 | 252 | 30.36 | 29281.3 |

| Tool pair | Read-set Jaccard | Per-gene Pearson | Per-barcode Pearson | Cell-type agreement | Mean type Jaccard | Cluster ARI |
|---|---:|---:|---:|---:|---:|---:|
| seqproc / splitcode | 0.9298 | 0.9953 | 0.9809 | 0.9431 | 0.7592 | 0.6482 |
| seqproc / matchbox | 0.6169 | 0.9856 | 0.9483 | 0.8720 | 0.5884 | 0.5760 |
| splitcode / matchbox | 0.6573 | 0.9885 | 0.9564 | 0.9005 | 0.6691 | 0.5976 |

Shared called cells: **211**.  All-tool cell-type agreement: **0.8673**.  Mean per-type Jaccard: **0.6722**.
