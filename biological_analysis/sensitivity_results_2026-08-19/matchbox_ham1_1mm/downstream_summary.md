# Final SPLiT-seq PE downstream summary

| Tool | Input pairs | Valid barcodes | Called cells | Inflection rank | STARsolo time (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| seqproc | 56,991,381 | 99.1886% | 225 | 257 | 47.57 | 29309.3 |
| splitcode | 53,495,595 | 99.9723% | 220 | 255 | 43.45 | 29292.8 |
| matchbox | 39,758,892 | 99.2976% | 215 | 254 | 34.71 | 29280.7 |

| Tool pair | Read-set Jaccard | Per-gene Pearson | Per-barcode Pearson | Cell-type agreement | Mean type Jaccard | Cluster ARI |
|---|---:|---:|---:|---:|---:|---:|
| seqproc / splitcode | 0.9298 | 0.9949 | 0.9722 | 0.9163 | 0.6819 | 0.6148 |
| seqproc / matchbox | 0.6976 | 0.9889 | 0.9608 | 0.9349 | 0.7609 | 0.5788 |
| splitcode / matchbox | 0.7432 | 0.9919 | 0.9632 | 0.9349 | 0.7736 | 0.5935 |

Shared called cells: **215**.  All-tool cell-type agreement: **0.9023**.  Mean per-type Jaccard: **0.7388**.
