# Final SPLiT-seq PE downstream summary

| Tool | Input pairs | Valid barcodes | Called cells | Inflection rank | STARsolo time (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| seqproc | 56,991,381 | 99.1886% | 225 | 257 | 47.57 | 29309.3 |
| splitcode | 53,495,595 | 99.9723% | 220 | 255 | 43.45 | 29292.8 |
| matchbox | 35,160,366 | 99.9747% | 211 | 252 | 30.21 | 29282.9 |

| Tool pair | Read-set Jaccard | Per-gene Pearson | Per-barcode Pearson | Cell-type agreement | Mean type Jaccard | Cluster ARI |
|---|---:|---:|---:|---:|---:|---:|
| seqproc / splitcode | 0.9298 | 0.9949 | 0.9722 | 0.9194 | 0.6989 | 0.6199 |
| seqproc / matchbox | 0.6169 | 0.9860 | 0.9538 | 0.8910 | 0.6309 | 0.6144 |
| splitcode / matchbox | 0.6573 | 0.9885 | 0.9564 | 0.9005 | 0.6691 | 0.5976 |

Shared called cells: **211**.  All-tool cell-type agreement: **0.8673**.  Mean per-type Jaccard: **0.6663**.
