# Final SPLiT-seq PE downstream summary

| Tool | Input pairs | Valid barcodes | Called cells | Inflection rank | STARsolo time (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| seqproc | 56,991,381 | 99.9731% | 225 | 258 | 46.54 | 29311.4 |
| splitcode | 53,495,595 | 99.9723% | 220 | 255 | 41.78 | 29294.2 |
| matchbox | 50,266,141 | 99.9739% | 221 | 251 | 45.17 | 29297.1 |

| Tool pair | Read-set Jaccard | Per-gene Pearson | Per-barcode Pearson | Cell-type agreement | Mean type Jaccard | Cluster ARI |
|---|---:|---:|---:|---:|---:|---:|
| seqproc / splitcode | 0.930 | 0.995 | 0.981 | 0.936 | 0.744 | 0.634 |
| seqproc / matchbox | 0.870 | 0.995 | 0.972 | 0.932 | 0.737 | 0.635 |
| splitcode / matchbox | 0.823 | 0.991 | 0.962 | 0.936 | 0.732 | 0.606 |

Shared called cells: **220**.  All-tool cell-type agreement: **0.909**.  Mean per-type Jaccard: **0.738**.
