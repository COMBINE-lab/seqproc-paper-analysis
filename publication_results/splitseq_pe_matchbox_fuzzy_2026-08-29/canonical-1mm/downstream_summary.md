# Final SPLiT-seq PE downstream summary

| Tool | Input pairs | Valid barcodes | Called cells | Inflection rank | STARsolo time (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| seqproc | 56,991,381 | 99.1886% | 225 | 257 | 47.57 | 29309.3 |
| splitcode | 53,495,595 | 99.9723% | 220 | 255 | 43.45 | 29292.8 |
| matchbox | 50,266,141 | 99.9739% | 221 | 251 | 43.37 | 29294.0 |

| Tool pair | Read-set Jaccard | Per-gene Pearson | Per-barcode Pearson | Cell-type agreement | Mean type Jaccard | Cluster ARI |
|---|---:|---:|---:|---:|---:|---:|
| seqproc / splitcode | 0.930 | 0.995 | 0.972 | 0.914 | 0.686 | 0.591 |
| seqproc / matchbox | 0.870 | 0.995 | 0.983 | 0.923 | 0.704 | 0.678 |
| splitcode / matchbox | 0.823 | 0.991 | 0.962 | 0.936 | 0.732 | 0.606 |

Shared called cells: **220**.  All-tool cell-type agreement: **0.891**.  Mean per-type Jaccard: **0.708**.
