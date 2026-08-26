# Final SPLiT-seq PE downstream summary

| Tool | Input pairs | Valid barcodes | Called cells | Inflection rank | STARsolo time (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| seqproc | 56,991,381 | 99.1886% | 225 | 257 | 47.57 | 29309.3 |
| splitcode | 53,495,595 | 99.9723% | 220 | 255 | 43.45 | 29292.8 |
| matchbox | 39,758,892 | 99.2976% | 215 | 254 | 34.71 | 29280.7 |

| Tool pair | Read-set Jaccard | Per-gene Pearson | Per-barcode Pearson | Cell-type agreement | Mean type Jaccard | Cluster ARI |
|---|---:|---:|---:|---:|---:|---:|
| seqproc / splitcode | 0.930 | 0.995 | 0.972 | 0.916 | 0.682 | 0.615 |
| seqproc / matchbox | 0.698 | 0.989 | 0.961 | 0.935 | 0.761 | 0.579 |
| splitcode / matchbox | 0.743 | 0.992 | 0.963 | 0.935 | 0.774 | 0.593 |

Shared called cells: **215**.  All-tool cell-type agreement: **0.902**.  Mean per-type Jaccard: **0.739**.
