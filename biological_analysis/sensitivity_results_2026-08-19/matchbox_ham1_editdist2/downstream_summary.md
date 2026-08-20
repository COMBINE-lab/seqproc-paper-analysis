# Final SPLiT-seq PE downstream summary

| Tool | Input pairs | Valid barcodes | Called cells | Inflection rank | STARsolo time (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| seqproc | 56,991,381 | 99.9731% | 225 | 258 | 46.54 | 29311.4 |
| splitcode | 53,495,595 | 99.9723% | 220 | 255 | 41.78 | 29294.2 |
| matchbox | 39,758,892 | 99.9744% | 217 | 254 | 33.55 | 29282.8 |

| Tool pair | Read-set Jaccard | Per-gene Pearson | Per-barcode Pearson | Cell-type agreement | Mean type Jaccard | Cluster ARI |
|---|---:|---:|---:|---:|---:|---:|
| seqproc / splitcode | 0.9298 | 0.9953 | 0.9809 | 0.9447 | 0.7708 | 0.6473 |
| seqproc / matchbox | 0.6976 | 0.9888 | 0.9576 | 0.9032 | 0.6752 | 0.5473 |
| splitcode / matchbox | 0.7432 | 0.9922 | 0.9676 | 0.9078 | 0.6987 | 0.6017 |

Shared called cells: **217**.  All-tool cell-type agreement: **0.8848**.  Mean per-type Jaccard: **0.7149**.
