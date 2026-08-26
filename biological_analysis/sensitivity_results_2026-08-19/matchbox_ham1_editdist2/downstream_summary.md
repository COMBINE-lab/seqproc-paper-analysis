# Final SPLiT-seq PE downstream summary

| Tool | Input pairs | Valid barcodes | Called cells | Inflection rank | STARsolo time (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| seqproc | 56,991,381 | 99.9731% | 225 | 258 | 46.54 | 29311.4 |
| splitcode | 53,495,595 | 99.9723% | 220 | 255 | 41.78 | 29294.2 |
| matchbox | 39,758,892 | 99.9744% | 217 | 254 | 33.55 | 29282.8 |

| Tool pair | Read-set Jaccard | Per-gene Pearson | Per-barcode Pearson | Cell-type agreement | Mean type Jaccard | Cluster ARI |
|---|---:|---:|---:|---:|---:|---:|
| seqproc / splitcode | 0.930 | 0.995 | 0.981 | 0.945 | 0.771 | 0.647 |
| seqproc / matchbox | 0.698 | 0.989 | 0.958 | 0.903 | 0.675 | 0.547 |
| splitcode / matchbox | 0.743 | 0.992 | 0.968 | 0.908 | 0.699 | 0.602 |

Shared called cells: **217**.  All-tool cell-type agreement: **0.885**.  Mean per-type Jaccard: **0.715**.
