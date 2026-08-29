# Paper data → source-code traceability

> **Revision notice (2026-08-29):** this document principally traces the
> original preprint and therefore contains stale benchmark and validity values.
> The current journal-revision sources are the frozen base
> [`performance`](publication_results/journal_performance_2026-08-17/) bundle,
> the guarded fuzzy-linker SPLiT-seq PE Matchbox
> [`performance and downstream`](publication_results/splitseq_pe_matchbox_fuzzy_2026-08-29/)
> override, the LR Matchbox
> [`anchor-first`](publication_results/lr_splitseq_matchbox_anchor_first_2026-08-29/)
> performance/accuracy override, and the corresponding combined one-run,
> 32-thread
> [`accuracy`](publication_results/journal_rerun_2026-08-29-matchbox-pe-lr-anchor-first/)
> bundle. The
> current SPLiT-seq PE vendor-set comparison is the full-data
> [`splitseq_pe_splitpipe_vendor_full.json`](publication_results/journal_rerun_2026-08-17/splitseq_pe_splitpipe_vendor_full.json),
> with the guarded fuzzy-linker Matchbox comparison recorded separately in the
> 2026-08-29 downstream bundle.
> Do not use the old values below to rebuild the revised manuscript.

Every table cell, figure quantity, and prose number in the original seqproc preprint, mapped to the script that generates it, the exact command, and the output file. Built by 10 tracer agents reading `/home/ubuntu/paper` against this repo, then reconciled by hand.

**How to read:** each section is one paper artifact, with a **Regenerate** command block and a table. The last column is `✅` (reproducible from this repo) or `🔴 NOT REPRODUCIBLE` (needs an external input — flagged below). Values are as they appear in the current preprint (`sections/*.tex`, `seqproc.tex`).

**Prerequisites** (same for every command): the three tool binaries (`SEQPROC_BIN`/`MATCHBOX_BIN`/`SPLITCODE_BIN`), the four SRA FASTQs in `$SEQPROC_DATA_DIR`, the 10x whitelists (`TENX_WHITELIST_V3`/`V2`), and — for downstream — STARsolo + a GRCm38 index. See `README.md` "Reproducing the paper from scratch". Commands assume repo root and `--threads 32`.

---

## ⚠️ Five reproduction gaps (read before trusting a re-run)

The tracers found that a subset of headline numbers **cannot be reproduced by running the repo scripts as-is** — they are correct in the paper but the repo's stored inputs/outputs are stale or external. This is the "if it can't be done easily, flag it" part:

1. **Input sizes are ENA-confirmed manual corrections.** Table 1 uses **77,621,181** (PE) and **234,382,218** (10x), but `scripts/data_config.py:SRA_INFO` and `results_final/benchmark_results.json` still hold the **old wrong** 86,820,578 and 56,514,800. `count_fastq_reads` on the FASTQs reproduces the *wrong* PE value. The **Emitted %** column is then derived by hand as `emit / ENA-input`, not the script's `recovery_rate`. → *A reviewer re-running the benchmark gets different Emitted % unless they use the ENA denominators.*
2. **Table 1 runtime/memory come from a snapshot, not the JSON.** The corrected full-dataset runtime/memory live in `results_final/fig3_summary_table.png`; `notebooks/regenerate_figures.ipynb` explicitly warns that `benchmark_results.json` holds stale 1M-subset numbers that "must NOT be used." → *Trust `fig3_summary_table`, not the JSON, for those cells.*
3. **The original-preprint split-pipe workflow was not reproducible from this repository.** The journal revision resolves this for licensed users with the checksummed container recipe and protocol configuration in `containers/`, then validates the full output with `scripts/splitpipe_full_concordance.py`. The recovered command exactly reproduces the archived 10-million-pair ID set before running the full input. Use remains subject to the Parse Biosciences license.
4. **Appendix B engineering numbers are from a different repo.** The allocator/SIMD/recycling percentages (B.2–B.8) come from the **seqproc Rust** `cargo bench` suite, not this analysis repo. Only the **LR recovery progression** (the `tab:lr_progression` cells) is reproducible here, via `scripts/lr_perf_rerun.py`.
5. **The two schematic figures are hand-drawn.** `fig:seqproc_workflow` (`Figures/seqproc_flow.pdf`) and `fig:illustrative_example` (`Figures/antisequence_graph.pdf`) have no generator — vector art committed in `paper/Figures/`.

**37 original-preprint data points are marked 🔴** below (input sizes, all Appendix B engineering values, the historical vendor concordance, and the SI throughput timings). The journal-revision vendor comparison is now separately reproducible as described above.

---

## Pairwise downstream concordance (`tab:pairwise_concordance`, results.tex:87–88)

*(Reconstructed by hand — the automated tracer for this table hit the schema-retry cap.)*

**Regenerate:**
```bash
bash biological_analysis/run_downstream.sh        # STARsolo + count_concordance.py + biological_analysis.py
```
| Datum | Value | Source script → output | Repro |
|---|---|---|---|
| Per-gene total-UMI Pearson (log1p), each pair | 0.994–0.999 | `count_concordance.py` (`per_gene_total_pearson_logspace`) → `full_run_results/count_concordance.json` | ✅ |
| Per-barcode total-UMI Pearson (log1p), each pair | 0.974–0.992 | `count_concordance.py` (`per_barcode_umi_pearson_logspace`) → `count_concordance.json` | ✅ |
| Per-gene Spearman (log1p), each pair | ≥0.993 | `count_concordance.py` → `count_concordance.json` | ✅ |
| Cell-type agreement (seqproc·splitcode / seqproc·matchbox / splitcode·matchbox) | 0.936 / 0.927 / 0.914 | `biological_analysis.py` (marker-based typing) → `full_run_results/biological_metrics.json` | ✅ |
| Pre-processed emit counts (seqproc / matchbox / splitcode) | 63,482,886 / 62,757,663 / 58,005,118 | `run_downstream.sh` STARsolo inputs (results.tex:91) | ✅ |

---

## Table 1 — main benchmark (`tab:updated_results_table`)

**Regenerate:**
```bash
python3 scripts/run_paper_benchmarks.py --reads full --threads 32 --replicates 3 --datasets lr_splitseq
python3 scripts/run_paper_benchmarks.py --reads full --threads 32 --replicates 3 --datasets 10x_short
python3 scripts/run_paper_benchmarks.py --reads full --threads 32 --replicates 3 --datasets sciseq
python3 scripts/run_paper_benchmarks.py --reads full --threads 32 --replicates 3 --datasets splitseq_pe
python3 scripts/run_paper_benchmarks.py --reads full --threads 32 --replicates 3 --datasets lr_splitseq ; python3 scripts/splitcode_lr_dual_validate.py --dataset full --threads 8 --outdir results/splitcode_lr_dual
python3 scripts/concordance_analysis.py --threads 32 --datasets splitseq_pe ; python3 scripts/validate_pe_recall_against_vtotal.py --outdir results/pe_recall_vtotal
python3 scripts/splitcode_lr_dual_validate.py --dataset full --threads 8 --outdir results/splitcode_lr_dual ; python3 scripts/validate_lr_recall_against_vtotal.py --dual-pass-fq results/splitcode_lr_dual/splitcode_dual_combined_out.fq --outdir results/lr_recall_vtotal
python3 scripts/concordance_analysis.py --threads 32 --datasets 10x_short ; python3 scripts/validate_10x_sci_recall_against_vtotal.py --chemistry 10x --outdir results/vtotal_10x
```
| Datum | Value | Paper loc | Source script → output | Repro |
|---|---|---|---|---|
| SPLiT-seq PE input size (row header): 77,621,181 | `77,621,181` | results.tex tab:updated_results_table, SP… | NOT REPRODUCIBLE from repo scripts as printed — ENA-confirmed manual value (per user memo… → none (manual); repo path would be results/paper_figure… | 🔴 NOT REPRODUCIBLE |
| LR-SPLiT-seq input size (row header): 5,764,421 | `5,764,421` | results.tex tab:updated_results_table, LR… | scripts/data_config.py (SRA_INFO['SRR13948564'].full_reads = 5_764_421); also hardcoded t… → results/paper_figures/benchmark_results.json (lr_split… | ✅ |
| 10x Chromium v2 input size (row header): 234,382,218 | `234,382,218` | results.tex tab:updated_results_table, 10… | scripts/run_paper_benchmarks.py count_fastq_reads on 10x R1 FASTQ (measures 234,382,218; … → results/paper_figures/benchmark_results.json (10x_shor… | ✅ |
| sci-RNA-seq3 input size (row header): 22,088,821 | `22,088,821` | results.tex tab:updated_results_table, sc… | scripts/data_config.py (SRA_INFO['SRR7827254'].full_reads = 22_088_821); measured by coun… → results/paper_figures/benchmark_results.json (sciseq.t… | ✅ |
| SPLiT-seq PE Emitted% (seqproc/matchbox/splitcode): 84.6 / 78.4 / 91.8 | `84.6 / 78.4 / 91.8` | results.tex tab:updated_results_table, PE… | scripts/run_paper_benchmarks.py (emit counts reads_out = 65,705,225 / 60,828,797 / 71,290… → results/paper_figures/benchmark_results.json (splitseq… | ✅ |
| LR-SPLiT-seq Emitted% (seqproc/matchbox/splitcode): 49.7 / 39.6 / 51.8 | `49.7 / 39.6 / 51.8` | results.tex tab:updated_results_table, LR… | seqproc & matchbox: scripts/run_paper_benchmarks.py (reads_out/5,764,421). splitcode: DUA… → results/paper_figures/benchmark_results.json (lr_split… | ✅ |
| 10x Chromium v2 Emitted% (all three tools): 100 / 100 / 100 | `100 / 100 / 100` | results.tex tab:updated_results_table, 10… | scripts/run_paper_benchmarks.py (reads_out = total_reads = 234,382,218 for all three -> 1… → results/paper_figures/benchmark_results.json (10x_shor… | ✅ |
| sci-RNA-seq3 Emitted% (seqproc/matchbox/splitcode): 90.0 / 90.4 / 89.2 | `90.0 / 90.4 / 89.2` | results.tex tab:updated_results_table, sc… | scripts/run_paper_benchmarks.py (reads_out 19,879,895 / 19,980,569 / 19,695,955 divided b… → results/paper_figures/benchmark_results.json (sciseq.t… | ✅ |
| SPLiT-seq PE Precision% (seqproc/matchbox/splitcode): 92.33 / 94.41 /… | `92.33 / 94.41 / 84.96` | results.tex tab:updated_results_table, PE… | DERIVED (no dedicated script): precision = intersection_with_v_total / emit, both from sc… → results/pe_recall_vtotal/pe_recall_vtotal_results.json… | ✅ |
| LR-SPLiT-seq Precision% (seqproc/matchbox/splitcode): 20.73 / 23.22 /… | `20.73 / 23.22 / 19.80` | results.tex tab:updated_results_table, LR… | DERIVED: precision = intersection/emit from scripts/validate_lr_recall_against_vtotal.py … → results/lr_recall_vtotal/lr_recall_vtotal_results.json… | ✅ |
| 10x Chromium v2 Precision% (all three tools): 100 / 100 / 100 | `100 / 100 / 100` | results.tex tab:updated_results_table, 10… | DERIVED: intersection/emit from scripts/validate_10x_sci_recall_against_vtotal.py --chemi… → results/vtotal_10x/vtotal_recall_10x.json (intersectio… | ✅ |
| sci-RNA-seq3 Precision% (seqproc/matchbox/splitcode): 99.44 / 98.94 /… | `99.44 / 98.94 / 99.99` | results.tex tab:updated_results_table, sc… | DERIVED: intersection/emit from scripts/validate_10x_sci_recall_against_vtotal.py --chemi… → results/vtotal_sci/vtotal_recall_sci.json (intersectio… | ✅ |
| SPLiT-seq PE Recall% (seqproc/matchbox/splitcode): 99.98 / 94.65 / 99… | `99.98 / 94.65 / 99.83` | results.tex tab:updated_results_table, PE… | scripts/validate_pe_recall_against_vtotal.py (tools.<tool>.recall_pct = 100*intersection/… → results/pe_recall_vtotal/pe_recall_vtotal_results.json… | ✅ |
| LR-SPLiT-seq Recall% (seqproc/matchbox/splitcode): 98.72 / 87.97 / 98… | `98.72 / 87.97 / 98.22` | results.tex tab:updated_results_table, LR… | scripts/validate_lr_recall_against_vtotal.py (tools.{seqproc,matchbox,splitcode_dual}.rec… → results/lr_recall_vtotal/lr_recall_vtotal_results.json… | ✅ |
| 10x Chromium v2 Recall% (all three tools): 100 / 100 / 100 | `100 / 100 / 100` | results.tex tab:updated_results_table, 10… | scripts/validate_10x_sci_recall_against_vtotal.py --chemistry 10x (tools.<tool>.recall_pc… → results/vtotal_10x/vtotal_recall_10x.json (tools.<tool… | ✅ |
| sci-RNA-seq3 Recall% (seqproc/matchbox/splitcode): 100 / 100 / 99.63 | `100 / 100 / 99.63` | results.tex tab:updated_results_table, sc… | scripts/validate_10x_sci_recall_against_vtotal.py --chemistry sci (tools.<tool>.recall_pc… → results/vtotal_sci/vtotal_recall_sci.json (tools.<tool… | ✅ |
| SPLiT-seq PE F1 (seqproc/matchbox/splitcode): 0.960 / 0.945 / 0.918 | `0.960 / 0.945 / 0.918` | results.tex tab:updated_results_table, PE… | DERIVED (no script): F1 = 2*P*R/(P+R) from the PE Precision/Recall cells (validate_pe_rec… → none (hand-computed from results/pe_recall_vtotal/pe_r… | ✅ |
| LR-SPLiT-seq F1 (seqproc/matchbox/splitcode): 0.343 / 0.367 / 0.330 | `0.343 / 0.367 / 0.330` | results.tex tab:updated_results_table, LR… | DERIVED: F1 = 2PR/(P+R) from LR Precision/Recall cells. seqproc: 2*0.2073*0.9872/(0.2073+… → none (hand-computed from results/lr_recall_vtotal/lr_r… | ✅ |
| 10x Chromium v2 F1 (all three tools): 1.000 | `1.000` | results.tex tab:updated_results_table, 10… | DERIVED: F1 from P=R=100 -> 1.000 → none (hand-computed) | ✅ |
| sci-RNA-seq3 F1 (seqproc/matchbox/splitcode): 0.997 / 0.995 / 0.998 | `0.997 / 0.995 / 0.998` | results.tex tab:updated_results_table, sc… | DERIVED: F1 = 2PR/(P+R) from sci Precision/Recall cells. seqproc: 2*0.9944*1.0/(1.9944)=0… → none (hand-computed from results/vtotal_sci/vtotal_rec… | ✅ |
| SPLiT-seq PE Runtime s (seqproc/matchbox/splitcode): 95.0 / 1,290.9 /… | `95.0 / 1,290.9 / 114.7` | results.tex tab:updated_results_table, PE… | scripts/run_paper_benchmarks.py (mean_runtime over 3 replicates, --reads full, 32 threads) → results/paper_figures/benchmark_results.json (splitseq… | ✅ |
| LR-SPLiT-seq Runtime s (seqproc/matchbox/splitcode): 11.1 / 24.8 / 54… | `11.1 / 24.8 / 54.0` | results.tex tab:updated_results_table, LR… | seqproc & matchbox: scripts/run_paper_benchmarks.py (mean_runtime). splitcode: scripts/sp… → results/paper_figures/benchmark_results.json (lr_split… | ✅ |
| 10x Chromium v2 Runtime s (seqproc/matchbox/splitcode): 429.5 / 635.0… | `429.5 / 635.0 / 416.7` | results.tex tab:updated_results_table, 10… | scripts/run_paper_benchmarks.py (mean_runtime, --reads full, 32 threads) → results/paper_figures/benchmark_results.json (10x_shor… | ✅ |
| sci-RNA-seq3 Runtime s (seqproc/matchbox/splitcode): 26.8 / 198.4 / 3… | `26.8 / 198.4 / 32.4` | results.tex tab:updated_results_table, sc… | scripts/run_paper_benchmarks.py (mean_runtime, --reads full, 32 threads) → results/paper_figures/benchmark_results.json (sciseq.t… | ✅ |
| SPLiT-seq PE Memory MB (seqproc/matchbox/splitcode): 176 / 23,049 / 1… | `176 / 23,049 / 1,293` | results.tex tab:updated_results_table, PE… | scripts/run_paper_benchmarks.py (mean_memory_mb = mean peak RSS over 3 replicates; peak c… → results/paper_figures/benchmark_results.json (splitseq… | ✅ |
| LR-SPLiT-seq Memory MB (seqproc/matchbox/splitcode): 249 / 1,054 / 531 | `249 / 1,054 / 531` | results.tex tab:updated_results_table, LR… | seqproc & matchbox: scripts/run_paper_benchmarks.py (mean_memory_mb). splitcode: scripts/… → results/paper_figures/benchmark_results.json (lr_split… | ✅ |
| 10x Chromium v2 Memory MB (seqproc/matchbox/splitcode): 140 / 42,585 … | `140 / 42,585 / 1,259` | results.tex tab:updated_results_table, 10… | scripts/run_paper_benchmarks.py (mean_memory_mb, --reads full, 32 threads) → results/paper_figures/benchmark_results.json (10x_shor… | ✅ |
| sci-RNA-seq3 Memory MB (seqproc/matchbox/splitcode): 138 / 3,198 / 1,… | `138 / 3,198 / 1,850` | results.tex tab:updated_results_table, sc… | scripts/run_paper_benchmarks.py (mean_memory_mb, --reads full, 32 threads) → results/paper_figures/benchmark_results.json (sciseq.t… | ✅ |
| Caption/prose: V_total_PE covers 78.17% of PE input; V_total_LR cover… | `78.17% (PE), 10.44% (LR)` | results.tex line 52 (caption context) and… | DERIVED: 78.17 = /V_total_PE/ 60,675,548 / PE input 77,621,181; 10.44 = /V_total_LR/ 601,… → results/pe_recall_vtotal/pe_recall_vtotal_results.json… | ✅ |
| Prose: V_total_PE = 60,675,548 genuine reads; seqproc 99.98% (60,666,… | `60,675,548; 60,666,173; 60,570,75…` | results.tex prose line 60 | scripts/validate_pe_recall_against_vtotal.py (v_total; tools.<tool>.intersection_with_v_t… → results/pe_recall_vtotal/pe_recall_vtotal_results.json | ✅ |
| Prose: seqproc recovers 95,422 more genuine reads than splitcode whil… | `95,422; ~5.6 million` | results.tex prose line 62 | DERIVED: 95,422 = seqproc∩V_total (60,666,173) - splitcode∩V_total (60,570,751), from val… → results/pe_recall_vtotal/pe_recall_vtotal_results.json… | ✅ |
| Prose: V_total_LR = 601,603 (10.44% of input); seqproc 98.72% (593,87… | `601,603; 593,873; 590,883; 529,24…` | results.tex prose line 68 | scripts/validate_lr_recall_against_vtotal.py (v_total; intersection_with_v_total per tool… → results/lr_recall_vtotal/lr_recall_vtotal_results.json… | ✅ |
| Prose: the validation script processed 5,764,421 LR reads in 162 s (~… | `162 s; ~35,600 reads/s; ~15x` | results.tex prose line 54 | scripts/validate_lr_recall_against_vtotal.py — the V_total validator runtime is PRINTED t… → none (stdout timing only); ratio depends on seqproc LR… | ✅ |
| Prose: runtime multipliers — 1.2x-2.2x faster than next-fastest, up t… | `1.2x-2.2x; 13.6x; 1/2-1/304; 42.6…` | results.tex prose line 57 | DERIVED from Table 1 runtime/memory cells (run_paper_benchmarks.py). 13.6x = matchbox PE … → results/paper_figures/benchmark_results.json (mean_run… | ✅ |
| Prose: full LR seqproc/matchbox pairwise Jaccard 0.796; 279 reads uni… | `0.796; 279; 2.28 million` | results.tex prose line 70 | scripts/concordance_analysis.py (lr_splitseq pairwise Jaccard between emitted-ID sets) → results/concordance/lr_splitseq/results.json (concorda… | ✅ |
| Prose: 99.2% of the 5,787,314 reads unique to splitcode are structura… | `99.2%; 5,787,314` | results.tex prose lines 60 and 66 | 5,787,314 unique-to-splitcode from scripts/concordance_analysis.py/discordant_analysis.py… → results/concordance/splitseq_pe/results.json (discorda… | ✅ |
| Prose (vendor 10M concordance): split-pipe accepted 7,539,920 (75.4%)… | `7,539,920 (75.4%); seqproc 8,367,…` | results.tex prose lines 64-66 | NOT REPRODUCIBLE without the proprietary Parse Biosciences split-pipe tool (its barcode_h… → Shipped reference: results_final/splitpipe_valid_ids_1… | 🔴 NOT REPRODUCIBLE |
| Experimental setup: 32 threads all tools, 3 replicates, full dataset,… | `32 threads; 3 replicates; CV < 5%` | results.tex Experimental setup, line 5 | scripts/run_paper_benchmarks.py (--threads 32 --replicates 3 --reads full; mean over 3 ru… → results/paper_figures/benchmark_results.json (mean_* a… | ✅ |

## Emitted-set UpSet plot (`fig:emitted_set_upset`)

**Regenerate:**

```bash
python3 scripts/generate_emitted_set_upset.py \
  --results-dir publication_results/journal_rerun_2026-08-29-matchbox-pe-lr-anchor-first \
  --output-prefix publication_results/journal_rerun_2026-08-29-matchbox-pe-lr-anchor-first/fig_emitted_set_upset
```

The generator reads each final `*_accuracy_metrics.json` artifact and its
canonical accession bitmaps, verifies the recorded SHA-256 checksums and
emitted-record totals, unions multi-product and multi-orientation outputs, and
computes all seven mutually exclusive three-tool intersections. It writes:

- `publication_results/journal_rerun_2026-08-29-matchbox-pe-lr-anchor-first/fig_emitted_set_upset.json`
  with full unrounded values and source provenance;
- the corresponding exact-count CSV; and
- publication-ready SVG, PDF, and PNG renderings.

The plotted percentages exactly reproduce the current accuracy artifacts:
SPLiT-seq PE `81.27/11.18/5.40/0.06/1.12/0.40/0.57`, LR-SPLiT-seq dual
`1.46/76.01/0.05/0.00/21.85/0.64/0.00`, 10x Chromium v2
`100.00/0/0/0/0/0/0`, and sci-RNA-seq3
`98.04/1.04/0.35/0.00/0.58/0.00/0.00` in matrix-column order. The panel union
denominators appear in the plot; exact intersection counts remain in the
machine-readable outputs.

## Table 4 — V_total sizes (`tab:v_total_splitseq`)

**Regenerate:**
```bash
python3 scripts/data_config.py
python3 scripts/edit_tolerant_validity.py $SEQPROC_DATA_DIR/SRR6750041_R2.fastq --chem pe --max-linker-edit 6 --out results/vtotal/v_total_pe_ids.txt
python3 scripts/edit_tolerant_validity.py $SEQPROC_DATA_DIR/SRR6750041_R2.fastq --chem pe --max-linker-edit 6
python3 scripts/edit_tolerant_validity.py $SEQPROC_DATA_DIR/SRR13948564_full.fastq --chem lr --max-linker-edit 6 --out results/vtotal/v_total_lr_ids.txt
python3 scripts/edit_tolerant_validity.py $SEQPROC_DATA_DIR/SRR13948564_full.fastq --chem lr --max-linker-edit 6
python3 scripts/validate_10x_sci_recall_against_vtotal.py --chemistry 10x --outdir results/vtotal_10x
python3 -c "print(round(100*234382218/234382218,2))"  # validator: python3 scripts/validate_10x_sci_recall_against_vtotal.py --chemistry 10x --outdir results/vtotal_10x
python3 scripts/validate_10x_sci_recall_against_vtotal.py --chemistry sci --outdir results/vtotal_sci
```
| Datum | Value | Paper loc | Source script → output | Repro |
|---|---|---|---|---|
| SPLiT-seq PE accession = SRR6750041 | `SRR6750041` | tab:v_total_splitseq, row 1 (SPLiT-seq PE… | scripts/data_config.py (SRA_INFO key 'SRR6750041', line 72); also README dataset table (l… → none (constant in data_config.py:SRA_INFO); echoed by … | ✅ |
| SPLiT-seq PE /V_total/ = 60,675,548 | `60675548` | tab:v_total_splitseq, row 1 (SPLiT-seq PE… | scripts/edit_tolerant_validity.py (--chem pe: edlib HW infix linker search, edit<=6, barc… → stdout JSON field 'valid' (also writes ID list to --ou… | ✅ |
| SPLiT-seq PE % of Total = 78.17% | `78.17%` | tab:v_total_splitseq, row 1 (SPLiT-seq PE… | scripts/edit_tolerant_validity.py (--chem pe) stdout field 'pct_of_scanned' = 100*valid/t… → stdout JSON field 'pct_of_scanned'. NOT pinned in any … | ✅ |
| LR-SPLiT-seq accession = SRR13948564 | `SRR13948564` | tab:v_total_splitseq, row 2 (LR-SPLiT-seq… | scripts/data_config.py (SRA_INFO key 'SRR13948564', line 77); also README dataset table (… → none (constant in data_config.py:SRA_INFO) | ✅ |
| LR-SPLiT-seq /V_total/ = 601,603 | `601603` | tab:v_total_splitseq, row 2 (LR-SPLiT-seq… | scripts/edit_tolerant_validity.py (--chem lr: edlib HW infix linker search edit<=6, BOTH … → stdout JSON field 'valid' (+ --out ID list). NOT pinne… | ✅ |
| LR-SPLiT-seq % of Total = 10.44% | `10.44%` | tab:v_total_splitseq, row 2 (LR-SPLiT-seq… | scripts/edit_tolerant_validity.py (--chem lr) stdout 'pct_of_scanned' = 100*601603/5,764,… → stdout JSON field 'pct_of_scanned'. NOT pinned in any … | ✅ |
| 10x Chromium v2 accession = SRR8315379 | `SRR8315379` | tab:v_total_splitseq, row 3 (10x Chromium… | scripts/data_config.py (SRA_INFO key 'SRR8315379', line 82); also README dataset table (l… → none (constant in data_config.py:SRA_INFO) | ✅ |
| 10x Chromium v2 /V_total/ = 234,382,218 | `234382218` | tab:v_total_splitseq, row 3 (10x Chromium… | scripts/validate_10x_sci_recall_against_vtotal.py --chemistry 10x (calls run_paper_benchm… → results/vtotal_10x/vtotal_recall_10x.json (field 'v_to… | ✅ |
| 10x Chromium v2 % of Total = 100.00% | `100.00%` | tab:v_total_splitseq, row 3 (10x Chromium… | Derived: /V_total//total = 234,382,218 / 234,382,218. Denominator is the ENA-corrected to… → computed by hand from vtotal_recall_10x.json 'v_total'… | ✅ |
| sci-RNA-seq3 accession = SRR7827254 | `SRR7827254` | tab:v_total_splitseq, row 4 (sci-RNA-seq3… | scripts/data_config.py (SRA_INFO key 'SRR7827254', line 87); also README dataset table (l… → none (constant in data_config.py:SRA_INFO) | ✅ |
| sci-RNA-seq3 /V_total/ = 19,767,975 | `19767975` | tab:v_total_splitseq, row 4 (sci-RNA-seq3… | scripts/validate_10x_sci_recall_against_vtotal.py --chemistry sci (calls run_paper_benchm… → results/vtotal_sci/vtotal_recall_sci.json (field 'v_to… | ✅ |
| sci-RNA-seq3 % of Total = 89.49% | `89.49%` | tab:v_total_splitseq, row 4 (sci-RNA-seq3… | Derived: 19,767,975 / 22,088,821. Denominator 22,088,821 = data_config.py:SRA_INFO['SRR78… → computed by hand; validate_10x_sci_recall_against_vtot… | ✅ |
| V_total anchored-linker search tolerance = edit distance 6 (edlib), b… | `6` | Supplementary Note S1, line 362 ('within … | scripts/edit_tolerant_validity.py (ap default --max-linker-edit 6, line 80; edlib.align m… → n/a (algorithm parameter, not an emitted datum) | ✅ |
| linker length = 30 bp; 10x R1 = 26 bp (16 bp barcode + 10 bp UMI); sc… | `30 bp / 26 bp (16+10) / Hamming<=1` | Supplementary Note S1, line 362 | 30bp linker: edit_tolerant_validity.py LINKERS (line 26-27) and run_paper_benchmarks.py L… → n/a (structural parameters embedded in analyzer code) | ✅ |
| Validator throughput ~35,600 reads/second; V_total_LR compute ~162 s … | `35,600 reads/s; 162 s; 36 min` | Supplementary Note S1, line 388 | NOT REPRODUCIBLE as exact figures — these are single-run timing/throughput estimates for … → none (printed to stdout only; not persisted). 35,600 r… | 🔴 NOT REPRODUCIBLE |

## Table 5 — per-cell-type Jaccard (`tab:jaccard`)

**Regenerate:**
```bash
biological_analysis/run_downstream.sh --r1 <R1.fastq> --r2 <R2.fastq> --genome <STAR_INDEX> --outdir downstream_out --threads 8 --min-umi 200   # (internally runs: python biological_analysis/scripts/biological_analysis.py downstream_out/analysis 200 seqproc:downstream_out/sp_Solo.out/Gene splitcode:downstream_out/sc_Solo.out/Gene matchbox:downstream_out/mb_Solo.out/Gene)
biological_analysis/run_downstream.sh --r1 <R1.fastq> --r2 <R2.fastq> --genome <STAR_INDEX> --outdir downstream_out --threads 8 --min-umi 200
python biological_analysis/scripts/jaccard_confusion.py downstream_out/analysis 200 seqproc:downstream_out/sp_Solo.out/Gene splitcode:downstream_out/sc_Solo.out/Gene matchbox:downstream_out/mb_Solo.out/Gene
fasterq-dump SRR6750041 -O $SEQPROC_DATA_DIR  # then biological_analysis/run_downstream.sh --r1 <R1> --r2 <R2> --genome <STAR_INDEX> --outdir downstream_out
```
| Datum | Value | Paper loc | Source script → output | Repro |
|---|---|---|---|---|
| Neuron per-type Jaccard = 0.959 | `0.959` | Table 5 (tab:jaccard), seqproc.tex:401, N… | biological_analysis/scripts/biological_analysis.py (per-cell-type Jaccard loop, lines 127… → biological_analysis/full_run_results/biological_metric… | ✅ |
| Astrocyte per-type Jaccard = 0.912 | `0.912` | Table 5 (tab:jaccard), seqproc.tex:402, A… | biological_analysis/scripts/biological_analysis.py (lines 127-138) → biological_analysis/full_run_results/biological_metric… | ✅ |
| Oligodendrocyte per-type Jaccard = 0.930 | `0.930` | Table 5 (tab:jaccard), seqproc.tex:403, O… | biological_analysis/scripts/biological_analysis.py (lines 127-138) → biological_analysis/full_run_results/biological_metric… | ✅ |
| OPC per-type Jaccard = 0.426 (bold, low-agreement type) | `0.426` | Table 5 (tab:jaccard), seqproc.tex:404, O… | biological_analysis/scripts/biological_analysis.py (lines 127-138) → biological_analysis/full_run_results/biological_metric… | ✅ |
| Microglia per-type Jaccard = 0.361 (bold, lowest-agreement type) | `0.361` | Table 5 (tab:jaccard), seqproc.tex:405, M… | biological_analysis/scripts/biological_analysis.py (lines 127-138) → biological_analysis/full_run_results/biological_metric… | ✅ |
| Endothelial per-type Jaccard = 1.000 (perfect agreement) | `1.000` | Table 5 (tab:jaccard), seqproc.tex:406, E… | biological_analysis/scripts/biological_analysis.py (lines 127-138) → biological_analysis/full_run_results/biological_metric… | ✅ |
| Mean per-type Jaccard = 0.764 | `0.764` | Table 5 (tab:jaccard) 'Mean' row, seqproc… | biological_analysis/scripts/biological_analysis.py (mean_ct_jac, line 138) → biological_analysis/full_run_results/biological_metric… | ✅ |
| Neuron shared-cell count (Cells column) = 142 | `142` | Table 5 (tab:jaccard), seqproc.tex:401, N… | biological_analysis/scripts/jaccard_confusion.py (consensus per-shared-cell type counts, … → downstream_out/analysis/jaccard_confusion.json (consensus… | ✅ |
| Astrocyte shared-cell count (Cells column) = 28 | `28` | Table 5 (tab:jaccard), seqproc.tex:402, A… | biological_analysis/scripts/jaccard_confusion.py (lines 36-40) → downstream_out/analysis/jaccard_confusion.json (consensus… | ✅ |
| Oligodendrocyte shared-cell count (Cells column) = 19 | `19` | Table 5 (tab:jaccard), seqproc.tex:403, O… | biological_analysis/scripts/jaccard_confusion.py (lines 36-40) → downstream_out/analysis/jaccard_confusion.json (consensus… | ✅ |
| OPC shared-cell count (Cells column) = 17 | `17` | Table 5 (tab:jaccard), seqproc.tex:404, O… | biological_analysis/scripts/jaccard_confusion.py (lines 36-40) → downstream_out/analysis/jaccard_confusion.json (consensus… | ✅ |
| Microglia shared-cell count (Cells column) = 9 | `9` | Table 5 (tab:jaccard), seqproc.tex:405, M… | biological_analysis/scripts/jaccard_confusion.py (lines 36-40) → downstream_out/analysis/jaccard_confusion.json (consensus… | ✅ |
| Endothelial shared-cell count (Cells column) = 5 | `5` | Table 5 (tab:jaccard), seqproc.tex:406, E… | biological_analysis/scripts/jaccard_confusion.py (lines 36-40) → downstream_out/analysis/jaccard_confusion.json (consensus… | ✅ |
| Table covers the 220 cells called by all three tools (caption) | `220` | Table 5 (tab:jaccard) caption, seqproc.te… | biological_analysis/scripts/biological_analysis.py (shared = set.intersection of per-tool… → biological_analysis/full_run_results/biological_metric… | ✅ |
| Dataset = full SPLiT-seq PE, SRR6750041 (caption) | `SRR6750041` | Table 5 (tab:jaccard) caption, seqproc.te… | biological_analysis/run_downstream.sh (pipeline driver; STARsolo Gene matrices feeding biolo… → N/A (input accession; downstream outputs land in biolo… | ✅ |
| 'every other type is at or above 0.91' (prose claim about the four no… | `>= 0.91 (min of {0.959, 0.912, 0.…` | Supplementary Note S2, seqproc.tex:393 ('… | biological_analysis/scripts/biological_analysis.py (per-type Jaccard, lines 127-138) → biological_analysis/full_run_results/biological_metric… | ✅ |
| Of shared cells any tool labels microglia, dissenting tools relabel 2… | `22 of 25` | Supplementary Note S2, seqproc.tex:393 | biological_analysis/scripts/jaccard_confusion.py (confusion_by_type loop, lines 43-52: ce… → downstream_out/analysis/jaccard_confusion.json (confusion… | ✅ |
| Of cells labeled OPC, 19 of 28 dissenting instances are microglia | `19 of 28` | Supplementary Note S2, seqproc.tex:393 | biological_analysis/scripts/jaccard_confusion.py (confusion_by_type loop, lines 43-52: di… → downstream_out/analysis/jaccard_confusion.json (confusion… | ✅ |

## Appendix B — LR optimization (`tab:lr_progression`)

**Regenerate:**
```bash
cd /home/ubuntu/combine-lab/seqproc && cargo bench
cat /home/ubuntu/seqproc-paper-analysis-clean/configs/seqproc/10x_v2.geom
n/a (design/profiling claim in the Rust engine)
cd /home/ubuntu/combine-lab/seqproc && samply record ./target/release/seqproc ... (profile, not scripted)
cd /home/ubuntu/combine-lab/seqproc && cargo bench -- 10x_N=1000000
cd /home/ubuntu/combine-lab/seqproc && cargo bench -- 10x_N=10000000
grep -rn 1024 /home/ubuntu/combine-lab/ANTISEQUENCE/src /home/ubuntu/combine-lab/seqproc/src
cd /home/ubuntu/combine-lab/seqproc && cargo bench (with chunk/thread env sweep)
```
| Datum | Value | Paper loc | Source script → output | Repro |
|---|---|---|---|---|
| Criterion collects 50--100 samples per benchmark (B.1) | `50--100 samples` | Section B.1, first paragraph | NOT REPRODUCIBLE from this repo — Criterion.rs sample_size is set in the Rust micro-bench… → combine-lab/seqproc/target/criterion/ (not in this rep… | 🔴 NOT REPRODUCIBLE |
| Benchmark scales: 1K, 10K, 1M, 10M synthetic reads (B.1) | `1,000 / 10,000 / 1,000,000 / 10,0…` | Section B.1, second paragraph | NOT REPRODUCIBLE from this repo — read counts are literals in the Rust bench (e.g. make_f… → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| 10x Chromium v2 trivial geometry string 1{b[16]u[10]}2{r:} (B.1) | `1{b[16]u[10]}2{r:}` | Section B.1, second paragraph | configs/seqproc/10x_v2.geom (this repo) and hard-coded in /home/ubuntu/combine-lab/seqpro… → configs/seqproc/10x_v2.geom | ✅ |
| 4--8 malloc/free cycles per read for a paired-end workload (B.2) | `4--8 malloc/free per read` | Section B.2, first paragraph | NOT REPRODUCIBLE — narrative characterization of the pre-optimization ANTISEQUENCE Read s… → none | 🔴 NOT REPRODUCIBLE |
| ~48% of CPU samples in condition-variable waits and allocator functio… | `~48%` | Section B.2, first paragraph | NOT REPRODUCIBLE from this repo — from a samply CPU sampling profile of ANTISEQUENCE visu… → none shipped | 🔴 NOT REPRODUCIBLE |
| Batched read recycling: 4.8x cumulative speedup, 1.21 s -> 252 ms (10… | `4.8x (1.21 s -> 252 ms)` | Section B.2, fourth paragraph (also resta… | NOT REPRODUCIBLE from this repo — Criterion micro-benchmark bench_10x_large (10x_N=100000… → combine-lab/seqproc/target/criterion/ (no shipped refe… | 🔴 NOT REPRODUCIBLE |
| Recycling pass alone yields 38% reduction at 10M reads (B.2) | `38% reduction` | Section B.2, fourth paragraph | NOT REPRODUCIBLE from this repo — Criterion bench_10x_10m (10x_N=10000000) in /home/ubunt… → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| Recycling batch / first-chunk size of 1,024 reads (B.2) | `1,024` | Section B.2 (para 2 and 4) | NOT REPRODUCIBLE from this repo — chunk-size constant in the ANTISEQUENCE/seqproc engine … → none (engine constant) | 🔴 NOT REPRODUCIBLE |
| SmallVec<[Mapping; 4]> (inline capacity 4) yields ~10% improvement (B… | `~10% (inline capacity 4)` | Section B.2, last paragraph (also Table t… | NOT REPRODUCIBLE from this repo — Criterion benchmark delta in combine-lab/seqproc/benche… → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| In-place SetOp + constant folding: ~5% improvement (B.2 / B.8) | `~5%` | Section B.2 last paragraph; Table tab:opt… | NOT REPRODUCIBLE from this repo — Criterion benchmark delta in combine-lab/seqproc/benche… → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| Original input-reader chunk size of 256 reads; parking_lot::Mutex swa… | `256 -> 1,024; ~10%` | Section B.3, second paragraph (also Table… | NOT REPRODUCIBLE from this repo — synchronization change measured by Criterion in combine… → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| Chunk-size x thread sweep: 1,024 chunk with 2--4 threads ~35% improve… | `~35% (chunk 1,024, 2--4 threads)` | Section B.3, third paragraph | NOT REPRODUCIBLE from this repo — manual Criterion sweep over chunk sizes {256,512,1024} … → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| Catastrophic 2x regression at chunk size 512 with 4 threads; ~125 us … | `2x regression; ~125 us/chunk` | Section B.3, third paragraph (regression … | NOT REPRODUCIBLE from this repo — Criterion/profiling result in combine-lab/seqproc/bench… → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| 8 threads no improvement over 4; ~244 chunks per thread (1,024-read c… | `8 threads: no gain; ~244 chunks/t…` | Section B.3, last paragraph (also Table t… | NOT REPRODUCIBLE from this repo (thread-scaling) — Criterion in combine-lab/seqproc. The … → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| FASTQ serialization/file output ~6% of total runtime (NullOutputOp is… | `~6% (=> 3% max total gain)` | Section B.4, first paragraph | NOT REPRODUCIBLE from this repo — measured via the NullOutputOp graph node / ANTISEQ_STUB… → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| Thread-local writer cache: 2--5% improvement (B.4) | `2--5%` | Section B.4, second paragraph (also Table… | NOT REPRODUCIBLE from this repo — Criterion delta in combine-lab/seqproc; FxHashMap write… → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| FASTQ record batching: 20--30% regression (failed experiment) (B.4) | `20--30% regression` | Section B.4, third paragraph (also Table … | NOT REPRODUCIBLE from this repo — negative Criterion result in combine-lab/seqproc; the c… → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| Vectored I/O (writev/IoSlice): no measurable improvement (B.4) | `no effect` | Section B.4, last paragraph (also Table t… | NOT REPRODUCIBLE from this repo — Criterion null result in combine-lab/seqproc. → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| Myers' bit-vector for patterns <=64 bp; Ukkonen DP for >64 bp; adapti… | `<=64 bp / >64 bp / <4 patterns / …` | Section B.5, itemized list and adaptive-s… | NOT REPRODUCIBLE from this repo — algorithm parameters/thresholds in the ANTISEQUENCE mat… → none (engine source constants) | 🔴 NOT REPRODUCIBLE |
| Seed-only edit distance 2.4x SLOWER than Hamming on LR-SPLiT-seq 1M: … | `2.4x slower (38.7 s vs 16.2 s)` | Section B.5, 'dramatic failure' paragraph | NOT REPRODUCIBLE from this repo — historical timing of the now-removed seed-only edit-dis… → none shipped | 🔴 NOT REPRODUCIBLE |
| Adaptive edit distance 2.3x FASTER than Hamming (2.0 s vs 4.7 s); ada… | `2.3x faster (2.0 s vs 4.7 s); 16.…` | Section B.5, 'dramatic failure' paragraph… | NOT REPRODUCIBLE as-is from this repo — single-linker forward-only edit-vs-Hamming timing… → combine-lab/seqproc/target/criterion/ (this repo: resu… | 🔴 NOT REPRODUCIBLE |
| Edit vs Hamming recovery gain: +15.8% LR-SPLiT-seq (B.5) | `+15.8%` | Section B.5, last paragraph | scripts/concordance_analysis.py — field edit_gain_pct (line 517: (edit_only - ham_only)/l… → results/concordance/lr_splitseq/results.json (hamming_… | ✅ |
| Edit vs Hamming recovery gain: +4.8% SPLiT-seq PE (B.5) | `+4.8%` | Section B.5, last paragraph | scripts/concordance_analysis.py — edit_gain_pct for splitseq_pe (splitseq_filter_edit.geo… → results/concordance/splitseq_pe/results.json (hamming_… | ✅ |
| Edit vs Hamming recovery gain: +0.9% sci-RNA-seq3 (B.5) | `+0.9%` | Section B.5, last paragraph | scripts/concordance_analysis.py — edit_gain_pct for sciseq3 (sciseq3_edit.geom vs sciseq3… → results/concordance/sciseq3/results.json (hamming_vs_e… | ✅ |
| Long reads arrive ~50/50 forward vs reverse-complement; forward-only … | `~50/50; ~half` | Section B.6, first paragraph | NOT REPRODUCIBLE (as a formal statistic) — data characteristic of PacBio/Nanopore; implic… → none | 🔴 NOT REPRODUCIBLE |
| Table tab:lr_progression row 1 — Forward-only, Hamming: 22.0% recover… | `22.0% / 4.5 s` | Table tab:lr_progression (Appendix B.6), … | NOT AUTO-PRODUCED by any shipped script — requires running seqproc with configs/seqproc/s… → no shipped JSON (config exists: configs/seqproc/splits… | ✅ |
| Table tab:lr_progression row 2 — Forward-only, edit distance: 25.8% r… | `25.8% / 2.8 s` | Table tab:lr_progression (Appendix B.6), … | NOT AUTO-PRODUCED by any shipped script — requires seqproc with configs/seqproc/splitseq_… → no shipped JSON (config: configs/seqproc/splitseq_sing… | ✅ |
| Table tab:lr_progression row 3 — Orientation-aware, Hamming: 43.1% re… | `43.1% / 8.4 s` | Table tab:lr_progression (Appendix B.6), … | NOT AUTO-PRODUCED by any shipped script — requires seqproc with configs/seqproc/splitseq_… → no shipped JSON (config: configs/seqproc/splitseq_sing… | ✅ |
| Table tab:lr_progression row 4 — Orientation-aware, edit distance: 49… | `49.9% / 5.3 s` | Table tab:lr_progression (Appendix B.6), … | scripts/lr_perf_rerun.py — seqproc block (lines 92-123): runs configs/seqproc/splitseq_si… → results/lr_perf/lr_splitseq_perf_results.json (tools.s… | ✅ |
| Table tab:lr_progression row 5 — matchbox, forward-only: 24.0% recove… | `24.0% / 3.5 s` | Table tab:lr_progression (Appendix B.6), … | NOT AUTO-PRODUCED by any shipped script — requires matchbox with configs/matchbox/splitse… → no shipped JSON (config: configs/matchbox/splitseq_sin… | ✅ |
| Table tab:lr_progression row 6 — matchbox, dual-orientation: 39.7% re… | `39.7% / 7.6 s` | Table tab:lr_progression (Appendix B.6), … | scripts/lr_perf_rerun.py — matchbox block (lines 125-164): runs configs/matchbox/splitseq… → results/lr_perf/lr_splitseq_perf_results.json (tools.m… | ✅ |
| Table tab:lr_progression caption params — 1M-read subset of SRR139485… | `1,000,000 reads / 4 threads / e=0…` | Table tab:lr_progression caption (Appendi… | scripts/lr_perf_rerun.py — LR_FASTQ = data/SRR13948564_1M.fastq (line 40), --threads defa… → results/lr_perf/lr_splitseq_perf_results.json (dataset… | ✅ |
| seqproc (both features) recovers 101K more reads than matchbox dual-o… | `+101K reads; 3 unique to matchbox` | Table tab:lr_progression caption (Appendi… | PARTIALLY REPRODUCIBLE — the emit-ID sets are produced by scripts/lr_perf_rerun.py (seqpr… → read-ID sets in results/lr_perf/seqproc_rep1_R1.fq and… | ✅ |
| Manual-RC investigation — union of forward + manual-RC geometries rec… | `406,446 / 430,807 / +24,361 (6%)` | Section B.6, last paragraph | NOT REPRODUCIBLE from this repo — the 'manually written reverse-complement geometry' is N… → none shipped | 🔴 NOT REPRODUCIBLE |
| #[match_ori(either)] annotation used by TryOrientationOp (B.6) | `#[match_ori(either)]` | Section B.6, last paragraph | configs/seqproc/splitseq_singleend_ann.geom and splitseq_singleend_edit_ann.geom (this re… → configs/seqproc/splitseq_singleend_edit_ann.geom, conf… | ✅ |
| Bit-shift overflow in HammingLookup for patterns longer than 8 bytes … | `8 bytes` | Section B.7, fifth bullet | NOT REPRODUCIBLE from this repo — bug/fix in the ANTISEQUENCE HammingLookup encoding (com… → none (engine code) | 🔴 NOT REPRODUCIBLE |
| Test suite grew from 8 tests (Sep 2025) to 402 tests (Mar 2026), 80.2… | `8 -> 402 tests; 80.2% coverage` | Section B.7, last paragraph | NOT REPRODUCIBLE from this repo — these are the seqproc/ANTISEQUENCE Rust cargo-test coun… → none shipped in this repo | 🔴 NOT REPRODUCIBLE |
| Table tab:optimization_summary impact cells — read recycling -30 to -… | `see per-row values` | Table tab:optimization_summary (Appendix … | NOT REPRODUCIBLE from this repo (except +25.8% and 2.3x which trace to lr_perf_rerun.py p… → combine-lab/seqproc/target/criterion/ (this repo: resu… | 🔴 NOT REPRODUCIBLE |
| Cumulative ~4.8x speedup on 10x geometry (1.21 s -> 252 ms, single th… | `4.8x (1.21 s -> 252 ms)` | Section B.8, penultimate paragraph | NOT REPRODUCIBLE from this repo — Criterion bench_10x_large (10x_N=1000000) in combine-la… → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| ~2.4x speedup on sci-RNA-seq3 geometry (521 ms -> 218 ms) — final sum… | `2.4x (521 ms -> 218 ms)` | Section B.8, penultimate paragraph | NOT REPRODUCIBLE from this repo — Criterion bench_sci3_large (sci3_N=1000000) in combine-… → combine-lab/seqproc/target/criterion/ | 🔴 NOT REPRODUCIBLE |
| With 4 threads on real data, seqproc processes 1M reads from any of t… | `< 6 s (1M reads, 4 threads)` | Section B.8, penultimate paragraph (refer… | scripts/run_paper_benchmarks.py — per-tool 'runtime' field (Table 3 runtime column). Cros… → results/paper_figures/benchmark_results.json (per-data… | ✅ |
| Peak memory footprint 23--35 MB — 6--87x less than alternatives (B.8) | `23--35 MB; 6--87x` | Section B.8, penultimate paragraph and la… | scripts/run_paper_benchmarks.py — 'peak_mem_mb' field (Table 3 Memory column); the 6--87x… → results/paper_figures/benchmark_results.json (per-tool… | ✅ |
| Only ~1,024 reads (one chunk) alive at any time; memory O(1) in input… | `~1,024; O(1)` | Section B.8, last paragraph | NOT REPRODUCIBLE as a measurement — architectural claim about the recycling design (combi… → results/paper_figures/benchmark_results.json (peak_mem… | 🔴 NOT REPRODUCIBLE |
| Seven-month optimization campaign, September 2025 -- March 2026 (B.1 … | `Sep 2025 -- Mar 2026 (7 months)` | Section B (intro paragraph) and B.7 | NOT REPRODUCIBLE — process/timeline statement; verifiable only from the git history of co… → none | 🔴 NOT REPRODUCIBLE |

## Figures

**Regenerate:**
```bash
(none) — illustrative diagram edited by hand; committed as a static PDF asset
biological_analysis/run_downstream.sh --r1 R1.fastq --r2 R2.fastq --genome STAR_INDEX --outdir downstream_out --threads 8 --min-umi 200   # figure step (run_downstream.sh:86): python3 biological_analysis/scripts/count_concordance.py downstream_out/analysis seqproc:downstream_out/sp_Solo.out/Gene splitcode:downstream_out/sc_Solo.out/Gene matchbox:downstream_out/mb_Solo.out/Gene
python3 biological_analysis/scripts/knee_barcoderanks.py seqproc:downstream_out/sp_Solo.out/Gene splitcode:downstream_out/sc_Solo.out/Gene matchbox:downstream_out/mb_Solo.out/Gene   # prints 'inflection: rank <r> (UMI>=...)' per tool
python3 biological_analysis/scripts/knee_barcoderanks.py seqproc:downstream_out/sp_Solo.out/Gene splitcode:downstream_out/sc_Solo.out/Gene matchbox:downstream_out/mb_Solo.out/Gene
python3 biological_analysis/scripts/count_concordance.py downstream_out/analysis seqproc:downstream_out/sp_Solo.out/Gene splitcode:downstream_out/sc_Solo.out/Gene matchbox:downstream_out/mb_Solo.out/Gene
```
| Datum | Value | Paper loc | Source script → output | Repro |
|---|---|---|---|---|
| Figure fig:seqproc_workflow — cartoon of the seqproc workflow (panels… | `Figures/seqproc_flow.pdf (schemat…` | sections/methods.tex:2-7, \includegraphic… | NOT REPRODUCIBLE — hand-authored vector schematic; no script in seqproc-paper-analysis-cl… → /home/ubuntu/paper/Figures/seqproc_flow.pdf (committed… | 🔴 NOT REPRODUCIBLE |
| Figure fig:illustrative_example — how the EFGDL sci-RNA-seq3 spec fro… | `Figures/antisequence_graph.pdf (s…` | sections/methods.tex:292-297, \includegra… | NOT REPRODUCIBLE — hand-authored vector schematic; no script in seqproc-paper-analysis-cl… → /home/ubuntu/paper/Figures/antisequence_graph.pdf (com… | 🔴 NOT REPRODUCIBLE |
| Figure fig:count_concordance — 3-panel quantification-concordance fig… | `fig_count_concordance.pdf (render…` | seqproc.tex:418-424, \includegraphics{fig… | biological_analysis/scripts/count_concordance.py (driven by biological_analysis/run_phase… → biological_analysis/full_run_results/count_concordance… | ✅ |
| Panel A — barcode-rank inflection ranks 257, 255, and 257 for seqproc… | `257 / 255 / 257 (range 255-257)` | seqproc.tex:421 (fig:count_concordance ca… | biological_analysis/scripts/knee_barcoderanks.py (DropletUtils barcodeRanks port, functio… → biological_analysis/full_run_results/count_concordance… | ✅ |
| Panel A — UMI threshold 101 to 109 (UMI value at the barcodeRanks inf… | `101 to 109` | seqproc.tex:416 ('(UMI threshold 101 to 1… | biological_analysis/scripts/knee_barcoderanks.py (barcode_ranks -> infl_umi) → stdout of knee_barcoderanks.py ('UMI>=<infl_umi>'); co… | ✅ |
| Panel A — nonzero/recovered-signal barcode counts 14,434 to 14,936 (n… | `14,434 to 14,936` | seqproc.tex:421 (fig:count_concordance ca… | biological_analysis/scripts/count_concordance.py (knees[n]['n_barcodes'], code line 69) → biological_analysis/full_run_results/count_concordance… | ✅ |
| Panel A — top-barcode UMI 18,847 to 19,644 (splitcode 18,847; seqproc… | `18,847 to 19,644` | results.tex:93 ('top-barcode UMI (18,847 … | biological_analysis/scripts/count_concordance.py (knees[n]['top_umi'], code line 69) → biological_analysis/full_run_results/count_concordance… | ✅ |
| Panel B — per-barcode total UMI Pearson on log1p scale, seqproc vs ot… | `0.978 (seqproc/splitcode) and 0.9…` | seqproc.tex:421 (fig:count_concordance ca… | biological_analysis/scripts/count_concordance.py (per_barcode_umi_pearson_logspace, code … → biological_analysis/full_run_results/count_concordance… | ✅ |
| Panel C — per-gene total UMI Pearson on log1p scale, seqproc vs other… | `0.994 (seqproc/splitcode) and 0.9…` | seqproc.tex:421 (fig:count_concordance ca… | biological_analysis/scripts/count_concordance.py (per_gene_total_pearson_logspace, code l… → biological_analysis/full_run_results/count_concordance… | ✅ |

## PE prose numbers (results.tex)

### Current full split-pipe vendor-set comparison

Build and run split-pipe 1.4.0 as documented in `containers/README.md`. Then
regenerate the full JSON, CSV, and compact vendor bitmap from the fresh
`barcode_head.fastq` and the three final accepted-ID bitmaps:

```bash
python3 scripts/splitpipe_full_concordance.py \
  --splitpipe-fastq /path/to/splitpipe/process/barcode_head.fastq \
  --input-records 77621181 \
  --archived-vendor-ids results_final/splitpipe_valid_ids_10M.txt.gz \
  --archived-records 10000000 \
  --input-r1 /path/to/SRR6750041_1.fastq.gz \
  --input-r2 /path/to/SRR6750041_2.fastq.gz \
  --campaign-input-r1 /path/to/SRR6750041_R1.fastq \
  --campaign-input-r2 /path/to/SRR6750041_R2.fastq \
  --splitpipe-run-def /path/to/splitpipe/process/run_proc_def.json \
  --splitpipe-log /path/to/splitpipe/split-pipe_v1_4_0.log \
  --splitpipe-config configs/split-pipe/splitseq_pe_v1.par
```

The output is `splitseq_pe_splitpipe_vendor_full.{json,csv,vendor.raw}` in the
current journal-results directory. It validates the full split-pipe FASTQ,
verifies every final tool bitmap against the accuracy artifact, records the
full provenance chain, and confirms that the first 10 million fresh calls have
zero symmetric difference from the archive. The historical table below records
the original-preprint workflow and values.

**Regenerate:**
```bash
seqkit head -n 10000000 $SEQPROC_DATA_DIR/SRR6750041_1.fastq > R1_10M.fastq; seqkit head -n 10000000 $SEQPROC_DATA_DIR/SRR6750041_2.fastq > R2_10M.fastq   # (10M subset; not scripted in-repo)
zcat results_final/splitpipe_valid_ids_10M.txt.gz / wc -l   # -> 7539920 (verify the shipped reference set size without split-pipe)
seqproc --geom configs/seqproc/splitseq_filter_edit.geom --file1 R1_10M.fastq --file2 R2_10M.fastq --out1 sp_R1.fq --out2 sp_R2.fq --threads 8 -a configs/seqproc/splitseq_bc3_seq2seq.tsv -a configs/seqproc/splitseq_bc2_seq2seq.tsv -a configs/seqproc/splitseq_bc1_seq2seq.tsv && python3 scripts/verify_vendor_concordance_pe.py --seqproc-r2 sp_R2.fq --splitcode-r2 sc_R2.fq
python3 scripts/verify_vendor_concordance_pe.py --seqproc-r2 sp_R2.fq --splitcode-r2 sc_R2.fq
python3 scripts/verify_vendor_concordance_pe.py --seqproc-r2 sp_R2.fq --splitcode-r2 sc_R2.fq   # (Jacc column); or python3 -c "i=7512336;e=8367328;g=7539920;print(round(i/(e+g-i),4))" -> 0.8949
splitcode -c configs/splitcode/splitseq_paper.config --assign -N 2 -t 8 -m mapping.txt -o sc_R1.fq,sc_R2.fq R1_10M.fastq R2_10M.fastq && python3 scripts/verify_vendor_concordance_pe.py --seqproc-r2 sp_R2.fq --splitcode-r2 sc_R2.fq
python3 -c "print(round(100-82.52,2),'%'); print(9120071-7526156,'reads')"   # -> 17.48 % ; 1593915 reads
python3 scripts/discordant_analysis.py
```
| Datum | Value | Paper loc | Source script → output | Repro |
|---|---|---|---|---|
| 10 million read subset of SRR6750041 (denominator for the 75.4% / 83.… | `10,000,000` | sections/results.tex:64, split-pipe vendo… | NOT REPRODUCIBLE from a committed script — it is a fixed input: the first 10,000,000 read… → n/a (input subset; the shipped ground-truth results_fi… | 🔴 NOT REPRODUCIBLE |
| split-pipe accepted N reads as having a valid barcode (vendor referen… | `7,539,920 (75.4%)` | sections/results.tex:64 ('split-pipe acce… | NOT REPRODUCIBLE — requires the proprietary Parse Biosciences split-pipe tool (its proces… → results_final/vendor_concordance_pe_10M.json (field "s… | 🔴 NOT REPRODUCIBLE |
| seqproc emitted reads on the 10M subset | `8,367,328 (83.7%)` | sections/results.tex:64 ('seqproc emitted… | scripts/verify_vendor_concordance_pe.py (counts seqproc's R2 output; open-source, reprodu… → results_final/vendor_concordance_pe_10M.json (tools.se… | ✅ |
| seqproc emitted reads that were in the split-pipe valid set (intersec… | `7,512,336` | sections/results.tex:64 ('of which 7{,}51… | scripts/verify_vendor_concordance_pe.py (intersects seqproc R2 IDs with shipped GT set). … → results_final/vendor_concordance_pe_10M.json (tools.se… | 🔴 NOT REPRODUCIBLE |
| seqproc precision vs split-pipe vendor reference | `89.8%` | sections/results.tex:64 ('giving a precis… | scripts/verify_vendor_concordance_pe.py (P = intersection/emit). Committed value in resul… → results_final/vendor_concordance_pe_10M.json (tools.se… | ✅ |
| seqproc recall vs split-pipe vendor reference | `99.6%` | sections/results.tex:64 ('a recall of 99.… | scripts/verify_vendor_concordance_pe.py (R = intersection//split-pipe valid/). Committed … → results_final/vendor_concordance_pe_10M.json (tools.se… | ✅ |
| seqproc F1 vs split-pipe vendor reference | `0.944` | sections/results.tex:64 ('an F1 score of … | scripts/verify_vendor_concordance_pe.py (F1 = 2PR/(P+R)). Committed value in results_fina… → results_final/vendor_concordance_pe_10M.json (tools.se… | ✅ |
| seqproc Jaccard index vs split-pipe vendor reference | `0.895` | sections/results.tex:64 ('and a Jaccard i… | scripts/verify_vendor_concordance_pe.py (J = intersection/union; printed by the script). … → NOT in JSON (derived); verify script prints Jacc=0.895… | ✅ |
| splitcode emitted reads on the 10M subset | `9,120,071 (91.2%)` | sections/results.tex:66 ('splitcode emitt… | scripts/verify_vendor_concordance_pe.py (counts splitcode R2 output; open-source, reprodu… → results_final/vendor_concordance_pe_10M.json (tools.sp… | ✅ |
| splitcode recall vs split-pipe vendor reference | `99.8%` | sections/results.tex:66 ('recovered a com… | scripts/verify_vendor_concordance_pe.py. Committed value in results_final/vendor_concorda… → results_final/vendor_concordance_pe_10M.json (tools.sp… | ✅ |
| splitcode precision vs split-pipe vendor reference | `82.5%` | sections/results.tex:66 ('lower precision… | scripts/verify_vendor_concordance_pe.py. Committed value in results_final/vendor_concorda… → results_final/vendor_concordance_pe_10M.json (tools.sp… | ✅ |
| splitcode F1 vs split-pipe vendor reference | `0.903` | sections/results.tex:66 ('F1 0.903') | scripts/verify_vendor_concordance_pe.py. Committed value in results_final/vendor_concorda… → results_final/vendor_concordance_pe_10M.json (tools.sp… | ✅ |
| fraction of splitcode output falling outside the split-pipe valid set | `About 17.5% (roughly 1.6 million …` | sections/results.tex:66 ('About 17.5\% of… | Derived from scripts/verify_vendor_concordance_pe.py outputs (100 - precision, and emit -… → results_final/vendor_concordance_pe_10M.json (derived … | ✅ |
| reads unique to splitcode on the full SPLiT-seq PE dataset (discordan… | `5,787,314` | sections/results.tex:60 ('Of the 5{,}787{… | scripts/discordant_analysis.py (writes splitcode_unique.total). Also surfaced by scripts/… → results/concordance/splitseq_pe/discordant_analysis.js… | ✅ |
| percentage of splitcode-unique reads that are structurally invalid | `99.2%` | sections/results.tex:60 ('99.2\% are stru… | scripts/discordant_analysis.py (structurally_invalid / total for splitcode_unique). Repro… → results/concordance/splitseq_pe/discordant_analysis.js… | ✅ |
| pairwise Jaccard index of seqproc and matchbox emitted read sets (SPL… | `0.926` | sections/results.tex:60 ('seqproc and mat… | scripts/concordance_analysis.py (compute_concordance -> pairwise jaccard). Reproducible. → results/concordance/concordance_results.json (splitseq… | ✅ |
| reads unique to matchbox (seqproc superset of matchbox output, SPLiT-… | `676` | sections/results.tex:60 ('only 676 reads … | scripts/concordance_analysis.py (compute_concordance -> matchbox_only for the matchbox/se… → results/concordance/concordance_results.json (splitseq… | ✅ |
| size of matchbox's SPLiT-seq PE output ('out of 60.8 million') | `60.8 million (60,828,797)` | sections/results.tex:60 ('out of 60.8 mil… | scripts/concordance_analysis.py (id-set size for matchbox = concordance.tools.matchbox / … → results/concordance/concordance_results.json (splitseq… | ✅ |

## LR / sci / 10x prose numbers (results.tex)

**Regenerate:**
```bash
python3 scripts/concordance_analysis.py --reads full --threads 32 --datasets lr_splitseq
python3 scripts/data_config.py --reads full  # prints resolved LR path + count
python3 scripts/edit_tolerant_validity.py $SEQPROC_DATA_DIR/SRR13948564_full.fastq --chem lr --max-linker-edit 6 --out results/lr_recall_vtotal/lr_valid_ids.txt
python3 scripts/edit_tolerant_validity.py $SEQPROC_DATA_DIR/SRR13948564_full.fastq --chem lr --out results/lr_recall_vtotal/lr_valid_ids.txt
n/a (hard-coded linker constant; length 30)
python3 scripts/edit_tolerant_validity.py $SEQPROC_DATA_DIR/SRR13948564_full.fastq --chem lr --out results/lr_recall_vtotal/lr_valid_ids.txt  # then: comm/set-intersect lr_valid_ids.txt with results/concordance/lr_splitseq/seqproc_edit_ids.txt, divide by 601603
python3 scripts/splitcode_lr_dual_validate.py --dataset full --threads 8 --outdir results/splitcode_lr_dual  # then intersect splitcode_dual_combined_out.fq IDs with edit_tolerant lr_valid_ids.txt / 601603
python3 scripts/edit_tolerant_validity.py $SEQPROC_DATA_DIR/SRR13948564_full.fastq --chem lr --out results/lr_recall_vtotal/lr_valid_ids.txt  # then intersect results/concordance/lr_splitseq/matchbox_ids.txt / 601603
```
| Datum | Value | Paper loc | Source script → output | Repro |
|---|---|---|---|---|
| LR-SPLiT-seq pairwise Jaccard index of 0.796 (seqproc vs matchbox) | `0.796` | sections/results.tex:70 (LR prose): "a pa… | scripts/concordance_analysis.py (compute_concordance(); jaccard() at line 116; pairwise e… → results/concordance/lr_splitseq/results.json  (pairwis… | ✅ |
| Only 279 reads unique to matchbox (out of 2.28 million) on full LR | `279` | sections/results.tex:70: "only 279 reads … | scripts/concordance_analysis.py (compute_concordance -> pairwise 'matchbox_only' = /match… → results/concordance/lr_splitseq/results.json (pairwise… | ✅ |
| matchbox recovers 2.28 million reads on full LR (denominator for the … | `2.28 million (~2,282,711)` | sections/results.tex:70: "out of 2.28 mil… | scripts/concordance_analysis.py (result['tools']['matchbox'] = len(matchbox emit ID set) … → results/concordance/lr_splitseq/results.json (tools.ma… | ✅ |
| Full LR-SPLiT-seq dataset = 5.76 million reads (5,764,421) | `5,764,421` | sections/results.tex:70 ("5.76 million re… | scripts/data_config.py (SRA_INFO['SRR13948564']['full_reads'] = 5_764_421); verified at r… → results/paper_figures/benchmark_results.json (lr_split… | ✅ |
| LR genuine reference set (V_total) = 601,603 reads | `601,603` | sections/results.tex:68: "a genuine refer… | scripts/edit_tolerant_validity.py (--chem lr: edlib infix linker match within edit distan… → stdout JSON {"valid": 601603, "pct_of_scanned": 10.44}… | ✅ |
| LR V_total is 10.44% of input | `10.44%` | sections/results.tex:68 ("which is 10.44%… | scripts/edit_tolerant_validity.py (--chem lr; prints pct_of_scanned = 100*valid/total) → stdout JSON pct_of_scanned=10.44 (= 601603/5764421) | ✅ |
| V_total requires the 30 bp linker located within a small edit distance | `30 bp` | sections/results.tex:68: "its 30 bp linke… | scripts/edit_tolerant_validity.py:27 (LINKERS['lr'][0] = 'GTGGCCGATGTTTCGCATCGGCGTACGACT'… → n/a (config constant, not an emitted metric) | ✅ |
| seqproc recovers 98.72% of LR V_total (593,873 reads) | `98.72% (593,873)` | sections/results.tex:68 ("\seqproc recove… | scripts/edit_tolerant_validity.py (produces the 601,603 valid_ids denominator) INTERSECTE… → NOT emitted as JSON by a committed script. Inputs pres… | ✅ |
| splitcode (dual-pass) recovers 98.22% of LR V_total (590,883 reads) | `98.22% (590,883)` | sections/results.tex:68 ("\splitcode 98.2… | scripts/edit_tolerant_validity.py (601,603 valid_ids) INTERSECTED with splitcode dual-pas… → NOT emitted as JSON by a committed script. Inputs: res… | ✅ |
| matchbox recovers 87.97% of LR V_total (529,242 reads) | `87.97% (529,242)` | sections/results.tex:68 ("\matchbox 87.97… | scripts/edit_tolerant_validity.py (601,603 valid_ids) INTERSECTED with matchbox LR emit I… → NOT emitted as JSON by a committed script. Inputs: res… | ✅ |
| seqproc emits 49.7% of LR input | `49.7%` | sections/results.tex:68 ("\seqproc emits … | scripts/run_paper_benchmarks.py (recovery_rate = mean reads_out / total_reads * 100, comp… → results/paper_figures/benchmark_results.json (lr_split… | ✅ |
| matchbox emits 39.6% of LR input | `39.6%` | sections/results.tex:68 ("\matchbox 39.6%… | scripts/run_paper_benchmarks.py (recovery_rate for matchbox on lr_splitseq) → results/paper_figures/benchmark_results.json (lr_split… | ✅ |
| splitcode emits 27.5% of LR input in a single forward pass | `27.5%` | sections/results.tex:68: "\splitcode ... … | scripts/splitcode_lr_dual_validate.py (forward.emitted / total_input_reads) → results/splitcode_lr_dual/splitcode_lr_dual_results.js… | ✅ |
| splitcode emits 51.8% of LR input as the union of forward + reverse-c… | `51.8%` | sections/results.tex:68 ("51.8% as the un… | scripts/splitcode_lr_dual_validate.py (dual_pass.emitted_union / total_input_reads; print… → results/splitcode_lr_dual/splitcode_lr_dual_results.js… | ✅ |
| Validation script processed full LR input in 162 seconds | `162 s` | sections/results.tex:54: "processed the f… | scripts/edit_tolerant_validity.py (or the V_total validity analyzer) wall-clock; validate… → not persisted to JSON (wall-clock timing printed to st… | ✅ |
| Validation throughput approximately 35,600 reads per second (LR) | `~35,600 reads/s` | sections/results.tex:54: "(approximately … | derived: 5,764,421 / 162 s (scripts/edit_tolerant_validity.py runtime) → not persisted (derived from timing) | ✅ |
| Validation script is roughly 15x slower than seqproc on LR | `~15x` | sections/results.tex:54: "roughly 15x slo… | derived: 162 s (validation script) / 11.1 s (seqproc LR runtime from run_paper_benchmarks… → results/paper_figures/benchmark_results.json (lr_split… | ✅ |
| seqproc LR-SPLiT-seq runtime = 11.1 s | `11.1 s` | sections/results.tex:54 ("\seqproc's 11.1… | scripts/run_paper_benchmarks.py (mean runtime over 3 replicates) → results/paper_figures/benchmark_results.json (lr_split… | ✅ |
| 10x Chromium v2: all three tools achieve 100% recovery | `100%` | sections/results.tex:60 ("all three tools… | scripts/run_paper_benchmarks.py (recovery_rate = 100 for 10x) + scripts/validate_10x_sci_… → results/paper_figures/benchmark_results.json (tenx rec… | ✅ |
| sci-RNA-seq3: pairwise Jaccard index exceeds 0.985 for all tool pairs | `>0.985` | sections/results.tex:60: "the pairwise Ja… | scripts/concordance_analysis.py (compute_concordance -> pairwise jaccard for all 3 tool p… → results/concordance/sciseq3/results.json (pairwise[].j… | ✅ |

## Headline speed/memory + SI Notes S1/S2

**Regenerate:**
```bash
python3 scripts/run_paper_benchmarks.py --threads 32 --replicates 3
python3 scripts/edit_tolerant_validity.py $SEQPROC_DATA_DIR/SRR6750041_R2.fastq --chem pe --out results/vtotal_pe/valid_ids.txt
python3 scripts/edit_tolerant_validity.py $SEQPROC_DATA_DIR/SRR13948564_full.fastq --chem lr --out results/vtotal_lr/valid_ids.txt
python3 scripts/validate_10x_sci_recall_against_vtotal.py --chemistry 10x --outdir results/vtotal_10x
python3 scripts/validate_10x_sci_recall_against_vtotal.py --chemistry sci --outdir results/vtotal_sci
time python3 scripts/edit_tolerant_validity.py $SEQPROC_DATA_DIR/SRR13948564_full.fastq --chem lr --out results/vtotal_lr/valid_ids.txt
time python3 scripts/edit_tolerant_validity.py $SEQPROC_DATA_DIR/SRR6750041_R2.fastq --chem pe --out results/vtotal_pe/valid_ids.txt
python3 scripts/run_paper_benchmarks.py --threads 32 --replicates 3   # for the 11.1 s LR runtime cell
```
| Datum | Value | Paper loc | Source script → output | Repro |
|---|---|---|---|---|
| "faster on most chemistries (up to $13.6\times$)" | `13.6x` | Conclusion, seqproc.tex:283; also Results… | scripts/run_paper_benchmarks.py (derived ratio of Table 3 runtime cells: matchbox PE 1290… → results/paper_figures/benchmark_results.json (canonica… | ✅ |
| "uses far less memory (up to $304\times$ less)" / "1/304 of the memor… | `304x` | Conclusion, seqproc.tex:283; Results pros… | scripts/run_paper_benchmarks.py (derived: matchbox 10x memory 42585.4MB / seqproc 10x 140… → results/paper_figures/benchmark_results.json (absent N… | ✅ |
| 10x Chromium v2 runtime: splitcode marginally faster than seqproc | `416.7 s vs 429.5 s` | Results prose, results.tex:57 | scripts/run_paper_benchmarks.py (Table 3 Runtime cells) → results/paper_figures/benchmark_results.json (absent);… | ✅ |
| "$1.2\times$ to $2.2\times$ faster than the next-fastest tool" (on PE… | `1.2x to 2.2x` | Results prose, results.tex:57 | scripts/run_paper_benchmarks.py (derived: sci splitcode 32.43/seqproc 26.76=1.21; LR matc… → results/paper_figures/benchmark_results.json (absent N… | ✅ |
| "requires 1/2 to 1/304 of the memory of the alternatives" | `1/2 to 1/304` | Results prose, results.tex:57 | scripts/run_paper_benchmarks.py (derived from Table 3 Memory: min ~2x on LR seqproc 249 v… → results/paper_figures/benchmark_results.json (absent N… | ✅ |
| 10x memory: matchbox 42.6 GB vs seqproc 140 MB | `42.6 GB vs 140 MB` | Results prose, results.tex:57 | scripts/run_paper_benchmarks.py (Table 3 Memory cells: matchbox 42,585 MB, seqproc 140 MB) → results/paper_figures/benchmark_results.json (absent);… | ✅ |
| V_total (SPLiT-seq PE), size and % of Total | `60,675,548 (78.17%)` | Note S1 table tab:v_total_splitseq (seqpr… | scripts/edit_tolerant_validity.py --chem pe (edlib HW infix, --max-linker-edit 6, Hamming… → NOT shipped; script prints JSON {total, valid, pct_of_… | ✅ |
| V_total (LR-SPLiT-seq), size and % of Total | `601,603 (10.44%)` | Note S1 table tab:v_total_splitseq (seqpr… | scripts/edit_tolerant_validity.py --chem lr (edlib HW infix, edit<=6, BOTH orientations, … → NOT shipped; stdout JSON {total, valid, pct_of_scanned… | ✅ |
| V_total (10x Chromium v2), size and % of Total | `234,382,218 (100.00%)` | Note S1 table tab:v_total_splitseq (seqpr… | scripts/validate_10x_sci_recall_against_vtotal.py --chemistry 10x (uses run_paper_benchma… → results/vtotal_10x/vtotal_recall_10x.json (key `v_tota… | ✅ |
| V_total (sci-RNA-seq3), size and % of Total | `19,767,975 (89.49%)` | Note S1 table tab:v_total_splitseq (seqpr… | scripts/validate_10x_sci_recall_against_vtotal.py --chemistry sci (uses run_paper_benchma… → results/vtotal_sci/vtotal_recall_sci.json (key `v_tota… | ✅ |
| Time to compute V_total_LR on the full LR dataset | `162 seconds` | Note S1, seqproc.tex:388; results.tex:54 | scripts/edit_tolerant_validity.py --chem lr (wall-clock runtime of the single-threaded Py… → NOT stored — wall-clock measurement (the script prints… | ✅ |
| Validator throughput (single-threaded Python) | `~35,600 reads/s` | Note S1, seqproc.tex:388; results.tex:54 | scripts/edit_tolerant_validity.py (derived: 5,764,421 reads / 162 s = 35,583) → NOT stored — derived from wall-clock time and input re… | ✅ |
| Time to compute V_total_PE on the full PE dataset | `about 36 minutes` | Note S1, seqproc.tex:388 | scripts/edit_tolerant_validity.py --chem pe (wall-clock on 77.6M R2 reads) → NOT stored — wall-clock measurement | ✅ |
| Validation script ~15x slower than seqproc on LR (seqproc's 11.1 s) | `~15x; 11.1 s` | results.tex:52,54 ("roughly 15x the runti… | 162 s (edit_tolerant_validity.py) / seqproc LR runtime 11.1 s (scripts/run_paper_benchmar… → 11.1 s from results/paper_figures/benchmark_results.js… | ✅ |
| Barcode-rank inflection rank per tool (seqproc/splitcode/matchbox) | `257, 255, 257 (UMI threshold 101-…` | Note S2, seqproc.tex:416; results.tex:92;… | biological_analysis/scripts/knee_barcoderanks.py (barcode_ranks -> infl_rank / infl_umi, … → biological_analysis/full_run_results/count_concordance… | ✅ |
| Microglia label confusion: dissenting relabels going to OPC | `22 of 25` | Note S2, seqproc.tex:393 | biological_analysis/scripts/jaccard_confusion.py (confusion_by_type: cells_any_tool + dis… → downstream_out/analysis/jaccard_confusion.json (confusion… | ✅ |
| OPC label confusion: dissenting instances that are microglia | `19 of 28` | Note S2, seqproc.tex:393 | biological_analysis/scripts/jaccard_confusion.py (confusion_by_type.OPC) → downstream_out/analysis/jaccard_confusion.json (confusion… | ✅ |
| Per-cell-type Jaccard table (tab:jaccard) values | `Neuron 0.959, Astrocyte 0.912, Ol…` | Note S2 table tab:jaccard (seqproc.tex:40… | biological_analysis/scripts/biological_analysis.py (celltype_jaccard_per_type + celltype_… → biological_analysis/full_run_results/biological_metric… | ✅ |
| tab:jaccard "Cells" column (per-type shared-cell counts) and shared-c… | `Neuron 142, Astrocyte 28, Oligode…` | Note S2 table tab:jaccard "Cells" column … | biological_analysis/scripts/jaccard_confusion.py (consensus_type_counts; each shared cell… → downstream_out/analysis/jaccard_confusion.json (consensus… | ✅ |
| Nonzero-barcode counts underlying fig:count_concordance panel A | `14,434 to 14,936 barcodes` | fig:count_concordance caption (seqproc.te… | biological_analysis/scripts/count_concordance.py (n_barcodes per tool, from STARsolo raw … → biological_analysis/full_run_results/count_concordance… | ✅ |
| Top-barcode UMI range | `18,847 to 19,644` | results.tex:93 ("near-identical top-barco… | biological_analysis/scripts/count_concordance.py (top_umi per tool) → biological_analysis/full_run_results/count_concordance… | ✅ |
| Per-barcode UMI Pearson (log1p), fig:count_concordance panel B | `0.978 and 0.992` | fig:count_concordance caption (seqproc.te… | biological_analysis/scripts/count_concordance.py (per_barcode_umi_pearson_logspace) → biological_analysis/full_run_results/count_concordance… | ✅ |
| Per-gene total UMI Pearson (log1p), fig:count_concordance panel C | `0.994 and 0.999` | fig:count_concordance caption (seqproc.te… | biological_analysis/scripts/count_concordance.py (per_gene_total_pearson_logspace) → biological_analysis/full_run_results/count_concordance… | ✅ |
| Synthetic-curve validation of the inflection detector ("recovers the … | `within 10% of the known cliff on …` | Note S2, seqproc.tex:416 ("On synthetic b… | biological_analysis/scripts/test_knee_point.py (test_barcoderanks_inflection_on_known_cli… → PASS/FAIL to stdout (no JSON) | ✅ |
