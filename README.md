# seqproc Paper Analysis

Reproducibility repository for the seqproc preprint. Every number, table cell,
and figure in the paper is produced by a script in this repository, applied to
the publicly available SRA datasets listed below. This branch (`paper-submission`)
contains only the scripts, configs, and tests needed to reproduce the paper.
Wider development history (exploratory scripts, abandoned chemistries,
intermediate findings docs) lives on the main development branch.

The paper benchmarks **seqproc** against **matchbox** and **splitcode** on four
single-cell RNA-seq chemistries, measuring runtime, memory, and a per-chemistry
recovery metric (see *Recovery metrics* below).

---

## Paper artifacts produced by this repository

| Paper artifact | Produced by (in this repo) | Key output files |
|---|---|---|
| **Table 3** (cross-chemistry benchmark — Total / Emitted / Validated / Recovery / Runtime / Memory) | `run_paper_benchmarks.py` + `concordance_analysis.py` + `splitcode_lr_dual_validate.py` + `validate_lr_recall_against_vtotal.py` + `validate_pe_recall_against_vtotal.py` | `results/paper_figures/benchmark_results.json`, `results/concordance/<chemistry>/results.json`, `results/lr_recall_vtotal/lr_recall_vtotal_results.json`, `results/pe_recall_vtotal/pe_recall_vtotal_results.json`, `results/splitcode_lr_dual/splitcode_lr_dual_results.json` |
| **Figure 3** (Recovery by Dataset and Tool — grouped bar chart) | `generate_figures.py` consuming the benchmark and V_total JSONs | `results/paper_figures/fig_recovery_comparison.{pdf,png}` |
| **Figure 4** (Hamming vs Edit Distance — relative recovery gain) | `concordance_analysis.py` (hamming_vs_edit subset) + `generate_figures.py` | `results/paper_figures/fig_hamming_vs_edit.{pdf,png}` |
| **Supplementary Figure S1** (Discordant read recovery breakdown) | `discordant_analysis.py` + `generate_figures.py` | `results/paper_figures/fig_discordant_summary.{pdf,png}`, `results/paper_figures/discordant_analysis.json` |
| **Supplementary Note S2** (V_total methodology) | `validate_lr_recall_against_vtotal.py` (LR) + `validate_pe_recall_against_vtotal.py` (PE) | JSON outputs listed in Table 3 row above |
| **Appendix B Table** (LR-SPLiT-seq optimization progression on 1M-read subset) | `lr_perf_rerun.py` | `results/lr_perf/lr_perf_rerun_results.json` |
| **`fig:count_concordance`** (downstream quantification concordance: knee, per-barcode, per-gene) | `biological_analysis/scripts/count_concordance.py` (driven by `run_phase2a.sh`) | `biological_analysis/full_run_results/count_concordance.{pdf,png,json}` |
| **`sec:downstream` metrics + `tab:pairwise_concordance`** (cell-type agreement, per-type + read-set Jaccard, cluster ARI, co-clustering) | `biological_analysis/scripts/biological_analysis.py` + `read_set_jaccard.py` | `biological_analysis/full_run_results/biological_metrics.json`, `read_set_jaccard.json` |
| **Merged `fig:downstream`** (knee, per-gene, scorecard, cell-type composition) | `biological_analysis/scripts/make_downstream_figure.py` | `biological_analysis/full_run_results/downstream_concordance{,_jaccard}.{pdf,png}` |
| **Supplementary `tab:jaccard`** (per-cell-type Jaccard) | `biological_analysis/scripts/biological_analysis.py` (emits `jaccard_supplement.md`) | `biological_analysis/full_run_results/jaccard_supplement.md` |
| **Tolerance sensitivity** (yield and concordance vs linker tolerance) | `biological_analysis/scripts/tolerance_sweep.sh` | `biological_analysis/full_run_results/tolerance_sensitivity.{pdf,png}`, `yield_by_tolerance.csv`, `concordance_by_tolerance.csv` |
| **Downstream pipeline driver** (STARsolo quantification + all downstream outputs, one command) | `biological_analysis/run_phase2a.sh` (see `biological_analysis/CLUSTER_RUN.md`) | `phase2a_out/analysis/` |

The corresponding LaTeX in the paper repo cites the figure and JSON paths
listed above; regenerating the JSONs and figures from this repo regenerates
those exact files.

---

## Recovery metrics — three layers, per chemistry

The paper reports a single column called *Recovery (%)* but the metric behind
it differs across chemistries because the tools' LR-SPLiT-seq output formats
diverge in ways that prevent a single validator from operating on all three.
The methodology in detail:

| Chemistry | Recovery metric | Numerator | Denominator |
|---|---|---|---|
| SPLiT-seq PE | **Recall vs V_total_PE** | `\|tool emit IDs ∩ V_total_PE\|` | `\|V_total_PE\|` (validity script on raw R2) |
| LR-SPLiT-seq | **Recall vs V_total_LR** | `\|tool emit IDs ∩ V_total_LR\|` | `\|V_total_LR\| = 489,305` (validity script on raw input) |
| 10x Chromium v2 | Valid Recovery (all 100% — trivial chemistry) | Validated reads in tool output | Total input reads |
| sci-RNA-seq3 | Valid Recovery | Validated reads in tool output | Total input reads |

`V_total` is the read-ID set the *independent* structural-validity script
accepts when run on the raw input FASTQ — i.e., the strict ground-truth set
of recoverable reads, defined entirely by the data and the validator's
criteria (exact `str.find` linker scan + Hamming ≤ 1 whitelist match).
A tool's Recovery is its recall against that set.

For more on why LR-SPLiT-seq requires the V_total framing (the three tools
emit incompatible output formats — full FASTQ, compact barcode-only FASTQ,
and TSV), see Supplementary Note S2 in the paper.

---

## Scripts: what each one does, when to run it, what it produces

All scripts live under `scripts/`. Run them from the repo root with the
Python environment from `requirements.txt`.

### Driver scripts (run these to reproduce the paper)

#### `run_paper_benchmarks.py` — the main benchmark runner

Runs all three tools on each of the four chemistries, three replicates each,
collecting wall-clock runtime, peak memory, emitted-read counts, and
structurally validated read counts. Produces the bulk of Table 3.

```bash
python3 scripts/run_paper_benchmarks.py --threads 32 --replicates 3
```

**Inputs:** tool binaries (set via `SEQPROC_BIN`, `MATCHBOX_BIN`,
`SPLITCODE_BIN`), raw FASTQs in `$SEQPROC_DATA_DIR`, configs in `configs/`.

**Output:** `results/paper_figures/benchmark_results.json` with one entry per
dataset, each containing `total_reads`, per-tool `reads_out`, `recovery_rate`,
`runtime`, `peak_mem_mb`, plus a sidecar `.benchmark_checkpoint.json` that
retains per-replicate `reads_valid` counts.

**Paper artifacts:** Table 3 (Total / Emitted / Validated / Runtime / Memory
columns for short-read chemistries).

#### `concordance_analysis.py` — pairwise tool concordance + Hamming-vs-edit

For each chemistry, runs all three tools and computes pairwise Jaccard indices
between their emitted-read ID sets. Also re-runs seqproc with Hamming-only
configs and compares against the edit-distance baseline.

```bash
python3 scripts/concordance_analysis.py --threads 32 --datasets all
# or restrict:
python3 scripts/concordance_analysis.py --datasets lr_splitseq splitseq_pe
```

**Output:** `results/concordance/<chemistry>/results.json` per chemistry plus
`results/concordance/<chemistry>/{seqproc_edit,matchbox,splitcode}_ids.txt`
holding the actual emitted-read ID sets. The latter are the inputs to the
V_total intersection scripts below.

**Paper artifacts:** the `hamming_vs_edit` subsection of each chemistry's JSON
feeds Figure 4 (Hamming vs Edit Distance). The ID-set caches are used by the
discordant-analysis and V_total scripts.

#### `discordant_analysis.py` — structural validation of tool-unique reads

Identifies reads recovered by only one of the three tools (the "discordant"
set), then runs an independent linker + whitelist validator on those records
to determine how many represent genuine signal vs. false-positive emit. On
SPLiT-seq PE this is where the 99.2% structurally-invalid figure for splitcode's
unique reads comes from.

```bash
python3 scripts/discordant_analysis.py
```

**Output:** `results/paper_figures/discordant_analysis.json` with per-tool,
per-chemistry breakdown of unique-emit validity rates.

**Paper artifacts:** Supplementary Figure S1, and the prose at
`results.tex:53` ("99.2% lacked valid linker sequences").

#### `splitcode_lr_dual_validate.py` — splitcode dual-pass on LR-SPLiT-seq + raw V_total_LR

Splitcode does not natively support orientation-aware matching, which on a
randomly-oriented chemistry like LR-SPLiT-seq leaves a fair-comparison gap.
This script runs splitcode twice (forward and reverse-complement), unions the
emitted-read sets, computes per-pass and combined validity, and also runs the
strict raw-input validator to produce V_total_LR (the denominator for the
LR-SPLiT-seq Recovery cells of Table 3).

```bash
python3 scripts/splitcode_lr_dual_validate.py \
    --dataset full --threads 8 \
    --outdir results/splitcode_lr_dual
```

**Output:** `results/splitcode_lr_dual/splitcode_lr_dual_results.json` with
the dual-pass aggregate counts, summed runtime, per-pass validity precision,
and `V_total_LR` (489,305 reads on this dataset). Also writes the deduplicated
combined FASTQ to the outdir for downstream intersection analysis.

**Paper artifacts:** the splitcode dual-pass row of Table 3 (LR-SPLiT-seq);
the V_total_LR denominator referenced throughout the LR row.

#### `validate_lr_recall_against_vtotal.py` — exact LR Validated cells

Takes the V_total_LR ID set from the strict raw-input validator and intersects
it with each tool's emitted-read ID set (loaded from
`results/concordance/lr_splitseq/*_ids.txt` and from the dual-pass FASTQ).
Outputs exact `|tool_emit ∩ V_total_LR|` counts per tool — these are the
Validated column values in the LR-SPLiT-seq rows of Table 3.

```bash
python3 scripts/validate_lr_recall_against_vtotal.py \
    --dual-pass-fq results/splitcode_lr_dual/splitcode_dual_combined_out.fq \
    --outdir results/lr_recall_vtotal
```

**Output:** `results/lr_recall_vtotal/lr_recall_vtotal_results.json` with the
four `|tool ∩ V_total|` counts and recall percentages (98.14%, 87.56%, 53.85%,
100.00% in the canonical run).

**Paper artifacts:** Table 3 LR-SPLiT-seq Validated cells (480,192 / 428,424 /
263,499 / 489,305).

#### `validate_pe_recall_against_vtotal.py` — exact PE Validated cells

The SPLiT-seq PE analogue of the LR script above. Computes `V_total_PE` by
running the paired-end validator on raw R1+R2, then intersects against each
tool's emit IDs in `results/concordance/splitseq_pe/`.

```bash
python3 scripts/validate_pe_recall_against_vtotal.py \
    --outdir results/pe_recall_vtotal
```

**Output:** `results/pe_recall_vtotal/pe_recall_vtotal_results.json`.

**Paper artifacts:** Table 3 SPLiT-seq PE Validated cells under the V_total
framing.

#### `lr_perf_rerun.py` — appendix LR optimization progression

Reproduces Appendix B's table showing the per-feature progression of seqproc
on LR-SPLiT-seq (forward-only Hamming → forward-only edit → orientation-aware
Hamming → orientation-aware edit). Used on a 1M-read subset for fast iteration.

```bash
python3 scripts/lr_perf_rerun.py --threads 4 --reps 3
```

**Output:** `results/lr_perf/lr_perf_rerun_results.json`.

**Paper artifacts:** Table 4 in `paper/sections/appendix_optimization.tex`
(LR-SPLiT-seq optimization progression).

#### `generate_figures.py` — render all paper figures

Reads the JSON outputs of the scripts above and produces the publication
figures.

```bash
python3 scripts/generate_figures.py
```

**Output:** under `results/paper_figures/`:
- `fig_recovery_comparison.{pdf,png}` (Figure 3)
- `fig_hamming_vs_edit.{pdf,png}` (Figure 4)
- `fig_discordant_summary.{pdf,png}` (Supplementary Figure S1)

### Supporting modules (do not run directly)

#### `data_config.py` — shared dataset / config / tool path resolution

Single source of truth for: SRA accessions, FASTQ file paths, per-chemistry
seqproc geometry files, matchbox configs, splitcode configs, and barcode
map paths. All driver scripts import from here. Run as `python3 scripts/data_config.py`
to print the resolved paths for the current environment.

### Shell wrappers

#### `setup_and_run.sh` — clean-machine bootstrap

End-to-end setup script: installs Python + micromamba env, clones and builds
tool binaries (seqproc, matchbox, splitcode), downloads the four SRA
datasets, sets the env vars, runs the full benchmark pipeline. Designed for
a fresh Linux machine with internet + gcc. Re-running is idempotent.

```bash
chmod +x scripts/setup_and_run.sh
./scripts/setup_and_run.sh
```

#### `dry_run_test.sh` — environment validation

Sanity-checks that every binary, config, FASTQ, and Python dependency
required by the driver scripts is present and resolvable. Exits non-zero if
anything is missing. Run this before any benchmark run on a new machine.

```bash
./scripts/dry_run_test.sh
```

#### `run_all.sh` — convenience wrapper

Runs the benchmark, concordance, discordant, and figure-generation scripts
in sequence with consistent threading and replicate count.

```bash
./scripts/run_all.sh --reads full --threads 32
```

---

## Reproducing the paper from scratch

If the environment is already configured (binaries on path, FASTQs in
`$SEQPROC_DATA_DIR`, Python env active), the full pipeline is:

```bash
# 1. Validate environment
./scripts/dry_run_test.sh

# 2. Main benchmark (Table 3 short-read rows + runtime/memory for all)
python3 scripts/run_paper_benchmarks.py --threads 32 --replicates 3

# 3. Pairwise concordance + Hamming-vs-edit data (Figure 4)
python3 scripts/concordance_analysis.py --threads 32 --datasets all

# 4. Discordant-read structural validation (Supp Fig S1)
python3 scripts/discordant_analysis.py

# 5. LR-SPLiT-seq splitcode dual-pass and raw V_total_LR
python3 scripts/splitcode_lr_dual_validate.py --dataset full --threads 8 \
    --outdir results/splitcode_lr_dual

# 6. LR-SPLiT-seq V_total intersection (Table 3 LR Validated cells)
python3 scripts/validate_lr_recall_against_vtotal.py \
    --dual-pass-fq results/splitcode_lr_dual/splitcode_dual_combined_out.fq \
    --outdir results/lr_recall_vtotal

# 7. SPLiT-seq PE V_total intersection (Table 3 PE Validated cells)
python3 scripts/validate_pe_recall_against_vtotal.py \
    --outdir results/pe_recall_vtotal

# 8. Optimization progression (Appendix B)
python3 scripts/lr_perf_rerun.py --threads 4 --reps 3

# 9. Render figures from collected JSONs
python3 scripts/generate_figures.py
```

For a clean machine without any of the above, use `./scripts/setup_and_run.sh`
which performs steps 1 through 9 plus all the environment provisioning.

---

## Datasets

| Chemistry | SRA accession | Total reads (full) | Mode |
|---|---|---:|---|
| SPLiT-seq PE | `SRR6750041` | 86,820,578 | paired-end |
| LR-SPLiT-seq | `SRR13948564` | 5,764,421 | single-end (PacBio HiFi) |
| 10x Chromium v2 | `SRR8315379` | 56,514,800 | paired-end |
| sci-RNA-seq3 | `SRR7827254` | 22,088,821 | paired-end |

Read counts are authoritative as of the canonical run; they are also
encoded in `scripts/data_config.py:SRA_INFO` and verified against
`count_fastq_reads()` at runtime.

FASTQ files are gitignored. Download with `fasterq-dump <accession>` into
`$SEQPROC_DATA_DIR` (default `data/`) before running anything.

---

## Configurations

`configs/` contains the geometry / config files actually used by the paper:

```
configs/
  seqproc/
    splitseq_filter_edit.geom         # SPLiT-seq PE, edit distance + whitelist
    splitseq_filter_hamming6.geom     # SPLiT-seq PE, Hamming baseline (Fig 4)
    splitseq_replacement_edit.geom    # SPLiT-seq PE, replacement-mode variant
    splitseq_singleend_edit_ann.geom  # LR-SPLiT-seq, orientation + edit (canonical)
    splitseq_singleend_ann.geom       # LR-SPLiT-seq, orientation + Hamming
    splitseq_singleend_edit.geom      # LR-SPLiT-seq, forward-only + edit
    splitseq_singleend.geom           # LR-SPLiT-seq, forward-only + Hamming
    splitseq_singleend_primer_edit.geom  # LR-SPLiT-seq, primer-based variant
    10x_v2.geom                       # 10x Chromium v2
    sciseq3.geom                      # sci-RNA-seq3, Hamming
    sciseq3_edit.geom                 # sci-RNA-seq3, edit distance
    splitseq_bc{1,2,3}_seq2seq.tsv    # SPLiT-seq barcode maps (seq → seq)
    splitseq_bc{1_whitelist_6bp,23_whitelist}.txt  # whitelists

  matchbox/
    splitseq_replacement.mb           # SPLiT-seq PE
    splitseq_singleend.mb             # LR-SPLiT-seq, forward-only
    splitseq_singleend_dual.mb        # LR-SPLiT-seq, dual-orientation
    10x_v2.mb                         # 10x Chromium v2
    sciseq3.mb                        # sci-RNA-seq3
    r2_r3.txt                         # matchbox support file for SPLiT-seq

  splitcode/
    splitseq_paper.config             # SPLiT-seq PE
    splitseq_singleend.config         # LR-SPLiT-seq
    10x_v2.config                     # 10x Chromium v2
    sciseq3.config                    # sci-RNA-seq3
```

The previous `10x_longread_*` configs (from an aborted long-read 10x
experiment that did not make the paper) have been removed in this branch.

---

## Tests

A regression test suite guards the benchmark pipeline against
silently-incorrect behavior:

```bash
python3 -m pytest tests/ -v
```

Coverage:
- Script syntax (every Python script compiles)
- Module imports (each script importable, with expected constants)
- Path resolution against `data_config.py`
- DATASETS dict consistency (every referenced file exists)
- Bug-fix regressions (subprocess exit code parsing, barcode-distance
  scoping, paired-end ID extraction)
- DNA reverse-complement primitive correctness
- Jaccard index correctness
- Recovery rates within expected bounds (when cluster artifacts are present)
- Banned-term guard (catches accidental development-iteration jargon
  leaking into the published code)
- Structural-validity analyzer correctness on synthetic and real reads
  (`test_se_validity.py`, `test_se_validity_real_data.py`)

A pre-commit hook in `.git/hooks/pre-commit` runs the suite automatically
before every commit; commits failing the suite are rejected.

---

## Repository layout

```
scripts/
  run_paper_benchmarks.py             # main benchmark runner (Table 3)
  concordance_analysis.py             # pairwise concordance + hamming-vs-edit (Fig 4)
  discordant_analysis.py              # unique-read structural validation (Supp Fig S1)
  generate_figures.py                 # render publication figures
  splitcode_lr_dual_validate.py       # splitcode dual-pass + V_total_LR (Table 3 LR rows)
  validate_lr_recall_against_vtotal.py  # exact LR Validated cells (Table 3)
  validate_pe_recall_against_vtotal.py  # exact PE Validated cells (Table 3)
  lr_perf_rerun.py                    # LR optimization progression (Appendix B)
  data_config.py                      # shared paths / datasets / configs module
  setup_and_run.sh                    # clean-machine bootstrap (binaries + data + run)
  dry_run_test.sh                     # environment validation
  run_all.sh                          # benchmark + concordance + discordant + figures

configs/
  seqproc/                            # EFGDL geometries (.geom) + whitelists
  matchbox/                           # matchbox pattern files (.mb)
  splitcode/                          # splitcode configs (.config)

tests/
  test_benchmark_pipeline.py          # bug-fix unit tests
  test_regression.py                  # whole-pipeline regression suite
  test_se_validity.py                 # SplitSeqSingleEndValidityAnalyzer unit tests
  test_se_validity_real_data.py       # validator behavior on actual FASTQ samples
  conftest.py                         # shared pytest fixtures

data/                                 # (gitignored) raw FASTQs from SRA
results/                              # (symlink) all JSON outputs and figures
README.md                             # this file
requirements.txt                      # Python dependencies (matplotlib, numpy, etc.)
pytest.ini                            # pytest configuration
.gitignore
```

---

## Environment

- Python 3.6+ (3.11 recommended; cluster runs use micromamba env `bench`)
- numpy, matplotlib (`pip install -r requirements.txt`)
- pytest for the test suite
- Tool binaries: seqproc, matchbox, splitcode — pinned versions installed by
  `setup_and_run.sh` from their respective repositories
- For SRA downloads: `sra-tools` (`fasterq-dump`), installed by `setup_and_run.sh`

Set the following before running scripts manually (`setup_and_run.sh` exports
them automatically):

```bash
export SEQPROC_PROJECT_ROOT=$(pwd)              # this repo
export SEQPROC_DATA_DIR=$WORKDIR/data           # FASTQ location
export SEQPROC_BIN=<path>/seqproc               # tool binaries
export MATCHBOX_BIN=<path>/matchbox
export SPLITCODE_BIN=<path>/splitcode
```

---

## Citation

If you use this analysis, please cite the seqproc paper.
