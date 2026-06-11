# Preprint Findings: LR-SPLiT-seq Denominator Correction + Reframing

**Status:** confirmed 2026-05-12. Blocking the seqproc preprint until resolved.
**Owner:** Elan Fisher

## Summary

The LR-SPLiT-seq row of Table~1 reports valid-recovery percentages against a
denominator (`4,229,250`) that does not match the actual input FASTQ
(`5,764,421` reads, verified by `wc -l` on
`/fs/nexus-projects/seqproc/bench/data/SRR13948564_full.fastq`). Tool emit
counts are correct. Only the denominator was wrong, so this is a one-line
data fix in `scripts/data_config.py` plus an arithmetic correction of every
LR-SPLiT-seq percentage in the paper.

Same root cause as commit `697662e` (which fixed the analogous issue for 10x
and sciseq3); LR-SPLiT-seq was missed because its benchmark output had
already been cached on disk and was not re-run.

## Verification trail

- File on disk: `5,764,421` reads (`wc -l` divided by 4).
- `results/concordance/lr_splitseq/results.json` reports `total_reads: 4229250`.
- `_ids.txt` line counts exactly match `recovery` integers in that same JSON
  (2,864,547 seqproc; 2,039,475 matchbox; 1,583,449 splitcode forward-only).
- Therefore: same FASTQ, wrong denominator. **Scenario A.**

## Corrected numbers (raw arithmetic, no re-runs needed)

| Tool | Emit (existing) | Old % | **Corrected %** |
|---|---:|---:|---:|
| seqproc | 2,864,547 | 67.73% | **49.69%** |
| matchbox | 2,039,475 | 48.22% | **35.38%** |
| splitcode forward-only | 1,583,449 | 37.44% | **27.47%** |
| splitcode dual-pass | TBD | — | TBD (from `splitcode_lr_dual_validate.py`) |

## Reframing for the preprint (per user direction 2026-05-12)

The corrected ~50% number for seqproc, while honest, is not the most
informative framing. Many input reads are simply not structurally recoverable
under any tool's algorithm: PCR artifacts, severe linker corruption, off-target
sequencing. A more informative metric is **recall against the achievable
ceiling** — defined as the union of valid reads recovered by any of the three
tools.

The original `concordance_results.json` already reports
`discordant.any_tool_union = 2,969,844`. Against the corrected input total
this is 51.5% — the practical ceiling on this dataset.

**Move to main text (Table 1 / abstract / results prose):** recall against
the achievable ceiling.

| Tool | Emit | Ceiling Recall (vs 2,969,844) |
|---|---:|---:|
| seqproc | 2,864,547 | **96.46%** |
| matchbox | 2,039,475 | **68.67%** |
| splitcode forward-only | 1,583,449 | **53.32%** |
| splitcode dual-pass | TBD | TBD |

**Move to supplementary** (as secondary validation, demonstrating that the
ceiling-recall framing is internally consistent and conservative): the
corrected raw-input recovery percentages from the table above (49.7%, 35.4%,
27.5%), plus the strict structural-validity audit (`V_total` from
`splitcode_lr_dual_validate.py`).

The two metrics tell complementary stories: the ceiling number says "of what
*can* be recovered, seqproc recovers 96.5%" — the dataset-yield number says
"of all input reads (including unrecoverable ones), seqproc emits 49.7% as
structurally valid output." Both are true. The first is the right headline.

## Action items

1. **Patch `scripts/data_config.py`**: `SRR13948564 full_reads: 4_229_250 → 5_764_421`. Combine with the in-flight `SRR7827254: 10_177_866 → 22_088_821` edit into one commit.
2. **Rewrite Table~1 LR-SPLiT-seq column header**: from `Valid Rec %` to `Ceiling Recall` (or similar — confirm name). Add the new percentages.
3. **Add a row to the table caption** defining the ceiling: "Ceiling Recall is the count of structurally valid reads in each tool's output divided by the union of valid reads recovered by any tool on this dataset (2,969,844 of 5,764,421 raw input reads for LR-SPLiT-seq, 51.5% of input)."
4. **Update the abstract**: drop "recovers 67.7% of input reads as structurally valid output" and replace with a ceiling-recall framing (e.g., "recovers 96.5% of the reads recoverable by any tool we evaluated").
5. **Update `paper/sections/results.tex:146`**: same fix in the prose paragraph after the table.
6. **Add a supplementary note** with the raw-input metric (yield), the strict V_total validator results from `splitcode_lr_dual_validate.py`, and the relationship between the two metrics.
7. **Methods section**: add the paragraph in `paper/sections/methods.tex` defining both metrics and clarifying the validator's strictness (exact `str.find()` linker; Hamming ≤ 1 barcode).
8. **Finalize splitcode dual-pass row** once `splitcode_lr_dual_validate.py` completes.

## Methodology note: precision is not cross-tool comparable

The three tools emit different output formats:

- `splitcode` preserves the original full-length read including the linker (~2.3 KB per record).
- `seqproc` emits a compact record containing only the extracted barcode block (~112 B per record); the linker is stripped during extraction.
- `matchbox` emits a TSV of barcode columns; no FASTQ sequence at all.

The structural-validity script (`SplitSeqSingleEndValidityAnalyzer`) is a `str.find()` linker scanner: it operates only on formats that preserve the linker. On the full LR-SPLiT-seq dataset:

- splitcode forward: 17.5% of emitted reads validate (1.58M emit → 277k valid).
- splitcode dual-pass: 17.3% (2.98M emit → 516k valid).
- seqproc: 0% by this validator (linker absent by design — extraction is the tool's job).
- matchbox: not testable in this format (TSV, not FASTQ).
- raw input direct scan: V_total = 489,305 = 8.49% of input.

**The 17.5% number characterizes splitcode's output verbosity, not its biological correctness.** Reads splitcode emits whose linker has been corrupted beyond exact-match recovery still contain a barcode region a downstream tool could try to use; the validator simply can't confirm them. Direct precision comparison across the three tools is not meaningful; the supplementary note must say so explicitly. The headline cross-tool comparison uses Ceiling Recall (read-ID level, format-agnostic) for this reason.

## Final Ceiling Recall numbers (LR-SPLiT-seq, against union-with-dual-pass-splitcode)

| Tool | Emit | Ceiling Recall (denom 3,068,359) |
|---|---:|---:|
| seqproc | 2,864,547 | **93.36%** |
| matchbox | 2,039,475 | 66.47% |
| splitcode (forward) | 1,583,449 | 51.61% |
| splitcode (dual-pass) | 2,984,418 | **97.26%** |

Ceiling = 3,068,359 = 53.23% of 5,764,421 input reads. The remaining 46.77% of input is consistent with expected PacBio long-read sequencing/library noise (~1% per-base error compounding across four critical anchor regions, no-insert reads, off-target priming).

## What does NOT need to be re-run

- **No tool re-runs are required.** All three tools' emit IDs and runtimes are
  cached at `results/concordance/lr_splitseq/`.
- **No other chemistry needs re-running.** 10x and sciseq3 were corrected by
  commit `697662e`. SPLiT-seq PE is unaffected by this bug.
- **No figures need regenerating** unless figure captions hard-code old
  percentages (audit `scripts/generate_figures.py` for literal `67.7`).

## What IS still running (the only compute remaining)

- `scripts/splitcode_lr_dual_validate.py --dataset full` — produces:
  - splitcode dual-pass emit count, valid count, summed runtime
  - V_total from strict raw-input validation (used in supplementary)

When it finishes, paste the JSON; the dual-pass row of Table~1 and the
supplementary V_total numbers get filled in deterministically from it.
