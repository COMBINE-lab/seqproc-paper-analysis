# Sample-barcode (`s`) feature validation

This directory holds the downstream validation that verifies the new
`s` (SampleBarcode) interval specifier added to EFGDL.

The specifier was motivated by multiplexed single-cell protocols — in
particular 10x Genomics Chromium Flex, where R1 carries the GEM cell
barcode and R2 carries a separate probe barcode that encodes the
sample of origin. Mechanically `s` behaves identically to `b`;
semantically it makes cell identity and sample identity distinguishable
in the output so downstream demultiplexers can route each barcode to
its correct use.

## What is validated

1. **Byte equivalence with `b` on real data.** Running seqproc with the
   same geometry but substituting `s[16]` for `b[16]` in the cell slot
   on a 100k-read subset of SRR8315379 (10x Chromium v2) produces
   byte-identical R1/R2 outputs under `--preserve-order`. This confirms
   that the new token is a pure labeling alias with no behavioural
   drift. (Procedure only; not automated here — see the `how to run`
   section below.)

2. **Correct extraction on synthetic Flex-like data.** `synth_flex.py`
   emits a 50,000-read paired-end dataset with 500 cells drawn from a
   nominal 40/30/20/10 sample prior and 4 fixed 8 bp probe barcodes.
   `analyze_flex.py` runs three assertions:
   - exactly the 4 known probe barcodes are recovered;
   - observed read-level proportions match the *realized* cell
     distribution within 1 percentage point;
   - every cell barcode appears with exactly one sample (no leakage),
     and no read disagrees with the ground-truth cell -> sample map.

   Tolerances are stated at 1 pp because the dominant source of
   variance is which cells are drawn from the 40/30/20/10 prior, not
   the read-level sampling itself. The script compares observed
   proportions to the per-sample realized cell fractions (computed
   from the truth file), which is the correct read-level expectation.

## How to run

From this directory:

```bash
# 1. generate synthetic FASTQ + truth
mkdir -p work
python3 synth_flex.py work/flex_R1.fastq work/flex_R2.fastq work/flex_truth.tsv

# 2. run seqproc with the Flex geometry
printf '1{b[16]u[12]}2{r[50]x[28]s[8]}' > work/geom_flex.efgdl
seqproc \
  -g work/geom_flex.efgdl \
  -1 work/flex_R1.fastq -2 work/flex_R2.fastq \
  -o work/out_flex_R1.fastq -w work/out_flex_R2.fastq \
  --preserve-order -s work/summary_flex.json

# 3. validate
python3 analyze_flex.py --work-dir work
```

Expected terminal output ends with
`[PASS] synthetic Flex recovery is clean: 4/4 probe barcodes,
proportions within 1%, no cell-sample leak.`

## Files

- `synth_flex.py` — deterministic (seed 42) generator for the Flex-like
  dataset. Arguments: `R1.fastq R2.fastq truth.tsv`.
- `analyze_flex.py` — three-check validator. `--work-dir` or
  `FLEX_WORK_DIR` selects the directory containing the truth file and
  the seqproc output.

## Reference result

On seed 42 with seqproc 0.1.0 and the geometry above:

| sample   | barcode   | nominal | realized | observed | \|obs−real\| |
|----------|-----------|--------:|---------:|---------:|-------------:|
| SAMPLE01 | ACGTACGT  |  0.4000 |   0.3600 |   0.3641 |       0.0041 |
| SAMPLE02 | TGCATGCA  |  0.3000 |   0.3080 |   0.3068 |       0.0012 |
| SAMPLE03 | GATTACAG  |  0.2000 |   0.2140 |   0.2092 |       0.0048 |
| SAMPLE04 | CCGGCCGG  |  0.1000 |   0.1180 |   0.1199 |       0.0019 |

All 500 cells appear with exactly one sample; 0 / 50,000 reads
disagree with the truth map.
