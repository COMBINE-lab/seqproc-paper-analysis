#!/usr/bin/env bash
# Phase 2A end-to-end: seqproc / splitcode / matchbox -> symmetric barcode reads -> STARsolo
# -> biological analysis (cell-calling, clustering, typing, joint-embedding concordance)
# -> runtime+RAM report -> bundle.
#
# Turnkey: validated on the sandbox box; on the cluster it is ONE command pointed at the larger
# dataset + the cluster's full-genome index.
#
#   run_phase2a.sh --r1 R1.fastq --r2 R2.fastq --genome STAR_INDEX --outdir OUT \
#                  [--threads 8] [--min-umi 200] [--python PYBIN]
#
# Env overrides (optional): SEQPROC_BIN, SPLITCODE_BIN, MATCHBOX_BIN.
set -euo pipefail

THREADS=8; MIN_UMI=200; PYBIN=""
R1=""; R2=""; GENOME=""; OUT=""
while [ $# -gt 0 ]; do case "$1" in
  --r1) R1=$2; shift 2;; --r2) R2=$2; shift 2;; --genome) GENOME=$2; shift 2;;
  --outdir) OUT=$2; shift 2;; --threads) THREADS=$2; shift 2;;
  --min-umi) MIN_UMI=$2; shift 2;; --python) PYBIN=$2; shift 2;;
  *) echo "unknown arg: $1" >&2; exit 1;; esac; done
[ -n "$R1$R2$GENOME$OUT" ] || { echo "usage: run_phase2a.sh --r1 --r2 --genome --outdir [--threads --min-umi --python]"; exit 1; }

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
SEQPROC=${SEQPROC_BIN:-$ROOT/../combine-lab/seqproc/target/release/seqproc}
SPLITCODE=${SPLITCODE_BIN:-$ROOT/../splitcode/build/src/splitcode}
MATCHBOX=${MATCHBOX_BIN:-$ROOT/../matchbox/target/release/matchbox}
WL=$HERE/configs/splitseq_bc_whitelist_96.txt
[ -n "$PYBIN" ] || PYBIN=$HERE/.venv_phase2a/bin/python
mkdir -p "$OUT"
# resolve inputs to absolute paths, then run from the repo root so the seqproc geom's
# repo-relative whitelist path (configs/seqproc/...) resolves no matter where this is invoked.
OUT=$(cd "$OUT" && pwd); R1=$(readlink -f "$R1"); R2=$(readlink -f "$R2"); GENOME=$(readlink -f "$GENOME")
cd "$ROOT"
echo "[config] threads=$THREADS min_umi=$MIN_UMI python=$PYBIN"

# ---- timing + peak-RAM instrumentation ----
RES="$OUT/resources.csv"; echo "step,tool,seconds,peak_ram_mb" > "$RES"
now()     { date +%s.%N; }
elapsed() { awk -v a="$1" -v b="$2" 'BEGIN{printf "%.2f", b-a}'; }
ram_mb()  { awk '/Maximum resident set size/{printf "%.1f", $NF/1024}' "$1"; }
record()  { echo "$1,$2,$3,$4" >> "$RES"; }   # step,tool,seconds,ram_mb

echo "[1/7] seqproc (observed 3-round geom)"
TF=$(mktemp); t0=$(now)
/usr/bin/time -v "$SEQPROC" -g "$HERE/configs/splitseq_quant_observed.geom" -1 "$R1" -2 "$R2" \
  -o "$OUT/sp_cdna.fq" -w "$OUT/sp_bc.fq" -t "$THREADS" 2>"$TF"
record barcode seqproc "$(elapsed "$t0" "$(now)")" "$(ram_mb "$TF")"; rm -f "$TF"

echo "[2/7] splitcode (tag-relative extraction, stitched)"
TF=$(mktemp); t0=$(now)
( /usr/bin/time -v "$SPLITCODE" -c "$HERE/configs/splitseq_extract.config" -N 2 -t "$THREADS" --x-only -p \
    -x '1:0<u[10]>,<b3[8]>{linker1},{linker1}<b2[8]>,{linker2}<b1[8]>' "$R1" "$R2" 2>"$TF" ) \
  | python3 "$HERE/scripts/splitcode_quant_extract.py" --r1 "$R1" --out-cdna "$OUT/sc_cdna.fq" --out-bc "$OUT/sc_bc.fq"
record barcode splitcode "$(elapsed "$t0" "$(now)")" "$(ram_mb "$TF")"; rm -f "$TF"

echo "[3/7] matchbox (positional extraction, stitched)"
TF=$(mktemp); t0=$(now)
( /usr/bin/time -v "$MATCHBOX" -s "$HERE/configs/splitseq_matchbox.mb" -t "$THREADS" "$R2" 2>"$TF" ) \
  | python3 "$HERE/scripts/matchbox_quant_extract.py" --r1 "$R1" --out-cdna "$OUT/mb_cdna.fq" --out-bc "$OUT/mb_bc.fq"
record barcode matchbox "$(elapsed "$t0" "$(now)")" "$(ram_mb "$TF")"; rm -f "$TF"

echo "[4/7] STARsolo (identical CB_UMI_Complex config for all three tools)"
declare -A TOOLNAME=( [sp]=seqproc [sc]=splitcode [mb]=matchbox )
for tool in sp sc mb; do
  TF=$(mktemp); t0=$(now)
  /usr/bin/time -v STAR --runThreadN "$THREADS" --genomeDir "$GENOME" \
    --soloType CB_UMI_Complex \
    --readFilesIn "$OUT/${tool}_cdna.fq" "$OUT/${tool}_bc.fq" \
    --soloCBwhitelist "$WL" "$WL" "$WL" \
    --soloCBposition 0_10_0_17 0_18_0_25 0_26_0_33 \
    --soloUMIposition 0_0_0_9 --soloCBmatchWLtype 1MM \
    --soloFeatures Gene --soloCellFilter None \
    --outSAMtype None --outFileNamePrefix "$OUT/${tool}_" 2>"$TF"
  record starsolo "${TOOLNAME[$tool]}" "$(elapsed "$t0" "$(now)")" "$(ram_mb "$TF")"; rm -f "$TF"
done

echo "[5/7] biological analysis + count concordance"
TF=$(mktemp); t0=$(now)
/usr/bin/time -v "$PYBIN" "$HERE/scripts/biological_analysis.py" "$OUT/analysis" "$MIN_UMI" \
  "seqproc:$OUT/sp_Solo.out/Gene" "splitcode:$OUT/sc_Solo.out/Gene" "matchbox:$OUT/mb_Solo.out/Gene" 2>"$TF"
record analysis biological "$(elapsed "$t0" "$(now)")" "$(ram_mb "$TF")"; rm -f "$TF"
TF=$(mktemp); t0=$(now)
/usr/bin/time -v "$PYBIN" "$HERE/scripts/count_concordance.py" "$OUT/analysis" \
  "seqproc:$OUT/sp_Solo.out/Gene" "splitcode:$OUT/sc_Solo.out/Gene" "matchbox:$OUT/mb_Solo.out/Gene" 2>"$TF"
record analysis count "$(elapsed "$t0" "$(now)")" "$(ram_mb "$TF")"; rm -f "$TF"

echo "[6/7] resource report (figure + table)"
"$PYBIN" "$HERE/scripts/resource_report.py" "$RES" "$OUT/analysis"

echo "[7/7] bundle results for transfer back"
BUNDLE="$OUT/phase2a_bundle.tar.gz"
tar -czf "$BUNDLE" -C "$OUT" \
  analysis resources.csv \
  sp_Solo.out/Gene/Summary.csv sc_Solo.out/Gene/Summary.csv mb_Solo.out/Gene/Summary.csv \
  sp_Solo.out/Barcodes.stats sc_Solo.out/Barcodes.stats mb_Solo.out/Barcodes.stats 2>/dev/null || true
echo "DONE. Figures+metrics+resources: $OUT/analysis/   Transferable bundle: $BUNDLE"
