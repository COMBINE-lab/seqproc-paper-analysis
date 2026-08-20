#!/usr/bin/env bash
# Quantify already-materialized final SPLiT-seq PE products with one pinned
# STARsolo configuration, then regenerate every downstream metric and figure.
set -euo pipefail

THREADS=32
MIN_UMI=200
CB_MATCH=1MM
PYBIN=""
GENOME=""
OUT=""
SEQPROC_R1=""; SEQPROC_R2=""
SPLITCODE_R1=""; SPLITCODE_R2=""
MATCHBOX_R1=""; MATCHBOX_R2=""
SEQPROC_BITMAP=""; SPLITCODE_BITMAP=""; MATCHBOX_BITMAP=""

usage() {
  cat <<'EOF'
usage: run_current_downstream.sh --genome STAR_INDEX --outdir OUT \
  --seqproc-r1 FILE --seqproc-r2 FILE \
  --splitcode-r1 FILE --splitcode-r2 FILE \
  --matchbox-r1 FILE --matchbox-r2 FILE \
  [--seqproc-bitmap FILE --splitcode-bitmap FILE --matchbox-bitmap FILE] \
  [--threads 32] [--min-umi 200] [--cb-match 1MM] [--python PYTHON]
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --genome) GENOME=$2; shift 2;;
    --outdir) OUT=$2; shift 2;;
    --seqproc-r1) SEQPROC_R1=$2; shift 2;;
    --seqproc-r2) SEQPROC_R2=$2; shift 2;;
    --splitcode-r1) SPLITCODE_R1=$2; shift 2;;
    --splitcode-r2) SPLITCODE_R2=$2; shift 2;;
    --matchbox-r1) MATCHBOX_R1=$2; shift 2;;
    --matchbox-r2) MATCHBOX_R2=$2; shift 2;;
    --seqproc-bitmap) SEQPROC_BITMAP=$2; shift 2;;
    --splitcode-bitmap) SPLITCODE_BITMAP=$2; shift 2;;
    --matchbox-bitmap) MATCHBOX_BITMAP=$2; shift 2;;
    --threads) THREADS=$2; shift 2;;
    --min-umi) MIN_UMI=$2; shift 2;;
    --cb-match) CB_MATCH=$2; shift 2;;
    --python) PYBIN=$2; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2;;
  esac
done

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
STARBIN=${STAR_BIN:-$HERE/tools/STAR-2.7.11b/bin/STAR}
WL=$ROOT/configs/seqproc/splitseq_bc8_whitelist.txt
[ -n "$PYBIN" ] || PYBIN=$HERE/.venv_downstream/bin/python

for var in GENOME OUT SEQPROC_R1 SEQPROC_R2 SPLITCODE_R1 SPLITCODE_R2 MATCHBOX_R1 MATCHBOX_R2; do
  [ -n "${!var}" ] || { echo "missing required argument: $var" >&2; usage >&2; exit 2; }
done
for path in "$GENOME" "$SEQPROC_R1" "$SEQPROC_R2" "$SPLITCODE_R1" "$SPLITCODE_R2" "$MATCHBOX_R1" "$MATCHBOX_R2" "$STARBIN" "$PYBIN" "$WL"; do
  [ -e "$path" ] || { echo "required path does not exist: $path" >&2; exit 2; }
done

mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
GENOME=$(readlink -f "$GENOME")
SEQPROC_R1=$(readlink -f "$SEQPROC_R1"); SEQPROC_R2=$(readlink -f "$SEQPROC_R2")
SPLITCODE_R1=$(readlink -f "$SPLITCODE_R1"); SPLITCODE_R2=$(readlink -f "$SPLITCODE_R2")
MATCHBOX_R1=$(readlink -f "$MATCHBOX_R1"); MATCHBOX_R2=$(readlink -f "$MATCHBOX_R2")

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONHASHSEED=0
export MPLCONFIGDIR="$OUT/.matplotlib"
mkdir -p "$MPLCONFIGDIR"

RES="$OUT/resources.csv"
echo "step,tool,seconds,peak_ram_mb" > "$RES"
now() { date +%s.%N; }
elapsed() { awk -v a="$1" -v b="$2" 'BEGIN{printf "%.2f", b-a}'; }
ram_mb() { awk '/Maximum resident set size/{printf "%.1f", $NF/1024}' "$1"; }
record() { echo "$1,$2,$3,$4" >> "$RES"; }

declare -A R1=( [seqproc]="$SEQPROC_R1" [splitcode]="$SPLITCODE_R1" [matchbox]="$MATCHBOX_R1" )
declare -A R2=( [seqproc]="$SEQPROC_R2" [splitcode]="$SPLITCODE_R2" [matchbox]="$MATCHBOX_R2" )
declare -A PREFIX=( [seqproc]=sp [splitcode]=sc [matchbox]=mb )

echo "[1/5] validate and checksum the three final FASTQ pairs"
for tool in seqproc splitcode matchbox; do
  "$PYBIN" "$HERE/scripts/validate_fastq_pairs.py" \
    --name "$tool" --r1 "${R1[$tool]}" --r2 "${R2[$tool]}" \
    --expected-r1-length 66 --expected-r2-length 34 \
    --output "$OUT/${tool}_fastq_provenance.json" > "$OUT/${tool}_fastq_validation.log"
done

echo "[2/5] STARsolo with one shared CB_UMI_Complex configuration"
"$STARBIN" --version > "$OUT/STAR.version.txt"
for tool in seqproc splitcode matchbox; do
  prefix=${PREFIX[$tool]}
  tf=$(mktemp); t0=$(now)
  /usr/bin/time -v "$STARBIN" --runThreadN "$THREADS" --genomeDir "$GENOME" \
    --soloType CB_UMI_Complex \
    --readFilesIn "${R1[$tool]}" "${R2[$tool]}" \
    --soloCBwhitelist "$WL" "$WL" "$WL" \
    --soloCBposition 0_10_0_17 0_18_0_25 0_26_0_33 \
    --soloUMIposition 0_0_0_9 \
    --soloCBmatchWLtype "$CB_MATCH" \
    --soloFeatures Gene --soloCellFilter None \
    --outSAMtype None --outFileNamePrefix "$OUT/${prefix}_" 2>"$tf"
  cp "$tf" "$OUT/${prefix}_STAR.time.txt"
  record starsolo "$tool" "$(elapsed "$t0" "$(now)")" "$(ram_mb "$tf")"
  rm -f "$tf"
done

ARGS=(
  "seqproc:$OUT/sp_Solo.out/Gene"
  "splitcode:$OUT/sc_Solo.out/Gene"
  "matchbox:$OUT/mb_Solo.out/Gene"
)

echo "[3/5] biological analysis and pairwise count concordance"
tf=$(mktemp); t0=$(now)
/usr/bin/time -v "$PYBIN" "$HERE/scripts/biological_analysis.py" \
  "$OUT/analysis" "$MIN_UMI" "${ARGS[@]}" 2>"$tf"
record analysis biological "$(elapsed "$t0" "$(now)")" "$(ram_mb "$tf")"
cp "$tf" "$OUT/biological_analysis.time.txt"; rm -f "$tf"

tf=$(mktemp); t0=$(now)
/usr/bin/time -v "$PYBIN" "$HERE/scripts/count_concordance.py" \
  "$OUT/analysis" "${ARGS[@]}" 2>"$tf"
record analysis count "$(elapsed "$t0" "$(now)")" "$(ram_mb "$tf")"
cp "$tf" "$OUT/count_concordance.time.txt"; rm -f "$tf"

echo "[4/5] confusion supplement, read-set concordance, and run provenance"
"$PYBIN" "$HERE/scripts/jaccard_confusion.py" \
  "$OUT/analysis" "$MIN_UMI" "${ARGS[@]}" > "$OUT/analysis/jaccard_confusion.log"

if [ -n "$SEQPROC_BITMAP$SPLITCODE_BITMAP$MATCHBOX_BITMAP" ]; then
  [ -n "$SEQPROC_BITMAP" ] && [ -n "$SPLITCODE_BITMAP" ] && [ -n "$MATCHBOX_BITMAP" ] || {
    echo "supply either all three bitmap arguments or none" >&2
    exit 2
  }
  "$PYBIN" "$HERE/scripts/read_set_jaccard_bitmaps.py" \
    --output "$OUT/analysis/read_set_jaccard.json" \
    "seqproc:$SEQPROC_BITMAP" "splitcode:$SPLITCODE_BITMAP" "matchbox:$MATCHBOX_BITMAP"
fi

"$PYBIN" "$HERE/scripts/downstream_provenance.py" \
  --output "$OUT/downstream_provenance.json" --results-root "$OUT" \
  --genome "$GENOME" --star "$STARBIN" --python "$PYBIN" --whitelist "$WL" \
  --threads "$THREADS" --min-umi "$MIN_UMI" \
  --cb-match-whitelist-type "$CB_MATCH" \
  --fastq-provenance "$OUT/seqproc_fastq_provenance.json" \
  --fastq-provenance "$OUT/splitcode_fastq_provenance.json" \
  --fastq-provenance "$OUT/matchbox_fastq_provenance.json" \
  --workflow-file "$HERE/run_current_downstream.sh" \
  --workflow-file "$HERE/requirements.lock" \
  --workflow-file "$HERE/scripts/biological_analysis.py" \
  --workflow-file "$HERE/scripts/count_concordance.py" \
  --workflow-file "$HERE/scripts/jaccard_confusion.py" \
  --workflow-file "$HERE/scripts/read_set_jaccard_bitmaps.py"

"$PYBIN" "$HERE/scripts/summarize_downstream_run.py" "$OUT" \
  --json "$OUT/analysis/downstream_summary.json" \
  --markdown "$OUT/analysis/downstream_summary.md"

# The summary is itself a result, so refresh output hashes after creating it.
"$PYBIN" "$HERE/scripts/downstream_provenance.py" \
  --output "$OUT/downstream_provenance.json" --results-root "$OUT" \
  --genome "$GENOME" --star "$STARBIN" --python "$PYBIN" --whitelist "$WL" \
  --threads "$THREADS" --min-umi "$MIN_UMI" \
  --cb-match-whitelist-type "$CB_MATCH" \
  --fastq-provenance "$OUT/seqproc_fastq_provenance.json" \
  --fastq-provenance "$OUT/splitcode_fastq_provenance.json" \
  --fastq-provenance "$OUT/matchbox_fastq_provenance.json" \
  --workflow-file "$HERE/run_current_downstream.sh" \
  --workflow-file "$HERE/requirements.lock" \
  --workflow-file "$HERE/scripts/biological_analysis.py" \
  --workflow-file "$HERE/scripts/count_concordance.py" \
  --workflow-file "$HERE/scripts/jaccard_confusion.py" \
  --workflow-file "$HERE/scripts/read_set_jaccard_bitmaps.py" \
  --workflow-file "$HERE/scripts/summarize_downstream_run.py"

echo "[5/5] bundle compact outputs"
tar -czf "$OUT/downstream_bundle.tar.gz" -C "$OUT" \
  analysis resources.csv downstream_provenance.json STAR.version.txt \
  seqproc_fastq_provenance.json splitcode_fastq_provenance.json matchbox_fastq_provenance.json \
  sp_Solo.out/Gene/Summary.csv sc_Solo.out/Gene/Summary.csv mb_Solo.out/Gene/Summary.csv \
  sp_Solo.out/Barcodes.stats sc_Solo.out/Barcodes.stats mb_Solo.out/Barcodes.stats
echo "DONE: $OUT"
