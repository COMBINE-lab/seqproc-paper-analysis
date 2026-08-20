#!/usr/bin/env bash
# Quantify a pre-materialized Matchbox radius-1-whitelist sensitivity result
# alongside the primary seqproc and splitcode STARsolo outputs.
set -euo pipefail

THREADS=32
MIN_UMI=200
CB_MATCH=1MM
INPUT_RECORDS=77621181
PYBIN=""
GENOME=""
PRIMARY=""
OUT=""
MATCHBOX_R1=""; MATCHBOX_R2=""
REFERENCE_BITMAP=""
SEQPROC_BITMAP=""; SPLITCODE_BITMAP=""; MATCHBOX_EXACT_BITMAP=""

usage() {
  cat <<'EOF'
usage: run_matchbox_ham1_sensitivity.sh --genome STAR_INDEX \
  --primary-run PRIMARY_DOWNSTREAM_RUN --outdir OUT \
  --matchbox-r1 FILE --matchbox-r2 FILE \
  --reference-bitmap FILE --seqproc-bitmap FILE --splitcode-bitmap FILE \
  --matchbox-exact-bitmap FILE [--input-records 77621181] \
  [--threads 32] [--min-umi 200] [--cb-match 1MM] [--python PYTHON]

This script deliberately reuses the primary seqproc and splitcode STARsolo
matrices. Only the Matchbox expanded-whitelist read set is remapped. The
expanded configuration remains a sensitivity analysis, not a replacement for
the benchmark configuration.
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --genome) GENOME=$2; shift 2;;
    --primary-run) PRIMARY=$2; shift 2;;
    --outdir) OUT=$2; shift 2;;
    --matchbox-r1) MATCHBOX_R1=$2; shift 2;;
    --matchbox-r2) MATCHBOX_R2=$2; shift 2;;
    --reference-bitmap) REFERENCE_BITMAP=$2; shift 2;;
    --seqproc-bitmap) SEQPROC_BITMAP=$2; shift 2;;
    --splitcode-bitmap) SPLITCODE_BITMAP=$2; shift 2;;
    --matchbox-exact-bitmap) MATCHBOX_EXACT_BITMAP=$2; shift 2;;
    --input-records) INPUT_RECORDS=$2; shift 2;;
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
AUDIT=$ROOT/tools/bin/fastq-numeric-audit
WL=$ROOT/configs/seqproc/splitseq_bc8_whitelist.txt
[ -n "$PYBIN" ] || PYBIN=$HERE/.venv_downstream/bin/python

for var in GENOME PRIMARY OUT MATCHBOX_R1 MATCHBOX_R2 REFERENCE_BITMAP \
  SEQPROC_BITMAP SPLITCODE_BITMAP MATCHBOX_EXACT_BITMAP; do
  [ -n "${!var}" ] || { echo "missing required argument: $var" >&2; usage >&2; exit 2; }
done
for path in "$GENOME" "$PRIMARY" "$MATCHBOX_R1" "$MATCHBOX_R2" \
  "$REFERENCE_BITMAP" "$SEQPROC_BITMAP" "$SPLITCODE_BITMAP" \
  "$MATCHBOX_EXACT_BITMAP" "$STARBIN" "$AUDIT" "$PYBIN" "$WL"; do
  [ -e "$path" ] || { echo "required path does not exist: $path" >&2; exit 2; }
done

primary_cb_match=$("$PYBIN" -c \
  'import json,sys; print(json.load(open(sys.argv[1]))["parameters"]["cb_match_whitelist_type"])' \
  "$PRIMARY/downstream_provenance.json")
[ "$primary_cb_match" = "$CB_MATCH" ] || {
  echo "primary run used $primary_cb_match but sensitivity requested $CB_MATCH" >&2
  echo "rerun the primary matrices with the same --cb-match setting first" >&2
  exit 2
}

mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
PRIMARY=$(readlink -f "$PRIMARY")
GENOME=$(readlink -f "$GENOME")
MATCHBOX_R1=$(readlink -f "$MATCHBOX_R1"); MATCHBOX_R2=$(readlink -f "$MATCHBOX_R2")

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONHASHSEED=0 MPLCONFIGDIR="$OUT/.matplotlib"
mkdir -p "$MPLCONFIGDIR" "$OUT/analysis"

echo "[1/6] validate expanded-whitelist Matchbox FASTQs"
"$PYBIN" "$HERE/scripts/validate_fastq_pairs.py" \
  --name matchbox --r1 "$MATCHBOX_R1" --r2 "$MATCHBOX_R2" \
  --expected-r1-length 66 --expected-r2-length 34 \
  --output "$OUT/matchbox_fastq_provenance.json" > "$OUT/matchbox_fastq_validation.log"
cp "$PRIMARY/seqproc_fastq_provenance.json" "$OUT/seqproc_fastq_provenance.json"
cp "$PRIMARY/splitcode_fastq_provenance.json" "$OUT/splitcode_fastq_provenance.json"

echo "[2/6] audit accessions and evaluate the conservative structural reference"
"$AUDIT" "$MATCHBOX_R1" 1 "$INPUT_RECORDS" "$OUT/matchbox_expanded.bitmap" \
  > "$OUT/matchbox_numeric_audit.json"
"$PYBIN" "$ROOT/scripts/evaluate_numeric_bitmaps.py" \
  --reference "$REFERENCE_BITMAP" --input-records "$INPUT_RECORDS" \
  --tool "seqproc=$SEQPROC_BITMAP" --tool "splitcode=$SPLITCODE_BITMAP" \
  --tool "matchbox_exact=$MATCHBOX_EXACT_BITMAP" \
  --tool "matchbox_expanded=$OUT/matchbox_expanded.bitmap" \
  --out "$OUT/accuracy.json" > "$OUT/accuracy.log"

echo "[3/6] STARsolo for expanded-whitelist Matchbox"
ln -sfn "$PRIMARY/sp_Solo.out" "$OUT/sp_Solo.out"
ln -sfn "$PRIMARY/sc_Solo.out" "$OUT/sc_Solo.out"
"$STARBIN" --version > "$OUT/STAR.version.txt"
tf=$(mktemp); t0=$(date +%s.%N)
/usr/bin/time -v "$STARBIN" --runThreadN "$THREADS" --genomeDir "$GENOME" \
  --soloType CB_UMI_Complex --readFilesIn "$MATCHBOX_R1" "$MATCHBOX_R2" \
  --soloCBwhitelist "$WL" "$WL" "$WL" \
  --soloCBposition 0_10_0_17 0_18_0_25 0_26_0_33 \
  --soloUMIposition 0_0_0_9 --soloCBmatchWLtype "$CB_MATCH" \
  --soloFeatures Gene --soloCellFilter None --outSAMtype None \
  --outFileNamePrefix "$OUT/mb_" 2>"$tf"
t1=$(date +%s.%N)
cp "$tf" "$OUT/mb_STAR.time.txt"
seconds=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}')
rss=$(awk '/Maximum resident set size/{printf "%.1f", $NF/1024}' "$tf")
rm -f "$tf"
awk -F, 'NR==1 || ($1=="starsolo" && ($2=="seqproc" || $2=="splitcode"))' \
  "$PRIMARY/resources.csv" > "$OUT/resources.csv"
echo "starsolo,matchbox,$seconds,$rss" >> "$OUT/resources.csv"

ARGS=(
  "seqproc:$OUT/sp_Solo.out/Gene"
  "splitcode:$OUT/sc_Solo.out/Gene"
  "matchbox:$OUT/mb_Solo.out/Gene"
)

echo "[4/6] regenerate biological and count-concordance analyses"
"$PYBIN" "$HERE/scripts/biological_analysis.py" "$OUT/analysis" "$MIN_UMI" "${ARGS[@]}"
"$PYBIN" "$HERE/scripts/count_concordance.py" "$OUT/analysis" "${ARGS[@]}"
"$PYBIN" "$HERE/scripts/jaccard_confusion.py" "$OUT/analysis" "$MIN_UMI" "${ARGS[@]}" \
  > "$OUT/analysis/jaccard_confusion.log"
"$PYBIN" "$HERE/scripts/read_set_jaccard_bitmaps.py" \
  --output "$OUT/analysis/read_set_jaccard.json" \
  "seqproc:$SEQPROC_BITMAP" "splitcode:$SPLITCODE_BITMAP" \
  "matchbox:$OUT/matchbox_expanded.bitmap"

echo "[5/6] provenance and summary"
"$PYBIN" "$HERE/scripts/downstream_provenance.py" \
  --output "$OUT/downstream_provenance.json" --results-root "$OUT" \
  --genome "$GENOME" --star "$STARBIN" --python "$PYBIN" --whitelist "$WL" \
  --threads "$THREADS" --min-umi "$MIN_UMI" \
  --cb-match-whitelist-type "$CB_MATCH" \
  --fastq-provenance "$OUT/seqproc_fastq_provenance.json" \
  --fastq-provenance "$OUT/splitcode_fastq_provenance.json" \
  --fastq-provenance "$OUT/matchbox_fastq_provenance.json" \
  --workflow-file "$HERE/run_matchbox_ham1_sensitivity.sh" \
  --workflow-file "$HERE/requirements.lock" \
  --workflow-file "$HERE/scripts/biological_analysis.py" \
  --workflow-file "$HERE/scripts/count_concordance.py" \
  --workflow-file "$HERE/scripts/jaccard_confusion.py" \
  --workflow-file "$HERE/scripts/read_set_jaccard_bitmaps.py"
"$PYBIN" "$HERE/scripts/summarize_downstream_run.py" "$OUT" \
  --json "$OUT/analysis/downstream_summary.json" \
  --markdown "$OUT/analysis/downstream_summary.md"

echo "[6/6] refresh provenance after summary creation and bundle compact outputs"
"$PYBIN" "$HERE/scripts/downstream_provenance.py" \
  --output "$OUT/downstream_provenance.json" --results-root "$OUT" \
  --genome "$GENOME" --star "$STARBIN" --python "$PYBIN" --whitelist "$WL" \
  --threads "$THREADS" --min-umi "$MIN_UMI" \
  --cb-match-whitelist-type "$CB_MATCH" \
  --fastq-provenance "$OUT/seqproc_fastq_provenance.json" \
  --fastq-provenance "$OUT/splitcode_fastq_provenance.json" \
  --fastq-provenance "$OUT/matchbox_fastq_provenance.json" \
  --workflow-file "$HERE/run_matchbox_ham1_sensitivity.sh" \
  --workflow-file "$HERE/requirements.lock"
tar -czf "$OUT/downstream_bundle.tar.gz" -C "$OUT" \
  analysis accuracy.json resources.csv downstream_provenance.json STAR.version.txt \
  seqproc_fastq_provenance.json splitcode_fastq_provenance.json matchbox_fastq_provenance.json \
  sp_Solo.out/Gene/Summary.csv sc_Solo.out/Gene/Summary.csv mb_Solo.out/Gene/Summary.csv \
  sp_Solo.out/Barcodes.stats sc_Solo.out/Barcodes.stats mb_Solo.out/Barcodes.stats
echo "DONE: $OUT"
