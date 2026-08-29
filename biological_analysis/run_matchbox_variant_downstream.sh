#!/usr/bin/env bash
# Recompute the Matchbox-dependent part of the SPLiT-seq PE downstream analysis
# while reusing unchanged seqproc and splitcode STARsolo matrices from a pinned
# base run. This supports both the primary canonical-barcode configuration and
# externally expanded-barcode sensitivity analyses without mislabeling either.
set -euo pipefail

THREADS=32
MIN_UMI=200
CB_MATCH=EditDist_2
INPUT_RECORDS=77621181
PYBIN=""
GENOME=""
BASE_RUN=""
OUT=""
VARIANT=""
MATCHBOX_R1=""
MATCHBOX_R2=""
MATCHBOX_PROVENANCE=""
REFERENCE_BITMAP=""
SEQPROC_BITMAP=""
SPLITCODE_BITMAP=""
MATCHBOX_BITMAP=""

usage() {
  cat <<'EOF'
usage: run_matchbox_variant_downstream.sh --genome STAR_INDEX \
  --base-run BASE_DOWNSTREAM_RUN --outdir OUT --variant NAME \
  --matchbox-r1 FILE --matchbox-r2 FILE \
  --matchbox-fastq-provenance FILE \
  --reference-bitmap FILE --seqproc-bitmap FILE --splitcode-bitmap FILE \
  --matchbox-bitmap FILE [--input-records 77621181] \
  [--threads 32] [--min-umi 200] [--cb-match EditDist_2] [--python PYTHON]

The base run must use the requested STARsolo barcode-correction mode. Its
seqproc and splitcode matrices and FASTQ provenance are reused; only Matchbox
is remapped and the joint downstream analysis is regenerated.
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --genome) GENOME=$2; shift 2;;
    --base-run) BASE_RUN=$2; shift 2;;
    --outdir) OUT=$2; shift 2;;
    --variant) VARIANT=$2; shift 2;;
    --matchbox-r1) MATCHBOX_R1=$2; shift 2;;
    --matchbox-r2) MATCHBOX_R2=$2; shift 2;;
    --matchbox-fastq-provenance) MATCHBOX_PROVENANCE=$2; shift 2;;
    --reference-bitmap) REFERENCE_BITMAP=$2; shift 2;;
    --seqproc-bitmap) SEQPROC_BITMAP=$2; shift 2;;
    --splitcode-bitmap) SPLITCODE_BITMAP=$2; shift 2;;
    --matchbox-bitmap) MATCHBOX_BITMAP=$2; shift 2;;
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
WL=$ROOT/configs/seqproc/splitseq_bc8_whitelist.txt
[ -n "$PYBIN" ] || PYBIN=$HERE/.venv_downstream/bin/python

for var in GENOME BASE_RUN OUT VARIANT MATCHBOX_R1 MATCHBOX_R2 \
  MATCHBOX_PROVENANCE REFERENCE_BITMAP SEQPROC_BITMAP SPLITCODE_BITMAP \
  MATCHBOX_BITMAP; do
  [ -n "${!var}" ] || { echo "missing required argument: $var" >&2; usage >&2; exit 2; }
done
for path in "$GENOME" "$BASE_RUN" "$MATCHBOX_R1" "$MATCHBOX_R2" \
  "$MATCHBOX_PROVENANCE" "$REFERENCE_BITMAP" "$SEQPROC_BITMAP" \
  "$SPLITCODE_BITMAP" "$MATCHBOX_BITMAP" "$STARBIN" "$PYBIN" "$WL"; do
  [ -e "$path" ] || { echo "required path does not exist: $path" >&2; exit 2; }
done
for path in "$BASE_RUN/downstream_provenance.json" \
  "$BASE_RUN/seqproc_fastq_provenance.json" \
  "$BASE_RUN/splitcode_fastq_provenance.json" \
  "$BASE_RUN/sp_Solo.out/Gene/raw/matrix.mtx" \
  "$BASE_RUN/sc_Solo.out/Gene/raw/matrix.mtx"; do
  [ -e "$path" ] || { echo "base-run artifact does not exist: $path" >&2; exit 2; }
done

base_cb_match=$("$PYBIN" -c \
  'import json,sys; print(json.load(open(sys.argv[1]))["parameters"]["cb_match_whitelist_type"])' \
  "$BASE_RUN/downstream_provenance.json")
[ "$base_cb_match" = "$CB_MATCH" ] || {
  echo "base run used $base_cb_match but variant requested $CB_MATCH" >&2
  exit 2
}

mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
BASE_RUN=$(readlink -f "$BASE_RUN")
GENOME=$(readlink -f "$GENOME")
MATCHBOX_R1=$(readlink -f "$MATCHBOX_R1")
MATCHBOX_R2=$(readlink -f "$MATCHBOX_R2")

"$PYBIN" - "$MATCHBOX_PROVENANCE" "$MATCHBOX_R1" "$MATCHBOX_R2" <<'PY'
import json
import os
import sys

record = json.load(open(sys.argv[1]))
expected = (os.path.realpath(sys.argv[2]), os.path.realpath(sys.argv[3]))
observed = (os.path.realpath(record["r1"]["path"]), os.path.realpath(record["r2"]["path"]))
if observed != expected:
    raise SystemExit(f"FASTQ provenance paths {observed!r} do not match requested inputs {expected!r}")
if record.get("records", 0) <= 0:
    raise SystemExit("FASTQ provenance has no positive record count")
PY

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONHASHSEED=0 MPLCONFIGDIR="$OUT/.matplotlib"
mkdir -p "$MPLCONFIGDIR" "$OUT/analysis"

echo "[1/5] stage frozen provenance and evaluate read-set agreement ($VARIANT)"
cp "$BASE_RUN/seqproc_fastq_provenance.json" "$OUT/seqproc_fastq_provenance.json"
cp "$BASE_RUN/splitcode_fastq_provenance.json" "$OUT/splitcode_fastq_provenance.json"
cp "$MATCHBOX_PROVENANCE" "$OUT/matchbox_fastq_provenance.json"
cp "$MATCHBOX_BITMAP" "$OUT/matchbox.bitmap"
"$PYBIN" "$ROOT/scripts/evaluate_numeric_bitmaps.py" \
  --reference "$REFERENCE_BITMAP" --input-records "$INPUT_RECORDS" \
  --tool "seqproc=$SEQPROC_BITMAP" --tool "splitcode=$SPLITCODE_BITMAP" \
  --tool "matchbox=$OUT/matchbox.bitmap" \
  --out "$OUT/accuracy.json" > "$OUT/accuracy.log"

echo "[2/5] STARsolo for Matchbox only ($CB_MATCH)"
ln -sfn "$BASE_RUN/sp_Solo.out" "$OUT/sp_Solo.out"
ln -sfn "$BASE_RUN/sc_Solo.out" "$OUT/sc_Solo.out"
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
  "$BASE_RUN/resources.csv" > "$OUT/resources.csv"
echo "starsolo,matchbox,$seconds,$rss" >> "$OUT/resources.csv"

ARGS=(
  "seqproc:$OUT/sp_Solo.out/Gene"
  "splitcode:$OUT/sc_Solo.out/Gene"
  "matchbox:$OUT/mb_Solo.out/Gene"
)

echo "[3/5] regenerate biological and count-concordance analyses"
"$PYBIN" "$HERE/scripts/biological_analysis.py" "$OUT/analysis" "$MIN_UMI" "${ARGS[@]}"
"$PYBIN" "$HERE/scripts/count_concordance.py" "$OUT/analysis" "${ARGS[@]}"
"$PYBIN" "$HERE/scripts/jaccard_confusion.py" "$OUT/analysis" "$MIN_UMI" "${ARGS[@]}" \
  > "$OUT/analysis/jaccard_confusion.log"
"$PYBIN" "$HERE/scripts/read_set_jaccard_bitmaps.py" \
  --output "$OUT/analysis/read_set_jaccard.json" \
  "seqproc:$SEQPROC_BITMAP" "splitcode:$SPLITCODE_BITMAP" \
  "matchbox:$OUT/matchbox.bitmap"

echo "[4/5] provenance and compact summary"
"$PYBIN" "$HERE/scripts/downstream_provenance.py" \
  --output "$OUT/downstream_provenance.json" --results-root "$OUT" \
  --genome "$GENOME" --star "$STARBIN" --python "$PYBIN" --whitelist "$WL" \
  --threads "$THREADS" --min-umi "$MIN_UMI" \
  --cb-match-whitelist-type "$CB_MATCH" \
  --fastq-provenance "$OUT/seqproc_fastq_provenance.json" \
  --fastq-provenance "$OUT/splitcode_fastq_provenance.json" \
  --fastq-provenance "$OUT/matchbox_fastq_provenance.json" \
  --workflow-file "$HERE/run_matchbox_variant_downstream.sh" \
  --workflow-file "$HERE/requirements.lock" \
  --workflow-file "$HERE/scripts/biological_analysis.py" \
  --workflow-file "$HERE/scripts/count_concordance.py" \
  --workflow-file "$HERE/scripts/jaccard_confusion.py" \
  --workflow-file "$HERE/scripts/read_set_jaccard_bitmaps.py"
"$PYBIN" "$HERE/scripts/summarize_downstream_run.py" "$OUT" \
  --json "$OUT/analysis/downstream_summary.json" \
  --markdown "$OUT/analysis/downstream_summary.md"

echo "[5/5] refresh provenance and bundle compact outputs"
"$PYBIN" "$HERE/scripts/downstream_provenance.py" \
  --output "$OUT/downstream_provenance.json" --results-root "$OUT" \
  --genome "$GENOME" --star "$STARBIN" --python "$PYBIN" --whitelist "$WL" \
  --threads "$THREADS" --min-umi "$MIN_UMI" \
  --cb-match-whitelist-type "$CB_MATCH" \
  --fastq-provenance "$OUT/seqproc_fastq_provenance.json" \
  --fastq-provenance "$OUT/splitcode_fastq_provenance.json" \
  --fastq-provenance "$OUT/matchbox_fastq_provenance.json" \
  --workflow-file "$HERE/run_matchbox_variant_downstream.sh" \
  --workflow-file "$HERE/requirements.lock" \
  --workflow-file "$HERE/scripts/summarize_downstream_run.py"
"$PYBIN" - "$OUT/run_metadata.json" "$VARIANT" "$CB_MATCH" <<'PY'
import json
import sys
from datetime import datetime, timezone

json.dump(
    {
        "schema_version": 1,
        "variant": sys.argv[2],
        "cb_match_whitelist_type": sys.argv[3],
        "created_utc": datetime.now(timezone.utc).isoformat(),
    },
    open(sys.argv[1], "w"),
    indent=2,
)
PY
tar -czf "$OUT/downstream_bundle.tar.gz" -C "$OUT" \
  analysis accuracy.json resources.csv run_metadata.json downstream_provenance.json STAR.version.txt \
  seqproc_fastq_provenance.json splitcode_fastq_provenance.json matchbox_fastq_provenance.json \
  sp_Solo.out/Gene/Summary.csv sc_Solo.out/Gene/Summary.csv mb_Solo.out/Gene/Summary.csv \
  sp_Solo.out/Barcodes.stats sc_Solo.out/Barcodes.stats mb_Solo.out/Barcodes.stats
echo "DONE: $OUT"
