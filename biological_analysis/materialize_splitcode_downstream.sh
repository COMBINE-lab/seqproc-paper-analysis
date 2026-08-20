#!/usr/bin/env bash
# Re-run the frozen publication splitcode configuration once while retaining
# the accepted cDNA mate that the accuracy campaign intentionally discarded.
set -euo pipefail

RAW_R1=""; RAW_R2=""; OUT=""; THREADS=32
usage() {
  echo "usage: materialize_splitcode_downstream.sh --r1 RAW_R1 --r2 RAW_R2 --outdir OUT [--threads 32]" >&2
}
while [ "$#" -gt 0 ]; do
  case "$1" in
    --r1) RAW_R1=$2; shift 2;;
    --r2) RAW_R2=$2; shift 2;;
    --outdir) OUT=$2; shift 2;;
    --threads) THREADS=$2; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "unknown argument: $1" >&2; usage; exit 2;;
  esac
done
[ -n "$RAW_R1" ] && [ -n "$RAW_R2" ] && [ -n "$OUT" ] || { usage; exit 2; }

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
SPLITCODE=${SPLITCODE_BIN:-$ROOT/../competitors/splitcode-v0.31.6/build/src/splitcode}
CONFIG=$ROOT/configs/splitcode/publication_splitseq_pe.config
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
RAW_R1=$(readlink -f "$RAW_R1"); RAW_R2=$(readlink -f "$RAW_R2")
CDNA=$OUT/sc_cdna.fastq
BC=$OUT/umi_bc3_bc2_bc1.fastq

if [ -e "$CDNA" ] || [ -e "$BC" ]; then
  echo "refusing to overwrite an existing materialized output in $OUT" >&2
  exit 2
fi

cd "$OUT"
python3 - "$OUT/command.json" "$SPLITCODE" "$CONFIG" "$RAW_R1" "$RAW_R2" "$THREADS" <<'PY'
import json, sys
out, binary, config, r1, r2, threads = sys.argv[1:]
argv = [binary, "--config", config, "--assign", "--mapping", "mapping.txt",
        "--no-outb", "--select", "0", "--nFastqs", "2", "--threads", threads,
        "--output", "sc_cdna.fastq,/dev/null", r1, r2]
open(out, "w").write(json.dumps({"argv": argv, "cwd": str(__import__('pathlib').Path.cwd())}, indent=2) + "\n")
PY
/usr/bin/time -v "$SPLITCODE" --config "$CONFIG" --assign \
  --mapping "$OUT/mapping.txt" --no-outb --select 0 --nFastqs 2 \
  --threads "$THREADS" --output "$CDNA,/dev/null" "$RAW_R1" "$RAW_R2" \
  > "$OUT/stdout.txt" 2> "$OUT/time-and-stderr.txt"

"$HERE/.venv_downstream/bin/python" "$HERE/scripts/validate_fastq_pairs.py" \
  --name splitcode --r1 "$CDNA" --r2 "$BC" \
  --expected-r1-length 66 --expected-r2-length 34 \
  --output "$OUT/fastq_provenance.json" > "$OUT/validation.log"
echo "splitcode downstream pair ready: $OUT"
