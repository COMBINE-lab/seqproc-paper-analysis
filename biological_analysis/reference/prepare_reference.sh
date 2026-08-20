#!/usr/bin/env bash
# Download, verify, and index the exact mouse reference used downstream.
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
BIO=$(cd "$HERE/.." && pwd)
DOWNLOADS=$HERE/downloads
SOURCE=$HERE/GRCm38_ensembl102
INDEX=$HERE/star_GRCm38_ensembl102
TOOLS=$BIO/tools

FA_GZ=Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
GTF_GZ=Mus_musculus.GRCm38.102.gtf.gz
STAR_TAR=STAR-2.7.11b.tar.gz
FA_URL=https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/$FA_GZ
GTF_URL=https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/$GTF_GZ
STAR_URL=https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz
FA_SHA=285bc481d583ab65b13d91853bf743acf950710afb3302264a4b4f116b6049c1
GTF_SHA=8321415404aaf788c7da79774488ff227ac006d09a57ce6c616573a510338f64
STAR_SHA=3f65305e4112bd154c7e22b333dcdaafc681f4a895048fa30fa7ae56cac408e7

mkdir -p "$DOWNLOADS" "$SOURCE" "$INDEX" "$TOOLS/downloads"
download() {
  local url=$1 output=$2
  [ -f "$output" ] || curl -L --fail --retry 5 --continue-at - --output "$output" "$url"
}
verify() {
  local expected=$1 path=$2 actual
  actual=$(sha256sum "$path" | awk '{print $1}')
  [ "$actual" = "$expected" ] || {
    echo "checksum mismatch for $path: $actual != $expected" >&2
    exit 1
  }
}

download "$FA_URL" "$DOWNLOADS/$FA_GZ"
download "$GTF_URL" "$DOWNLOADS/$GTF_GZ"
download "$STAR_URL" "$TOOLS/downloads/$STAR_TAR"
verify "$FA_SHA" "$DOWNLOADS/$FA_GZ"
verify "$GTF_SHA" "$DOWNLOADS/$GTF_GZ"
verify "$STAR_SHA" "$TOOLS/downloads/$STAR_TAR"

if [ ! -x "$TOOLS/STAR-2.7.11b/bin/STAR" ]; then
  [ -d "$TOOLS/STAR-2.7.11b" ] || tar -xzf "$TOOLS/downloads/$STAR_TAR" -C "$TOOLS"
  make -C "$TOOLS/STAR-2.7.11b/source" STAR -j "${BUILD_THREADS:-32}"
  mkdir -p "$TOOLS/STAR-2.7.11b/bin"
  cp "$TOOLS/STAR-2.7.11b/source/STAR" "$TOOLS/STAR-2.7.11b/bin/STAR"
fi
STARBIN=$TOOLS/STAR-2.7.11b/bin/STAR

FA=${FA_GZ%.gz}
GTF=${GTF_GZ%.gz}
if [ ! -f "$SOURCE/$FA" ]; then
  cp "$DOWNLOADS/$FA_GZ" "$SOURCE/$FA_GZ"
  gzip -dk "$SOURCE/$FA_GZ"
fi
if [ ! -f "$SOURCE/$GTF" ]; then
  cp "$DOWNLOADS/$GTF_GZ" "$SOURCE/$GTF_GZ"
  gzip -dk "$SOURCE/$GTF_GZ"
fi

if [ ! -s "$INDEX/Genome" ]; then
  "$STARBIN" --runMode genomeGenerate --runThreadN "${INDEX_THREADS:-32}" \
    --genomeDir "$INDEX" --genomeFastaFiles "$SOURCE/$FA" \
    --sjdbGTFfile "$SOURCE/$GTF" --sjdbOverhang 65
fi

python3 "$HERE/write_reference_manifest.py" \
  --output "$INDEX/seqproc_reference_manifest.json" \
  --star "$STARBIN" --index "$INDEX" \
  --fasta "$SOURCE/$FA" --fasta-url "$FA_URL" --fasta-gzip-sha256 "$FA_SHA" \
  --gtf "$SOURCE/$GTF" --gtf-url "$GTF_URL" --gtf-gzip-sha256 "$GTF_SHA" \
  --star-source-url "$STAR_URL" --star-source-sha256 "$STAR_SHA" \
  --sjdb-overhang 65
echo "reference ready: $INDEX"
