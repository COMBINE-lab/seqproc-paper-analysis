#!/usr/bin/env bash
# Pairwise concordance by linker tolerance: run the full pipeline at tolerances 1-4 (seqproc edit,
# matchbox edit, splitcode Hamming) and record per-gene Pearson(log) between tool pairs.
set -uo pipefail
cd /home/ubuntu/seqproc-paper-analysis-clean
SP=/home/ubuntu/combine-lab/seqproc/target/release/seqproc
MB=/home/ubuntu/matchbox/target/release/matchbox
SC=/home/ubuntu/splitcode/build/src/splitcode
STAR=/usr/lib/rna-star/bin/STAR-avx2
PY=/home/ubuntu/.venv_nb/bin/python
IDX=biological_analysis/refs/star_GRCm38
WL=biological_analysis/configs/splitseq_bc_whitelist_96.txt
R1=data/SRR6750041_1M_R1.fastq; R2=data/SRR6750041_1M_R2.fastq
G=biological_analysis/configs/splitseq_quant_observed.geom
MBC=biological_analysis/configs/splitseq_matchbox.mb
SCC=biological_analysis/configs/splitseq_extract.config
XPAT='1:0<u[10]>,<b3[8]>{linker1},{linker1}<b2[8]>,{linker2}<b1[8]>'
rm -rf /tmp/consw/T*; mkdir -p /tmp/consw
echo "tol,pair,per_gene_pearson_log,per_gene_spearman" > /tmp/consw/conc.csv
for T in 1 2 3 4; do
  O=/tmp/consw/T$T; mkdir -p "$O"
  sed "s/edit([0-9]*)/edit($T)/g" "$G" > "$O/g.geom"
  "$SP" -g "$O/g.geom" -1 "$R1" -2 "$R2" -o "$O/sp_cdna.fq" -w "$O/sp_bc.fq" -t 8 >/dev/null 2>&1
  err=$(awk -v n=$T 'BEGIN{printf "%.4f",(n+0.3)/30}')
  sed "s/linker1~[0-9.]*/linker1~$err/; s/linker2~[0-9.]*/linker2~$err/" "$MBC" > "$O/mb.mb"
  "$MB" -s "$O/mb.mb" -t 8 "$R2" 2>/dev/null | "$PY" biological_analysis/scripts/matchbox_quant_extract.py --r1 "$R1" --out-cdna "$O/mb_cdna.fq" --out-bc "$O/mb_bc.fq" >/dev/null 2>&1
  awk -F'\t' -v d=$T 'BEGIN{OFS="\t"} /^linker[12]\t/{$4=d} {print}' "$SCC" > "$O/sc.config"
  "$SC" -c "$O/sc.config" -N 2 -t 8 --x-only -p -x "$XPAT" "$R1" "$R2" 2>/dev/null | "$PY" biological_analysis/scripts/splitcode_quant_extract.py --r1 "$R1" --out-cdna "$O/sc_cdna.fq" --out-bc "$O/sc_bc.fq" >/dev/null 2>&1
  for t in sp sc mb; do
    "$STAR" --runThreadN 8 --genomeDir "$IDX" --soloType CB_UMI_Complex \
      --readFilesIn "$O/${t}_cdna.fq" "$O/${t}_bc.fq" \
      --soloCBwhitelist "$WL" "$WL" "$WL" --soloCBposition 0_10_0_17 0_18_0_25 0_26_0_33 \
      --soloUMIposition 0_0_0_9 --soloCBmatchWLtype 1MM --soloFeatures Gene --soloCellFilter None \
      --outSAMtype None --outFileNamePrefix "$O/${t}_" >/dev/null 2>&1
  done
  "$PY" biological_analysis/scripts/count_concordance.py "$O/analysis" \
    "seqproc:$O/sp_Solo.out/Gene" "splitcode:$O/sc_Solo.out/Gene" "matchbox:$O/mb_Solo.out/Gene" >/dev/null 2>&1
  "$PY" -c "
import json
d=json.load(open('$O/analysis/count_concordance.json'))
pg=d['per_gene_total_pearson_logspace']; sp=d['per_gene_total_spearman']
for k in pg: print(f'$T,{k},{pg[k]},{sp[k]}')
" >> /tmp/consw/conc.csv
  echo "tolerance $T done"
done
echo "SWEEP DONE"
cat /tmp/consw/conc.csv
