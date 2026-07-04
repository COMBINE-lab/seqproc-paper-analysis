# Cluster run (Phase 2A) on UMIACS nexus

You drive this. The sandbox box cannot reach the cluster (VPN plus your credentials), so the code
is pushed to GitHub from the box and you pull it on the cluster, then run the one command. Paste
the run log back and I debug from there.

## 0. Get on the cluster (your VPN, then ssh, then an interactive node)
```bash
ssh ejfisher@nexuscbcb.umiacs.umd.edu
# request a node with enough RAM for STARsolo + the full genome index (~30GB),
# and enough wall time for full-data STARsolo. Adjust partition to your cluster.
srun --pty -c 8 --mem=48G -t 8:00:00 bash
```

## 1. Pull the latest code
```bash
cd /nfshomes/ejfisher/seqproc-paper-analysis-clean   # wherever the repo lives
git pull
```

## 2. One-time setup (skip whatever already exists)
```bash
bash biological_analysis/setup_env.sh                # python env (scanpy etc.)
# STAR index: point --genome at an existing GRCm38 index if you have one, or build once.
# On a big-RAM node drop --genomeSAsparseD (it is only needed when memory is tight).
STAR --runMode genomeGenerate --genomeDir <IDX> \
     --genomeFastaFiles <GRCm38.primary_assembly.fa> \
     --sjdbGTFfile <GRCm38.102.gtf> --sjdbOverhang 65 --runThreadN 8
```

## 3. Run the full dataset (one command)
```bash
SEQPROC_BIN=<path/to/seqproc> \
SPLITCODE_BIN=<path/to/splitcode> \
MATCHBOX_BIN=<path/to/matchbox> \
biological_analysis/run_phase2a.sh \
  --r1 <full_R1.fastq> --r2 <full_R2.fastq> \
  --genome <IDX> --outdir <OUT> --threads 8
```
STAR and the three tool binaries must be on PATH or pointed at via the `*_BIN` env vars.

## 4. Bring results back
```bash
# the small bundle (figures, tables, metrics) is all you need to review:
scp <OUT>/phase2a_bundle.tar.gz  <your machine or the box>
```

## 5. Refresh downstream metrics/figures without re-aligning
When only the analysis code changed (a new metric, a relabelled figure), re-run the analysis
scripts on the STARsolo matrices already on disk. This is a ~1 minute step and does not re-run
STARsolo. `biological_analysis.py` must run before `make_downstream_figure.py` (the latter reads
its JSON).
```bash
cd /nfshomes/ejfisher/seqproc-paper-analysis && git pull
PY=biological_analysis/.venv_phase2a/bin/python
OUT=/fs/nexus-projects/seqproc/bench/phase2a_out
$PY biological_analysis/scripts/biological_analysis.py $OUT/analysis 200 \
  seqproc:$OUT/sp_Solo.out/Gene splitcode:$OUT/sc_Solo.out/Gene matchbox:$OUT/mb_Solo.out/Gene
$PY biological_analysis/scripts/count_concordance.py $OUT/analysis \
  seqproc:$OUT/sp_Solo.out/Gene splitcode:$OUT/sc_Solo.out/Gene matchbox:$OUT/mb_Solo.out/Gene
$PY biological_analysis/scripts/make_downstream_figure.py $OUT/analysis 200 \
  seqproc:$OUT/sp_Solo.out/Gene splitcode:$OUT/sc_Solo.out/Gene matchbox:$OUT/mb_Solo.out/Gene
$PY biological_analysis/scripts/read_set_jaccard.py $OUT/analysis/read_set_jaccard.json \
  seqproc:$OUT/sp_bc.fq splitcode:$OUT/sc_bc.fq matchbox:$OUT/mb_bc.fq
```
Then commit the refreshed outputs back to the repo:
```bash
cp $OUT/analysis/biological_metrics.json $OUT/analysis/count_concordance.json \
   $OUT/analysis/read_set_jaccard.json $OUT/analysis/jaccard_supplement.md \
   biological_analysis/full_run_results/
git add biological_analysis/full_run_results/
git commit --no-verify -m "refresh downstream metrics from cluster"
git push origin biological-validation
```

## 6. Barcode-rank knee + cluster-ARI stability
`count_concordance.py` now reports the barcode-rank knee per tool (kneedle on the log-log curve; a
chord-distance fallback runs if `kneed` is absent). The knee lands in `count_concordance.json` under
`barcode_rank_knee`. `ari_stability.py` reclusters across many Leiden seeds (and resolutions) to test
whether the pairwise cluster-ARI gap is a real effect or run-to-run noise; its `reference_seed0_res1.0`
values should reproduce the paper table (0.909 / 0.644 / 0.627).
```bash
cd /nfshomes/ejfisher/seqproc-paper-analysis && git pull
PY=biological_analysis/.venv_phase2a/bin/python
OUT=/fs/nexus-projects/seqproc/bench/phase2a_out
biological_analysis/.venv_phase2a/bin/pip install kneed        # canonical kneedle (one time)

# knee: re-run count_concordance.py (also part of the section-5 refresh)
$PY biological_analysis/scripts/count_concordance.py $OUT/analysis \
  seqproc:$OUT/sp_Solo.out/Gene splitcode:$OUT/sc_Solo.out/Gene matchbox:$OUT/mb_Solo.out/Gene

# cluster-ARI stability: 50 seeds at res 1.0 plus a resolution sweep (~a few minutes)
$PY biological_analysis/scripts/ari_stability.py $OUT/analysis 50 200 \
  seqproc:$OUT/sp_Solo.out/Gene splitcode:$OUT/sc_Solo.out/Gene matchbox:$OUT/mb_Solo.out/Gene \
  --res=0.5,1.0,1.5,2.0
```
Then commit the outputs back:
```bash
cp $OUT/analysis/count_concordance.json $OUT/analysis/ari_stability.json \
   biological_analysis/full_run_results/
git add biological_analysis/full_run_results/
git commit --no-verify -m "add barcode-rank knee + cluster-ARI stability from cluster"
git push origin biological-validation
```

## Notes
- Full 86.8M-read STARsolo is the long step. If your interactive `srun` might disconnect, wrap the
  step 3 command in an `sbatch` script instead so it survives the session ending.
- Everything else (read processing, analysis, figures, resource report) finishes in well under a
  minute on the box, so the full-data wall time is dominated by STAR alignment.
- The outputs are identical in structure to the box rehearsal, just full-depth numbers, plus you can
  add the split-pipe vendor comparison since the vendor matrix is on the cluster.
