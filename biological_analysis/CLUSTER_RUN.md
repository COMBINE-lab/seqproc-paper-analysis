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

## Notes
- Full 86.8M-read STARsolo is the long step. If your interactive `srun` might disconnect, wrap the
  step 3 command in an `sbatch` script instead so it survives the session ending.
- Everything else (read processing, analysis, figures, resource report) finishes in well under a
  minute on the box, so the full-data wall time is dominated by STAR alignment.
- The outputs are identical in structure to the box rehearsal, just full-depth numbers, plus you can
  add the split-pipe vendor comparison since the vendor matrix is on the cluster.
