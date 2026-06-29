# Biological validation (Phase 2)

Downstream, biology-level validation of seqproc. See
`paper/biological_validation_prd.md` for the full plan.

Workflow: develop each phase here on a small sample with the chr19 reference,
prove the pipeline runs, then run the real full-data analysis on the cluster
with the full GRCm38 reference. One phase at a time.

## Phase 2A. SPLiT-seq, are splitcode's extra reads real or junk

Run seqproc and splitcode output through one STARsolo and compare the count
matrices, plus an isolated-extras analysis on the reads only splitcode keeps.

Status: scaffolding. chr19 dev index built. Tool configs in progress.

## Layout
- `configs/`  tool and STARsolo configs (committed)
- `scripts/`  pipeline scripts (committed)
- `refs/`     STAR indices (gitignored, regenerable)
- `results/`  outputs (gitignored)
- `notebooks/` final figures (committed once stable)

A full reproducibility README is written at the end of the leg.
