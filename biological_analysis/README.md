# Biological validation

Downstream, biology-level validation of seqproc. See
`paper/biological_validation_prd.md` for the full plan.

Workflow: the analysis runs on the cluster against the full GRCm38 reference.
(An earlier chr19 dev-index rehearsal step was removed once the full-data run
was in place.)

## Downstream validation: SPLiT-seq, are splitcode's extra reads real or junk

Run seqproc and splitcode output through one STARsolo and compare the count
matrices, plus an isolated-extras analysis on the reads only splitcode keeps.

Status: complete. Full-data results in `full_run_results/` feed the paper's `sec:downstream`, `fig:count_concordance`, and `tab:jaccard`.

## Layout
- `configs/`  tool and STARsolo configs (committed)
- `scripts/`  pipeline scripts (committed)
- `refs/`     STAR indices (gitignored, regenerable)
- `results/`  outputs (gitignored)
- `notebooks/` final figures (committed once stable)

A full reproducibility README is written at the end of the leg.
