# Forward-only supplementary counterpart of publication_lr_splitseq_dual_core.
prefix = r:
umi = u[10]
#[ambig_policy = accept]
bc3 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc8_whitelist.txt", 1)
#[search(relative)] #[edit(3)] linker_a = f[GTGGCCGATGTTTCGCATCGGCGTACGACT]
#[ambig_policy = accept]
bc2 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc8_whitelist.txt", 1)
#[edit(3)] linker_b = f[ATCCACGTGCTTGAGACTGTGG]
#[ambig_policy = accept]
bc1 = filter_within_dist(b[6], "configs/seqproc/splitseq_bc1_whitelist_6bp.txt", 1)
rest = r:

1{<prefix><umi><bc3><linker_a><bc2><linker_b><bc1><rest>}

-> 1{<umi><bc3><bc2><bc1>}
