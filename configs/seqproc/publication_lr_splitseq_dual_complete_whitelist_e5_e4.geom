# Optional high-specificity LR-SPLiT-seq mode.
#
# In addition to requiring complete fixed-length components, each observed
# barcode must be within Hamming distance one of the corresponding whitelist.
# `accept` is explicit because the structural reference accepts a sequence
# that is close to more than one distinct whitelist entry; this filter does
# not select or substitute a canonical barcode.

prefix = r:
umi = u[10]
#[ambig_policy = accept]
bc3 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc23_whitelist.txt", 1)
#[search(relative)] #[edit(5)] linker_a = f[GTGGCCGATGTTTCGCATCGGCGTACGACT]
between_linkers = r:
#[ambig_policy = accept]
bc2 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc23_whitelist.txt", 1)
#[search(relative)] #[edit(4)] linker_b = f[ATCCACGTGCTTGAGACTGTGG]
#[ambig_policy = accept]
bc1 = filter_within_dist(b[6], "configs/seqproc/splitseq_bc1_whitelist_6bp.txt", 1)
rest = r:

#[match_ori(either)]
1{<prefix><umi><bc3><linker_a><between_linkers><bc2><linker_b><bc1><rest>}

-> 1{<umi><bc3><bc2><bc1>}
