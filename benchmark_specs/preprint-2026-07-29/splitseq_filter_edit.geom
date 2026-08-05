# SPLiT-seq paired-end geometry from the analysis-repository snapshot used by
# the July 29, 2026 preprint.

read1 = r:
skip2 = x[2]
umi = u[10]

bc3 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc23_whitelist.txt", 1)
#[search(relative)] #[edit(6)] l1 = f[GTGGCCGCTGTTTCGCATCGGCGTACGACT]

bc2 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc23_whitelist.txt", 1)
#[search(relative)] #[edit(6)] l2 = f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG]

bc1 = filter_within_dist(b[6], "configs/seqproc/splitseq_bc1_whitelist_6bp.txt", 1)
rest = r:

1{<read1>}
2{<skip2><umi><bc3><l1><bc2><l2><bc1><rest>}

-> 1{<read1>} 2{<umi><bc3><bc2><bc1>}
