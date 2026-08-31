# Corrected cross-tool SPLiT-seq PE geometry.  The physical R2 layout is
# UMI(10)+BC3(8)+L1(30)+BC2(8)+L2(30)+BC1(8).  The preprint configuration
# accidentally modeled a 2-nt skip + 8-nt UMI, transcribed one L1 base, and
# truncated BC1 to the 6-nt form used by the long-read protocol.
read1 = r:
umi = u[10]
bc3 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc8_whitelist.txt", 1)
#[search(relative)] #[edit(3)] l1 = f[GTGGCCGATGTTTCGCATCGGCGTACGACT]
bc2 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc8_whitelist.txt", 1)
#[search(relative)] #[edit(3)] l2 = f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG]
bc1 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc8_whitelist.txt", 1)
rest = r:

1{<read1>}
2{<umi><bc3><l1><bc2><l2><bc1><rest>}
-> 1{<read1>} 2{<umi><bc3><bc2><bc1>}
