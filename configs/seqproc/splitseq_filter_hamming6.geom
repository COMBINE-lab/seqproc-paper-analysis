# SPLiT-seq Round 2 - Approach 1: Anchors + Filter using HAMMING DISTANCE (tolerance 6)
# Uses anchor_relative for linkers with hamming distance and filter_within_dist for barcodes
# Tolerance 6 = ~0.2 fraction of 30bp linker, matching edit(6) and matchbox -e 0.2

read1 = r:
skip2 = x[2]
umi = u[10]

# BC3: 8bp, filter against whitelist (dist 1)
bc3 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc23_whitelist.txt", 1)

# Linker1 - 30bp - Anchor Relative (Hamming distance 6, ~0.2 fraction)
#[search(relative)] #[hamming(6)] l1 = f[GTGGCCGCTGTTTCGCATCGGCGTACGACT]

# BC2: 8bp, filter against whitelist (dist 1)
bc2 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc23_whitelist.txt", 1)

# Linker2 - 30bp - Anchor Relative (Hamming distance 6, ~0.2 fraction)
#[search(relative)] #[hamming(6)] l2 = f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG]

# BC1: 6bp (truncated in 94bp read with 10bp UMI), filter against whitelist (dist 1)
bc1 = filter_within_dist(b[6], "configs/seqproc/splitseq_bc1_whitelist_6bp.txt", 1)
rest = r:

1{<read1>}
2{<skip2><umi><bc3><l1><bc2><l2><bc1><rest>}

-> 1{<read1>} 2{<umi><bc3><bc2><bc1>}
