# SPLiT-seq Round 2 - Approach 1: Anchors + Filter
# Uses anchor_relative for linkers and filter_within_dist for barcodes

read1 = r:
skip2 = x[2]
umi = u[8]

# BC3: 8bp, filter against whitelist (dist 1)
bc3 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc23_whitelist.txt", 1)

# Linker1 - 30bp - Anchor Relative (Hamming 3)
l1 = anchor_relative(hamming(f[GTGGCCGCTGTTTCGCATCGGCGTACGACT], 3))

# BC2: 8bp, filter against whitelist (dist 1)
bc2 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc23_whitelist.txt", 1)

# Linker2 - 30bp - Anchor Relative (Hamming 3)
l2 = anchor_relative(hamming(f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG], 3))

# BC1: 8bp (Full), filter against whitelist (dist 1)
bc1 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc1_whitelist.txt", 1)
rest = r:

1{<read1>}
2{<skip2><umi><bc3><l1><bc2><l2><bc1><rest>}

-> 1{<read1>} 2{<umi><bc3><bc2><bc1>}
