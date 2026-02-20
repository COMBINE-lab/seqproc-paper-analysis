# SPLiT-seq Round 2 - Approach 2: Positional + Filter
# Uses x[30] to skip linkers (assuming fixed position) and filter_within_dist for barcodes

read1 = r:
skip2 = x[2]
umi = u[10]

# BC3: 8bp, filter against whitelist (dist 1)
bc3 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc23_whitelist.txt", 1)

# Linker1 - 30bp - Skip (Positional assumption)
skip_l1 = x[30]

# BC2: 8bp, filter against whitelist (dist 1)
bc2 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc23_whitelist.txt", 1)

# Linker2 - 30bp - Skip (Positional assumption)
skip_l2 = x[30]

# BC1: 6bp (truncated), filter against whitelist (dist 1)
bc1 = filter_within_dist(b[6], "configs/seqproc/splitseq_bc1_whitelist_6bp.txt", 1)
rest = r:

1{<read1>}
2{<skip2><umi><bc3><skip_l1><bc2><skip_l2><bc1><rest>}

-> 1{<read1>} 2{<umi><bc3><bc2><bc1>}
