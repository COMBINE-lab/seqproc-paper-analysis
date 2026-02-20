# SPLiT-seq Round 2 geometry for REAL data (SRR6750041) - FIXED for 30bp L2
# Structure on R2: [NN:2][UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:30][BC1:6][rest]

read1 = r:
skip2 = x[2]
umi = u[10]
bc3 = b[8]

# Linker1 - 30bp - STRICT POSITIONAL MATCH
l1 = hamming(f[GTGGCCGCTGTTTCGCATCGGCGTACGACT], 5)

bc2 = b[8]

# Linker2 - 30bp (Full length for SRR6750041) - STRICT POSITIONAL MATCH
# ATCCACGTGCTTGAGAGGCCAGAGCATTCG
l2 = hamming(f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG], 5)

bc1 = b[6]
rest = r:

1{<read1>}
2{<skip2><umi><bc3><l1><bc2><l2><bc1><rest>}

-> 1{<read1>} 2{<umi><bc3><bc2><bc1>}
