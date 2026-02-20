# SPLiT-seq Round 2 geometry for REAL data (SRR6750041, SRR13948564)
# Structure on R2: [NN:2][UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:30][BC1:6][rest]
#
# Real data has linker variants - using higher hamming distances for tolerance
# Skip first 2 bases to match matchbox _:|2| pattern

read1 = r:
skip2 = x[2]
umi = u[10]
bc3 = b[8]

# Linker1 - 30bp with hamming distance 9 (30% tolerance for variants)
l1 = anchor_relative(hamming(f[GTGGCCGCTGTTTCGCATCGGCGTACGACT], 9))

bc2 = b[8]

# Linker2 - 30bp with hamming distance 5
l2 = anchor_relative(hamming(f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG], 5))

bc1 = b[6]
rest = r:

1{<read1>}
2{<skip2><umi><bc3><l1><bc2><l2><bc1><rest>}

-> 1{<read1>} 2{<umi><bc3><bc2><bc1>}
