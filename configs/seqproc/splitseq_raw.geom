# SPLiT-seq RAW extraction (no barcode replacement) for comparison with matchbox
# For REAL data (SRR6750041) - GCT variant linkers
# CORRECTED Structure on R2: [NN:2][UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:30][BC1:6][rest]
# Note: UMI is 10bp (positions 2-11, per Rosenberg 2018 and STARsolo config)

read1 = r:
skip2 = x[2]
umi = u[10]
bc3 = b[8]

# Linker1 - 30bp with hamming distance 3 (validated)
l1 = anchor_relative(hamming(f[GTGGCCGCTGTTTCGCATCGGCGTACGACT], 3))

bc2 = b[8]

# Linker2 - 30bp with hamming distance 3 (validated)
l2 = anchor_relative(hamming(f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG], 3))

bc1 = b[6]
rest = r:

1{<read1>}
2{<skip2><umi><bc3><l1><bc2><l2><bc1><rest>}

-> 1{<read1>} 2{<umi><bc3><bc2><bc1>}
