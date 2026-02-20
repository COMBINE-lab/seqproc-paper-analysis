# SPLiT-seq with FAST barcode replacement using pre-computed mismatch variants
# Uses exact hash lookup instead of hamming distance matching
# For REAL data (SRR6750041) - GCT variant linkers
# Structure on R2: [NN:2][UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:30][BC1:6][rest]

read1 = r:
skip2 = x[2]
umi = u[10]

# BC3 with exact lookup (mismatch variants pre-computed in expanded TSV)
bc3_def = b[8]

# Linker1 - 30bp with hamming distance 6 (~20% error rate)
l1 = anchor_relative(hamming(f[GTGGCCGCTGTTTCGCATCGGCGTACGACT], 6))

# BC2 with exact lookup
bc2_def = b[8]

# Linker2 - 30bp with hamming distance 3
l2 = anchor_relative(hamming(f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG], 3))

# BC1 with exact lookup (6bp truncated in 94bp read with 10bp UMI)
bc1_def = b[6]
rest = r:

1{<read1>}
2{
    <skip2>
    <umi>
    map(<bc3_def>, $0, self)
    <l1>
    map(<bc2_def>, $1, self)
    <l2>
    map(<bc1_def>, $2, self)
    <rest>
}

-> 1{<read1>} 2{<umi><bc3_def><bc2_def><bc1_def>}
