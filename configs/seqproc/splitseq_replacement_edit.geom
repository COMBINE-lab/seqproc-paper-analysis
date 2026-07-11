# SPLiT-seq with barcode REPLACEMENT using map_with_edit
# For REAL data (SRR6750041) - GCT variant linkers
# CORRECTED Structure on R2: [NN:2][UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:30][BC1:6][rest]

read1 = r:
skip2 = x[2]
umi = u[10]

# BC3 with whitelist replacement (1 edit tolerance)
bc3_def = b[8]

# Linker1 - 30bp with edit distance 6 (~0.2 fraction)
#[search(relative)] #[edit(6)] l1 = f[GTGGCCGCTGTTTCGCATCGGCGTACGACT]

# BC2 with whitelist replacement
bc2_def = b[8]

# Linker2 - 30bp with edit distance 6 (~0.2 fraction)
#[search(relative)] #[edit(6)] l2 = f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG]

# BC1 with whitelist replacement (6bp truncated in 94bp read with 10bp UMI)
bc1_def = b[6]
rest = r:

1{<read1>}
2{
    <skip2>
    <umi>
    map_with_edit(<bc3_def>, $0, self, 1)
    <l1>
    map_with_edit(<bc2_def>, $1, self, 1)
    <l2>
    map_with_edit(<bc1_def>, $2, self, 1)
    <rest>
}

-> 1{<read1>} 2{<umi><bc3_def><bc2_def><bc1_def>}
