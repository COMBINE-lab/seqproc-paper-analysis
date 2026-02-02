# SPLiT-seq single-end with PRIMER ANCHORING using EDIT DISTANCE
# Structure: [skip_start][UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:22][BC1:8][rest]

# Anchor on Linker1 with edit distance 5
linker1 = anchor_relative(edit(f[GTGGCCGATGTTTCGCATCGGCGTACGACT], 5))

# Elements
skip_start = r:
umi = u[10]
bc3 = b[8]
bc2 = b[8]
linker2 = x[22]
bc1 = b[8]
rest = r:

# Structure
1{<skip_start><umi><bc3><linker1><bc2><linker2><bc1><rest>}

# Output
-> 1{<umi><bc3><bc2><bc1>}
