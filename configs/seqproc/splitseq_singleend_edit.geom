# SPLiT-seq geometry for SINGLE-END long-read data (SRR13948564) using EDIT DISTANCE
# Structure: [UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:22][BC1:8][rest]
#
# TRUE single-read mode: only one input file, one read definition
# Uses anchor_relative to find linkers with fuzzy matching (edit distance)

umi = u[10]
bc3 = b[8]
linker_a = anchor_relative(edit(f[GTGGCCGATGTTTCGCATCGGCGTACGACT], 9))
bc2 = b[8]
linker_b = anchor_relative(edit(f[ATCCACGTGCTTGAGACTGTGG], 7))
bc1 = b[8]
rest = r:

# Single read with barcode structure
1{<umi><bc3><linker_a><bc2><linker_b><bc1><rest>}

# Output: extracted UMI+barcodes
-> 1{<umi><bc3><bc2><bc1>}
