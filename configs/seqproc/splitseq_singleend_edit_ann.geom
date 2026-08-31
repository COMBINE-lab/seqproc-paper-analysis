# SPLiT-seq geometry for SINGLE-END long-read data (SRR13948564) using EDIT DISTANCE
# Structure: [UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:22][BC1:6][rest]
#
# Orientation-aware: #[match_ori(either)] tries forward, then RC on failure
# Uses anchor_relative to find linkers with fuzzy matching (edit distance)

umi = u[10]
bc3 = b[8]
#[search(relative)] #[edit(6)] linker_a = f[GTGGCCGATGTTTCGCATCGGCGTACGACT]
bc2 = b[8]
#[search(relative)] #[edit(4)] linker_b = f[ATCCACGTGCTTGAGACTGTGG]
bc1 = b[6]
rest = r:

# Single read with barcode structure (try both orientations)
#[match_ori(either)]
1{<umi><bc3><linker_a><bc2><linker_b><bc1><rest>}

# Output: extracted UMI+barcodes
-> 1{<umi><bc3><bc2><bc1>}
