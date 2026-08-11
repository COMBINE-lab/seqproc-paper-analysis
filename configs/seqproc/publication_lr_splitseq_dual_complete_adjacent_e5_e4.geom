# LR-SPLiT-seq, dual orientation, complete components and strict cassette
# adjacency. Linker A may occur anywhere in the long read. Once found, BC2
# must occupy the next eight bases and linker B must begin immediately after
# BC2 (subject to edit-distance alignment).

prefix = r:
umi = u[10]
bc3 = b[8]
#[search(relative)] #[edit(5)] linker_a = f[GTGGCCGATGTTTCGCATCGGCGTACGACT]
bc2 = b[8]
#[edit(4)] linker_b = f[ATCCACGTGCTTGAGACTGTGG]
bc1 = b[6]
rest = r:

#[match_ori(either)]
1{<prefix><umi><bc3><linker_a><bc2><linker_b><bc1><rest>}

-> 1{<umi><bc3><bc2><bc1>}
