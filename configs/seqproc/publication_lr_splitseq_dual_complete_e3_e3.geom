# LR-SPLiT-seq, dual orientation, complete components, strict linker sweep.

prefix = r:
umi = u[10]
bc3 = b[8]
#[search(relative)] #[edit(3)] linker_a = f[GTGGCCGATGTTTCGCATCGGCGTACGACT]
between_linkers = r:
bc2 = b[8]
#[search(relative)] #[edit(3)] linker_b = f[ATCCACGTGCTTGAGACTGTGG]
bc1 = b[6]
rest = r:

#[match_ori(either)]
1{<prefix><umi><bc3><linker_a><between_linkers><bc2><linker_b><bc1><rest>}

-> 1{<umi><bc3><bc2><bc1>}
