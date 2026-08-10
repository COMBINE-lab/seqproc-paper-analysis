# Cross-tool forward-orientation LR-SPLiT-seq benchmark geometry. Thresholds
# match the practical splitcode configuration (edit distance 3 per linker).
umi = u[10]
bc3 = b[8]
#[search(relative)] #[edit(3)] linker_a = f[GTGGCCGATGTTTCGCATCGGCGTACGACT]
bc2 = b[8]
#[search(relative)] #[edit(3)] linker_b = f[ATCCACGTGCTTGAGACTGTGG]
bc1 = b[6]
rest = r:

1{<umi><bc3><linker_a><bc2><linker_b><bc1><rest>}
-> 1{<umi><bc3><bc2><bc1>}
