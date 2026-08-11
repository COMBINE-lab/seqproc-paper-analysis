# LR-SPLiT-seq, dual orientation, with complete fixed-length components.
#
# The leading unbounded interval absorbs sequence before the embedded barcode
# structure.  This makes UMI, BC3, and BC2 true fixed-length intervals relative
# to their following anchors; reads truncated within one of those components
# are rejected instead of emitting a sub-32-nt projection.

# Keep the prefix mapped during matching; the final projection omits it. Using
# an unbounded discard interval would remove the sequence before the anchor-
# relative fixed fields have been sliced.
prefix = r:
umi = u[10]
bc3 = b[8]
#[search(relative)] #[edit(6)] linker_a = f[GTGGCCGATGTTTCGCATCGGCGTACGACT]
# As with the leading prefix, retain any sequence between the anchors only long
# enough to slice the complete 8-nt BC2 immediately preceding linker B.
between_linkers = r:
bc2 = b[8]
#[search(relative)] #[edit(4)] linker_b = f[ATCCACGTGCTTGAGACTGTGG]
bc1 = b[6]
rest = r:

#[match_ori(either)]
1{<prefix><umi><bc3><linker_a><between_linkers><bc2><linker_b><bc1><rest>}

-> 1{<umi><bc3><bc2><bc1>}
