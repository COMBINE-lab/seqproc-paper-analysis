# SPLiT-seq single-end REVERSE COMPLEMENT geometry
# Original: [skip_start][UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:22][BC1:8][rest]
# RC:       [rest_rc][BC1_RC:8][Linker2_RC:22][BC2_RC:8][Linker1_RC:30][BC3_RC:8][UMI_RC:10][skip_end]
#
# Both linkers use matched 0.2 error rate tolerances:
#   L1_RC (30bp): hamming(6)  -- same as forward config
#   L2_RC (22bp): hamming(4)  -- same as forward config

# RC Linker 1: RC of GTGGCCGATGTTTCGCATCGGCGTACGACT
l1_rc = anchor_relative(hamming(f[AGTCGTACGCCGATGCGAAACATCGGCCAC], 6))
# RC Linker 2: RC of ATCCACGTGCTTGAGACTGTGG
l2_rc = anchor_relative(hamming(f[CCACAGTCTCAAGCACGTGGAT], 4))

# Elements
genomic = r:
bc1_rc = b[8]
# l2_rc is defined above
bc2_rc = b[8]
# l1_rc is defined above
bc3_rc = b[8]
umi_rc = u[10]
skip_end = r:

# Structure
1{<genomic><bc1_rc><l2_rc><bc2_rc><l1_rc><bc3_rc><umi_rc><skip_end>}

# Output (keep RC sequences for ID matching)
-> 1{<umi_rc><bc3_rc><bc2_rc><bc1_rc>}
