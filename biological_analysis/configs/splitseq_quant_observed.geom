# SPLiT-seq quantification geometry -- OBSERVED barcodes (no canonicalization).
# Recovery is via filter_within_dist (edit<=1, same as the benchmark), but the EMITTED
# barcode is the observed read sequence, not the whitelist sequence. STARsolo does all
# barcode correction downstream, identically for both tools. Symmetric with the splitcode
# -x extraction config. Layout: UMI 10 + bc3 8 + bc2 8 + bc1 8 = 34bp.
# All three rounds share the same 8bp whitelist (splitseq_bc23_whitelist.txt).
# Linkers are matched within edit distance 3 with indels allowed, the common budget used for
# all three tools in the downstream comparison. seqproc can raise this to 6 at negligible cost
# (Myers bit-vector), splitcode cannot (its cost grows near-exponentially), so 3 is the fair
# shared setting; the residual recovery gap then reflects indel handling alone.

read1 = r:
umi = u[10]
bc3 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc23_whitelist.txt", 1)
#[search(relative)] #[edit(3)] l1 = f[GTGGCCGCTGTTTCGCATCGGCGTACGACT]
bc2 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc23_whitelist.txt", 1)
#[search(relative)] #[edit(3)] l2 = f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG]
bc1 = filter_within_dist(b[8], "configs/seqproc/splitseq_bc23_whitelist.txt", 1)
rest = r:

1{<read1>}
2{<umi><bc3><l1><bc2><l2><bc1><rest>}

-> 1{<read1>} 2{<umi><bc3><bc2><bc1>}
