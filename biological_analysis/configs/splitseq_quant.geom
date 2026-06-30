# SPLiT-seq quantification geometry (downstream/biological analysis)
# Corrected to match the empirically decoded ground-truth R2 structure and to be
# SYMMETRIC with the splitcode @extract config (UMI 10 + bc3 8 + bc2 8 + bc1 8).
#
# Ground-truth R2 (94bp), verified by linker-anchored decode of SRR6750041:
#   [UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:30][BC1:8]
#   positions: umi 0:10, bc3 10:18, L1 18:48, bc2 48:56, L2 56:86, bc1 86:94
#
# Differences from configs/seqproc/splitseq_replacement_edit.geom (the benchmark
# recovery geom), and WHY:
#   - NO skip2. The first 2 bases vary across reads (not constant N); they are
#     real UMI bases. The old skip2=x[2] discarded 2 bits of UMI information.
#   - bc1 = b[8], not b[6]. bc1 is a full 8bp barcode at 86:94; the old b[6]
#     truncated it, losing 2bp of cell-barcode resolution.
# Output: 1{cDNA}  2{<umi:10><bc3:8><bc2:8><bc1:8>} = 34bp barcode read.

read1 = r:
umi = u[10]
bc3_def = b[8]
#[search(relative)] #[edit(6)] l1 = f[GTGGCCGCTGTTTCGCATCGGCGTACGACT]
bc2_def = b[8]
#[search(relative)] #[edit(6)] l2 = f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG]
bc1_def = b[8]
rest = r:

1{<read1>}
2{
    <umi>
    map_with_edit(<bc3_def>, $0, self, 1)
    <l1>
    map_with_edit(<bc2_def>, $1, self, 1)
    <l2>
    map_with_edit(<bc1_def>, $2, self, 1)
    <rest>
}

-> 1{<read1>} 2{<umi><bc3_def><bc2_def><bc1_def>}
