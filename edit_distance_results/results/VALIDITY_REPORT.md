# SPLiT-seq Read Validity Analysis Report

## Executive Summary
We investigated the discrepancy in read recovery between `seqproc` and `splitcode` using the SPLiT-seq dataset (SRR6750041).
After identifying and correcting mismatches in the expected read structure (UMI length and Linker2 length), we found that both tools recover a nearly identical number of high-quality reads, but differ significantly in their handling of low-quality data.

## Key Findings

1. **Structural Correction**:
   - The actual data structure differs from the initial configuration.
   - **UMI**: 8bp (was configured as 10bp)
   - **Linker2**: 30bp (was configured as 16bp)
   - **BC3 position**: Immediately follows UMI (at 8bp), not at 10bp.

2. **Read Recovery**:
   - **seqproc output**: 741,006 reads
   - **splitcode output**: 911,217 reads
   - **Difference**: 170,211 reads

3. **Validity Analysis**:
   - **High-Quality Yield (all barcodes within distance 1 of whitelist)**:
     - **seqproc**: 672,406 reads (90.7% of output)
     - **splitcode**: 681,982 reads (74.8% of output)
     - The effective yield is nearly identical (difference < 1.5%).

4. **Noise Profile**:
   - **splitcode** outputs ~170k additional reads that generally fail validation:
     - Only ~5.6% of these extra reads have valid barcodes (d≤1).
     - Most fail structural checks (missing linkers) or have heavily mutated barcodes.
   - **seqproc** acts as a stricter filter, outputting fewer total reads but maintaining much higher purity (>90% valid vs ~75% valid).

## Conclusion
The apparent "lower sensitivity" of `seqproc` is actually **higher precision**. It correctly filters out reads that lack valid SPLiT-seq structure or contain uncorrectable barcodes. When filtered for validity (which is required for downstream single-cell analysis), both tools provide equivalent sensitivity.
