import pysam
import sys
import os
import glob
import gzip
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

def load_seqproc_barcodes(r2_fastq):
    """
    Seqproc R2 output format (from geometry):
    UMI(10) + BC3(8) + BC2(8) + BC1(6)
    Total length: 32bp
    
    Returns: dict {read_id: (bc_seq, umi)}
    """
    print(f"Loading seqproc results from {r2_fastq}...")
    valid_reads = {}
    
    # Check if file exists
    if not os.path.exists(r2_fastq):
        print(f"Error: {r2_fastq} not found")
        return valid_reads

    with open(r2_fastq, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header: break
            seq = f.readline().strip()
            f.readline() # +
            f.readline() # qual
            
            # Parse ID (remove @ and length info)
            read_id = header.split()[0].replace('@', '')
            
            # Parse Barcode (last 22bp = BC3+BC2+BC1) and UMI (first 10bp)
            # Geometry: UMI(10) BC3(8) BC2(8) BC1(6)
            if len(seq) >= 32:
                umi = seq[:10]
                bc = seq[10:32] # BC3+BC2+BC1
                valid_reads[read_id] = (bc, umi)
                
    print(f"  Found {len(valid_reads)} valid seqproc reads")
    return valid_reads

def load_starsolo_barcodes(bam_file):
    """
    Load CB and UB tags from STARsolo BAM
    """
    print(f"Loading STARsolo results from {bam_file}...")
    valid_reads = {}
    
    if not os.path.exists(bam_file):
        print(f"Error: {bam_file} not found")
        return valid_reads
        
    samfile = pysam.AlignmentFile(bam_file, "rb")
    count = 0
    for read in samfile:
        if read.has_tag("CB"):
            cb = read.get_tag("CB")
            ub = read.get_tag("UB")
            # STARsolo CB is usually formatted as BC1_BC2_BC3 or similar?
            # We used whitelist_v1_6bp.txt whitelist_v1.txt whitelist_v1.txt
            # So CB should be 6+8+8 = 22bp?
            # Or does STARsolo output it as a single string?
            # Usually it matches the whitelist.
            
            # Let's clean the read ID
            read_id = read.query_name
            valid_reads[read_id] = (cb, ub)
            count += 1
            
    print(f"  Found {len(valid_reads)} valid STARsolo reads")
    return valid_reads

def load_splitpipe_barcodes(bam_file):
    """
    Load XC (cell bc) and XM (umi) tags from split-pipe BAM
    """
    print(f"Loading split-pipe results from {bam_file}...")
    valid_reads = {}
    
    if not os.path.exists(bam_file):
        print(f"Error: {bam_file} not found")
        return valid_reads
        
    samfile = pysam.AlignmentFile(bam_file, "rb")
    for read in samfile:
        # split-pipe tags: XC = cell barcode, XM = UMI
        if read.has_tag("XC"):
            cb = read.get_tag("XC")
            umi = read.get_tag("XM") if read.has_tag("XM") else ""
            read_id = read.query_name
            valid_reads[read_id] = (cb, umi)
            
    print(f"  Found {len(valid_reads)} valid split-pipe reads")
    return valid_reads

def main():
    # File paths
    seqproc_fq = "seqproc_fixed_R2.fq"
    starsolo_bam = "STARsolo_results/Aligned.out.bam"
    
    # Find split-pipe BAM
    splitpipe_bams = glob.glob("splitpipe_results/**/*.bam", recursive=True)
    splitpipe_bam = splitpipe_bams[0] if splitpipe_bams else None
    
    # Load data
    seqproc_data = load_seqproc_barcodes(seqproc_fq)
    starsolo_data = load_starsolo_barcodes(starsolo_bam)
    
    splitpipe_data = {}
    if splitpipe_bam:
        splitpipe_data = load_splitpipe_barcodes(splitpipe_bam)
    else:
        print("Warning: split-pipe BAM not found yet (maybe running?)")

    # Compare Sets of Read IDs
    s_ids = set(seqproc_data.keys())
    st_ids = set(starsolo_data.keys())
    sp_ids = set(splitpipe_data.keys())
    
    print("\n--- Concordance (Read IDs) ---")
    print(f"Seqproc: {len(s_ids)}")
    print(f"STARsolo: {len(st_ids)}")
    print(f"Split-pipe: {len(sp_ids)}")
    
    common_all = s_ids & st_ids & sp_ids
    print(f"Intersection (All 3): {len(common_all)}")
    
    if sp_ids:
        print(f"Seqproc vs Split-pipe Jaccard: {len(s_ids & sp_ids) / len(s_ids | sp_ids):.4f}")
        print(f"STARsolo vs Split-pipe Jaccard: {len(st_ids & sp_ids) / len(st_ids | sp_ids):.4f}")
    
    print(f"Seqproc vs STARsolo Jaccard: {len(s_ids & st_ids) / len(s_ids | st_ids):.4f}")

    # Venn Diagram
    if sp_ids:
        plt.figure(figsize=(10, 10))
        venn3([s_ids, st_ids, sp_ids], ('Seqproc', 'STARsolo', 'Split-pipe'))
        plt.title("Valid Read ID Overlap")
        plt.savefig("splitseq_concordance_venn.png")
        print("Venn diagram saved to splitseq_concordance_venn.png")

if __name__ == "__main__":
    main()
