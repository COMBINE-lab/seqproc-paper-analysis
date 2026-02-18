import sys
import os
import gzip
import argparse

try:
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2
except ImportError:
    plt = None
    venn2 = None

def open_file(filepath):
    # Check for gzip signature if extension is ambiguous, or just try/except
    # But simple extension check is usually enough
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')

def load_seqproc_reads(r2_fastq):
    """
    Seqproc output R2 contains valid reads.
    We just need the Read IDs to compare validity.
    """
    print(f"Loading seqproc results from {r2_fastq}...")
    valid_ids = set()
    
    if not os.path.exists(r2_fastq):
        print(f"Error: {r2_fastq} not found")
        return valid_ids

    with open_file(r2_fastq) as f:
        while True:
            header = f.readline().strip()
            if not header: break
            f.readline() # seq
            f.readline() # +
            f.readline() # qual
            
            # Parse ID: @SRR6750041.1 1 length=66 -> SRR6750041.1
            if not header.startswith('@'):
                continue
                
            read_id = header.split()[0].replace('@', '')
            valid_ids.add(read_id)
                
    print(f"  Found {len(valid_ids)} valid seqproc reads")
    return valid_ids

def load_splitpipe_reads(fastq_file):
    """
    Load valid Read IDs from split-pipe's barcode_head.fastq.
    """
    if fastq_file is None:
        fastq_file = "/home/ubuntu/splitpipe_results/process/barcode_head.fastq"
        
    print(f"Loading split-pipe results from {fastq_file}...")
    valid_ids = set()
    
    if not os.path.exists(fastq_file):
        print(f"Error: {fastq_file} not found")
        return valid_ids
        
    with open_file(fastq_file) as f:
        while True:
            header = f.readline().strip()
            if not header: break
            f.readline()
            f.readline()
            f.readline()
            
            # Header likely contains the original ID
            # split-pipe modifies header but keeps original after _OH_
            # Example: @04...__OH_@SRR6750041.71 1 length=66
            if "_OH_@" in header:
                read_id = header.split("_OH_@")[1].split()[0]
            else:
                # Fallback or standard format
                read_id = header.split()[0].replace('@', '')
            
            valid_ids.add(read_id)
            
    print(f"  Found {len(valid_ids)} valid split-pipe reads")
    return valid_ids

def main():
    parser = argparse.ArgumentParser(description="Compare Seqproc and Split-pipe results")
    parser.add_argument("-s", "--seqproc", required=True, help="Seqproc R2 FASTQ output")
    parser.add_argument("-p", "--splitpipe", required=True, help="Split-pipe barcode_head.fastq output")
    parser.add_argument("-n", "--total_reads", type=int, default=1000000, help="Total input reads (default: 1,000,000)")
    args = parser.parse_args()

    seqproc_file = args.seqproc
    splitpipe_file = args.splitpipe
    total_reads = args.total_reads
    
    s_ids = load_seqproc_reads(seqproc_file)
    p_ids = load_splitpipe_reads(splitpipe_file)
    
    intersection = s_ids & p_ids
    union = s_ids | p_ids
    
    print("\n--- Concordance Analysis (Valid Read Recovery) ---")
    print(f"Total Input Reads: {total_reads:,} (Subsampled)")
    print(f"Seqproc Recovered:   {len(s_ids):,} ({len(s_ids)/total_reads:.1%})")
    print(f"Split-pipe Recovered:{len(p_ids):,} ({len(p_ids)/total_reads:.1%})")
    print(f"Intersection:        {len(intersection):,}")
    if len(union) > 0:
        print(f"Jaccard Index:       {len(intersection)/len(union):.4f}")
    
    tp = len(intersection)
    fp = len(s_ids - p_ids)
    fn = len(p_ids - s_ids)
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    print(f"\n--- Metrics (Ref = Split-pipe) ---")
    print(f"Precision (Seqproc vs Ref): {precision:.4f}")
    print(f"Recall (Seqproc vs Ref):    {recall:.4f}")
    print(f"F1 Score:                   {f1:.4f}")
    
    print(f"\n--- Discrepancy Analysis ---")
    print(f"Reads in Seqproc ONLY:   {fp:,}")
    print(f"Reads in Split-pipe ONLY:{fn:,}", flush=True)
    
    if plt and venn2:
        plt.figure(figsize=(8, 8))
        venn2([s_ids, p_ids], ('Seqproc', 'Split-pipe'))
        plt.title("Valid Read Recovery Overlap")
        plt.savefig("splitseq_concordance_venn.png")
        print("\nPlot saved to splitseq_concordance_venn.png", flush=True)
    else:
        print("\nSkipping Venn diagram (matplotlib or matplotlib_venn not installed)", flush=True)

if __name__ == "__main__":
    main()
