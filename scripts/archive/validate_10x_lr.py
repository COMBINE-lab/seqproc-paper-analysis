import sys
import os
import gzip
from collections import Counter

def load_whitelist(whitelist_path):
    print(f"Loading whitelist from {whitelist_path}...")
    whitelist = set()
    with open(whitelist_path, 'r') as f:
        for line in f:
            # Assume one barcode per line, possibly CSV
            bc = line.strip().split(',')[0]
            if bc:
                whitelist.add(bc)
    print(f"Loaded {len(whitelist)} whitelisted barcodes.")
    return whitelist

def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        return 99
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def validate_barcodes(fastq_path, whitelist, max_dist=1):
    print(f"Validating barcodes in {fastq_path}...")
    total_reads = 0
    valid_reads = 0
    
    # Check if file exists
    if not os.path.exists(fastq_path):
        print(f"Error: File {fastq_path} not found.")
        return 0, 0

    try:
        if fastq_path.endswith('.gz'):
            f = gzip.open(fastq_path, 'rt')
        else:
            f = open(fastq_path, 'r')
            
        with f:
            while True:
                header = f.readline()
                if not header: break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                
                total_reads += 1
                
                # Extract Barcode (Assuming R1 structure or Seqproc output structure)
                # For 10x v2, BC is 16bp. 
                # Seqproc output for 10x usually outputs: [BC][UMI]...
                # Let's assume the first 16bp is the barcode for validation
                if len(seq) >= 16:
                    bc = seq[:16]
                    
                    if bc in whitelist:
                        valid_reads += 1
                    else:
                        # Re-implementing efficient check:
                        # Generate all 1-hamming neighbors of the observed barcode
                        # Check if any neighbor is in whitelist
                        # Alphabet: A, C, G, T, N
                        found = False
                        alphabet = ['A', 'C', 'G', 'T', 'N']
                        for i in range(len(bc)):
                            orig_char = bc[i]
                            for char in alphabet:
                                if char == orig_char: continue
                                neighbor = bc[:i] + char + bc[i+1:]
                                if neighbor in whitelist:
                                    valid_reads += 1
                                    found = True
                                    break
                            if found: break
                            
                if total_reads % 10000 == 0:
                    print(f"  Processed {total_reads} reads...", end='\r')
                    
    except Exception as e:
        print(f"Error processing {fastq_path}: {e}")
        
    print(f"  Processed {total_reads} reads.")
    return total_reads, valid_reads

def main():
    # Configuration
    # Whitelist path (10x v2)
    whitelist_path = "/home/ubuntu/combine-lab/seqproc_paper_analysis/results/processed/10x_barcode_whitelist.csv"
    
    # Datasets to validate (10x LR)
    # Output files from seqproc/matchbox/splitcode for 10x LR dataset
    # We need to identify where these are located.
    # Based on previous runs, they might be in results/paper_figures/10x_promethion/ or similar?
    # Or in the data/ folder if they were just extracted?
    # For now, let's allow passing files as arguments
    
    if len(sys.argv) < 2:
        print("Usage: python validate_10x_lr.py <output_fastq_1> <output_fastq_2> ...")
        sys.exit(1)
        
    whitelist = load_whitelist(whitelist_path)
    
    print(f"{'File':<50} | {'Total':<10} | {'Valid':<10} | {'% Valid':<10}")
    print("-" * 90)
    
    for fastq_file in sys.argv[1:]:
        total, valid = validate_barcodes(fastq_file, whitelist)
        pct = (valid / total * 100) if total > 0 else 0
        print(f"{os.path.basename(fastq_file):<50} | {total:<10} | {valid:<10} | {pct:<10.2f}")

if __name__ == "__main__":
    main()
