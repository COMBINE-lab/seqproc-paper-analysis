#!/usr/bin/env python3
"""
Calculate Jaccard Index for FULL dataset (seqproc vs matchbox).
"""

from collections import defaultdict
from pathlib import Path

# Paths to full dataset outputs
SEQPROC_FQ = "/home/ubuntu/seqproc_full_fwd.fq"
MATCHBOX_TSV = "/home/ubuntu/matchbox_full.tsv"

def parse_seqproc_ids(fq_path):
    print(f"Parsing Seqproc output: {fq_path} ...")
    ids = set()
    count = 0
    with open(fq_path, 'r') as f:
        while True:
            header = f.readline()
            if not header: break
            f.readline() # seq
            f.readline() # +
            f.readline() # qual
            
            # Read ID: @SRR... 123/1 -> SRR...
            read_id = header.strip().split()[0].replace('@', '')
            ids.add(read_id)
            count += 1
            if count % 1000000 == 0:
                print(f"  Parsed {count} reads...")
    print(f"Total Seqproc valid reads: {len(ids):,}")
    return ids

def parse_matchbox_ids(tsv_path):
    print(f"Parsing Matchbox output: {tsv_path} ...")
    ids = set()
    count = 0
    with open(tsv_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 1:
                read_id = parts[0]
                ids.add(read_id)
                count += 1
                if count % 1000000 == 0:
                    print(f"  Parsed {count} reads...")
    print(f"Total Matchbox valid reads: {len(ids):,}")
    return ids

def main():
    sp_ids = parse_seqproc_ids(SEQPROC_FQ)
    mb_ids = parse_matchbox_ids(MATCHBOX_TSV)
    
    intersection = sp_ids.intersection(mb_ids)
    union = sp_ids.union(mb_ids)
    
    jaccard = len(intersection) / len(union) if union else 0
    
    print("\n" + "="*40)
    print("RESULTS (Seqproc vs Matchbox - Full Dataset FWD)")
    print("="*40)
    print(f"Seqproc Count:  {len(sp_ids):,}")
    print(f"Matchbox Count: {len(mb_ids):,}")
    print(f"Intersection:   {len(intersection):,}")
    print(f"Union:          {len(union):,}")
    print(f"Jaccard Index:  {jaccard:.4f}")
    print("="*40)

if __name__ == "__main__":
    main()
