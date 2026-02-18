#!/usr/bin/env python3
"""
Calculate Jaccard Index for FULL dataset (Seqproc vs Matchbox) - COMBINED FWD + RC.
"""

from pathlib import Path

# Paths to full dataset outputs
SEQPROC_FWD_FQ = "/home/ubuntu/seqproc_full_fwd.fq"
SEQPROC_RC_FQ = "/home/ubuntu/seqproc_full_rc.fq"
MATCHBOX_FWD_TSV = "/home/ubuntu/matchbox_full.tsv"
MATCHBOX_RC_TSV = "/home/ubuntu/combine-lab/seqproc-paper-analysis-clean/matchbox_full_rc.tsv"

def parse_seqproc_ids(fq_path, label):
    print(f"Parsing Seqproc {label} output: {fq_path} ...")
    ids = set()
    count = 0
    try:
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
    except FileNotFoundError:
        print(f"Error: File not found {fq_path}")
    print(f"Total Seqproc {label} valid reads: {len(ids):,}")
    return ids

def parse_matchbox_ids(tsv_path, label):
    print(f"Parsing Matchbox {label} output: {tsv_path} ...")
    ids = set()
    count = 0
    try:
        with open(tsv_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 1:
                    read_id = parts[0]
                    ids.add(read_id)
                    count += 1
                    if count % 1000000 == 0:
                        print(f"  Parsed {count} reads...")
    except FileNotFoundError:
        print(f"Error: File not found {tsv_path}")
    print(f"Total Matchbox {label} valid reads: {len(ids):,}")
    return ids

def main():
    # Parse Seqproc
    sp_fwd = parse_seqproc_ids(SEQPROC_FWD_FQ, "FWD")
    sp_rc = parse_seqproc_ids(SEQPROC_RC_FQ, "RC")
    sp_union = sp_fwd.union(sp_rc)
    print(f"Seqproc Combined (Union) Valid Reads: {len(sp_union):,}")

    # Parse Matchbox
    mb_fwd = parse_matchbox_ids(MATCHBOX_FWD_TSV, "FWD")
    mb_rc = parse_matchbox_ids(MATCHBOX_RC_TSV, "RC")
    mb_union = mb_fwd.union(mb_rc)
    print(f"Matchbox Combined (Union) Valid Reads: {len(mb_union):,}")

    # Calculate Jaccard
    intersection = sp_union.intersection(mb_union)
    total_union = sp_union.union(mb_union)
    
    jaccard = len(intersection) / len(total_union) if total_union else 0
    
    print("\n" + "="*50)
    print("RESULTS (Seqproc vs Matchbox - Full Dataset COMBINED)")
    print("="*50)
    print(f"Seqproc Total Unique Valid:  {len(sp_union):,}")
    print(f"Matchbox Total Unique Valid: {len(mb_union):,}")
    print(f"Intersection:                {len(intersection):,}")
    print(f"Union of both tools:         {len(total_union):,}")
    print(f"Combined Jaccard Index:      {jaccard:.4f}")
    print("="*50)

if __name__ == "__main__":
    main()
