#!/usr/bin/env python3
"""
Calculate Intersection of Forward and RC passes for Seqproc and Matchbox.
"""

from pathlib import Path

# Paths to full dataset outputs
SEQPROC_FWD_FQ = "/home/ubuntu/seqproc_full_fwd.fq"
SEQPROC_RC_FQ = "/home/ubuntu/seqproc_full_rc.fq"
MATCHBOX_FWD_TSV = "/home/ubuntu/matchbox_full.tsv"
MATCHBOX_RC_TSV = "/home/ubuntu/combine-lab/seqproc-paper-analysis-clean/matchbox_full_rc.tsv"

def parse_ids(path, format_type):
    ids = set()
    try:
        with open(path, 'r') as f:
            if format_type == 'fastq':
                while True:
                    header = f.readline()
                    if not header: break
                    f.readline(); f.readline(); f.readline()
                    ids.add(header.strip().split()[0].replace('@', ''))
            elif format_type == 'tsv':
                for line in f:
                    parts = line.strip().split('\t')
                    if parts: ids.add(parts[0])
    except FileNotFoundError:
        print(f"File not found: {path}")
    return ids

def analyze_tool(fwd_path, rc_path, format_type, name, total_reads=5764421):
    print(f"\nAnalyzing {name}...")
    fwd = parse_ids(fwd_path, format_type)
    rc = parse_ids(rc_path, format_type)
    
    intersection = fwd.intersection(rc)
    union = fwd.union(rc)
    
    fwd_pct = len(fwd) / total_reads * 100
    rc_pct = len(rc) / total_reads * 100
    union_pct = len(union) / total_reads * 100
    intersect_pct = len(intersection) / total_reads * 100
    overlap_pct_of_union = len(intersection) / len(union) * 100 if union else 0
    
    print(f"  Fwd: {len(fwd):,} ({fwd_pct:.2f}%)")
    print(f"  RC:  {len(rc):,} ({rc_pct:.2f}%)")
    print(f"  Union: {len(union):,} ({union_pct:.2f}%)")
    print(f"  Intersection (Ambiguous): {len(intersection):,} ({intersect_pct:.2f}%)")
    print(f"  Intersection as % of Union: {overlap_pct_of_union:.2f}%")

def main():
    analyze_tool(SEQPROC_FWD_FQ, SEQPROC_RC_FQ, 'fastq', "Seqproc")
    analyze_tool(MATCHBOX_FWD_TSV, MATCHBOX_RC_TSV, 'tsv', "Matchbox")

if __name__ == "__main__":
    main()
