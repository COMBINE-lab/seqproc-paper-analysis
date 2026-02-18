#!/usr/bin/env python3
import os
import subprocess
import sys
import gzip
from collections import Counter

# Configuration
PROJECT_ROOT = "/home/ubuntu/combine-lab/seqproc-paper-analysis-clean"
SEQPROC_BIN = f"{PROJECT_ROOT}/../seqproc/target/release/seqproc"
MATCHBOX_BIN = f"{PROJECT_ROOT}/../matchbox/target/release/matchbox"
SPLITCODE_BIN = f"{PROJECT_ROOT}/../splitcode/build/src/splitcode"

DATASETS = {
    "splitseq_pe": {
        "r1": f"{PROJECT_ROOT}/data/SRR6750041_1M_R1.fastq",
        "r2": f"{PROJECT_ROOT}/data/SRR6750041_1M_R2.fastq",
        "mode": "paired",
        "seqproc_geom": f"{PROJECT_ROOT}/configs/seqproc/splitseq_filter.geom",
        "seqproc_maps": [
            f"{PROJECT_ROOT}/configs/seqproc/splitseq_bc3_seq2seq.tsv",
            f"{PROJECT_ROOT}/configs/seqproc/splitseq_bc2_seq2seq.tsv",
            f"{PROJECT_ROOT}/configs/seqproc/splitseq_bc1_seq2seq.tsv"
        ],
        "matchbox_config": f"{PROJECT_ROOT}/configs/matchbox/splitseq_replacement.mb",
        "splitcode_config": f"{PROJECT_ROOT}/configs/splitcode/splitseq_paper.config"
    },
    "splitseq_se": {
        "r1": f"{PROJECT_ROOT}/data/SRR13948564_1M.fastq",
        "mode": "single",
        "seqproc_geom": f"{PROJECT_ROOT}/configs/seqproc/splitseq_singleend_primer.geom",
        "seqproc_maps": [
            f"{PROJECT_ROOT}/configs/seqproc/splitseq_bc3_seq2seq.tsv",
            f"{PROJECT_ROOT}/configs/seqproc/splitseq_bc2_seq2seq.tsv",
            f"{PROJECT_ROOT}/configs/seqproc/splitseq_bc1_seq2seq.tsv"
        ],
        "matchbox_config": f"{PROJECT_ROOT}/configs/matchbox/splitseq_singleend.mb",
        "splitcode_config": f"{PROJECT_ROOT}/configs/splitcode/splitseq_singleend.config"
    },
    "10x_promethion": {
        "r1": f"{PROJECT_ROOT}/data/10x/ERR9958135_1M.fastq",
        "mode": "single",
        "seqproc_geom": f"{PROJECT_ROOT}/configs/seqproc/10x_longread_fwd.geom",
        "seqproc_geom_rev": f"{PROJECT_ROOT}/configs/seqproc/10x_longread_rev.geom",
        "matchbox_config": f"{PROJECT_ROOT}/configs/matchbox/10x_longread.mb",
        "splitcode_config": f"{PROJECT_ROOT}/configs/splitcode/10x_longread.config",
        "whitelist": "/home/ubuntu/3M-february-2018.txt.gz"
    }
}

def run_command(cmd, cwd=None):
    print(f"Running: {cmd}")
    # Don't suppress stderr so we can see errors
    subprocess.run(cmd, shell=True, check=True, cwd=cwd, stdout=subprocess.DEVNULL)

def get_read_ids(fastq_path):
    ids = set()
    if not os.path.exists(fastq_path):
        print(f"Warning: {fastq_path} not found.")
        return ids
    
    opener = gzip.open if fastq_path.endswith('.gz') else open
    with opener(fastq_path, 'rt') as f:
        while True:
            header = f.readline()
            if not header: break
            f.readline() # seq
            f.readline() # plus
            f.readline() # qual
            ids.add(header.split()[0].strip().replace('@', ''))
    return ids

def parse_matchbox_tsv(tsv_path):
    ids = set()
    barcodes = {}
    if not os.path.exists(tsv_path):
        print(f"Warning: {tsv_path} not found.")
        return ids, barcodes
        
    with open(tsv_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 1:
                rid = parts[0]
                ids.add(rid)
                if len(parts) >= 2:
                    barcodes[rid] = parts[1]
    return ids, barcodes

def analyze_splitseq(dataset_name, config):
    print(f"\n{'='*60}")
    print(f"Analyzing {dataset_name}...")
    print(f"{'='*60}")
    
    # Clean previous results
    os.system(f"rm -f {dataset_name}_seqproc* {dataset_name}_matchbox* {dataset_name}_splitcode*")
    
    # Seqproc
    print("  Running Seqproc...")
    maps = " ".join([f"-a {m}" for m in config.get('seqproc_maps', [])])
    if config['mode'] == 'paired':
        cmd = f"{SEQPROC_BIN} --geom {config['seqproc_geom']} {maps} --file1 {config['r1']} --file2 {config['r2']} --out1 {dataset_name}_seqproc_R1.fq --out2 {dataset_name}_seqproc_R2.fq --threads 8"
    else:
        cmd = f"{SEQPROC_BIN} --geom {config['seqproc_geom']} {maps} --file1 {config['r1']} --out1 {dataset_name}_seqproc_R1.fq --threads 8"
    run_command(cmd)
    
    # Matchbox
    print("  Running Matchbox...")
    mb_fq_mode = False
    
    # Determine if Matchbox config outputs FASTQ or TSV (stdout)
    # splitseq_pe (replacement) -> FASTQ
    # splitseq_se (singleend) -> TSV
    
    if dataset_name == "splitseq_pe":
        # FASTQ output mode (defined in .mb config)
        mb_fq_mode = True
        cmd = f"{MATCHBOX_BIN} -e 0.2 -t 8 -s {config['matchbox_config']} {config['r1']} -p {config['r2']} > /dev/null"
        run_command(cmd, cwd=PROJECT_ROOT)
        
        # Move generated FASTQs
        os.system(f"mv {PROJECT_ROOT}/mb_r*.fq . 2>/dev/null")
        if os.path.exists("mb_r2.fq"): os.rename("mb_r2.fq", f"{dataset_name}_matchbox_R2.fq")
        if os.path.exists("mb_r1.fq"): os.rename("mb_r1.fq", f"{dataset_name}_matchbox_R1.fq")
        
    else:
        # TSV output mode (stdout)
        mb_fq_mode = False
        tsv_out = f"{dataset_name}_matchbox.tsv"
        # We must use the absolute path for tsv_out because we run with cwd=PROJECT_ROOT
        # But we want the file in the current directory.
        abs_tsv_out = os.path.abspath(tsv_out)
        
        if config['mode'] == 'paired':
             cmd = f"{MATCHBOX_BIN} -e 0.2 -t 8 -s {config['matchbox_config']} {config['r1']} -p {config['r2']} > {abs_tsv_out}"
        else:
             cmd = f"{MATCHBOX_BIN} -e 0.2 -t 8 -s {config['matchbox_config']} {config['r1']} > {abs_tsv_out}"
        
        run_command(cmd, cwd=PROJECT_ROOT)

    # Splitcode (only if config exists)
    if 'splitcode_config' in config:
        print("  Running Splitcode...")
        # Need mapping file
        with open(f"{dataset_name}_mapping.txt", 'w') as f: f.write("")
        if config['mode'] == 'paired':
            cmd = f"{SPLITCODE_BIN} -c {config['splitcode_config']} --assign -N 2 -t 8 -m {dataset_name}_mapping.txt -o {dataset_name}_splitcode_R1.fq,{dataset_name}_splitcode_R2.fq {config['r1']} {config['r2']}"
        else:
            cmd = f"{SPLITCODE_BIN} -c {config['splitcode_config']} --assign -N 1 -t 8 -m {dataset_name}_mapping.txt -o {dataset_name}_splitcode_R1.fq {config['r1']}"
        run_command(cmd)

    # Calculate Jaccard
    print("  Calculating concordance...")
    
    target_read = "R2" if config['mode'] == 'paired' else "R1"
    
    sp_ids = get_read_ids(f"{dataset_name}_seqproc_{target_read}.fq")
    
    if mb_fq_mode:
        mb_ids = get_read_ids(f"{dataset_name}_matchbox_{target_read}.fq")
    else:
        mb_ids, _ = parse_matchbox_tsv(f"{dataset_name}_matchbox.tsv")
        
    sc_ids = set()
    if 'splitcode_config' in config:
        if config['mode'] == 'paired':
             sc_ids = get_read_ids(f"{dataset_name}_splitcode_{target_read}.fq")
        else:
             sc_ids = get_read_ids(f"{dataset_name}_splitcode_{target_read}.fq")
        
    print(f"    Seqproc: {len(sp_ids)}")
    print(f"    Matchbox: {len(mb_ids)}")
    if sc_ids: print(f"    Splitcode: {len(sc_ids)}")
    
    # Jaccard Seqproc vs Matchbox
    intersect_sm = len(sp_ids.intersection(mb_ids))
    union_sm = len(sp_ids.union(mb_ids))
    jaccard_sm = intersect_sm / union_sm if union_sm > 0 else 0
    print(f"    Jaccard (Seqproc vs Matchbox): {jaccard_sm:.4f} (Intersection: {intersect_sm})")
    
    if sc_ids:
        # Jaccard Seqproc vs Splitcode
        intersect_ss = len(sp_ids.intersection(sc_ids))
        union_ss = len(sp_ids.union(sc_ids))
        jaccard_ss = intersect_ss / union_ss if union_ss > 0 else 0
        print(f"    Jaccard (Seqproc vs Splitcode): {jaccard_ss:.4f} (Intersection: {intersect_ss})")

def analyze_10x(dataset_name, config):
    print(f"\n{'='*60}")
    print(f"Analyzing {dataset_name}...")
    print(f"{'='*60}")
    
    # Clean
    os.system(f"rm -f {dataset_name}_seqproc* {dataset_name}_matchbox*")
    
    # Seqproc (Dual pass)
    print("  Running Seqproc (Dual Pass)...")
    cmd_fwd = f"{SEQPROC_BIN} --geom {config['seqproc_geom']} --file1 {config['r1']} --out1 {dataset_name}_seqproc_fwd.fq --threads 8"
    run_command(cmd_fwd)
    cmd_rev = f"{SEQPROC_BIN} --geom {config['seqproc_geom_rev']} --file1 {config['r1']} --out1 {dataset_name}_seqproc_rev.fq --threads 8"
    run_command(cmd_rev)
    # Merge ids
    sp_fwd_ids = get_read_ids(f"{dataset_name}_seqproc_fwd.fq")
    sp_rev_ids = get_read_ids(f"{dataset_name}_seqproc_rev.fq")
    sp_ids = sp_fwd_ids.union(sp_rev_ids)
    
    # Matchbox
    print("  Running Matchbox...")
    # This config prints to stdout, capture to TSV
    # Use absolute path for output to avoid CWD confusion
    abs_tsv = os.path.abspath(f"{dataset_name}_matchbox.tsv")
    cmd = f"{MATCHBOX_BIN} -e 0.2 -t 8 -s {config['matchbox_config']} {config['r1']} > {abs_tsv}"
    run_command(cmd, cwd=PROJECT_ROOT)
    
    # Extract IDs from TSV
    mb_ids, mb_barcodes = parse_matchbox_tsv(f"{dataset_name}_matchbox.tsv")
    
    print(f"    Seqproc: {len(sp_ids)}")
    print(f"    Matchbox: {len(mb_ids)}")
    
    # Jaccard
    intersect = len(sp_ids.intersection(mb_ids))
    union = len(sp_ids.union(mb_ids))
    jaccard = intersect / union if union > 0 else 0
    print(f"    Jaccard (Seqproc vs Matchbox): {jaccard:.4f}")
    
    # Discrepancy
    mb_only = mb_ids - sp_ids
    print(f"    Matchbox ONLY: {len(mb_only)}")
    
    if len(mb_only) > 0:
        print("    Investigating Matchbox-only reads...")
        # Load Whitelist
        print("      Loading whitelist...")
        whitelist = set()
        with gzip.open(config['whitelist'], 'rt') as f:
            for line in f:
                whitelist.add(line.strip().split()[0])
        print(f"      Loaded {len(whitelist)} barcodes.")
        
        # Check MB-only reads
        valid_mb_only = 0
        total_checked = 0
        
        for rid in mb_only:
            bc = mb_barcodes.get(rid)
            if bc:
                total_checked += 1
                if bc in whitelist:
                    valid_mb_only += 1
        
        print(f"      Analyzed {total_checked} Matchbox-only reads.")
        print(f"      Valid Barcodes (Exact Whitelist Match): {valid_mb_only} ({valid_mb_only/total_checked*100:.1f}%)")
        
        if valid_mb_only / total_checked > 0.5:
            print("      CONCLUSION: Matchbox is finding valid reads that Seqproc misses.")
        else:
            print("      CONCLUSION: Matchbox is outputting invalid barcodes (junk) or reversed sequences that Seqproc rejects.")

# Main
if __name__ == "__main__":
    analyze_splitseq("splitseq_pe", DATASETS["splitseq_pe"])
    analyze_splitseq("splitseq_se", DATASETS["splitseq_se"])
    analyze_10x("10x_promethion", DATASETS["10x_promethion"])
