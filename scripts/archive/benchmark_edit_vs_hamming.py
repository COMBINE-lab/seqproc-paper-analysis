#!/usr/bin/env python3
import subprocess
import time
import os
from pathlib import Path

# Paths
SEQPROC_BIN = "/home/ubuntu/combine-lab/seqproc/target/release/seqproc"
DATA_DIR = Path("/home/ubuntu/seqproc-paper-analysis-clean/data")
CONFIG_DIR = Path("/home/ubuntu/seqproc-paper-analysis-clean/configs/seqproc")
INPUT_FILE = DATA_DIR / "SRR13948564_1M.fastq"

# output files
HAMMING_OUT = "hamming_out.fq"
EDIT_OUT = "edit_out.fq"

def run_benchmark(tool_name, config_file, output_file):
    print(f"Benchmarking {tool_name}...")
    cmd = [
        SEQPROC_BIN,
        "-g", str(config_file),
        "-1", str(INPUT_FILE),
        "-o", output_file,
        "-t", "4"
    ]
    
    start_time = time.time()
    try:
        subprocess.run(cmd, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error running {tool_name}:")
        print(e.stderr.decode())
        return None, 0

    end_time = time.time()
    runtime = end_time - start_time
    
    # Count output reads
    num_lines = 0
    try:
        # Use wc -l to count lines quickly
        wc_proc = subprocess.run(["wc", "-l", output_file], capture_output=True, text=True)
        num_lines = int(wc_proc.stdout.split()[0])
    except Exception as e:
        print(f"Error counting lines: {e}")
        
    num_reads = num_lines // 4
    
    return runtime, num_reads

def main():
    if not INPUT_FILE.exists():
        print(f"Input file not found: {INPUT_FILE}")
        return

    # 1. Run Hamming Benchmark
    hamming_config = CONFIG_DIR / "splitseq_singleend.geom"
    hamming_time, hamming_reads = run_benchmark("Hamming", hamming_config, HAMMING_OUT)
    
    # 2. Run Edit Benchmark
    edit_config = CONFIG_DIR / "splitseq_singleend_edit.geom"
    edit_time, edit_reads = run_benchmark("Edit", edit_config, EDIT_OUT)
    
    print("\n" + "="*40)
    print(f"{'Metric':<15} | {'Hamming':<10} | {'Edit':<10}")
    print("-" * 40)
    print(f"{'Runtime (s)':<15} | {hamming_time:<10.4f} | {edit_time:<10.4f}")
    print(f"{'Reads Recovered':<15} | {hamming_reads:<10} | {edit_reads:<10}")
    print("-" * 40)
    
    if hamming_time > 0:
        speedup = hamming_time / edit_time if edit_time > 0 else 0
        print(f"Speedup (Hamming/Edit): {speedup:.2f}x")
        
    # Cleanup
    if os.path.exists(HAMMING_OUT):
        os.remove(HAMMING_OUT)
    if os.path.exists(EDIT_OUT):
        os.remove(EDIT_OUT)

if __name__ == "__main__":
    main()
