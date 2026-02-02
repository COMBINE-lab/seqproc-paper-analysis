#!/usr/bin/env python3

import json
import sys

# Load results
with open('/home/ubuntu/seqproc-paper-analysis-clean/results/paper_figures/benchmark_results.json') as f:
    results = json.load(f)

print("=" * 80)
print("EDIT DISTANCE BENCHMARK RESULTS")
print("=" * 80)
print()

# Print summary table
print("{:<25} {:<10} {:<12} {:<10} {:<12} {:<10}".format(
    "Dataset", "Tool", "Runtime (s)", "Mem (MB)", "Reads Out", "Recovery"))
print("-" * 80)

for dataset_key, dataset in results.items():
    dataset_name = dataset['name'][:24]  # Truncate if too long
    total_reads = dataset['total_reads']
    
    for tool_name, tool in dataset['tools'].items():
        print("{:<25} {:<10} {:<12.2f} {:<10.1f} {:<12,} {:<10.2f}%".format(
            dataset_name,
            tool_name,
            tool['mean_runtime'],
            tool['mean_memory_mb'],
            tool['mean_reads_out'],
            tool['recovery_rate']
        ))
    print()

print("\n" + "=" * 80)
print("KEY OBSERVATIONS:")
print("=" * 80)
print()

# Observations
observations = []

# SPLiT-seq PE
seqproc_time = results['splitseq_pe_raw']['tools']['seqproc']['mean_runtime']
matchbox_time = results['splitseq_pe_raw']['tools']['matchbox']['mean_runtime']
observations.append(f"SPLiT-seq PE: Seqproc is {matchbox_time/seqproc_time:.1f}x faster than Matchbox")

# 10x long reads
seqproc_gridion = results['10x_gridion']['tools']['seqproc']['mean_runtime']
matchbox_gridion = results['10x_gridion']['tools']['matchbox']['mean_runtime']
observations.append(f"10x GridION: Seqproc is {matchbox_gridion/seqproc_gridion:.1f}x faster")

# Memory usage
seqproc_mem = results['splitseq_pe_raw']['tools']['seqproc']['mean_memory_mb']
matchbox_mem = results['splitseq_pe_raw']['tools']['matchbox']['mean_memory_mb']
observations.append(f"Memory: Seqproc uses {matchbox_mem/seqproc_mem:.1f}x less memory")

for obs in observations:
    print(f"• {obs}")
print()

print("All benchmarks now use EDIT DISTANCE (Levenshtein) for consistent comparison.")
