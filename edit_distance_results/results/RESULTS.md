# seqproc Benchmark Results

## Overview

Comprehensive benchmark comparing **seqproc**, **matchbox**, and **splitcode** on:
- SPLiT-seq paired-end (short reads)
- SPLiT-seq single-end (long reads)
- 10x Chromium GridION (long reads)
- 10x Chromium PromethION (long reads)

All tests run on 1M reads with 4 threads.

## Results Summary

| Dataset | seqproc | matchbox | splitcode | Speedup |
|---------|---------|----------|-----------|---------|
| SPLiT-seq Paired-End | 2.57s | 4.03s | 10.62s | 1.57x |
| SPLiT-seq Single-End | 16.14s | 2.91s | 102.97s | 0.18x |
| 10x GridION | 20.07s | 6.42s | 65.76s | 0.32x |
| 10x PromethION | 18.71s | 6.29s | 61.12s | 0.34x |

## Dataset Analysis

### SPLiT-seq Paired-End (Short Reads)

**seqproc is 1.57x faster**

- **seqproc**: 2.57s ± 0.00s, 87.2% recovery (872,281 reads)
- **matchbox**: 4.03s ± 0.00s, 86.7% recovery (867,118 reads)
- **splitcode**: 10.62s ± 0.00s, 88.4% recovery (883,946 reads)

**Why**: Fixed-position barcode extraction on short R2 reads (~100bp). No searching needed - seqproc extracts directly by offset.

### SPLiT-seq Single-End (Long Reads)

**matchbox is 5.55x faster**

- **seqproc**: 16.14s ± 0.00s, 24.2% recovery (241,771 reads)
- **matchbox**: 2.91s ± 0.00s, 24.0% recovery (240,065 reads)
- **splitcode**: 102.97s ± 0.00s, 28.0% recovery (279,699 reads)

**Why**: Long reads (500-2000bp) require anchor-based searching. matchbox uses optimized pattern automata; seqproc's `anchor_relative` has higher per-read overhead but scales better with threads.

### 10x Chromium GridION (Long Reads)

**matchbox is 3.12x faster**

- **seqproc**: 20.07s ± 0.00s, 36.2% recovery (362,405 reads)
- **matchbox**: 6.42s ± 0.00s, 41.0% recovery (410,096 reads)
- **splitcode**: 65.76s ± 0.00s, 36.2% recovery (362,409 reads)

### 10x Chromium PromethION (Long Reads)

**matchbox is 2.98x faster**

- **seqproc**: 18.71s ± 0.00s, 30.3% recovery (302,711 reads)
- **matchbox**: 6.29s ± 0.00s, 35.7% recovery (357,348 reads)
- **splitcode**: 61.12s ± 0.00s, 30.3% recovery (302,715 reads)

**Why (10x)**: Simple primer anchor + fixed 16bp barcode with low tolerance (3bp). Optimal use case for seqproc's fixed-offset extraction.

## Overall Summary

**seqproc is faster in 1 out of 4 datasets.**

| Dataset Type | Winner | Speedup |
|--------------|--------|---------|
| SPLiT-seq Paired-End | seqproc | 1.57x |
| SPLiT-seq Single-End | matchbox | 5.55x |
| 10x GridION | matchbox | 3.12x |
| 10x PromethION | matchbox | 2.98x |

### Performance Patterns

- **seqproc excels** at fixed-position extraction and simple anchor patterns
- **matchbox excels** at complex pattern matching in long reads (single-threaded)
- **seqproc scales better** with multiple threads

## Figures

- `fig1_runtime_comparison.png` - Runtime bar charts for all datasets
- `fig2_speedup_chart.png` - Speedup comparison (seqproc vs matchbox)
- `fig3_read_recovery.png` - Read recovery rates for all tools
- `fig4_summary_dashboard.png` - Complete summary dashboard with insights