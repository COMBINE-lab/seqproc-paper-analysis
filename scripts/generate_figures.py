#!/usr/bin/env python3
"""
Generate all paper figures from concordance results
and existing performance benchmark results.

Produces:
  1. Concordance heatmaps (Jaccard index per dataset)
  2. Recovery comparison bar chart (all tools, all datasets)
  3. Hamming vs Edit distance comparison
  4. Discordant read summary
  5. Performance table (runtime + memory from benchmark results)
  6. Updated combined benchmark_results.json

Usage:
    python3 scripts/generate_figures.py
"""

import json
import os
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent

# Publication-quality defaults
plt.rcParams.update({
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})
CONCORDANCE_DIR = PROJECT_ROOT / "results" / "concordance"
PERF_JSON = PROJECT_ROOT / "results" / "paper_figures" / "benchmark_results.json"
PERF_BACKUP = PROJECT_ROOT / "results" / "paper_figures" / "benchmark_results_perf.json"
OUTPUT_DIR = PROJECT_ROOT / "results" / "paper_figures"

COLORS = {
    'seqproc': '#2E86AB',
    'matchbox': '#E94F37',
    'splitcode': '#7B2D8E',
}

TOOL_ORDER = ['seqproc', 'matchbox', 'splitcode']

# Dataset display order and labels
DS_ORDER = ['splitseq_pe', 'lr_splitseq', '10x_short', 'sciseq']
DS_LABELS = {
    'splitseq_pe': 'SPLiT-seq PE',
    'lr_splitseq': 'LR-SPLiT-seq',
    '10x_short': '10x Short',
    'sciseq': 'sci-RNA-seq3',
}

# Map performance benchmark keys to concordance keys
PERF_TO_CONCORDANCE = {
    'splitseq_pe_raw': 'splitseq_pe',
    'splitseq_se_raw': 'lr_splitseq',
    '10x_short': '10x_short',
    'sciseq': 'sciseq',
}


def load_data():
    """Load concordance results and performance benchmark results.

    Performance data is backed up on first run to prevent loss when this script
    overwrites benchmark_results.json with merged concordance output.
    """
    concordance = {}
    concordance_path = CONCORDANCE_DIR / "concordance_results.json"
    if concordance_path.exists():
        with open(concordance_path) as f:
            concordance = json.load(f)

    # Back up original performance data before it gets overwritten
    if not PERF_BACKUP.exists() and PERF_JSON.exists():
        shutil.copy2(PERF_JSON, PERF_BACKUP)
        print(f"Backed up performance data to: {PERF_BACKUP.name}")

    # Read performance from backup (original data), fall back to main JSON
    perf = {}
    perf_source = PERF_BACKUP if PERF_BACKUP.exists() else PERF_JSON
    if perf_source.exists():
        with open(perf_source) as f:
            perf = json.load(f)

    return concordance, perf


def fig_concordance_heatmaps(concordance, output_dir):
    """Generate concordance heatmap grid (one subplot per dataset)."""
    datasets = [k for k in DS_ORDER if k in concordance]
    if not datasets:
        print("  [SKIP] No concordance data")
        return

    n = len(datasets)
    fig, axes = plt.subplots(1, n, figsize=(7.0, 1.8), squeeze=False)

    for idx, ds_key in enumerate(datasets):
        ax = axes[0][idx]
        res = concordance[ds_key]
        conc = res.get("concordance", {})
        pairwise = conc.get("pairwise", [])

        # Build Jaccard matrix
        tools = TOOL_ORDER
        matrix = np.ones((3, 3))  # diagonal = 1.0
        for pair in pairwise:
            ta, tb = pair["tool_a"], pair["tool_b"]
            j = pair["jaccard"]
            if ta in tools and tb in tools:
                i1, i2 = tools.index(ta), tools.index(tb)
                matrix[i1][i2] = j
                matrix[i2][i1] = j

        im = ax.imshow(matrix, cmap='YlOrRd', vmin=0.0, vmax=1.0, aspect='equal')

        # Annotate cells
        for i in range(3):
            for j in range(3):
                val = matrix[i][j]
                color = 'white' if val < 0.5 else 'black'
                ax.text(j, i, f'{val:.3f}', ha='center', va='center',
                        fontsize=7, fontweight='bold', color=color)

        ax.set_xticks(range(3))
        ax.set_yticks(range(3))
        ax.set_xticklabels([t.capitalize() for t in tools])
        ax.set_yticklabels([t.capitalize() for t in tools])
        ax.set_title(DS_LABELS.get(ds_key, ds_key), fontweight='bold')

    fig.suptitle('Pairwise Concordance (Jaccard Index)', fontsize=10, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / 'fig_concordance_heatmaps.pdf', bbox_inches='tight')
    plt.savefig(output_dir / 'fig_concordance_heatmaps.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: fig_concordance_heatmaps.png")


def fig_recovery_comparison(concordance, output_dir):
    """Generate grouped bar chart of recovery rates across all datasets."""
    datasets = [k for k in DS_ORDER if k in concordance]
    if not datasets:
        return

    fig, ax = plt.subplots(figsize=(7.0, 3.0))

    x = np.arange(len(datasets))
    width = 0.25

    for i, tool in enumerate(TOOL_ORDER):
        rates = []
        for ds_key in datasets:
            pct = concordance[ds_key].get("recovery_pct", {}).get(tool, 0)
            rates.append(pct)
        bars = ax.bar(x + i * width, rates, width, label=tool.capitalize(),
                      color=COLORS[tool], edgecolor='white', linewidth=0.5)
        # Add value labels
        for bar, val in zip(bars, rates):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.8,
                        f'{val:.1f}%', ha='center', va='bottom', fontsize=6)

    ax.set_ylabel('Recovery Rate (%)')
    ax.set_title('Read Recovery by Dataset and Tool', fontweight='bold')
    ax.set_xticks(x + width)
    ax.set_xticklabels([DS_LABELS.get(k, k) for k in datasets])
    ax.legend()
    ax.set_ylim(0, 110)
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_dir / 'fig_recovery_comparison.pdf', bbox_inches='tight')
    plt.savefig(output_dir / 'fig_recovery_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: fig_recovery_comparison.png")


def fig_hamming_vs_edit(concordance, output_dir):
    """Generate hamming vs edit distance comparison chart."""
    datasets = []
    for ds_key in DS_ORDER:
        if ds_key in concordance and concordance[ds_key].get("hamming_vs_edit"):
            datasets.append(ds_key)

    if not datasets:
        print("  [SKIP] No hamming vs edit data")
        return

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8))

    # Left: Paired bar chart (hamming vs edit reads)
    ax = axes[0]
    x = np.arange(len(datasets))
    width = 0.35

    ham_reads = [concordance[k]["hamming_vs_edit"]["hamming_reads"] for k in datasets]
    edit_reads = [concordance[k]["hamming_vs_edit"]["edit_reads"] for k in datasets]

    ax.bar(x - width/2, ham_reads, width, label='Hamming', color='#95C8D8', edgecolor='white')
    ax.bar(x + width/2, edit_reads, width, label='Edit', color='#2E86AB', edgecolor='white')

    ax.set_ylabel('Recovered Reads')
    ax.set_title('seqproc: Hamming vs Edit Distance', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([DS_LABELS.get(k, k) for k in datasets])
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Format y-axis with commas
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))

    # Right: Edit distance gain percentage
    ax = axes[1]
    gains = [concordance[k]["hamming_vs_edit"]["edit_gain_pct"] for k in datasets]
    colors = ['#2E86AB' if g > 0 else '#E94F37' for g in gains]

    bars = ax.bar(x, gains, 0.5, color=colors, edgecolor='white')
    for bar, val in zip(bars, gains):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                f'+{val:.1f}%', ha='center', va='bottom', fontsize=7, fontweight='bold')

    ax.set_ylabel('Edit Distance Gain (%)')
    ax.set_title('Additional Reads from Edit Distance', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([DS_LABELS.get(k, k) for k in datasets])
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_dir / 'fig_hamming_vs_edit.pdf', bbox_inches='tight')
    plt.savefig(output_dir / 'fig_hamming_vs_edit.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: fig_hamming_vs_edit.png")


def fig_discordant_summary(concordance, output_dir):
    """Generate stacked bar chart showing concordance breakdown."""
    datasets = [k for k in DS_ORDER if k in concordance]
    if not datasets:
        return

    fig, ax = plt.subplots(figsize=(7.0, 3.5))

    x = np.arange(len(datasets))

    # For each dataset, show: consensus | seqproc-only | matchbox-only | splitcode-only
    consensus = []
    sp_only = []
    mb_only = []
    sc_only = []
    total_reads = []

    for ds_key in datasets:
        disc = concordance[ds_key].get("discordant", {})
        con = disc.get("all_tools_consensus", 0)
        union = disc.get("any_tool_union", 0)
        consensus.append(con)
        sp_only.append(disc.get("seqproc_unique", 0))
        mb_only.append(disc.get("matchbox_unique", 0))
        sc_only.append(disc.get("splitcode_unique", 0))
        total_reads.append(concordance[ds_key].get("total_reads", 0))

    # The "shared but not consensus" region
    shared_not_all = []
    for i, ds_key in enumerate(datasets):
        union_val = concordance[ds_key].get("discordant", {}).get("any_tool_union", 0)
        remaining = union_val - consensus[i] - sp_only[i] - mb_only[i] - sc_only[i]
        shared_not_all.append(max(0, remaining))

    bar_width = 0.6

    # Stack: consensus, shared_not_all, seqproc-only, matchbox-only, splitcode-only
    b1 = ax.bar(x, consensus, bar_width, label='All tools agree', color='#4CAF50')
    bottom = np.array(consensus, dtype=float)

    b2 = ax.bar(x, shared_not_all, bar_width, bottom=bottom, label='Shared (2 tools)', color='#8BC34A')
    bottom += np.array(shared_not_all, dtype=float)

    b3 = ax.bar(x, sp_only, bar_width, bottom=bottom, label='seqproc only', color=COLORS['seqproc'])
    bottom += np.array(sp_only, dtype=float)

    b4 = ax.bar(x, mb_only, bar_width, bottom=bottom, label='matchbox only', color=COLORS['matchbox'])
    bottom += np.array(mb_only, dtype=float)

    b5 = ax.bar(x, sc_only, bar_width, bottom=bottom, label='splitcode only', color=COLORS['splitcode'])

    # Add total reads line
    for i in range(len(datasets)):
        ax.axhline(y=total_reads[i], color='gray', linestyle='--', alpha=0.3, xmin=0, xmax=1)

    ax.set_ylabel('Reads')
    ax.set_title('Read Recovery Breakdown by Concordance', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([DS_LABELS.get(k, k) for k in datasets])
    ax.legend(loc='upper right')
    ax.grid(axis='y', alpha=0.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))

    plt.tight_layout()
    plt.savefig(output_dir / 'fig_discordant_summary.pdf', bbox_inches='tight')
    plt.savefig(output_dir / 'fig_discordant_summary.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: fig_discordant_summary.png")


def fig_performance_table(perf, concordance, output_dir):
    """Generate performance summary table figure combining perf benchmarks + concordance recovery."""
    # Use perf for runtime/memory (3-replicate means), concordance for recovery
    # Map perf keys to display order

    columns = ['Dataset', 'Tool', 'Recovery %', 'Runtime (s)', 'Memory (MB)']
    cell_text = []
    row_colors = []

    for ds_key in DS_ORDER:
        # Find performance key
        perf_key = None
        for k, v in PERF_TO_CONCORDANCE.items():
            if v == ds_key:
                perf_key = k
                break

        label = DS_LABELS.get(ds_key, ds_key)

        for tool in TOOL_ORDER:
            # Recovery from concordance
            rec_pct = concordance.get(ds_key, {}).get("recovery_pct", {}).get(tool, "N/A")
            if isinstance(rec_pct, (int, float)):
                rec_str = f"{rec_pct:.1f}%"
            else:
                rec_str = str(rec_pct)

            # Runtime and memory from perf benchmarks
            perf_data = perf.get(perf_key, perf.get(ds_key, {})) if perf_key else perf.get(ds_key, {})
            if tool in perf_data.get("tools", {}):
                perf_tool = perf_data["tools"][tool]
                rt_mean = perf_tool.get("mean_runtime", 0)
                rt_std = perf_tool.get("std_runtime", 0)
                mem = perf_tool.get("mean_memory_mb", 0)
                rt_str = f"{rt_mean:.1f} +/- {rt_std:.1f}"
                mem_str = f"{mem:.0f}"
            else:
                rt_str = "N/A"
                mem_str = "N/A"

            cell_text.append([label, tool.capitalize(), rec_str, rt_str, mem_str])
            row_colors.append(COLORS.get(tool, '#FFFFFF'))

    if not cell_text:
        print("  [SKIP] No performance data")
        return

    fig, ax = plt.subplots(figsize=(14, len(cell_text) * 0.4 + 2))
    ax.axis('off')

    table = ax.table(cellText=cell_text, colLabels=columns, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)

    # Style
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight='bold')
            cell.set_facecolor('#e0e0e0')
        elif col == 1:
            # Color tool name cells
            tool_name = cell_text[row - 1][1].lower()
            cell.set_facecolor(COLORS.get(tool_name, 'white'))
            cell.set_text_props(color='white', weight='bold')
        cell.set_edgecolor('#cccccc')

    plt.title("Benchmark Summary: Recovery, Runtime, and Memory",
              fontweight='bold', fontsize=13, pad=15)
    plt.tight_layout()
    plt.savefig(output_dir / 'fig_performance_table.pdf', bbox_inches='tight')
    plt.savefig(output_dir / 'fig_performance_table.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: fig_performance_table.png")


def update_benchmark_json(concordance, perf, output_dir):
    """Update benchmark_results.json with concordance data merged with performance benchmarks."""
    combined = {}

    for ds_key in DS_ORDER:
        perf_key = None
        for k, v in PERF_TO_CONCORDANCE.items():
            if v == ds_key:
                perf_key = k
                break

        conc = concordance.get(ds_key, {})
        label = DS_LABELS.get(ds_key, ds_key)

        entry = {
            "name": label,
            "total_reads": conc.get("total_reads", 0),
            "tools": {},
            "concordance": conc.get("concordance", {}),
            "discordant": conc.get("discordant", {}),
            "hamming_vs_edit": conc.get("hamming_vs_edit", {}),
        }

        for tool in TOOL_ORDER:
            tool_data = {}

            # Concordance recovery
            rec = conc.get("recovery", {}).get(tool, 0)
            rec_pct = conc.get("recovery_pct", {}).get(tool, 0)
            tool_data["reads_out"] = rec
            tool_data["recovery_rate"] = rec_pct

            # Performance benchmarks
            perf_data = perf.get(perf_key, perf.get(ds_key, {})) if perf_key else perf.get(ds_key, {})
            if tool in perf_data.get("tools", {}):
                perf_tool = perf_data["tools"][tool]
                tool_data["mean_runtime"] = perf_tool.get("mean_runtime", 0)
                tool_data["std_runtime"] = perf_tool.get("std_runtime", 0)
                tool_data["mean_memory_mb"] = perf_tool.get("mean_memory_mb", 0)

            if tool_data:
                entry["tools"][tool] = tool_data

        combined[ds_key] = entry

    out_path = output_dir / "benchmark_results.json"
    with open(out_path, 'w') as f:
        json.dump(combined, f, indent=2)
    print(f"  Saved: {out_path}")

    return combined


def main():
    print("=" * 70)
    print("FIGURE GENERATION")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    concordance, perf = load_data()
    if not concordance:
        print("[ERROR] No concordance results found. Run concordance_analysis.py first.")
        return

    print(f"Concordance datasets: {list(concordance.keys())}")
    print(f"Performance datasets: {list(perf.keys())}")

    print("\nGenerating figures...")

    # NOTE: Initial LR-SPLiT-seq benchmark used forward-only config (23.8% recovery)
    # Current config uses annotation+edit (49.9% recovery) -- this is the updated number
    # The initial runtime for LR-SPLiT-seq is stale (forward-only was 2.1s,
    # annotation+edit is ~5.1s). Current perf is in concordance_results.json.

    fig_concordance_heatmaps(concordance, OUTPUT_DIR)
    fig_recovery_comparison(concordance, OUTPUT_DIR)
    fig_hamming_vs_edit(concordance, OUTPUT_DIR)
    fig_discordant_summary(concordance, OUTPUT_DIR)
    fig_performance_table(perf, concordance, OUTPUT_DIR)

    print("\nUpdating benchmark_results.json...")
    combined = update_benchmark_json(concordance, perf, OUTPUT_DIR)

    print(f"\n{'='*70}")
    print("ALL FIGURES GENERATED")
    print(f"{'='*70}")
    print(f"Output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
