#!/usr/bin/env python3
"""Regenerate precision-recall figure with corrected labels."""

import json
import matplotlib.pyplot as plt
from pathlib import Path

# Load existing results
results_file = Path('results/precision_recall/precision_recall_results.json')
with open(results_file) as f:
    data = json.load(f)

tolerance_levels = data['tolerance_levels']
results = data['results']

COLORS = {
    'seqproc': '#4DBBD5',
    'matchbox': '#E64B35',
    'splitcode': '#3C5488',
}
MARKERS = {0: 'o', 1: 's', 2: '^', 3: 'D'}

fig, ax = plt.subplots(figsize=(8, 6))

for tool in ['seqproc', 'matchbox']:
    if tool not in results:
        continue
    frac_correct = results[tool]['frac_correct']
    frac_incorrect = results[tool]['frac_incorrect']
    
    # Plot line connecting points
    ax.plot(frac_incorrect, frac_correct, 
            color=COLORS[tool], linewidth=2, label=tool, alpha=0.8)
    
    # Plot points with different markers for tolerance levels
    for i, (fi, fc) in enumerate(zip(frac_incorrect, frac_correct)):
        ax.scatter(fi, fc, color=COLORS[tool], marker=MARKERS[i], 
                  s=100, edgecolor='black', linewidth=1, zorder=5)

# Add legend for tolerance levels (CORRECTED labels)
for i, tol in enumerate(tolerance_levels):
    ax.scatter([], [], color='gray', marker=MARKERS[i], s=80, 
              label=f'Tolerance {tol}', edgecolor='black')

ax.set_xlabel('Fraction of reads with incorrect barcodes', fontsize=12)
ax.set_ylabel('Fraction of reads with correct barcodes', fontsize=12)
ax.set_title('Precision-Recall Analysis\n(Synthetic SPLiT-seq Data)', fontsize=14, fontweight='bold')
ax.legend(loc='lower right', fontsize=10)
ax.grid(True, alpha=0.3)

# Set axis limits
tools = list(results.keys())
max_incorrect = max(max(results[t]['frac_incorrect']) for t in tools)
ax.set_xlim(-0.005, max(0.1, max_incorrect * 1.2))
ax.set_ylim(0, 1.05)

# Add explanatory note
ax.text(0.02, 0.02, 
        'Tolerance = max allowed mismatches for whitelist matching\n'
        'Higher tolerance → more reads recovered but more errors',
        transform=ax.transAxes, fontsize=9, verticalalignment='bottom',
        style='italic', color='gray', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()

output_dir = Path('results/paper_figures_fixed')
output_dir.mkdir(parents=True, exist_ok=True)
plt.savefig(output_dir / 'fig3_precision_recall.png', dpi=300, bbox_inches='tight')
plt.savefig(output_dir / 'fig3_precision_recall.pdf', bbox_inches='tight')
plt.close()

print(f"Fixed precision-recall figure saved to {output_dir}")
