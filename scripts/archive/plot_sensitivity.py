import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn import venn2

def plot_sensitivity_analysis():
    # Data from final run (Approach 1: Anchors + Filter)
    ref_total = 759571
    seqproc_total = 766584
    intersection = 720915
    
    seqproc_unique = seqproc_total - intersection
    ref_unique = ref_total - intersection
    
    # Metrics
    union = seqproc_total + ref_unique
    jaccard = intersection / union
    recall = intersection / ref_total
    precision = intersection / seqproc_total
    f1 = 2 * (precision * recall) / (precision + recall)
    
    # 1. Stacked Bar: Read Composition
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Bar 1: Split-pipe
    p1 = ax.bar('Split-pipe (Ref)', intersection, color='#2E86AB', label='Shared (Valid)')
    p2 = ax.bar('Split-pipe (Ref)', ref_unique, bottom=intersection, color='#A23B72', label='Unique to Tool')
    
    # Bar 2: Seqproc
    p3 = ax.bar('Seqproc', intersection, color='#2E86AB')
    p4 = ax.bar('Seqproc', seqproc_unique, bottom=intersection, color='#F18F01', label='Seqproc Extra')
    
    ax.set_ylabel('Number of Reads')
    ax.set_title('Read Recovery Comparison')
    ax.legend()
    
    # Add text labels
    for c in ax.containers:
        ax.bar_label(c, label_type='center', color='white', weight='bold')
        
    plt.tight_layout()
    plt.savefig('figure_read_fate.png', dpi=300)
    print("Saved figure_read_fate.png")
    
    # 2. Metrics Bar Chart
    fig, ax = plt.subplots(figsize=(8, 6))
    metrics = ['Recall (Sensitivity)', 'Precision', 'F1 Score', 'Jaccard Index']
    values = [recall, precision, f1, jaccard]
    colors = ['#2E86AB', '#F18F01', '#C73E1D', '#3B1F2B']
    
    bars = ax.bar(metrics, values, color=colors)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel('Score')
    ax.set_title('Concordance Metrics (Seqproc vs Split-pipe)')
    
    ax.bar_label(bars, fmt='%.4f', padding=3)
    
    plt.tight_layout()
    plt.savefig('figure_metrics.png', dpi=300)
    print("Saved figure_metrics.png")
    
    # 3. Venn Diagram
    plt.figure(figsize=(8, 6))
    venn = venn2(subsets=(ref_unique, seqproc_unique, intersection), 
                 set_labels=('Split-pipe', 'Seqproc'))
    
    venn.get_patch_by_id('10').set_color('#A23B72')
    venn.get_patch_by_id('01').set_color('#F18F01')
    venn.get_patch_by_id('11').set_color('#2E86AB')
    
    plt.title('Exact Read ID Overlap')
    plt.savefig('figure_venn.png', dpi=300)
    print("Saved figure_venn.png")

if __name__ == "__main__":
    plot_sensitivity_analysis()
