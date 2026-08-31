"""Shared publication style for Phase 2A figures.

Fixed colorblind-safe (Okabe-Ito) colors per tool so every figure is consistent, clean
spines/fonts sized for a journal column, and a save() that writes both PNG (300 dpi) and a
vector PDF with editable text. Import and call set_paper_style() at the top of a plotting script.
"""
import matplotlib as mpl

TOOL_COLORS = {
    "seqproc":   "#0072B2",   # blue
    "splitcode": "#D55E00",   # vermillion
    "matchbox":  "#009E73",   # bluish green
}
_FALLBACK = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]

# Pairwise comparisons are composite categories, so they deliberately do not
# reuse the colors assigned to individual tools.  The first three colors are
# optimized for the three-tool figures in the paper and remain distinguishable
# in grayscale; the remaining colors support analyses with more tool pairs.
PAIR_COLORS = ["#6F4C9B", "#4D6A6D", "#B88900", "#A05195", "#665191", "#2F4B7C"]

def tool_color(name, i=0):
    return TOOL_COLORS.get(name, _FALLBACK[i % len(_FALLBACK)])

def pair_color(i=0):
    return PAIR_COLORS[i % len(PAIR_COLORS)]

def set_paper_style():
    mpl.rcParams.update({
        "savefig.dpi": 300, "figure.dpi": 110, "savefig.bbox": "tight",
        "pdf.fonttype": 42, "ps.fonttype": 42,            # editable text in vector editors
        "font.family": "sans-serif", "figure.facecolor": "white",
        "font.size": 11, "axes.titlesize": 12, "axes.titleweight": "regular",
        "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "legend.fontsize": 9.5, "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
    })

def panel(ax, letter, dx=-0.14, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=14,
            fontweight="bold", va="bottom", ha="right")

def save(fig, path_noext):
    fig.savefig(path_noext + ".png")
    fig.savefig(path_noext + ".pdf")
    return path_noext + ".png"
