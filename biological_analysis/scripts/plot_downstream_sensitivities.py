#!/usr/bin/env python3
"""Plot the four frozen SPLiT-seq PE downstream sensitivity configurations.

The primary Matchbox variant uses canonical barcode lists, fuzzy linkers, and
captured-component length guards; ``expanded`` changes only the barcode lists.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from paper_style import panel, save, set_paper_style, tool_color


TOOLS = ("seqproc", "splitcode", "matchbox")
PAIR_ORDER = ("seqproc|splitcode", "seqproc|matchbox", "splitcode|matchbox")
CONFIG_COLORS = ("#0072B2", "#56B4E9", "#D55E00", "#E69F00")


def load_summary(path: Path) -> dict:
    if path.is_dir():
        path = path / "downstream_summary.json"
    return json.loads(path.read_text())


def annotate_bars(ax, bars, fmt, *, rotation=0):
    for bar in bars:
        value = bar.get_height()
        ax.annotate(
            fmt(value),
            (bar.get_x() + bar.get_width() / 2, value),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=rotation,
        )


def annotate_bars_inside(ax, bars, fmt):
    span = ax.get_ylim()[1] - ax.get_ylim()[0]
    for bar in bars:
        value = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value - 0.025 * span,
            fmt(value),
            ha="center",
            va="top",
            fontsize=7.5,
            rotation=90,
            color="white",
            fontweight="bold",
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-editdist2", type=Path, required=True)
    parser.add_argument("--exact-1mm", type=Path, required=True)
    parser.add_argument("--expanded-editdist2", type=Path, required=True)
    parser.add_argument("--expanded-1mm", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    labels = (
        "Canonical\nEditDist_2",
        "Canonical\n1MM",
        "Expanded\nEditDist_2",
        "Expanded\n1MM",
    )
    summaries = [
        load_summary(args.exact_editdist2),
        load_summary(args.exact_1mm),
        load_summary(args.expanded_editdist2),
        load_summary(args.expanded_1mm),
    ]
    x = np.arange(len(labels))
    set_paper_style()
    # Reserve a narrow band between each title and plotting region for its
    # legend. Keeping legends outside the axes avoids covering bars or labels
    # without introducing a large gap between the two panel rows.
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.4))

    # A: STAR-valid barcode fraction. The narrow, explicitly labeled axis is
    # intentional: all four upstream products have already passed their own
    # structural filters, so the remaining downstream correction loss is small.
    ax = axes[0, 0]
    ax.set_ylim(98.7, 100.08)
    width = 0.24
    for index, tool in enumerate(TOOLS):
        values = [100 * s["tools"][tool]["STARsolo"]["Reads With Valid Barcodes"] for s in summaries]
        bars = ax.bar(
            x + (index - 1) * width,
            values,
            width,
            label=tool,
            color=tool_color(tool, index),
        )
        annotate_bars_inside(ax, bars, lambda value: f"{value:.2f}")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Reads with valid barcodes (%)")
    ax.set_title("STARsolo barcode validity (expanded y-axis)", pad=36)
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=3,
        frameon=False,
        borderaxespad=0,
    )
    panel(ax, "A")

    # B: thresholded cell calls. Preserve a zero baseline for counts.
    ax = axes[0, 1]
    for index, tool in enumerate(TOOLS):
        values = [s["tools"][tool]["called_cells_min_umi"] for s in summaries]
        bars = ax.bar(
            x + (index - 1) * width,
            values,
            width,
            label=tool,
            color=tool_color(tool, index),
        )
        annotate_bars(ax, bars, lambda value: f"{value:.0f}")
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 245)
    ax.set_ylabel("Called cells (at least 200 UMI)")
    ax.set_title("Cell-call stability", pad=36)
    panel(ax, "B")

    # C: conservative summary of count concordance: the minimum correlation
    # among all three tool pairs, so no favorable pair is selected post hoc.
    ax = axes[1, 0]
    metrics = (
        ("Per-gene Pearson", "per_gene_pearson_log1p", "#0072B2"),
        ("Per-barcode Pearson", "per_barcode_pearson_log1p", "#D55E00"),
    )
    width2 = 0.34
    for index, (name, key, color) in enumerate(metrics):
        values = [min(s["pairs"][pair][key] for pair in PAIR_ORDER) for s in summaries]
        bars = ax.bar(x + (index - 0.5) * width2, values, width2, label=name, color=color)
        annotate_bars(ax, bars, lambda value: f"{value:.3f}")
    ax.set_xticks(x, labels)
    ax.set_ylim(0.74, 1.015)
    ax.set_ylabel("Minimum pairwise correlation")
    ax.set_title("Worst-pair count concordance", pad=36)
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=2,
        frameon=False,
        borderaxespad=0,
    )
    panel(ax, "C")

    # D: metrics calculated on shared called cells. These are deliberately kept
    # separate from count correlations because low-signal marker labels and
    # independent Leiden clusterings are less stable.
    ax = axes[1, 1]
    biological = (
        ("All-tool type agreement", "all_tool_celltype_agreement", "#009E73"),
        ("Mean type Jaccard", "celltype_jaccard_mean", "#CC79A7"),
        ("Mean cluster ARI", "cluster_ari_mean", "#E69F00"),
    )
    width3 = 0.24
    for index, (name, key, color) in enumerate(biological):
        values = [s[key] for s in summaries]
        bars = ax.bar(x + (index - 1) * width3, values, width3, label=name, color=color)
        annotate_bars(ax, bars, lambda value: f"{value:.3f}")
    ax.set_xticks(x, labels)
    ax.set_ylim(0.5, 1.01)
    ax.set_ylabel("Agreement score")
    ax.set_title("Shared-cell biological concordance", pad=36)
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=3,
        frameon=False,
        borderaxespad=0,
    )
    panel(ax, "D")

    fig.tight_layout(w_pad=2.2, h_pad=3.0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    print(save(fig, str(args.output)))


if __name__ == "__main__":
    main()
