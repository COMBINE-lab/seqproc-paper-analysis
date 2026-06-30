#!/usr/bin/env python3
"""Runtime + peak-RAM report for the Phase 2A pipeline.

  resource_report.py <resources.csv> <outdir>

resources.csv (written by run_phase2a.sh) has columns: step,tool,seconds,peak_ram_mb.
Produces a paper-ready figure (read-processing runtime and peak RAM per tool) and markdown +
CSV tables (full per-step breakdown) for the supplement. RAM for a tool step is that tool
executable's peak RSS; STARsolo RAM is dominated by the shared genome index.
"""
import sys, os, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from paper_style import set_paper_style, tool_color, panel, save

TOOLS = ["seqproc", "splitcode", "matchbox"]

def main():
    res_csv, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    set_paper_style()
    data = {}                                   # data[step][tool] = (seconds, ram_mb)
    for r in csv.DictReader(open(res_csv)):
        data.setdefault(r["step"], {})[r["tool"]] = (float(r["seconds"]), float(r["peak_ram_mb"]))

    bc_t = [data["barcode"][t][0] for t in TOOLS]
    bc_r = [data["barcode"][t][1] for t in TOOLS]
    ss_t = [data["starsolo"][t][0] for t in TOOLS]
    colors = [tool_color(t, i) for i, t in enumerate(TOOLS)]

    # ---- figure: read-processing runtime, peak RAM, and total pipeline runtime ----
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    x = np.arange(len(TOOLS))
    b0 = ax[0].bar(x, bc_t, 0.6, color=colors)
    ax[0].set_xticks(x); ax[0].set_xticklabels(TOOLS); ax[0].set_ylabel("Runtime (s)")
    ax[0].set_title("Read processing runtime"); ax[0].bar_label(b0, fmt="%.1f", padding=2)
    ax[0].set_ylim(0, max(bc_t) * 1.18); panel(ax[0], "A")

    b1 = ax[1].bar(x, bc_r, 0.6, color=colors)
    ax[1].set_xticks(x); ax[1].set_xticklabels(TOOLS); ax[1].set_ylabel("Peak RAM (MB)")
    ax[1].set_title("Read processing peak memory"); ax[1].bar_label(b1, fmt="%.0f", padding=2)
    ax[1].set_ylim(0, max(bc_r) * 1.18); panel(ax[1], "B")

    # total pipeline runtime per tool (read processing + STARsolo), stacked.
    # color by COMPONENT (not tool; the tool is on the x-axis) so it matches the legend.
    ax[2].bar(x, bc_t, 0.6, color="#4C72B0", label="read processing")
    ax[2].bar(x, ss_t, 0.6, bottom=bc_t, color="0.75", label="STARsolo")
    tot = [bc_t[i] + ss_t[i] for i in range(len(TOOLS))]
    for xi, tv in zip(x, tot):
        ax[2].text(xi, tv + max(tot) * 0.012, f"{tv:.1f}", ha="center", va="bottom")
    ax[2].set_xticks(x); ax[2].set_xticklabels(TOOLS); ax[2].set_ylabel("Runtime (s)")
    ax[2].set_title("Total runtime per tool"); ax[2].set_ylim(0, max(tot) * 1.2)
    ax[2].legend(loc="upper center"); panel(ax[2], "C")
    fig.tight_layout()
    print("saved", save(fig, os.path.join(outdir, "resource_usage")), "(+ .pdf)")

    # ---- tables (markdown + csv) ----
    def cell(step, tool, idx):
        return data.get(step, {}).get(tool, (float("nan"), float("nan")))[idx]
    def fmt(v):
        return f"{v:.1f}" if v == v else "n/a"

    md = ["| Step | seqproc | splitcode | matchbox |", "|---|---|---|---|"]
    for step, idx, label in [("barcode", 0, "Read processing (s)"),
                             ("starsolo", 0, "STARsolo (s)"),
                             ("barcode", 1, "Read processing RAM (MB)"),
                             ("starsolo", 1, "STARsolo RAM (MB)")]:
        md.append(f"| {label} | " + " | ".join(fmt(cell(step, t, idx)) for t in TOOLS) + " |")
    # shared downstream analysis, split into its two steps
    bio = cell("analysis", "biological", 0); cnt = cell("analysis", "count", 0)
    md.append(f"| Biological analysis (s, shared) | {fmt(bio)} (all tools) |||")
    md.append(f"| Count concordance (s, shared) | {fmt(cnt)} (all tools) |||")
    # total per tool (read processing + STARsolo)
    tot = [fmt(cell('barcode', t, 0) + cell('starsolo', t, 0)) for t in TOOLS]
    md.append("| **Total per tool (s)** | " + " | ".join(f"**{v}**" for v in tot) + " |")

    open(os.path.join(outdir, "resource_table.md"), "w").write("\n".join(md) + "\n")
    print("saved", os.path.join(outdir, "resource_table.md"))
    print("\n".join(md))

if __name__ == "__main__":
    main()
