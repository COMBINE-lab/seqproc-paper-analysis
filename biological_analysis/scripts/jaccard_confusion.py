#!/usr/bin/env python3
"""Cell-type-label confusion across tools on the shared cells.

Answers two open items for the paper:
  (1) consensus per-type counts (each shared cell -> one type) that sum to |shared|,
      fixing the tab:jaccard "Cells" column (currently sums to 224, not 220).
  (2) the confusion structure: for OPC and microglia, of the cells any tool calls
      that type, how many are unanimous and where the dissenting labels go. This
      tests whether the low Jaccard is OPC<->oligodendrocyte lineage boundary
      (evidence) rather than generic noise (assertion).

Reuses biological_analysis.load/process/type_cells so labels are IDENTICAL to the
table. Run with the SAME args you passed to biological_analysis.py:

  python biological_analysis/scripts/jaccard_confusion.py \
      <outdir> <min_umi> seqproc:<gene_dir> splitcode:<gene_dir> matchbox:<gene_dir>
"""
import sys, os, json
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import biological_analysis as ba


def main():
    outdir, min_umi = sys.argv[1], int(sys.argv[2])
    tools = [s.split(":", 1) for s in sys.argv[3:]]
    names = [n for n, _ in tools]
    os.makedirs(outdir, exist_ok=True)

    proc = {n: ba.process(ba.load(d), min_umi=min_umi) for n, d in tools}
    called = {n: set(proc[n].obs_names) for n in names}
    shared = sorted(set.intersection(*called.values()))
    lab = {n: dict(zip(proc[n].obs_names, proc[n].obs["cell_type"])) for n in names}
    triples = {c: tuple(lab[n][c] for n in names) for c in shared}   # label per tool, fixed order

    # (1) consensus type per shared cell (majority; tie -> first tool) -> counts sum to |shared|
    def consensus(t):
        top = Counter(t).most_common()
        return top[0][0] if (len(top) == 1 or top[0][1] > top[1][1]) else t[0]
    consensus_counts = Counter(consensus(triples[c]) for c in shared)

    # (2) confusion by type: among cells ANY tool labels ct, unanimity and dissent destinations
    confusion = {}
    for ct in ba.MARKERS:
        cells_ct = [c for c in shared if ct in triples[c]]
        unanimous = sum(1 for c in cells_ct if all(x == ct for x in triples[c]))
        dissent = Counter(x for c in cells_ct for x in triples[c] if x != ct)
        confusion[ct] = {
            "cells_any_tool": len(cells_ct),
            "unanimous_all_three": unanimous,
            "dissent_goes_to": dict(dissent.most_common()),
        }

    out = {
        "tool_order": names,
        "n_shared": len(shared),
        "consensus_type_counts": dict(consensus_counts.most_common()),
        "consensus_sum": sum(consensus_counts.values()),   # should equal n_shared
        "confusion_by_type": confusion,
    }
    print(json.dumps(out, indent=2))
    dst = os.path.join(outdir, "jaccard_confusion.json")
    json.dump(out, open(dst, "w"), indent=2)
    print("saved", dst)


if __name__ == "__main__":
    main()
