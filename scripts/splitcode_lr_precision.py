#!/usr/bin/env python3
"""splitcode LR precision (DUAL-pass union) vs the short-read real-cell set.

splitcode encodes each read's assignment as a fixed-length code (the prefix of the output
read sequence); the mapping file translates code -> comma-separated tags. We decode:
  bc1 = first 6 bp of the bc1 tag's sequence (the round-1 identity; the 8bp tag's last 2 bp
        are RT-primer bases, not cell identity -- verified: all 20 bc1 tags' first 6 bp are
        in the round-1 pool)
  bc2 = the bc2 tag's sequence
  bc3 = the bc3 tag's sequence if splitcode tagged it, else extracted from the read
        (linker-anchored, snapped to the bc23 whitelist)
A read is CORRECT if EITHER pass yields a full CB (bc1_bc2_bc3) that is a real cell.
Emitted = union of read ids seen in either pass. precision = correct / emitted.

  splitcode_lr_precision.py real_cells.txt config \
     --pass out_fwd.fq map_fwd.txt --pass out_rev.fq map_rev.txt
"""
import sys, os, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import edit_tolerant_validity as V
from short_read_cb_caller import snap_dict, load


def load_tagseqs(config_path):
    m = {}
    for line in open(config_path):
        s = line.strip()
        if not s or s[0] in "#@":
            continue
        p = s.split("\t")
        if len(p) >= 2 and p[0] not in ("ID", "group") and set(p[1]) <= set("ACGTN"):
            m[p[0]] = p[1]
    return m


def load_mapping(path):
    m = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            m[p[0]] = p[1]
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("real_cells")
    ap.add_argument("config")
    ap.add_argument("--pass", dest="passes", nargs=2, action="append",
                    metavar=("OUT_FQ", "MAPPING"), required=True)
    ap.add_argument("--chem", default="lr")
    ap.add_argument("--code-len", type=int, default=16)
    ap.add_argument("--max-linker-edit", type=int, default=4)
    a = ap.parse_args()

    L1, _ = V.LINKERS[a.chem]
    wl = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "configs", "seqproc")
    wl23 = load(os.path.join(wl, "splitseq_bc23_whitelist.txt")); n23 = len(wl23[0])
    m23 = snap_dict(wl23)
    tagseq = load_tagseqs(a.config)
    real = set(l.strip() for l in open(a.real_cells) if l.strip())

    emitted = {}   # read_id -> correct (bool); union across passes
    full = 0       # reads (union) that formed a complete CB at all
    for out_fq, mapping_f in a.passes:
        mp = load_mapping(mapping_f)
        with open(out_fq) as f:
            while True:
                h = f.readline()
                if not h:
                    break
                seq = f.readline().rstrip("\n"); f.readline(); f.readline()
                rid = h[1:].split()[0]
                emitted.setdefault(rid, False)
                tags = mp.get(seq[:a.code_len], "")
                b1 = b2 = b3 = None
                for t in tags.split(","):
                    if t.startswith("bc1_") and t in tagseq:
                        b1 = tagseq[t][:6]
                    elif t.startswith("bc2_") and t in tagseq:
                        b2 = tagseq[t]
                    elif t.startswith("bc3_") and t in tagseq:
                        b3 = tagseq[t]
                if b1 and b2 and not b3:                      # bc3 not a tag -> extract
                    region = seq[a.code_len:]
                    l1 = V.find(L1, region, a.max_linker_edit)
                    if l1:
                        b3 = m23.get(region[l1[0] - n23:l1[0]])
                if b1 and b2 and b3:
                    full += 0 if emitted[rid] else 1         # count a read's first full CB once
                    if f"{b1}_{b2}_{b3}" in real:
                        emitted[rid] = True

    n = len(emitted); correct = sum(emitted.values())
    P = 100 * correct / n if n else 0.0
    print(f"emitted(union)={n}  formed_full_CB>={full}  correct={correct}  precision={P:.3f}%")


if __name__ == "__main__":
    main()
