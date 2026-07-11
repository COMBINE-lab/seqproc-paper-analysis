#!/usr/bin/env python3
"""Read-level concordance: pairwise Jaccard of the recovered read-ID sets. A direct,
embedding-free measure of how much the tools agree on which reads to keep.

  read_set_jaccard.py <out.json> <name1>:<bc1.fq> <name2>:<bc2.fq> [<name3>:<bc3.fq> ...]
"""
import sys, json, itertools

def ids(path):
    s = set()
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                s.add(line.split()[0])
    return s

def main():
    out = sys.argv[1]
    tools = [a.split(":", 1) for a in sys.argv[2:]]
    S = {n: ids(f) for n, f in tools}
    res = {}
    for a, b in itertools.combinations([n for n, _ in tools], 2):
        u = S[a] | S[b]
        res[f"{a}|{b}"] = round(len(S[a] & S[b]) / len(u), 4) if u else None
    json.dump(res, open(out, "w"), indent=2)
    print("read-set Jaccard:", res)

if __name__ == "__main__":
    main()
