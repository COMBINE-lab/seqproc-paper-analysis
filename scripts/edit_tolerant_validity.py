#!/usr/bin/env python3
"""Build a conservative, edit-tolerant SPLiT-seq structural reference.

Distinct plausible Linker1 placements are enumerated rather than accepting the
first globally best hit. A placement is valid only when the complete 10-nt UMI
and 8-nt BC3 precede Linker1, a complete 8-nt BC2 separates the linkers,
Linker2 begins at the expected post-BC2 offset, and the complete 6-nt BC1
follows Linker2. The three barcodes must each be within Hamming distance 1 of
their whitelist. This is a conservative structural reference, not experimental
ground truth.

Usage:
  edit_tolerant_validity.py R2.fastq --out valid_ids.txt [--chem pe|lr]
"""

import argparse
import json
import os
from collections import Counter

import edlib


LINKERS = {
    "pe": ("GTGGCCGCTGTTTCGCATCGGCGTACGACT", "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"),
    "lr": ("GTGGCCGATGTTTCGCATCGGCGTACGACT", "ATCCACGTGCTTGAGACTGTGG"),
}
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def ham1_set(path):
    barcodes = [line.strip() for line in open(path) if line.strip()]
    expanded = set(barcodes)
    for barcode in barcodes:
        for index in range(len(barcode)):
            for base in "ACGT":
                if base != barcode[index]:
                    expanded.add(barcode[:index] + base + barcode[index + 1:])
    return expanded, (len(barcodes[0]) if barcodes else 0)


def find_linker_candidates(query, target, max_edit, max_candidates=8):
    """Enumerate distinct non-overlapping linker placements within max_edit.

    Edlib reports only globally optimal locations. After recording a hit, this
    searches the sequence on each side of it, allowing a slightly worse genuine
    cassette to be considered when a better decoy linker occurs elsewhere.
    """
    minimum_span = max(1, len(query) - max_edit)
    pending = [(0, len(target))]
    candidates = []
    while pending and len(candidates) < max_candidates:
        begin, end = pending.pop()
        if end - begin < minimum_span:
            continue
        result = edlib.align(
            query,
            target[begin:end],
            mode="HW",
            task="locations",
            k=max_edit,
        )
        if result["editDistance"] < 0 or not result["locations"]:
            continue
        start, stop = result["locations"][0]
        start += begin
        stop += begin
        candidates.append((result["editDistance"], start, stop))
        if start - begin >= minimum_span:
            pending.append((begin, start))
        if end - (stop + 1) >= minimum_span:
            pending.append((stop + 1, end))
    return sorted(set(candidates))


def prefix_linker_matches(query, target, max_edit):
    """Return optimal linker matches constrained to target offset zero."""
    window = target[:len(query) + max_edit]
    if len(window) < len(query) - max_edit:
        return []
    result = edlib.align(query, window, mode="SHW", task="locations", k=max_edit)
    if result["editDistance"] < 0:
        return []
    return [
        (result["editDistance"], start, stop)
        for start, stop in result["locations"]
        if start == 0
    ]


def validate_orientation(
    seq,
    linker1,
    linker2,
    bc23_whitelist,
    bc1_whitelist,
    bc23_length,
    bc1_length,
    max_linker1_edit,
    max_linker2_edit,
    umi_length=10,
    max_candidates=8,
):
    """Validate one orientation and return a reason-coded result."""
    candidates = find_linker_candidates(
        linker1, seq, max_linker1_edit, max_candidates
    )
    if not candidates:
        return {"accepted": False, "reason": "no_linker1"}

    furthest_reason = "incomplete_umi_bc3_prefix"
    for linker1_edit, linker1_start, linker1_stop in candidates:
        if linker1_start < umi_length + bc23_length:
            continue
        bc3 = seq[linker1_start - bc23_length:linker1_start]
        if bc3 not in bc23_whitelist:
            furthest_reason = "invalid_bc3"
            continue

        bc2_start = linker1_stop + 1
        bc2_stop = bc2_start + bc23_length
        bc2 = seq[bc2_start:bc2_stop]
        if len(bc2) != bc23_length:
            furthest_reason = "incomplete_bc2"
            continue
        if bc2 not in bc23_whitelist:
            furthest_reason = "invalid_bc2"
            continue

        linker2_matches = prefix_linker_matches(
            linker2, seq[bc2_stop:], max_linker2_edit
        )
        if not linker2_matches:
            furthest_reason = "no_linker2_at_expected_offset"
            continue
        for linker2_edit, _, linker2_stop in linker2_matches:
            bc1_start = bc2_stop + linker2_stop + 1
            bc1 = seq[bc1_start:bc1_start + bc1_length]
            if len(bc1) != bc1_length:
                furthest_reason = "incomplete_bc1"
                continue
            if bc1 not in bc1_whitelist:
                furthest_reason = "invalid_bc1"
                continue
            return {
                "accepted": True,
                "reason": "accepted",
                "linker1_edit": linker1_edit,
                "linker2_edit": linker2_edit,
                "linker1_start": linker1_start,
            }
    return {"accepted": False, "reason": furthest_reason}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq")
    parser.add_argument(
        "--chem",
        choices=["pe", "lr"],
        default="pe",
        help="pe: check forward only; lr: check both orientations",
    )
    parser.add_argument("--out", default=None)
    parser.add_argument("--summary-json", default=None)
    parser.add_argument("--sample", type=int, default=0)
    parser.add_argument(
        "--max-linker-edit",
        type=int,
        default=6,
        help="backward-compatible default applied to both linkers",
    )
    parser.add_argument("--max-linker1-edit", type=int, default=None)
    parser.add_argument("--max-linker2-edit", type=int, default=None)
    parser.add_argument("--max-linker1-candidates", type=int, default=8)
    args = parser.parse_args()

    linker1, linker2 = LINKERS[args.chem]
    whitelist_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "configs", "seqproc"
    )
    bc23_whitelist, bc23_length = ham1_set(
        os.path.join(whitelist_dir, "splitseq_bc23_whitelist.txt")
    )
    bc1_whitelist, bc1_length = ham1_set(
        os.path.join(whitelist_dir, "splitseq_bc1_whitelist_6bp.txt")
    )
    max_linker1_edit = (
        args.max_linker1_edit
        if args.max_linker1_edit is not None
        else args.max_linker_edit
    )
    max_linker2_edit = (
        args.max_linker2_edit
        if args.max_linker2_edit is not None
        else args.max_linker_edit
    )

    valid = set()
    total = 0
    reasons = Counter()
    forward_failure_reasons = Counter()
    orientations = Counter()
    linker1_distances = Counter()
    linker2_distances = Counter()
    with open(args.fastq) as fastq:
        while True:
            header = fastq.readline()
            if not header:
                break
            seq = fastq.readline().strip()
            fastq.readline()
            fastq.readline()
            if args.sample and total >= args.sample:
                break
            total += 1
            read_id = header[1:].split()[0]
            result = validate_orientation(
                seq,
                linker1,
                linker2,
                bc23_whitelist,
                bc1_whitelist,
                bc23_length,
                bc1_length,
                max_linker1_edit,
                max_linker2_edit,
                max_candidates=args.max_linker1_candidates,
            )
            orientation = "forward"
            if not result["accepted"] and args.chem == "lr":
                forward_failure_reasons[result["reason"]] += 1
                result = validate_orientation(
                    revcomp(seq),
                    linker1,
                    linker2,
                    bc23_whitelist,
                    bc1_whitelist,
                    bc23_length,
                    bc1_length,
                    max_linker1_edit,
                    max_linker2_edit,
                    max_candidates=args.max_linker1_candidates,
                )
                orientation = "reverse"
            reasons[result["reason"]] += 1
            if result["accepted"]:
                valid.add(read_id)
                orientations[orientation] += 1
                linker1_distances[result["linker1_edit"]] += 1
                linker2_distances[result["linker2_edit"]] += 1

    summary = {
        "schema_version": "2.0.0",
        "fastq": os.path.basename(args.fastq),
        "chem": args.chem,
        "total": total,
        "valid": len(valid),
        "pct_of_scanned": round(100 * len(valid) / total, 2) if total else 0.0,
        "criteria": {
            "umi_length": 10,
            "bc23_length": bc23_length,
            "bc1_length": bc1_length,
            "max_linker1_edit": max_linker1_edit,
            "max_linker2_edit": max_linker2_edit,
            "max_linker1_candidates": args.max_linker1_candidates,
            "linker2_expected_immediately_after_bc2": True,
        },
        "accepted_orientation": dict(sorted(orientations.items())),
        "linker1_edit_histogram": dict(sorted(linker1_distances.items())),
        "linker2_edit_histogram": dict(sorted(linker2_distances.items())),
        "outcome_counts": dict(sorted(reasons.items())),
        "forward_failure_counts": dict(sorted(forward_failure_reasons.items())),
    }
    print(json.dumps(summary, sort_keys=True))
    if args.summary_json:
        with open(args.summary_json, "w") as output:
            json.dump(summary, output, indent=2, sort_keys=True)
            output.write("\n")
    if args.out:
        with open(args.out, "w") as output:
            output.write("\n".join(sorted(valid)))
        print("wrote", args.out)


if __name__ == "__main__":
    main()
