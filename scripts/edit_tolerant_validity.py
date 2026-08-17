#!/usr/bin/env python3
"""Build reproducible conservative structural-reference read sets.

This program is intentionally independent of every benchmarked tool.  It
recognizes the four paper workloads, checks FASTQ integrity (and mate identity
when a second FASTQ is supplied), emits accepted read IDs in input order, and
writes a reason-coded JSON report.  The accepted set is a conservative
*structural reference*, not biological ground truth.

The implementation has two performance properties that matter on the paper's
large inputs:

* common exact structures take a bytes-only fast path before calling edlib;
* ``--threads`` validates ordered batches in worker processes, while accepted
  IDs are streamed instead of accumulating tens of millions of Python strings.

Examples
--------
  # Paired-end SPLiT-seq (R2 contains the barcode cassette)
  edit_tolerant_validity.py R2.fastq --r2 R1.fastq --chem pe --threads 32 \
      --out pe.valid_ids.txt --summary-json pe.summary.json

  # Dual-orientation long-read SPLiT-seq core definition
  edit_tolerant_validity.py reads.fastq --chem lr --threads 32 \
      --max-linker1-edit 3 --max-linker2-edit 3 --out lr.valid_ids.txt

  # 10x v2 and sci-RNA-seq3 paired structural checks
  edit_tolerant_validity.py R1.fastq --r2 R2.fastq --chem tenx_v2
  edit_tolerant_validity.py R1.fastq --r2 R2.fastq --chem scirnaseq3
"""

from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing as mp
import os
import resource
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterable, Iterator

import edlib


# The PE and LR protocols use different Linker2 sequences.  Linker1 is the
# documented A-containing sequence in both; the C-containing PE sequence in the
# original benchmark was a transcription error.
LINKERS = {
    "pe": ("GTGGCCGATGTTTCGCATCGGCGTACGACT", "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"),
    "lr": ("GTGGCCGATGTTTCGCATCGGCGTACGACT", "ATCCACGTGCTTGAGACTGTGG"),
}
SCI_ANCHOR = "CAGAGC"
_DNA = b"ACGT"
_COMP_BYTES = bytes.maketrans(
    b"ACGTNacgtn", b"TGCANtgcan"
)
_COMP_TEXT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_SCRIPT_DIR = Path(__file__).resolve().parent
_CONFIG_DIR = _SCRIPT_DIR.parent / "configs" / "seqproc"


def revcomp(seq: str | bytes) -> str | bytes:
    """Reverse-complement text or bytes without allocating a lookup dict."""
    if isinstance(seq, bytes):
        return seq.translate(_COMP_BYTES)[::-1]
    return seq.translate(_COMP_TEXT)[::-1]


def _read_whitelist(path: str | os.PathLike[str]) -> tuple[bytes, ...]:
    entries: list[bytes] = []
    with open(path, "rb") as handle:
        for raw in handle:
            raw = raw.strip()
            if not raw or raw.startswith(b"#"):
                continue
            # Accept one-column lists, seq-to-seq TSVs, and simple CSVs.
            fields = raw.replace(b",", b"\t").split(b"\t")
            candidate = fields[-1].strip().upper()
            if candidate and all(base in _DNA for base in candidate):
                entries.append(candidate)
    if not entries:
        raise ValueError(f"whitelist is empty: {path}")
    lengths = {len(entry) for entry in entries}
    if len(lengths) != 1:
        raise ValueError(f"whitelist contains mixed lengths {sorted(lengths)}: {path}")
    return tuple(dict.fromkeys(entries))


def _hamming_owner_map(
    barcodes: Iterable[bytes], max_distance: int = 1
) -> dict[bytes, tuple[bytes, ...]]:
    """Map every Hamming-neighbor to its canonical owner(s).

    The paper uses distance one, for which direct expansion is much faster than
    scanning a whitelist per read.  Retaining owners (rather than only set
    membership) lets the report distinguish unique and ambiguous corrections.
    """
    if max_distance != 1:
        raise ValueError("only Hamming distance one is supported")
    owners: dict[bytes, set[bytes]] = defaultdict(set)
    for barcode in barcodes:
        owners[barcode].add(barcode)
        mutable = bytearray(barcode)
        for index, original in enumerate(barcode):
            for base in _DNA:
                if base != original:
                    mutable[index] = base
                    owners[bytes(mutable)].add(barcode)
            mutable[index] = original
    return {observed: tuple(sorted(values)) for observed, values in owners.items()}


def ham1_set(path):
    """Backward-compatible helper used by earlier tests and notebooks."""
    barcodes = _read_whitelist(path)
    owners = _hamming_owner_map(barcodes)
    return {item.decode() for item in owners}, len(barcodes[0])


def find_linker_candidates(
    query: str | bytes,
    target: str | bytes,
    max_edit: int,
    max_candidates: int | None = None,
):
    """Enumerate distinct non-overlapping linker placements within ``max_edit``.

    edlib returns globally optimal locations.  Recursively searching both sides
    prevents a better decoy elsewhere in a long read from hiding a slightly
    worse genuine cassette.  The default has no semantic candidate cap.  A
    caller-supplied cap is retained only for sensitivity/compatibility runs.
    """
    minimum_span = max(1, len(query) - max_edit)
    pending = [(0, len(target))]
    candidates = []
    while pending and (max_candidates is None or len(candidates) < max_candidates):
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
        # Record every tied optimum returned for this interval, then subdivide
        # around their union.  Genuine cassettes cannot contain overlapping L1s.
        locations = sorted(set(result["locations"]))
        leftmost = end
        rightmost = begin - 1
        for start, stop in locations:
            absolute_start = start + begin
            absolute_stop = stop + begin
            candidates.append((result["editDistance"], absolute_start, absolute_stop))
            leftmost = min(leftmost, absolute_start)
            rightmost = max(rightmost, absolute_stop)
            if max_candidates is not None and len(candidates) >= max_candidates:
                break
        if leftmost - begin >= minimum_span:
            pending.append((begin, leftmost))
        if end - (rightmost + 1) >= minimum_span:
            pending.append((rightmost + 1, end))
    return sorted(set(candidates))


def prefix_linker_matches(query, target, max_edit):
    """Return optimal linker matches constrained to target offset zero."""
    window = target[: len(query) + max_edit]
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


def _barcode_status(observed, whitelist) -> tuple[bool, bool]:
    owners = whitelist.get(observed)
    return owners is not None, owners is not None and len(owners) > 1


def _as_owner_map(values) -> dict:
    """Accept both the new owner map and the old expanded-set API."""
    if hasattr(values, "get"):
        return values
    return {value: (value,) for value in values}


def _accepted_split_result(
    linker1_edit: int,
    linker2_edit: int,
    linker1_start: int,
    ambiguous_barcodes: int,
    fast_path: bool,
):
    return {
        "accepted": True,
        "reason": "accepted_ambiguous_barcode" if ambiguous_barcodes else "accepted",
        "linker1_edit": linker1_edit,
        "linker2_edit": linker2_edit,
        "linker1_start": linker1_start,
        "ambiguous_barcodes": ambiguous_barcodes,
        "fast_path": fast_path,
    }


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
    max_candidates=None,
):
    """Validate one SPLiT-seq orientation and return a reason-coded result.

    The exact-layout branch handles the overwhelming majority of PE reads and a
    useful fraction of HiFi reads without entering the recursive edlib search.
    It is semantic-preserving: any failure falls through to the general path.
    """
    bc23_whitelist = _as_owner_map(bc23_whitelist)
    bc1_whitelist = _as_owner_map(bc1_whitelist)
    text_mode = isinstance(seq, str)
    if text_mode:
        seq = seq.encode()
        linker1 = linker1.encode()
        linker2 = linker2.encode()
        bc23_whitelist = {
            (key.encode() if isinstance(key, str) else key): value
            for key, value in bc23_whitelist.items()
        }
        bc1_whitelist = {
            (key.encode() if isinstance(key, str) else key): value
            for key, value in bc1_whitelist.items()
        }
    elif isinstance(linker1, str):
        linker1, linker2 = linker1.encode(), linker2.encode()

    expected_l1_start = umi_length + bc23_length
    expected_bc2_start = expected_l1_start + len(linker1)
    expected_l2_start = expected_bc2_start + bc23_length
    expected_bc1_start = expected_l2_start + len(linker2)
    if (
        len(seq) >= expected_bc1_start + bc1_length
        and seq[expected_l1_start:expected_bc2_start] == linker1
        and seq[expected_l2_start:expected_bc1_start] == linker2
    ):
        bc3_ok, bc3_ambig = _barcode_status(
            seq[umi_length:expected_l1_start], bc23_whitelist
        )
        bc2_ok, bc2_ambig = _barcode_status(
            seq[expected_bc2_start:expected_l2_start], bc23_whitelist
        )
        bc1_ok, bc1_ambig = _barcode_status(
            seq[expected_bc1_start:expected_bc1_start + bc1_length], bc1_whitelist
        )
        if bc3_ok and bc2_ok and bc1_ok:
            return _accepted_split_result(
                0, 0, expected_l1_start, sum((bc3_ambig, bc2_ambig, bc1_ambig)), True
            )

    candidates = find_linker_candidates(
        linker1, seq, max_linker1_edit, max_candidates
    )
    if not candidates:
        return {"accepted": False, "reason": "no_linker1", "fast_path": False}

    furthest_reason = "incomplete_umi_bc3_prefix"
    accepted = []
    for linker1_edit, linker1_start, linker1_stop in candidates:
        if linker1_start < umi_length + bc23_length:
            continue
        bc3 = seq[linker1_start - bc23_length:linker1_start]
        bc3_ok, bc3_ambig = _barcode_status(bc3, bc23_whitelist)
        if not bc3_ok:
            furthest_reason = "invalid_bc3"
            continue

        bc2_start = linker1_stop + 1
        bc2_stop = bc2_start + bc23_length
        bc2 = seq[bc2_start:bc2_stop]
        if len(bc2) != bc23_length:
            furthest_reason = "incomplete_bc2"
            continue
        bc2_ok, bc2_ambig = _barcode_status(bc2, bc23_whitelist)
        if not bc2_ok:
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
            bc1_ok, bc1_ambig = _barcode_status(bc1, bc1_whitelist)
            if not bc1_ok:
                furthest_reason = "invalid_bc1"
                continue
            accepted.append(
                _accepted_split_result(
                    linker1_edit,
                    linker2_edit,
                    linker1_start,
                    sum((bc3_ambig, bc2_ambig, bc1_ambig)),
                    False,
                )
            )
    if not accepted:
        return {"accepted": False, "reason": furthest_reason, "fast_path": False}
    accepted.sort(
        key=lambda result: (
            result["linker1_edit"] + result["linker2_edit"],
            result["ambiguous_barcodes"],
            result["linker1_start"],
        )
    )
    best = accepted[0]
    best_key = (
        best["linker1_edit"] + best["linker2_edit"],
        best["ambiguous_barcodes"],
    )
    tied = sum(
        (
            item["linker1_edit"] + item["linker2_edit"],
            item["ambiguous_barcodes"],
        ) == best_key
        for item in accepted
    )
    best["cassette_tie_count"] = tied
    if tied > 1:
        best["reason"] = "accepted_ambiguous_cassette"
    return best


@dataclass(frozen=True)
class Policy:
    chem: str
    linker1: bytes | None = None
    linker2: bytes | None = None
    bc23: dict[bytes, tuple[bytes, ...]] | None = None
    bc1: dict[bytes, tuple[bytes, ...]] | None = None
    bc23_length: int = 0
    bc1_length: int = 0
    max_linker1_edit: int = 0
    max_linker2_edit: int = 0
    max_candidates: int | None = None
    orientation: str = "forward"


def _splitseq_result(seq: bytes, policy: Policy) -> dict:
    kwargs = (
        policy.linker1,
        policy.linker2,
        policy.bc23,
        policy.bc1,
        policy.bc23_length,
        policy.bc1_length,
        policy.max_linker1_edit,
        policy.max_linker2_edit,
    )
    forward = validate_orientation(
        seq, *kwargs, max_candidates=policy.max_candidates
    )
    if policy.orientation == "forward":
        if forward["accepted"]:
            forward["orientation"] = "forward"
        return forward
    reverse = validate_orientation(
        revcomp(seq), *kwargs, max_candidates=policy.max_candidates
    )
    if forward["accepted"] and reverse["accepted"]:
        chosen = min(
            (forward, reverse),
            key=lambda item: (
                item["linker1_edit"] + item["linker2_edit"],
                item.get("ambiguous_barcodes", 0),
            ),
        )
        chosen = dict(chosen)
        chosen["orientation"] = "both"
        chosen["reason"] = "accepted_both_orientations"
        return chosen
    if forward["accepted"]:
        forward["orientation"] = "forward"
        return forward
    if reverse["accepted"]:
        reverse["orientation"] = "reverse"
        reverse["forward_failure_reason"] = forward["reason"]
        return reverse
    return {
        "accepted": False,
        "reason": reverse["reason"],
        "forward_failure_reason": forward["reason"],
        "fast_path": False,
    }


def validate_scirnaseq3(seq: bytes) -> dict:
    """Validate the two documented BC1 lengths without a global-best decoy bias."""
    candidates = []
    for start in (9, 10):
        suffix = seq[start:start + len(SCI_ANCHOR) + 1]
        if len(suffix) < len(SCI_ANCHOR) - 1:
            continue
        # Exact comparison is both common and substantially cheaper than edlib.
        if seq[start:start + len(SCI_ANCHOR)] == SCI_ANCHOR.encode():
            candidates.append((0, start, start + len(SCI_ANCHOR) - 1, True))
            continue
        for edit, match_start, stop in prefix_linker_matches(
            SCI_ANCHOR.encode(), suffix, 1
        ):
            if match_start == 0:
                candidates.append((edit, start, start + stop, False))
    if not candidates:
        return {"accepted": False, "reason": "no_anchor_at_allowed_offset"}
    complete = [item for item in candidates if item[2] + 1 + 18 <= len(seq)]
    if not complete:
        return {"accepted": False, "reason": "incomplete_umi_bc2"}
    best_edit = min(item[0] for item in complete)
    best = sorted({item for item in complete if item[0] == best_edit})
    best_starts = {item[1] for item in best}
    if len(best_starts) != 1:
        return {
            "accepted": False,
            "reason": "ambiguous_anchor_offset",
            "anchor_edit": best_edit,
        }
    item = best[0]
    return {
        "accepted": True,
        "reason": "accepted",
        "anchor_edit": item[0],
        "anchor_start": item[1],
        "fast_path": item[3],
    }


def _validate_record(seq: bytes, policy: Policy) -> dict:
    if policy.chem in ("pe", "lr"):
        return _splitseq_result(seq, policy)
    if policy.chem == "scirnaseq3":
        return validate_scirnaseq3(seq)
    if policy.chem == "tenx_v2":
        if len(seq) < 26:
            return {"accepted": False, "reason": "r1_shorter_than_26"}
        return {"accepted": True, "reason": "accepted", "fast_path": True}
    raise ValueError(f"unsupported chemistry: {policy.chem}")


def _normalize_read_id(header: bytes) -> bytes:
    read_id = header[1:].split(None, 1)[0]
    if read_id.endswith((b"/1", b"/2")):
        return read_id[:-2]
    return read_id


def _open_fastq(path: str) -> BinaryIO:
    if path == "-":
        return sys.stdin.buffer
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb", buffering=1024 * 1024)


def _fastq_records(handle: BinaryIO, label: str) -> Iterator[tuple[bytes, bytes]]:
    record = 0
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        record += 1
        if not sequence or not plus or not quality:
            raise ValueError(f"truncated FASTQ record {record} in {label}")
        header = header.rstrip(b"\r\n")
        sequence = sequence.rstrip(b"\r\n").upper()
        plus = plus.rstrip(b"\r\n")
        quality = quality.rstrip(b"\r\n")
        if not header.startswith(b"@") or not plus.startswith(b"+"):
            raise ValueError(f"malformed FASTQ record {record} in {label}")
        if len(sequence) != len(quality):
            raise ValueError(
                f"sequence/quality length mismatch at record {record} in {label}"
            )
        yield _normalize_read_id(header), sequence


def _paired_records(
    fastq: str, mate_fastq: str | None
) -> Iterator[tuple[bytes, bytes]]:
    with _open_fastq(fastq) as first:
        first_records = _fastq_records(first, fastq)
        if mate_fastq is None:
            yield from first_records
            return
        with _open_fastq(mate_fastq) as second:
            second_records = _fastq_records(second, mate_fastq)
            index = 0
            while True:
                left = next(first_records, None)
                right = next(second_records, None)
                if left is None and right is None:
                    return
                index += 1
                if left is None or right is None:
                    raise ValueError(
                        f"paired FASTQs contain different record counts near record {index}"
                    )
                if left[0] != right[0]:
                    raise ValueError(
                        "paired FASTQ ID mismatch at record "
                        f"{index}: {left[0].decode(errors='replace')} != "
                        f"{right[0].decode(errors='replace')}"
                    )
                yield left


_WORKER_POLICY: Policy | None = None
_WORKER_RETURN_IDS = False


def _init_worker(policy: Policy, return_ids: bool):
    global _WORKER_POLICY, _WORKER_RETURN_IDS
    _WORKER_POLICY = policy
    _WORKER_RETURN_IDS = return_ids


def _summarize_result(counters: Counter, result: dict):
    counters[f"outcome:{result['reason']}"] += 1
    if result.get("accepted"):
        counters["accepted"] += 1
        counters[f"orientation:{result.get('orientation', 'forward')}"] += 1
        counters[f"fast_path:{bool(result.get('fast_path'))}"] += 1
        if "linker1_edit" in result:
            counters[f"linker1_edit:{result['linker1_edit']}"] += 1
            counters[f"linker2_edit:{result['linker2_edit']}"] += 1
            counters[f"ambiguous_barcodes:{result.get('ambiguous_barcodes', 0)}"] += 1
        if "anchor_edit" in result:
            counters[f"anchor_edit:{result['anchor_edit']}"] += 1
            counters[f"anchor_start:{result['anchor_start']}"] += 1
    if "forward_failure_reason" in result:
        counters[f"forward_failure:{result['forward_failure_reason']}"] += 1


def _validate_batch(batch: list[tuple[bytes, bytes]]):
    assert _WORKER_POLICY is not None
    counters = Counter()
    accepted_ids = []
    for read_id, sequence in batch:
        result = _validate_record(sequence, _WORKER_POLICY)
        _summarize_result(counters, result)
        if result["accepted"] and _WORKER_RETURN_IDS:
            accepted_ids.append(read_id)
    return len(batch), accepted_ids, counters


def _batches(records: Iterable[tuple[bytes, bytes]], size: int, limit: int):
    batch = []
    seen = 0
    for record in records:
        if limit and seen >= limit:
            break
        batch.append(record)
        seen += 1
        if len(batch) == size:
            yield batch
            batch = []
    if batch:
        yield batch


def _counter_group(counter: Counter, prefix: str) -> dict[str, int]:
    marker = prefix + ":"
    return {
        key[len(marker):]: value
        for key, value in sorted(counter.items())
        if key.startswith(marker)
    }


def _make_policy(args: argparse.Namespace) -> Policy:
    if args.chem not in ("pe", "lr"):
        return Policy(chem=args.chem)
    linker1, linker2 = (item.encode() for item in LINKERS[args.chem])
    bc23_path = args.bc23_whitelist or str(
        _CONFIG_DIR / "splitseq_bc8_whitelist.txt"
    )
    if args.bc1_whitelist:
        bc1_path = args.bc1_whitelist
    elif args.chem == "pe":
        # PE round-one barcodes are the full 8-nt sequences.  The original
        # paper accidentally used their 6-nt LR truncations.
        bc1_path = bc23_path
    else:
        bc1_path = str(_CONFIG_DIR / "splitseq_bc1_whitelist_6bp.txt")
    bc23_values = _read_whitelist(bc23_path)
    bc1_values = _read_whitelist(bc1_path)
    max1 = args.max_linker1_edit
    max2 = args.max_linker2_edit
    if max1 is None:
        max1 = args.max_linker_edit if args.max_linker_edit is not None else 3
    if max2 is None:
        max2 = args.max_linker_edit if args.max_linker_edit is not None else 3
    orientation = args.orientation
    if orientation is None:
        orientation = "both" if args.chem == "lr" else "forward"
    return Policy(
        chem=args.chem,
        linker1=linker1,
        linker2=linker2,
        bc23=_hamming_owner_map(bc23_values),
        bc1=_hamming_owner_map(bc1_values),
        bc23_length=len(bc23_values[0]),
        bc1_length=len(bc1_values[0]),
        max_linker1_edit=max1,
        max_linker2_edit=max2,
        max_candidates=args.max_linker1_candidates,
        orientation=orientation,
    )


def _criteria(policy: Policy, paired: bool) -> dict:
    common = {
        "paired_fastq_integrity_checked": paired,
        "reference_kind": "conservative_structural_reference_not_ground_truth",
    }
    if policy.chem in ("pe", "lr"):
        common.update(
            {
                "umi_length": 10,
                "bc23_length": policy.bc23_length,
                "bc1_length": policy.bc1_length,
                "barcode_hamming_distance": 1,
                "ambiguous_barcode_membership": "accepted_and_reported",
                "max_linker1_edit": policy.max_linker1_edit,
                "max_linker2_edit": policy.max_linker2_edit,
                "max_linker1_candidates": policy.max_candidates,
                "linker2_expected_immediately_after_bc2": True,
                "orientation": policy.orientation,
            }
        )
    elif policy.chem == "scirnaseq3":
        common.update(
            {
                "allowed_anchor_offsets": [9, 10],
                "anchor": SCI_ANCHOR,
                "max_anchor_edit": 1,
                "minimum_trailing_umi_bc2_length": 18,
                "ambiguous_best_anchor_offset": "rejected_and_reported",
            }
        )
    else:
        common.update({"minimum_r1_length": 26})
    return common


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq", help="FASTQ containing the validated read")
    parser.add_argument(
        "--r2",
        help=(
            "optional mate FASTQ used for record-count and read-ID integrity checks; "
            "the first FASTQ remains the validated read"
        ),
    )
    parser.add_argument(
        "--chem",
        choices=["pe", "lr", "tenx_v2", "scirnaseq3"],
        default="pe",
    )
    parser.add_argument("--out")
    parser.add_argument("--summary-json")
    parser.add_argument("--sample", type=int, default=0)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument(
        "--max-linker-edit",
        type=int,
        default=None,
        help="compatibility option applied to both linkers (new default: 3)",
    )
    parser.add_argument("--max-linker1-edit", type=int)
    parser.add_argument("--max-linker2-edit", type=int)
    parser.add_argument(
        "--max-linker1-candidates",
        type=int,
        default=None,
        help="optional sensitivity cap; the primary definition is uncapped",
    )
    parser.add_argument("--bc23-whitelist")
    parser.add_argument("--bc1-whitelist")
    parser.add_argument("--orientation", choices=["forward", "both"])
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    if args.batch_size < 1:
        parser.error("--batch-size must be at least 1")
    if args.fastq == "-" and args.r2:
        parser.error("paired validation does not support stdin")
    return args


def main(argv=None):
    args = parse_args(argv)
    policy = _make_policy(args)
    started = time.perf_counter()
    counters = Counter()
    total = 0
    out_handle = open(args.out, "wb", buffering=1024 * 1024) if args.out else None
    try:
        batches = _batches(
            _paired_records(args.fastq, args.r2), args.batch_size, args.sample
        )
        if args.threads == 1:
            _init_worker(policy, out_handle is not None)
            results = map(_validate_batch, batches)
            for count, accepted_ids, batch_counts in results:
                total += count
                counters.update(batch_counts)
                if out_handle and accepted_ids:
                    out_handle.write(b"\n".join(accepted_ids) + b"\n")
        else:
            # imap preserves input order, so the streamed ID file is stable
            # across thread counts without retaining the whole reference set.
            with mp.Pool(
                args.threads,
                initializer=_init_worker,
                initargs=(policy, out_handle is not None),
            ) as pool:
                for count, accepted_ids, batch_counts in pool.imap(
                    _validate_batch, batches, chunksize=1
                ):
                    total += count
                    counters.update(batch_counts)
                    if out_handle and accepted_ids:
                        out_handle.write(b"\n".join(accepted_ids) + b"\n")
    finally:
        if out_handle:
            out_handle.close()

    elapsed = time.perf_counter() - started
    own_usage = resource.getrusage(resource.RUSAGE_SELF)
    child_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    valid = counters["accepted"]
    summary = {
        "schema_version": "3.0.0",
        "fastq": os.path.abspath(args.fastq) if args.fastq != "-" else "stdin",
        "mate_fastq": os.path.abspath(args.r2) if args.r2 else None,
        "chem": args.chem,
        "total": total,
        "valid": valid,
        "invalid": total - valid,
        "pct_of_scanned": round(100 * valid / total, 4) if total else 0.0,
        "elapsed_seconds": round(elapsed, 6),
        "reads_per_second": round(total / elapsed, 2) if elapsed else None,
        "user_cpu_seconds": round(own_usage.ru_utime + child_usage.ru_utime, 6),
        "system_cpu_seconds": round(own_usage.ru_stime + child_usage.ru_stime, 6),
        "controller_peak_rss_kib": own_usage.ru_maxrss,
        "maximum_worker_peak_rss_kib": child_usage.ru_maxrss,
        "threads": args.threads,
        "batch_size": args.batch_size,
        "criteria": _criteria(policy, args.r2 is not None),
        "accepted_orientation": _counter_group(counters, "orientation"),
        "linker1_edit_histogram": _counter_group(counters, "linker1_edit"),
        "linker2_edit_histogram": _counter_group(counters, "linker2_edit"),
        "anchor_edit_histogram": _counter_group(counters, "anchor_edit"),
        "anchor_start_histogram": _counter_group(counters, "anchor_start"),
        "barcode_ambiguity_histogram": _counter_group(
            counters, "ambiguous_barcodes"
        ),
        "fast_path_counts": _counter_group(counters, "fast_path"),
        "outcome_counts": _counter_group(counters, "outcome"),
        "forward_failure_counts": _counter_group(counters, "forward_failure"),
    }
    rendered = json.dumps(summary, indent=2, sort_keys=True)
    print(rendered)
    if args.summary_json:
        with open(args.summary_json, "w") as output:
            output.write(rendered + "\n")
    return summary


if __name__ == "__main__":
    main()
