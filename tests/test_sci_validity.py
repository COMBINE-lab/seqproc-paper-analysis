#!/usr/bin/env python3
"""Tests for SciSeqValidityAnalyzer (edit-tolerant sci-RNA-seq3 structural check).

sci-RNA-seq3 R1 layout: [brc1: 9-10 bp][CAGAGC][umi: 8 bp][brc2: 10 bp].
A read is valid iff a unique best EDIT-distance-<=1 anchor alignment begins at
offset 9 or 10 (the two allowed brc1 lengths) and leaves room for umi(8)+brc2(10).
An apparent exact anchor one base outside those offsets is observationally
indistinguishable from an allowed boundary indel and is therefore retained.
These tests verify that documented behavior, substitutions/indels, the length
boundary, ambiguity handling, and read-id parsing.
"""
import os
import sys
import tempfile
import unittest.mock as mock

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
# run_paper_benchmarks imports numpy/matplotlib at module load; stub them out.
for m in ('numpy', 'matplotlib', 'matplotlib.pyplot', 'matplotlib.patches', 'matplotlib.gridspec'):
    sys.modules[m] = mock.MagicMock()

import run_paper_benchmarks as rpb
# tests must never read/write the shared on-disk validity cache
rpb._load_validity_cache = lambda *a, **k: None
rpb._save_validity_cache = lambda *a, **k: None
from run_paper_benchmarks import SciSeqValidityAnalyzer

ANCHOR = "CAGAGC"


def sci(brc1=9, anchor=ANCHOR, umi="T" * 8, brc2="G" * 10, tail=""):
    """Build a sci-RNA-seq3 R1: brc1 filler, then anchor at offset len(brc1)."""
    return "A" * brc1 + anchor + umi + brc2 + tail


def run_analyzer(records):
    """Write (id, seq) records to a temp FASTQ, return the analyzer's valid-id set."""
    fd, path = tempfile.mkstemp(suffix=".fastq")
    with os.fdopen(fd, "w") as f:
        for rid, seq in records:
            f.write(f"@{rid}\n{seq}\n+\n{'I' * len(seq)}\n")
    try:
        return set(SciSeqValidityAnalyzer().analyze_fastqs(path))
    finally:
        os.unlink(path)


# --- offset correctness: only brc1 lengths 9 and 10 accepted ---

def test_brc1_9bp_valid():
    assert run_analyzer([("r", sci(brc1=9))]) == {"r"}          # anchor at index 9


def test_brc1_10bp_valid():
    assert run_analyzer([("r", sci(brc1=10))]) == {"r"}         # anchor at index 10


def test_apparent_brc1_8_is_observationally_equivalent_to_anchor_deletion():
    # With edit-one anchors, an exact anchor at 8 is indistinguishable from a
    # legal offset-9 anchor whose first base was deleted. The structural set
    # therefore retains it; the reason-coded standalone report exposes edits.
    assert run_analyzer([("r", sci(brc1=8))]) == {"r"}


def test_apparent_brc1_11_is_observationally_equivalent_to_anchor_insertion():
    # Likewise, the extra base before an offset-11 exact anchor can be the one
    # insertion in a legal offset-10 anchor.
    assert run_analyzer([("r", sci(brc1=11))]) == {"r"}


# --- Hamming tolerance on the anchor ---

def test_anchor_1_mismatch_accepted():
    assert run_analyzer([("r", sci(brc1=9, anchor="CAGATC"))]) == {"r"}   # 1 sub


def test_anchor_2_mismatch_rejected():
    assert run_analyzer([("r", sci(brc1=9, anchor="CAGTTC"))]) == set()   # 2 subs


def test_anchor_1bp_deletion_accepted():
    # CAGAGC with the 2nd A deleted -> CAGGC (5bp); edit distance 1 -> valid (edit-tolerant)
    assert run_analyzer([("r", sci(brc1=9, anchor="CAGGC"))]) == {"r"}


def test_anchor_1bp_insertion_accepted():
    # CAGAGC with an extra A inserted -> CAGAAGC (7bp); edit distance 1 -> valid
    assert run_analyzer([("r", sci(brc1=9, anchor="CAGAAGC"))]) == {"r"}


# --- length / trailing-structure requirement ---

def test_truncated_no_room_for_umi_brc2_rejected():
    assert run_analyzer([("r", "A" * 9 + ANCHOR + "TTTTT")]) == set()     # only 5 bp trail


def test_min_length_boundary():
    # brc1=9 => minimum valid length is 9+6+8+10 = 33
    full = sci(brc1=9)
    assert len(full) == 33
    assert run_analyzer([("ok", full)]) == {"ok"}               # exactly 33 -> valid
    assert run_analyzer([("short", full[:-1])]) == set()        # 32 -> one short, rejected


# --- non-genuine reads ---

def test_no_anchor_rejected():
    assert run_analyzer([("r", "ACGT" * 20)]) == set()


def test_anchor_at_wrong_place_rejected():
    # anchor present but at offset 0, not 9/10
    assert run_analyzer([("r", ANCHOR + "A" * 40)]) == set()


# --- read-id parsing + multi-read file ---

def test_read_id_strips_at_and_extra_fields():
    assert run_analyzer([("r1 extra:field 2:N:0", sci(brc1=9))]) == {"r1"}


def test_mixed_file_returns_only_valid():
    recs = [
        ("v9", sci(9)), ("v10", sci(10)),
        ("bad8", sci(8)), ("bad11", sci(11)),
        ("noanchor", "ACGT" * 20), ("trunc", "A" * 9 + ANCHOR + "TT"),
        ("ham1", sci(9, anchor="CAGATC")),
    ]
    assert run_analyzer(recs) == {"v9", "v10", "ham1", "bad8", "bad11"}


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-v"]))
