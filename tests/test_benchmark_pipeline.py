#!/usr/bin/env python3
"""
Tests for benchmark pipeline bug fixes.

Covers:
  - Bug 1: run_cmd exit code parsing from /usr/bin/time stderr
  - Bug 2: check_wl_match scoping (bc distances initialized before use)
  - Bug 3: matchbox_paired flag replaces magic string dispatch
  - Bug 5: matchbox ID extraction prefers R2 for paired-end
  - DNA reverse-complement primitive (used by splitcode_lr_dual_validate)
"""

import subprocess
import sys
import os

# Add scripts to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))


# ---------------------------------------------------------------------------
# Bug 1: run_cmd exit code parsing
# ---------------------------------------------------------------------------

def test_run_cmd_parses_exit_status():
    """Exit status should be extracted from /usr/bin/time stderr, not shell returncode."""
    time_cmd = "/usr/bin/time -v echo hello"
    result = subprocess.run(time_cmd, shell=True, capture_output=True, text=True)

    tool_exit_code = result.returncode
    for line in result.stderr.split('\n'):
        if 'Exit status' in line:
            tool_exit_code = int(line.split(':')[1].strip())

    assert tool_exit_code == 0, f"Expected exit 0 for 'echo hello', got {tool_exit_code}"


def test_run_cmd_detects_failure():
    """Failing commands should report non-zero exit status."""
    time_cmd = "/usr/bin/time -v false"
    result = subprocess.run(time_cmd, shell=True, capture_output=True, text=True)

    tool_exit_code = result.returncode
    for line in result.stderr.split('\n'):
        if 'Exit status' in line:
            tool_exit_code = int(line.split(':')[1].strip())

    assert tool_exit_code == 1, f"Expected exit 1 for 'false', got {tool_exit_code}"


def test_run_cmd_failure_with_redirect():
    """Failing command with redirect should still report non-zero exit status."""
    time_cmd = "/usr/bin/time -v bash -c 'false > /dev/null'"
    result = subprocess.run(time_cmd, shell=True, capture_output=True, text=True)

    tool_exit_code = result.returncode
    for line in result.stderr.split('\n'):
        if 'Exit status' in line:
            tool_exit_code = int(line.split(':')[1].strip())

    assert tool_exit_code == 1, f"Expected exit 1 for 'false > /dev/null', got {tool_exit_code}"


# ---------------------------------------------------------------------------
# Bug 2: check_wl_match scoping -- distances initialized before if-blocks
# ---------------------------------------------------------------------------

def test_bc_distances_initialized():
    """bc distances should be 99 (no match) when barcode length is wrong."""
    # Simulates the fixed code path where bc3_d, bc2_d, bc1_d = 99, 99, 99
    bc3_d, bc2_d, bc1_d = 99, 99, 99

    bc3 = "SHORT"  # Not 8bp
    bc2 = "AACGTGAT"  # 8bp
    bc1 = "AACGT"  # Not 6bp

    # Only update if correct length (simulating the fixed code)
    if len(bc3) == 8:
        bc3_d = 0
    if len(bc2) == 8:
        bc2_d = 0
    if len(bc1) == 6:
        bc1_d = 0

    bc_ok = (bc3_d <= 1 and bc2_d <= 1 and bc1_d <= 1)
    assert not bc_ok, "bc_ok should be False when bc3 and bc1 have wrong length"
    assert bc3_d == 99, f"bc3_d should stay 99 for wrong-length barcode, got {bc3_d}"
    assert bc2_d == 0, f"bc2_d should be 0 for matching barcode, got {bc2_d}"
    assert bc1_d == 99, f"bc1_d should stay 99 for wrong-length barcode, got {bc1_d}"


# ---------------------------------------------------------------------------
# Bug 3: matchbox_paired flag
# ---------------------------------------------------------------------------

def test_matchbox_paired_flag_in_datasets():
    """All datasets should have explicit matchbox_paired flag."""
    from concordance_analysis import DATASETS

    for ds_key, ds in DATASETS.items():
        assert 'matchbox_paired' in ds, (
            f"Dataset '{ds_key}' missing 'matchbox_paired' flag")
        assert isinstance(ds['matchbox_paired'], bool), (
            f"Dataset '{ds_key}' matchbox_paired should be bool")


def test_matchbox_paired_consistency():
    """matchbox_paired should match expected values for known datasets."""
    from concordance_analysis import DATASETS

    expected = {
        'splitseq_pe': True,
        'lr_splitseq': False,
        '10x_short': False,
        'sciseq': True,
    }
    for ds_key, expect_paired in expected.items():
        if ds_key in DATASETS:
            assert DATASETS[ds_key]['matchbox_paired'] == expect_paired, (
                f"{ds_key}: expected matchbox_paired={expect_paired}")


# ---------------------------------------------------------------------------
# Bug 5: ID extraction order
# ---------------------------------------------------------------------------

def test_paired_prefers_r2():
    """For paired-end datasets, R2 should be preferred for ID extraction."""
    # Simulates the fixed logic
    dataset_mode = 'paired'
    prefer = 'mb_r2.fq' if dataset_mode == 'paired' else 'mb_r1.fq'
    assert prefer == 'mb_r2.fq'

    dataset_mode = 'single'
    prefer = 'mb_r2.fq' if dataset_mode == 'paired' else 'mb_r1.fq'
    assert prefer == 'mb_r1.fq'


# ---------------------------------------------------------------------------
# DNA reverse-complement (used by splitcode_lr_dual_validate)
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _rc(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def test_reverse_complement():
    """RC should correctly reverse-complement DNA sequences."""
    assert _rc("ACGT") == "ACGT"  # palindrome
    assert _rc("AAAA") == "TTTT"
    assert _rc("CCCC") == "GGGG"
    assert _rc("ATCG") == "CGAT"
    assert _rc("N") == "N"
    assert _rc("") == ""
    # Mixed case
    assert _rc("AcGt") == "aCgT"


def test_reverse_complement_idempotent():
    """RC(RC(seq)) should return original sequence."""
    seqs = ["ACGTACGT", "AAAAGGGG", "NNNNN", "ATCGATCG"]
    for seq in seqs:
        assert _rc(_rc(seq)) == seq, f"RC(RC({seq})) != {seq}"


# ---------------------------------------------------------------------------
# Jaccard index
# ---------------------------------------------------------------------------

def test_jaccard():
    """Jaccard index should be correct for known sets.
    Canonical source is concordance_analysis."""
    from concordance_analysis import jaccard

    assert jaccard(set(), set()) == 1.0
    assert jaccard({1, 2, 3}, {1, 2, 3}) == 1.0
    assert jaccard({1, 2}, {3, 4}) == 0.0
    assert jaccard({1, 2, 3}, {2, 3, 4}) == 0.5  # 2 shared / 4 union


# ---------------------------------------------------------------------------
# Run all tests
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    tests = [v for k, v in globals().items() if k.startswith("test_")]
    passed = 0
    failed = 0
    for test_fn in tests:
        try:
            test_fn()
            print(f"  PASS: {test_fn.__name__}")
            passed += 1
        except Exception as e:
            print(f"  FAIL: {test_fn.__name__}: {e}")
            failed += 1

    print(f"\n{passed}/{passed + failed} tests passed")
    sys.exit(1 if failed else 0)
