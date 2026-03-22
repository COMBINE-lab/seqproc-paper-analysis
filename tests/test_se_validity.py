#!/usr/bin/env python3
"""Tests for SplitSeqSingleEndValidityAnalyzer.

Verifies that the analyzer:
1. Searches the entire read for the linker (not just positions 10-50)
2. Handles both forward and reverse-complement orientations
3. Correctly extracts barcodes relative to the linker position
4. Validates barcodes against whitelists with Hamming distance <= 1
"""

import os
import sys
import tempfile

# Add scripts dir to path so we can import the analyzer
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

# Patch out matplotlib/numpy so we can import run_paper_benchmarks without them
# (they may not be installed locally)
import unittest.mock as mock
sys.modules['numpy'] = mock.MagicMock()
sys.modules['matplotlib'] = mock.MagicMock()
sys.modules['matplotlib.pyplot'] = mock.MagicMock()
sys.modules['matplotlib.patches'] = mock.MagicMock()
sys.modules['matplotlib.gridspec'] = mock.MagicMock()

from run_paper_benchmarks import SplitSeqSingleEndValidityAnalyzer


# --- Constants from the analyzer ---
LINKER1 = "GTGGCCGATGTTTCGCATCGGCGTACGACT"  # 30bp

# Known valid barcodes from the whitelist files
# BC3 and BC2 are 8bp, BC1 is 8bp in seq2seq.tsv
BC3_VALID = "AACGTGAT"  # first entry in bc3_seq2seq.tsv
BC2_VALID = "AACGTGAT"  # first entry in bc2_seq2seq.tsv
# BC1 entries are 6bp in the tsv -- analyzer must handle this
BC1_VALID_6BP = "AAACAT"  # first entry in bc1_seq2seq.tsv

COMP = str.maketrans('ACGTacgt', 'TGCAtgca')


def rc(seq):
    return seq.translate(COMP)[::-1]


def random_seq(n):
    """Generate a deterministic pseudo-random sequence of length n."""
    import hashlib
    h = hashlib.sha256(str(n).encode()).hexdigest()
    bases = 'ACGT'
    result = []
    for i in range(n):
        idx = int(h[(i * 2) % len(h)], 16) % 4
        result.append(bases[idx])
    return ''.join(result)


def make_barcode_region_8bp_bc1():
    """Build the barcode region with 8bp BC1 (for testing purposes)."""
    umi = "ACGTACGTAC"  # 10bp
    linker2 = "ATCCACGTGCTTGAGACTGTGG"  # 22bp (from geom)
    return umi + BC3_VALID + LINKER1 + BC2_VALID + linker2 + BC3_VALID  # reuse BC3 as 8bp BC1


def make_barcode_region_6bp_bc1():
    """Build the barcode region with 6bp BC1 (actual whitelist)."""
    umi = "ACGTACGTAC"  # 10bp
    linker2 = "ATCCACGTGCTTGAGACTGTGG"  # 22bp
    return umi + BC3_VALID + LINKER1 + BC2_VALID + linker2 + BC1_VALID_6BP


def write_fastq(filepath, records):
    """Write FASTQ records as list of (read_id, seq) tuples."""
    with open(filepath, 'w') as f:
        for read_id, seq in records:
            f.write(f"@{read_id}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write("I" * len(seq) + "\n")


# --- Whitelist file helpers ---
def write_whitelist_8bp(filepath, entries):
    """Write a TSV whitelist with col1=col2 (8bp entries)."""
    with open(filepath, 'w') as f:
        for e in entries:
            f.write(f"{e}\t{e}\n")


def write_whitelist_6bp(filepath, entries):
    """Write a TSV whitelist with 6bp entries."""
    with open(filepath, 'w') as f:
        for e in entries:
            f.write(f"{e}\t{e}\n")


# ============================================================================
# Tests
# ============================================================================

def test_forward_read_linker_at_start():
    """Barcode region near the start of the read (position 18 for linker)."""
    with tempfile.TemporaryDirectory() as td:
        bc1_f = os.path.join(td, "bc1.tsv")
        bc2_f = os.path.join(td, "bc2.tsv")
        bc3_f = os.path.join(td, "bc3.tsv")
        fq = os.path.join(td, "test.fq")

        write_whitelist_8bp(bc3_f, [BC3_VALID, "AAACATCG", "ATGCCTAA"])
        write_whitelist_8bp(bc2_f, [BC2_VALID, "AAACATCG", "ATGCCTAA"])
        # BC1 whitelist has 8bp entries -- use BC3_VALID as a stand-in
        write_whitelist_8bp(bc1_f, [BC3_VALID, "AAACATCG", "ATGCCTAA"])

        barcode_region = make_barcode_region_8bp_bc1()
        # Barcode at the very start (small prefix)
        seq = barcode_region + random_seq(500)

        write_fastq(fq, [("read_fw_start", seq)])

        analyzer = SplitSeqSingleEndValidityAnalyzer(bc1_f, bc2_f, bc3_f)
        result = analyzer.analyze_fastqs(fq)

        assert "read_fw_start" in result, (
            f"Forward read with linker near start not found as valid. "
            f"Linker at pos {seq.find(LINKER1)}, read len {len(seq)}"
        )
        print("[PASS] test_forward_read_linker_at_start")


def test_forward_read_linker_deep():
    """Barcode region deep inside a long read (simulates skip_start = r:)."""
    with tempfile.TemporaryDirectory() as td:
        bc1_f = os.path.join(td, "bc1.tsv")
        bc2_f = os.path.join(td, "bc2.tsv")
        bc3_f = os.path.join(td, "bc3.tsv")
        fq = os.path.join(td, "test.fq")

        write_whitelist_8bp(bc3_f, [BC3_VALID])
        write_whitelist_8bp(bc2_f, [BC2_VALID])
        write_whitelist_8bp(bc1_f, [BC3_VALID])

        barcode_region = make_barcode_region_8bp_bc1()
        # 500bp prefix (simulates cDNA before barcode region)
        seq = random_seq(500) + barcode_region + random_seq(200)

        write_fastq(fq, [("read_fw_deep", seq)])

        analyzer = SplitSeqSingleEndValidityAnalyzer(bc1_f, bc2_f, bc3_f)
        result = analyzer.analyze_fastqs(fq)

        l1_pos = seq.find(LINKER1)
        assert "read_fw_deep" in result, (
            f"Forward read with linker deep in read not found. "
            f"Linker at pos {l1_pos} (read len {len(seq)}). "
            f"Old analyzer would have searched only positions 10-50."
        )
        print("[PASS] test_forward_read_linker_deep")


def test_rc_read():
    """Barcode region in reverse complement orientation."""
    with tempfile.TemporaryDirectory() as td:
        bc1_f = os.path.join(td, "bc1.tsv")
        bc2_f = os.path.join(td, "bc2.tsv")
        bc3_f = os.path.join(td, "bc3.tsv")
        fq = os.path.join(td, "test.fq")

        write_whitelist_8bp(bc3_f, [BC3_VALID])
        write_whitelist_8bp(bc2_f, [BC2_VALID])
        write_whitelist_8bp(bc1_f, [BC3_VALID])

        barcode_region = make_barcode_region_8bp_bc1()
        # Forward read with barcode, then RC the whole thing
        fw_seq = random_seq(300) + barcode_region + random_seq(400)
        rc_seq = rc(fw_seq)

        # Verify linker is NOT in forward but IS findable via RC
        assert rc_seq.find(LINKER1) < 0, "Linker should not be in forward of RC read"
        assert rc(rc_seq).find(LINKER1) >= 0, "Linker should be in RC of RC read"

        write_fastq(fq, [("read_rc", rc_seq)])

        analyzer = SplitSeqSingleEndValidityAnalyzer(bc1_f, bc2_f, bc3_f)
        result = analyzer.analyze_fastqs(fq)

        assert "read_rc" in result, (
            f"RC read not found as valid. "
            f"Old analyzer had no RC support."
        )
        print("[PASS] test_rc_read")


def test_no_linker_read():
    """Read with no linker should NOT be marked valid."""
    with tempfile.TemporaryDirectory() as td:
        bc1_f = os.path.join(td, "bc1.tsv")
        bc2_f = os.path.join(td, "bc2.tsv")
        bc3_f = os.path.join(td, "bc3.tsv")
        fq = os.path.join(td, "test.fq")

        write_whitelist_8bp(bc3_f, [BC3_VALID])
        write_whitelist_8bp(bc2_f, [BC2_VALID])
        write_whitelist_8bp(bc1_f, [BC3_VALID])

        # Pure random sequence -- no linker
        seq = random_seq(1000)
        assert seq.find(LINKER1) < 0, "Random seq should not contain linker"

        write_fastq(fq, [("read_no_linker", seq)])

        analyzer = SplitSeqSingleEndValidityAnalyzer(bc1_f, bc2_f, bc3_f)
        result = analyzer.analyze_fastqs(fq)

        assert "read_no_linker" not in result, "Read with no linker should not be valid"
        print("[PASS] test_no_linker_read")


def test_bad_barcodes():
    """Read with linker but invalid barcodes should NOT be valid."""
    with tempfile.TemporaryDirectory() as td:
        bc1_f = os.path.join(td, "bc1.tsv")
        bc2_f = os.path.join(td, "bc2.tsv")
        bc3_f = os.path.join(td, "bc3.tsv")
        fq = os.path.join(td, "test.fq")

        write_whitelist_8bp(bc3_f, [BC3_VALID])
        write_whitelist_8bp(bc2_f, [BC2_VALID])
        write_whitelist_8bp(bc1_f, [BC3_VALID])

        # Build a read with linker but completely wrong barcodes (>d1)
        umi = "ACGTACGTAC"
        bad_bc = "TTTTTTTT"  # Not in any whitelist and >d1 from all entries
        linker2 = "ATCCACGTGCTTGAGACTGTGG"
        barcode_region = umi + bad_bc + LINKER1 + bad_bc + linker2 + bad_bc
        seq = random_seq(200) + barcode_region + random_seq(200)

        write_fastq(fq, [("read_bad_bc", seq)])

        analyzer = SplitSeqSingleEndValidityAnalyzer(bc1_f, bc2_f, bc3_f)
        result = analyzer.analyze_fastqs(fq)

        assert "read_bad_bc" not in result, "Read with bad barcodes should not be valid"
        print("[PASS] test_bad_barcodes")


def test_mixed_reads():
    """Mix of valid fw, valid rc, and invalid reads."""
    with tempfile.TemporaryDirectory() as td:
        bc1_f = os.path.join(td, "bc1.tsv")
        bc2_f = os.path.join(td, "bc2.tsv")
        bc3_f = os.path.join(td, "bc3.tsv")
        fq = os.path.join(td, "test.fq")

        write_whitelist_8bp(bc3_f, [BC3_VALID])
        write_whitelist_8bp(bc2_f, [BC2_VALID])
        write_whitelist_8bp(bc1_f, [BC3_VALID])

        barcode_region = make_barcode_region_8bp_bc1()

        records = [
            # Valid forward, linker at various positions
            ("fw_pos50", random_seq(50) + barcode_region + random_seq(300)),
            ("fw_pos500", random_seq(500) + barcode_region + random_seq(100)),
            ("fw_pos0", barcode_region + random_seq(400)),
            # Valid RC
            ("rc_pos200", rc(random_seq(200) + barcode_region + random_seq(500))),
            ("rc_pos800", rc(random_seq(800) + barcode_region + random_seq(50))),
            # Invalid -- no linker
            ("invalid1", random_seq(800)),
            ("invalid2", random_seq(1200)),
        ]

        write_fastq(fq, records)

        analyzer = SplitSeqSingleEndValidityAnalyzer(bc1_f, bc2_f, bc3_f)
        result = analyzer.analyze_fastqs(fq)

        expected_valid = {"fw_pos50", "fw_pos500", "fw_pos0", "rc_pos200", "rc_pos800"}
        expected_invalid = {"invalid1", "invalid2"}

        for rid in expected_valid:
            assert rid in result, f"Expected {rid} to be valid but it wasn't"
        for rid in expected_invalid:
            assert rid not in result, f"Expected {rid} to be invalid but it was marked valid"

        assert len(result) == len(expected_valid), (
            f"Expected {len(expected_valid)} valid reads, got {len(result)}: {result}"
        )
        print(f"[PASS] test_mixed_reads ({len(result)} valid, {len(expected_invalid)} invalid)")


def test_bc1_length_mismatch():
    """BC1 whitelist has 6bp entries but read has 8bp at BC1 position.

    The analyzer must handle this correctly -- either by extracting only
    6bp for BC1 or by handling the length mismatch in _check_wl.
    This test documents the current behavior.
    """
    with tempfile.TemporaryDirectory() as td:
        bc1_f = os.path.join(td, "bc1.tsv")
        bc2_f = os.path.join(td, "bc2.tsv")
        bc3_f = os.path.join(td, "bc3.tsv")
        fq = os.path.join(td, "test.fq")

        write_whitelist_8bp(bc3_f, [BC3_VALID])
        write_whitelist_8bp(bc2_f, [BC2_VALID])
        # BC1 with 6bp entries (actual data)
        write_whitelist_6bp(bc1_f, [BC1_VALID_6BP])

        barcode_region = make_barcode_region_6bp_bc1()
        seq = random_seq(100) + barcode_region + random_seq(300)

        write_fastq(fq, [("read_6bp_bc1", seq)])

        analyzer = SplitSeqSingleEndValidityAnalyzer(bc1_f, bc2_f, bc3_f)

        # Check what the whitelist actually loaded
        print(f"    BC1 whitelist entries: {analyzer.bc1_wl}")
        print(f"    BC1 entry lengths: {[len(e) for e in analyzer.bc1_wl]}")

        result = analyzer.analyze_fastqs(fq)

        # This will currently FAIL because we extract 8bp but whitelist has 6bp
        # _hamming(8bp, 6bp) returns 99
        if "read_6bp_bc1" in result:
            print("[PASS] test_bc1_length_mismatch -- analyzer handles 6bp BC1 correctly")
        else:
            print("[FAIL] test_bc1_length_mismatch -- BC1 length mismatch!")
            print("       Analyzer extracts 8bp for BC1 but whitelist has 6bp entries.")
            print("       _hamming(8bp, 6bp) returns 99, so no read can ever be valid.")
            print("       THIS IS THE ROOT CAUSE of 0 valid reads.")
            return False
    return True


if __name__ == "__main__":
    print("=" * 60)
    print("SplitSeqSingleEndValidityAnalyzer Tests")
    print("=" * 60)

    # First run the BC1 length mismatch test -- this is the critical one
    print("\n--- BC1 Length Mismatch Test (root cause diagnostic) ---")
    bc1_ok = test_bc1_length_mismatch()

    if not bc1_ok:
        print("\n[CRITICAL] BC1 length mismatch must be fixed before other tests can pass.")
        print("Fix: extract 6bp for BC1 instead of 8bp, or adapt _check_wl.")
        sys.exit(1)

    print("\n--- Core Functionality Tests ---")
    test_forward_read_linker_at_start()
    test_forward_read_linker_deep()
    test_rc_read()
    test_no_linker_read()
    test_bad_barcodes()
    test_mixed_reads()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED")
    print("=" * 60)
