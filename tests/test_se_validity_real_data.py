#!/usr/bin/env python3
"""Real-data validation of SplitSeqSingleEndValidityAnalyzer fix.

Uses actual LR-SPLiT-seq reads (SRR13948564) to confirm:
1. The OLD buggy analyzer produces 0 valid reads (reproduces the bug)
2. The FIXED analyzer produces a non-trivial number of valid reads
3. The valid-read count is consistent with Phase 3 empirical data
   (Phase 3 showed ~43% barcoded, ~50% forward, ~74% exact linker match
    => expected ~16% valid on exact-match forward-only,
    with RC support => ~32%)

/bugfix methodology: test written to FAIL on old code, PASS on new code.
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

# Patch out heavy deps
import unittest.mock as mock
for mod in ['numpy', 'matplotlib', 'matplotlib.pyplot',
            'matplotlib.patches', 'matplotlib.gridspec']:
    sys.modules[mod] = mock.MagicMock()

from run_paper_benchmarks import SplitSeqSingleEndValidityAnalyzer

PROJECT_ROOT = os.path.join(os.path.dirname(__file__), '..')
REAL_FQ = os.path.join(PROJECT_ROOT, 'data', 'SRR13948564_1M.fastq')
BC1_MAP = os.path.join(PROJECT_ROOT, 'configs', 'seqproc', 'splitseq_bc1_seq2seq.tsv')
BC2_MAP = os.path.join(PROJECT_ROOT, 'configs', 'seqproc', 'splitseq_bc2_seq2seq.tsv')
BC3_MAP = os.path.join(PROJECT_ROOT, 'configs', 'seqproc', 'splitseq_bc3_seq2seq.tsv')

# 10K read subset for fast testing
SAMPLE_FQ = '/tmp/lr_splitseq_10k.fq'
SAMPLE_READS = 10_000

LINKER1 = "GTGGCCGATGTTTCGCATCGGCGTACGACT"  # 30bp


def ensure_sample():
    """Create 10K read sample if it doesn't exist."""
    if os.path.exists(SAMPLE_FQ):
        return
    if not os.path.exists(REAL_FQ):
        raise FileNotFoundError(
            f"Real data not found: {REAL_FQ}\n"
            f"Download SRR13948564 1M subset to run this test."
        )
    with open(REAL_FQ) as fin, open(SAMPLE_FQ, 'w') as fout:
        for i in range(SAMPLE_READS * 4):
            line = fin.readline()
            if not line:
                break
            fout.write(line)


# ============================================================================
# OLD BUGGY ANALYZER (embedded inline to reproduce the bug)
# ============================================================================

class OldBuggySEAnalyzer:
    """Exact copy of the old broken analyzer for regression testing.

    Bug 1: _find_linker searches only positions start..start+40 (40bp window).
            PacBio reads are ~1kb; barcode structure can be anywhere.
    Bug 2: No reverse complement support. ~50% of PacBio reads are RC.
    Bug 3: Extracts 8bp for BC1 but whitelist has 6bp entries.
            _hamming(8bp, 6bp) always returns 99 => every read fails BC1 check.
    """

    LINKER1 = "GTGGCCGATGTTTCGCATCGGCGTACGACT"

    def __init__(self, bc1_map, bc2_map, bc3_map):
        self.bc1_wl = self._load_whitelist(bc1_map)
        self.bc2_wl = self._load_whitelist(bc2_map)
        self.bc3_wl = self._load_whitelist(bc3_map)

    def _load_whitelist(self, path):
        wl = set()
        with open(path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    wl.add(parts[1])
        return wl

    @staticmethod
    def _hamming(s1, s2):
        if len(s1) != len(s2): return 99
        return sum(a != b for a, b in zip(s1, s2))

    def _check_wl(self, bc, wl):
        if bc in wl: return True
        for cand in wl:
            if self._hamming(bc, cand) <= 1: return True
        return False

    def _find_linker(self, read, linker, start=0):
        """BUG: Only searches positions start..start+40."""
        best_pos, best_dist = -1, 100
        search_end = min(len(read) - len(linker) + 1, start + 40)
        for i in range(start, search_end):
            dist = self._hamming(read[i:i+len(linker)], linker)
            if dist < best_dist:
                best_dist = dist
                best_pos = i
                if dist <= 1: break
        return best_pos, best_dist

    def analyze(self, fq_path):
        """Run old buggy analysis. Returns set of valid read IDs."""
        valid_ids = set()
        with open(fq_path) as f:
            while True:
                header = f.readline()
                if not header: break
                seq = f.readline().strip()
                f.readline(); f.readline()
                read_id = header.strip().split()[0].replace('@', '')

                # BUG: starts search at 10, window is only 40bp
                l1_pos, l1_dist = self._find_linker(seq, self.LINKER1, 10)
                if l1_dist > 3: continue
                if l1_pos < 18: continue

                bc3 = seq[l1_pos-8:l1_pos]
                bc2 = seq[l1_pos+30:l1_pos+38]
                # BUG: extracts 8bp for BC1 but whitelist has 6bp entries
                bc1 = seq[l1_pos+30+8+22:l1_pos+30+8+22+8]

                if len(bc3) != 8 or len(bc2) != 8 or len(bc1) != 8:
                    continue

                if (self._check_wl(bc3, self.bc3_wl) and
                    self._check_wl(bc2, self.bc2_wl) and
                    self._check_wl(bc1, self.bc1_wl)):
                    valid_ids.add(read_id)
        return valid_ids


# ============================================================================
# Diagnostic helpers
# ============================================================================

def count_linker_positions(fq_path, max_reads=10000):
    """Count where linker appears in reads, both orientations."""
    COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
    fw_found = 0
    rc_found = 0
    neither = 0
    deep_fw = 0  # linker at position > 50 (outside old search window)

    with open(fq_path) as f:
        count = 0
        while count < max_reads:
            header = f.readline()
            if not header: break
            seq = f.readline().strip()
            f.readline(); f.readline()
            count += 1

            fw_pos = seq.find(LINKER1)
            rc_seq = seq.translate(COMP)[::-1]
            rc_pos = rc_seq.find(LINKER1)

            if fw_pos >= 0:
                fw_found += 1
                if fw_pos > 50:
                    deep_fw += 1
            elif rc_pos >= 0:
                rc_found += 1
            else:
                neither += 1

    return {
        'total': count,
        'fw_found': fw_found,
        'rc_found': rc_found,
        'neither': neither,
        'deep_fw': deep_fw,
    }


# ============================================================================
# Tests
# ============================================================================

def test_old_analyzer_returns_zero_on_real_data():
    """BUG REPRODUCTION: The old analyzer must return 0 valid reads.

    Root cause: BC1 whitelist has 6bp entries, but _try_extract extracts
    8bp for BC1. _hamming(8bp, 6bp) always returns 99, so the BC1
    whitelist check always fails, yielding 0 valid reads for every read.
    This is independent of the search window and RC bugs.
    """
    ensure_sample()
    old = OldBuggySEAnalyzer(BC1_MAP, BC2_MAP, BC3_MAP)
    result = old.analyze(SAMPLE_FQ)

    # Verify the BC1 length mismatch is real
    bc1_lengths = {len(e) for e in old.bc1_wl}
    assert bc1_lengths == {6}, (
        f"Expected BC1 whitelist to have 6bp entries, got lengths: {bc1_lengths}"
    )

    assert len(result) == 0, (
        f"OLD analyzer should find 0 valid reads due to BC1 length mismatch, "
        f"but found {len(result)}"
    )
    print(f"[PASS] test_old_analyzer_returns_zero_on_real_data: "
          f"0 valid reads (confirmed BC1 6bp vs 8bp extraction bug)")


def test_fixed_analyzer_finds_valid_reads_on_real_data():
    """FIXED: The new analyzer must find a non-trivial number of valid reads.

    Phase 3 empirical data (1M reads) showed:
    - ~43% of reads have barcodes (linker detectable in either orientation)
    - ~50% of barcoded reads are forward, ~50% RC
    - HiFi reads have ~1% error rate => ~74% have error-free 30bp linker
    Expected on 10K sample: ~43% * 74% exact-match = ~32% with linker found
    With RC support, slightly higher. Conservatively expect >5% valid.
    """
    ensure_sample()

    with tempfile.TemporaryDirectory() as td:
        # Copy sample to temp to avoid cache interference
        import shutil
        tmp_fq = os.path.join(td, 'test.fq')
        shutil.copy(SAMPLE_FQ, tmp_fq)

        analyzer = SplitSeqSingleEndValidityAnalyzer(BC1_MAP, BC2_MAP, BC3_MAP)
        result = analyzer.analyze_fastqs(tmp_fq)

    valid_count = len(result)
    valid_pct = valid_count / SAMPLE_READS * 100

    # Sanity: must find SOME valid reads (old analyzer found 0)
    assert valid_count > 0, (
        f"FIXED analyzer found 0 valid reads -- fix did not work!"
    )

    # Phase 3 showed ~20% forward-only exact recovery on 1M reads.
    # With RC support, expect higher. Be conservative: >5% of 10K.
    min_expected = SAMPLE_READS * 0.05  # 500 reads
    assert valid_count > min_expected, (
        f"FIXED analyzer found only {valid_count} valid reads ({valid_pct:.1f}%). "
        f"Expected >{min_expected:.0f} (>5%) based on Phase 3 data."
    )

    # Upper bound sanity: should not exceed ~50% (the barcoded population)
    max_expected = SAMPLE_READS * 0.50
    assert valid_count < max_expected, (
        f"FIXED analyzer found {valid_count} valid reads ({valid_pct:.1f}%). "
        f"This exceeds the ~43% barcoded population -- something is wrong."
    )

    print(f"[PASS] test_fixed_analyzer_finds_valid_reads_on_real_data: "
          f"{valid_count:,} valid reads ({valid_pct:.1f}%) in {SAMPLE_READS:,} sample")


def test_linker_positions_confirm_deep_and_rc():
    """Diagnostic: confirm that real reads have linkers at deep positions
    and in RC orientation, validating that the old 40bp window and
    forward-only search were insufficient.
    """
    ensure_sample()
    stats = count_linker_positions(SAMPLE_FQ)

    print(f"  Linker position stats on {stats['total']:,} real reads:")
    print(f"    Forward linker found: {stats['fw_found']:,} "
          f"({stats['fw_found']/stats['total']*100:.1f}%)")
    print(f"    RC linker found:      {stats['rc_found']:,} "
          f"({stats['rc_found']/stats['total']*100:.1f}%)")
    print(f"    No linker:            {stats['neither']:,} "
          f"({stats['neither']/stats['total']*100:.1f}%)")
    print(f"    FW at position > 50:  {stats['deep_fw']:,} "
          f"({stats['deep_fw']/max(stats['fw_found'],1)*100:.1f}% of FW)")

    # Must have reads in RC orientation (confirms bug #3)
    assert stats['rc_found'] > 0, (
        f"Expected RC linker reads in real PacBio data, found 0"
    )

    # Phase 3 showed ~20% forward, ~20% RC, ~57% no linker
    # On 10K sample, expect at least 1000 in each orientation
    assert stats['fw_found'] > 500, (
        f"Expected >500 forward-linker reads, got {stats['fw_found']}"
    )
    assert stats['rc_found'] > 500, (
        f"Expected >500 RC-linker reads, got {stats['rc_found']}"
    )

    # Confirm deep linker positions exist (validates bug #2)
    assert stats['deep_fw'] > 0, (
        f"Expected some reads with linker at position > 50, found 0. "
        f"This would mean the old 40bp window was actually sufficient."
    )

    deep_pct = stats['deep_fw'] / max(stats['fw_found'], 1) * 100
    print(f"\n  [PASS] test_linker_positions_confirm_deep_and_rc")
    print(f"  RC reads confirm bug #3 (no RC support).")
    print(f"  {stats['deep_fw']} reads ({deep_pct:.1f}% of FW) have linker "
          f"at pos>50, confirming bug #2 (40bp window).")


def test_bc1_whitelist_length_is_6bp():
    """Confirm the real BC1 whitelist has 6bp entries.
    This is the prerequisite for the BC1 length mismatch bug.
    """
    wl = set()
    with open(BC1_MAP) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                wl.add(parts[1])

    lengths = {len(e) for e in wl}
    assert lengths == {6}, (
        f"Expected BC1 whitelist to have only 6bp entries, "
        f"got lengths: {lengths}"
    )
    print(f"[PASS] test_bc1_whitelist_length_is_6bp: "
          f"{len(wl)} entries, all 6bp (confirms bug #1 prerequisite)")


if __name__ == '__main__':
    print("=" * 70)
    print("SE Validity Analyzer: Real Data Validation (/bugfix)")
    print(f"Data: {REAL_FQ}")
    print(f"Sample: {SAMPLE_READS:,} reads")
    print("=" * 70)

    if not os.path.exists(REAL_FQ):
        print(f"\n[SKIP] Real data not available: {REAL_FQ}")
        print("This test requires SRR13948564_1M.fastq in data/")
        sys.exit(0)

    print("\n--- Step 1: Confirm bug prerequisites ---")
    test_bc1_whitelist_length_is_6bp()

    print("\n--- Step 2: Confirm linker positions require full-read + RC ---")
    test_linker_positions_confirm_deep_and_rc()

    print("\n--- Step 3: OLD analyzer returns 0 (bug reproduction) ---")
    test_old_analyzer_returns_zero_on_real_data()

    print("\n--- Step 4: FIXED analyzer finds valid reads ---")
    test_fixed_analyzer_finds_valid_reads_on_real_data()

    print("\n" + "=" * 70)
    print("ALL REAL-DATA TESTS PASSED")
    print("=" * 70)
