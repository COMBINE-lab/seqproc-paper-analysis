import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from generate_emitted_set_upset import exclusive_counts, percent_label


def test_exclusive_counts_partition_three_sets():
    # Four-element universe with memberships: all, A+B, A+C, and C only.
    a = 0b0111
    b = 0b0011
    c = 0b1101
    counts = exclusive_counts((a, b, c), 4)
    assert counts == {1: 0, 2: 0, 3: 1, 4: 1, 5: 1, 6: 0, 7: 1}
    assert sum(counts.values()) == 4


def test_percent_label_preserves_small_nonzero_intersections():
    assert percent_label(0, 0.0) == ""
    assert percent_label(10, 0.0018) == "<0.01%"
    assert percent_label(100, 12.345) == "12.35%"
