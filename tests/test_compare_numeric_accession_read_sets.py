import sys
from pathlib import Path

import pytest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from compare_numeric_accession_read_sets import intersection_size, read_bitset


def write_fastq(path: Path, indices: list[int]) -> None:
    path.write_text(
        "".join(f"@SRR1.{index} comment\nA\n+\nI\n" for index in indices)
    )


def test_numeric_accession_bitsets_support_exact_pairwise_counts(tmp_path):
    left = tmp_path / "left.fastq"
    right = tmp_path / "right.fastq"
    write_fastq(left, [1, 3, 5])
    write_fastq(right, [2, 3, 5])

    left_bits, left_count, left_digest = read_bitset(left, 5)
    right_bits, right_count, right_digest = read_bitset(right, 5)

    assert left_count == right_count == 3
    assert intersection_size(left_bits, right_bits) == 2
    assert len(left_digest) == len(right_digest) == 64


def test_duplicate_numeric_accession_id_is_rejected(tmp_path):
    path = tmp_path / "duplicate.fastq"
    write_fastq(path, [1, 1])

    with pytest.raises(ValueError, match="duplicate read IDs"):
        read_bitset(path, 2)
