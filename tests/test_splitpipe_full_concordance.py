import gzip
import importlib.util
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location(
    "splitpipe_full_concordance",
    SCRIPTS / "splitpipe_full_concordance.py",
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def write_fastq(path: Path, numeric_ids: list[int]) -> None:
    with path.open("wb") as handle:
        for numeric_id in numeric_ids:
            handle.write(
                b"@01_02_03__OH_@ACC."
                + str(numeric_id).encode()
                + b" "
                + str(numeric_id).encode()
                + b"/1\nACGT\n+\nIIII\n"
            )


def test_splitpipe_fastq_to_bitmap(tmp_path):
    path = tmp_path / "barcode_head.fastq"
    write_fastq(path, [1, 3, 8])
    bitmap, provenance = MODULE.splitpipe_fastq_bitmap(path, 8, "ACC")
    assert MODULE.common.popcount(bitmap) == 3
    assert bitmap == bytes([0b10000101])
    assert provenance["accepted_records"] == 3
    assert provenance["minimum_numeric_id"] == 1
    assert provenance["maximum_numeric_id"] == 8
    assert provenance["fastq_validated"] is True


def test_splitpipe_fastq_rejects_duplicate_ids(tmp_path):
    path = tmp_path / "barcode_head.fastq"
    write_fastq(path, [2, 2])
    with pytest.raises(ValueError, match="duplicate split-pipe read ID"):
        MODULE.splitpipe_fastq_bitmap(path, 8, "ACC")


def test_splitpipe_fastq_rejects_bad_quality_length(tmp_path):
    path = tmp_path / "barcode_head.fastq"
    path.write_bytes(b"@x__OH_@ACC.1 1/1\nACGT\n+\nIII\n")
    with pytest.raises(ValueError, match="sequence/quality length mismatch"):
        MODULE.splitpipe_fastq_bitmap(path, 8, "ACC")


def test_compressed_input_matches_campaign_file(tmp_path):
    plain = tmp_path / "input.fastq"
    compressed = tmp_path / "input.fastq.gz"
    plain.write_bytes(b"@ACC.1 1/1\nACGT\n+\nIIII\n")
    with gzip.open(compressed, "wb") as handle:
        handle.write(plain.read_bytes())
    provenance = MODULE.compressed_input_provenance(compressed, plain)
    assert provenance["decompressed_payload"]["byte_identical_to_campaign_input"]
    assert (
        provenance["decompressed_payload"]["sha256"]
        == provenance["campaign_uncompressed"]["sha256"]
    )
