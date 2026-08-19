import gzip
import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parent.parent
SPEC = importlib.util.spec_from_file_location(
    "verify_vendor_concordance_pe",
    ROOT / "scripts" / "verify_vendor_concordance_pe.py",
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def canonical_bitmap(path: Path, ids: list[int], numeric_id_max: int = 16) -> None:
    payload = bytearray((numeric_id_max + 7) // 8)
    for numeric_id in ids:
        byte, bit = divmod(numeric_id - 1, 8)
        payload[byte] |= 1 << bit
    path.write_bytes(
        MODULE.MAGIC
        + b"\0"
        + b"1\0"
        + str(numeric_id_max).encode()
        + b"\0ACC\0"
        + payload
    )


def test_vendor_bitmap_and_metrics(tmp_path):
    vendor_path = tmp_path / "vendor.txt.gz"
    with gzip.open(vendor_path, "wt") as handle:
        handle.write("ACC.1\nACC.3\nACC.8\n")
    vendor, provenance = MODULE.load_vendor_bitmap(vendor_path, 8)
    assert provenance["records"] == 3
    assert provenance["accession_prefix"] == "ACC"

    bitmap_path = tmp_path / "tool.bitmap"
    canonical_bitmap(bitmap_path, [1, 2, 3, 8])
    metadata, full = MODULE.load_canonical_bitmap(bitmap_path)
    assert metadata["numeric_id_max"] == 16
    observed = MODULE.metrics(MODULE.restrict_bitmap(full, 8), vendor, 8)
    assert observed["emitted_records"] == 4
    assert observed["intersection_records"] == 3
    assert observed["precision"] == pytest.approx(0.75)
    assert observed["recall"] == pytest.approx(1.0)
    assert observed["f1"] == pytest.approx(6 / 7)
    assert observed["jaccard"] == pytest.approx(0.75)


def test_vendor_bitmap_rejects_duplicates(tmp_path):
    vendor_path = tmp_path / "vendor.txt"
    vendor_path.write_text("ACC.1\nACC.1\n")
    with pytest.raises(ValueError, match="duplicate vendor read ID"):
        MODULE.load_vendor_bitmap(vendor_path, 8)


def test_prefix_comparison(tmp_path):
    subset = tmp_path / "subset.fastq"
    full = tmp_path / "full.fastq"
    subset.write_bytes(b"prefix")
    full.write_bytes(b"prefix-and-more")
    assert MODULE.files_are_prefix_equal(subset, full)
    full.write_bytes(b"prefiy-and-more")
    assert not MODULE.files_are_prefix_equal(subset, full)
