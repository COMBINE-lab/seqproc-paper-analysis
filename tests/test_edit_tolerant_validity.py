import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from edit_tolerant_validity import (
    LINKERS,
    Policy,
    _hamming_owner_map,
    _paired_records,
    _splitseq_result,
    revcomp,
    validate_orientation,
    validate_scirnaseq3,
)


L1, L2 = LINKERS["lr"]
UMI = "ACGTACGTAC"
BC3 = "AACGTGAT"
BC2 = "AAACATCG"
BC1 = "AAACAT"
BC23 = {BC3, BC2}
BC1S = {BC1}


def validate(seq):
    return validate_orientation(seq, L1, L2, BC23, BC1S, 8, 6, 6, 6)


def cassette(linker1=L1, separator=""):
    return UMI + BC3 + linker1 + BC2 + separator + L2 + BC1 + "TTTT"


def test_complete_exact_cassette_is_accepted():
    result = validate("GGGG" + cassette())
    assert result["accepted"]
    assert result["linker1_edit"] == 0
    assert result["linker2_edit"] == 0


def test_truncated_umi_prefix_is_rejected():
    result = validate(UMI[1:] + BC3 + L1 + BC2 + L2 + BC1)
    assert not result["accepted"]
    assert result["reason"] == "incomplete_umi_bc3_prefix"


def test_excessive_inter_linker_gap_is_rejected():
    # Up to six inserted bases may be explained as linker edits; seven may not.
    result = validate(cassette(separator="G" * 7))
    assert not result["accepted"]
    assert result["reason"] == "no_linker2_at_expected_offset"


def test_worse_valid_linker_is_found_after_better_decoy():
    invalid_decoy = "T" * 18 + L1 + "T" * 80
    one_edit_l1 = "A" + L1[1:]
    result = validate(invalid_decoy + cassette(linker1=one_edit_l1))
    assert result["accepted"]
    assert result["linker1_edit"] == 1


def test_reverse_orientation_can_be_validated_after_reverse_complement():
    reverse_read = revcomp("GGGG" + cassette())
    assert not validate(reverse_read)["accepted"]
    assert validate(revcomp(reverse_read))["accepted"]


def test_documented_pe_linker_and_full_eight_base_bc1_are_used():
    pe_l1, _ = LINKERS["pe"]
    assert pe_l1 == "GTGGCCGATGTTTCGCATCGGCGTACGACT"

    full_bc1 = "AAACATCG"
    result = validate_orientation(
        UMI + BC3 + pe_l1 + BC2 + LINKERS["pe"][1] + full_bc1,
        *LINKERS["pe"],
        {BC3, BC2},
        {full_bc1},
        8,
        8,
        3,
        3,
    )
    assert result["accepted"]
    assert result["fast_path"]


def test_barcode_membership_ambiguity_is_accepted_but_reported():
    # AAAT is one substitution from both AAAA and AATT.
    bc1_owners = _hamming_owner_map((b"AAAA", b"AATT"))
    result = validate_orientation(
        b"A" * 10 + b"AACGTGAT" + L1.encode() + b"AAACATCG"
        + L2.encode() + b"AAAT",
        L1.encode(),
        L2.encode(),
        _hamming_owner_map((b"AACGTGAT", b"AAACATCG")),
        bc1_owners,
        8,
        4,
        3,
        3,
    )
    assert result["accepted"]
    assert result["reason"] == "accepted_ambiguous_barcode"
    assert result["ambiguous_barcodes"] == 1


def test_lr_reports_reads_accepted_in_both_orientations():
    # Construct an artificial palindrome of two opposing valid cassettes.
    forward = cassette()
    sequence = forward + "G" * 12 + revcomp(forward)
    policy = Policy(
        chem="lr",
        linker1=L1.encode(),
        linker2=L2.encode(),
        bc23=_hamming_owner_map((BC3.encode(), BC2.encode())),
        bc1=_hamming_owner_map((BC1.encode(),)),
        bc23_length=8,
        bc1_length=6,
        max_linker1_edit=3,
        max_linker2_edit=3,
        orientation="both",
    )
    result = _splitseq_result(sequence.encode(), policy)
    assert result["accepted"]
    assert result["orientation"] == "both"


def test_scirnaseq3_checks_allowed_offsets_instead_of_global_best_decoy():
    # A perfect decoy at offset zero must not hide the one-edit legal anchor.
    sequence = b"CAGAGCTTT" + b"CGGAGC" + b"A" * 18
    result = validate_scirnaseq3(sequence)
    assert result["accepted"]
    assert result["anchor_start"] == 9
    assert result["anchor_edit"] == 1


def test_scirnaseq3_rejects_equal_best_offset_ambiguity():
    # This sequence has one-edit prefix matches beginning at both offsets.
    result = validate_scirnaseq3(b"A" * 9 + b"TAGAGC" + b"A" * 18)
    assert not result["accepted"]
    assert result["reason"] == "ambiguous_anchor_offset"


def _write_fastq(path, records):
    with open(path, "wb") as handle:
        for read_id, sequence in records:
            handle.write(
                b"@" + read_id + b"\n" + sequence + b"\n+\n" + b"I" * len(sequence) + b"\n"
            )


def test_paired_reader_checks_ids_and_accepts_slash_mate_suffixes(tmp_path):
    first = tmp_path / "r1.fastq"
    second = tmp_path / "r2.fastq"
    _write_fastq(first, [(b"read/1", b"A" * 26)])
    _write_fastq(second, [(b"read/2", b"C" * 30)])
    assert list(_paired_records(str(first), str(second))) == [(b"read", b"A" * 26)]


def test_paired_reader_rejects_mismatched_ids(tmp_path):
    first = tmp_path / "r1.fastq"
    second = tmp_path / "r2.fastq"
    _write_fastq(first, [(b"left", b"A" * 26)])
    _write_fastq(second, [(b"right", b"C" * 30)])
    try:
        list(_paired_records(str(first), str(second)))
    except ValueError as error:
        assert "ID mismatch" in str(error)
    else:
        raise AssertionError("mismatched read IDs must fail validation")
