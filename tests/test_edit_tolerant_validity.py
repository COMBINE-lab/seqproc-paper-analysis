import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from edit_tolerant_validity import LINKERS, revcomp, validate_orientation


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
