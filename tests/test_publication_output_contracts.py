from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def config(relative: str) -> str:
    return (ROOT / relative).read_text()


def test_splitseq_pe_splitcode_emits_complete_compact_projection():
    text = config("configs/splitcode/publication_splitseq_pe.config")

    assert "@extract 1:2<umi_bc3_bc2_bc1[8]>" in text
    for group in ("bc3", "bc2", "bc1"):
        assert f"@extract <umi_bc3_bc2_bc1{{{{{group}}}}}>" in text
def test_lr_publication_configs_project_six_base_bc1_and_require_both_linkers():
    splitcode = config("configs/splitcode/publication_lr_splitseq.config")

    assert "MINFINDSG" in splitcode
    assert "@no-chain" in splitcode
    assert "<bc1[6]>" in splitcode
    assert "<bc1[8]>" not in splitcode
    assert "bc1 = b[6]" in config(
        "configs/seqproc/splitseq_singleend_edit_ann.geom"
    )
    assert "bc1 = b[6]" in config("configs/seqproc/publication_lr_splitseq.geom")
    assert "bc1:|6|" in config(
        "configs/matchbox/publication_lr_splitseq.mb"
    )
    assert "bc1:|8|" not in config(
        "configs/matchbox/publication_lr_splitseq.mb"
    )


def test_10x_splitcode_is_length_filter_only():
    text = config("configs/splitcode/publication_10x_v2_filter.config")

    assert "@filter-len 26,0" in text
    assert "@extract" not in text
