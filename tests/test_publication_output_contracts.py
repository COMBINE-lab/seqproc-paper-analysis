from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def config(relative: str) -> str:
    return (ROOT / relative).read_text()


def test_splitseq_pe_splitcode_emits_complete_compact_projection():
    text = config("configs/splitcode/publication_splitseq_pe.config")

    assert "@extract 1:0<umi_bc3_bc2_bc1[10]>" in text
    assert "GTGGCCGATGTTTCGCATCGGCGTACGACT" in text
    assert "AAACATCG\tbc1\t1" in text
    for group in ("bc3", "bc2", "bc1"):
        assert f"@extract <umi_bc3_bc2_bc1{{{{{group}}}}}>" in text


def test_splitseq_pe_configs_model_the_documented_10_8_30_8_30_8_layout():
    seqproc = config("configs/seqproc/publication_splitseq_pe.geom")
    matchbox = config("configs/matchbox/publication_splitseq_pe.mb")

    assert "umi = u[10]" in seqproc
    assert "skip2" not in seqproc
    assert "bc1 = filter_within_dist(b[8]" in seqproc
    assert "GTGGCCGATGTTTCGCATCGGCGTACGACT" in seqproc
    assert "umi:|10|" in matchbox
    assert "bc1.round_23" in matchbox
    assert "rt_6bp.csv" not in matchbox
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


def test_scirnaseq3_configs_preserve_natural_bc1_length_and_guard_captures():
    seqproc = config("configs/seqproc/sciseq3_edit.geom")
    matchbox = config("configs/matchbox/publication_sciseq3.mb")

    assert "brc1 = b[9-10]" in seqproc
    assert "norm(" not in seqproc
    assert "bc1.seq.len() == 9" in matchbox
    assert "bc1.seq.len() == 10" in matchbox


def test_lr_core_configs_require_whitelists_and_complete_components():
    dual = config("configs/seqproc/publication_lr_splitseq_dual_core.geom")
    forward = config("configs/seqproc/publication_lr_splitseq_forward_core.geom")
    matchbox = config("configs/matchbox/publication_lr_splitseq_dual_core.mb")
    splitcode = config("configs/splitcode/publication_lr_splitseq_core.config")

    assert dual.count("filter_within_dist") == 3
    assert "#[match_ori(either)]" in dual
    assert "#[edit(3)] linker_a" in dual
    assert "#[edit(3)] linker_b" in dual
    assert "#[match_ori" not in forward
    assert "bc3_hit.seq.len() == 8" in matchbox
    assert "bcs1_seqs.contains" in matchbox
    assert "bc3:|8| linker1~0.0" in matchbox
    assert "bc3_hit:bc3.round_23~0.0 linker1" not in matchbox
    assert "@extract\t<prefix[18]>{linker1}" in splitcode
