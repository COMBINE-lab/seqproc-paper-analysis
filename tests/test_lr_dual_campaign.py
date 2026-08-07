import json
import stat
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from make_publication_core_manifest import DATASETS, DEFAULT_REPLICATES, DEFAULT_THREADS
from run_splitcode_dual_pass import main as dual_pass_main
from stage_lr_reverse_complement import stage


def test_publication_campaign_defaults_and_lr_roles():
    datasets = {item["name"]: item for item in DATASETS}

    assert DEFAULT_THREADS == (1, 4, 16, 32)
    assert DEFAULT_REPLICATES == 3
    assert datasets["lr_splitseq_dual"]["analysis_role"] == "primary"
    assert datasets["lr_splitseq_dual"]["splitcode_dual_pass"] is True
    assert datasets["lr_splitseq_dual"]["r1_reverse_complement"] is True
    assert datasets["lr_splitseq_forward"]["analysis_role"] == "supplementary"
    assert datasets["lr_splitseq_forward"]["splitcode_dual_pass"] is False


def test_reverse_complement_staging_preserves_ids_and_reverses_quality(tmp_path):
    source = tmp_path / "input.fastq"
    output = tmp_path / "input.rc.fastq"
    provenance = tmp_path / "input.rc.provenance.json"
    source.write_bytes(b"@SRR1.1 description\nACGTRYKN\n+comment\nABCDEFGH\n")

    payload = stage(source, output, provenance)

    assert output.read_bytes() == b"@SRR1.1 description\nNMRYACGT\n+comment\nHGFEDCBA\n"
    assert payload["records"] == 1
    assert payload["header_preserved"] is True
    assert (
        json.loads(provenance.read_text())["output"]["sha256"]
        == payload["output"]["sha256"]
    )


def test_splitcode_dual_wrapper_runs_two_passes_without_reconciliation(tmp_path):
    fake = tmp_path / "fake-splitcode"
    fake.write_text(
        "#!/usr/bin/env python3\n"
        "import pathlib, shutil, sys\n"
        "args = sys.argv[1:]\n"
        "out = pathlib.Path(args[args.index('--output') + 1])\n"
        "mapping = pathlib.Path(args[args.index('--mapping') + 1])\n"
        "shutil.copyfile(pathlib.Path(args[-1]), out)\n"
        "mapping.write_text('mapping\\n')\n"
    )
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    config = tmp_path / "config.txt"
    config.write_text("config\n")
    forward = tmp_path / "forward.fastq"
    reverse = tmp_path / "reverse.fastq"
    forward.write_text("@SRR1.1\nAC\n+\nII\n")
    reverse.write_text("@SRR1.1\nGT\n+\nII\n")
    forward_output = tmp_path / "forward.out.fastq"
    reverse_output = tmp_path / "reverse.out.fastq"
    report = tmp_path / "report.json"

    assert (
        dual_pass_main(
            [
                "--binary",
                str(fake),
                "--config",
                str(config),
                "--threads",
                "4",
                "--forward-input",
                str(forward),
                "--reverse-input",
                str(reverse),
                "--forward-output",
                str(forward_output),
                "--reverse-output",
                str(reverse_output),
                "--report",
                str(report),
            ]
        )
        == 0
    )

    payload = json.loads(report.read_text())
    assert [item["label"] for item in payload["passes"]] == [
        "forward",
        "reverse-complement",
    ]
    assert payload["duplicate_reconciliation_performed"] is False
    assert payload["duplicate_reconciliation_included_in_timing"] is False
    assert forward_output.read_text() == forward.read_text()
    assert reverse_output.read_text() == reverse.read_text()
    assert not (tmp_path / "splitcode-forward.mapping.txt").exists()
    assert not (tmp_path / "splitcode-reverse.mapping.txt").exists()
