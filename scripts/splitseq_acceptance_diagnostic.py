#!/usr/bin/env python3
"""Localize SPLiT-seq accepted-set changes to individual geometry stages.

This diagnostic constructs cumulative geometries ending after each barcode or
linker operation, runs each supplied seqproc binary on an identical paired
FASTQ prefix, and records output read names.  Pairwise set differences identify
the first stage at which two implementations disagree.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


STAGES = (
    ("bc3", 'bc3 = filter_within_dist(b[8], "{bc23}", 1)'),
    ("linker1", "#[search(relative)] #[edit(6)] l1 = f[GTGGCCGCTGTTTCGCATCGGCGTACGACT]"),
    ("bc2", 'bc2 = filter_within_dist(b[8], "{bc23}", 1)'),
    ("linker2", "#[search(relative)] #[edit(6)] l2 = f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG]"),
    ("bc1", 'bc1 = filter_within_dist(b[6], "{bc1}", 1)'),
)


def labeled_path(value: str) -> tuple[str, Path]:
    label, separator, path = value.partition("=")
    if not separator or not label or not path:
        raise argparse.ArgumentTypeError("expected LABEL=PATH")
    return label, Path(path).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def geometry_for(stage_index: int, bc23: Path, bc1: Path) -> str:
    definitions = ["read1 = r:", "skip2 = x[2]", "umi = u[8]"]
    pieces = ["<skip2>", "<umi>"]
    for name, definition in STAGES[: stage_index + 1]:
        definitions.append(definition.format(bc23=bc23, bc1=bc1))
        pieces.append(f"<{name if not name.startswith('linker') else 'l' + name[-1]}>")
    definitions.append("rest = r:")
    pieces.append("<rest>")
    output_labels = ["<umi>"] + [
        f"<{name}>" for name, _ in STAGES[: stage_index + 1] if name.startswith("bc")
    ]
    return "\n".join(
        definitions
        + [
            "",
            "1{<read1>}",
            f"2{{{''.join(pieces)}}}",
            "",
            f"-> 1{{<read1>}} 2{{{''.join(output_labels)}}}",
            "",
        ]
    )


def read_names(path: Path) -> set[str]:
    names: set[str] = set()
    with path.open("rb") as handle:
        for line_number, line in enumerate(handle):
            if line_number % 4 == 0:
                names.add(line[1:].split(None, 1)[0].decode("utf-8"))
    return names


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", action="append", type=labeled_path, required=True)
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r2", type=Path, required=True)
    parser.add_argument("--bc23", type=Path, required=True)
    parser.add_argument("--bc1", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    binaries = dict(args.binary)
    if len(binaries) != 2:
        parser.error("supply exactly two uniquely labeled binaries")
    args.outdir.mkdir(parents=True, exist_ok=True)
    bc23 = args.bc23.resolve()
    bc1 = args.bc1.resolve()
    results: dict[str, object] = {
        "inputs": {
            "r1": {"path": str(args.r1.resolve()), "sha256": sha256(args.r1)},
            "r2": {"path": str(args.r2.resolve()), "sha256": sha256(args.r2)},
            "bc23": {"path": str(bc23), "sha256": sha256(bc23)},
            "bc1": {"path": str(bc1), "sha256": sha256(bc1)},
        },
        "binaries": {
            label: {"path": str(path), "sha256": sha256(path)}
            for label, path in binaries.items()
        },
        "stages": [],
    }

    shared_prior = read_names(args.r1)
    for stage_index, (stage_name, _) in enumerate(STAGES):
        geometry = args.outdir / f"stage-{stage_index + 1}-{stage_name}.geom"
        geometry.write_text(geometry_for(stage_index, bc23, bc1))
        names_by_version: dict[str, set[str]] = {}
        stage_result: dict[str, object] = {
            "stage": stage_name,
            "geometry": str(geometry),
            "geometry_sha256": sha256(geometry),
            "versions": {},
        }
        for label, binary in binaries.items():
            output1 = args.outdir / f"{label}-{stage_name}-r1.fastq"
            output2 = args.outdir / f"{label}-{stage_name}-r2.fastq"
            command = [
                str(binary), "--geom", str(geometry),
                "--file1", str(args.r1.resolve()), "--file2", str(args.r2.resolve()),
                "--out1", str(output1), "--out2", str(output2), "--threads", "1",
            ]
            completed = subprocess.run(command, text=True, capture_output=True)
            if completed.returncode:
                raise RuntimeError(f"failed: {command}\n{completed.stderr}")
            names = read_names(output1)
            names_by_version[label] = names
            stage_result["versions"][label] = {
                "accepted": len(names),
                "output1_sha256": sha256(output1),
                "output2_sha256": sha256(output2),
                "stderr": completed.stderr,
            }
            output1.unlink()
            output2.unlink()

        first, second = binaries
        first_only = sorted(names_by_version[first] - names_by_version[second])
        second_only = sorted(names_by_version[second] - names_by_version[first])
        stage_result["comparison"] = {
            "intersection": len(names_by_version[first] & names_by_version[second]),
            f"{first}_only": len(first_only),
            f"{second}_only": len(second_only),
            f"{first}_only_examples": first_only[:20],
            f"{second}_only_examples": second_only[:20],
        }
        first_conditional = names_by_version[first] & shared_prior
        second_conditional = names_by_version[second] & shared_prior
        stage_result["conditional_on_shared_prior"] = {
            "shared_prior": len(shared_prior),
            f"{first}_accepted": len(first_conditional),
            f"{second}_accepted": len(second_conditional),
            f"{first}_only": len(first_conditional - second_conditional),
            f"{second}_only": len(second_conditional - first_conditional),
            f"{first}_only_examples": sorted(first_conditional - second_conditional)[:20],
            f"{second}_only_examples": sorted(second_conditional - first_conditional)[:20],
        }
        shared_prior = names_by_version[first] & names_by_version[second]
        results["stages"].append(stage_result)

    destination = args.outdir / "splitseq-acceptance-diagnostic.json"
    destination.write_text(json.dumps(results, indent=2) + "\n")
    print(destination)


if __name__ == "__main__":
    main()
