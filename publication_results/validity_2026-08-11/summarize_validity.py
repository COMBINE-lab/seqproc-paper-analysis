#!/usr/bin/env python3
"""Summarize the frozen 32-thread materialized-output validity campaign.

The script discovers successful attempts from their identity/run records, asks
the campaign's pinned Rust auditor to encode each emitted accession set as a
bitmap, and computes overlap against the structural reference sets.  It also
reports pairwise emitted-set concordance.  Long-read splitcode dual-orientation
output is the set union of its forward and reverse-complement passes.
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "runs"
OUT = ROOT / "aggregates" / "validity"
BITMAPS = OUT / "bitmaps"
AUDITOR = Path(
    "/scratch1/seqproc-ecosystem/seqproc-paper-analysis/tools/bin/fastq-numeric-audit"
)
V_TOTAL_FILES = {
    "splitseq_pe": ROOT / "vtotal" / "splitseq_pe.ids.txt",
    "lr_splitseq_dual": ROOT / "vtotal" / "lr_splitseq.ids.txt",
    "scirnaseq3": ROOT / "vtotal" / "scirnaseq3" / "v_total_ids.txt",
}
PRIMARY_DATASETS = ("splitseq_pe", "lr_splitseq_dual", "tenx_v2", "scirnaseq3")
TOOLS = ("seqproc", "matchbox", "splitcode")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


KNOWN_CONTRACT_EXCEPTIONS = {
    ("scirnaseq3", "matchbox"): "projected barcode length may extend to 30 under indel matching",
    ("scirnaseq3", "splitcode"): "--no-outb intentionally leaves accepted-read output empty; extracted projection is the product",
}


def successful_attempts() -> dict[tuple[str, str], tuple[Path, dict, list[str]]]:
    found: dict[tuple[str, str], tuple[Path, dict, list[str]]] = {}
    for identity_path in RUNS.glob("*/attempt-*/identity.json"):
        attempt = identity_path.parent
        run_path = attempt / "run.json"
        if not run_path.exists():
            continue
        run = json.loads(run_path.read_text())
        if run.get("exit_code") != 0 or run.get("missing_outputs"):
            continue
        identity = json.loads(identity_path.read_text())
        metadata = identity["metadata"]
        key = (metadata["dataset"], metadata["tool"])
        invalid_outputs = run.get("invalid_outputs", [])
        if invalid_outputs and key not in KNOWN_CONTRACT_EXCEPTIONS:
            continue
        if key in found:
            raise RuntimeError(f"multiple successful attempts for {key}")
        found[key] = (attempt, identity, invalid_outputs)
    return found


def parse_canonical(path: Path) -> tuple[dict[str, object], bytes]:
    data = path.read_bytes()
    fields = data.split(b"\0", 4)
    if len(fields) != 5 or fields[0] != b"fastq_numeric_accession_set_v1":
        raise RuntimeError(f"unexpected canonical bitmap format: {path}")
    metadata = {
        "domain": fields[0].decode(),
        "mate": int(fields[1]),
        "numeric_id_max": int(fields[2]),
        "accession_prefix": fields[3].decode(),
    }
    return metadata, fields[4]


def materialized_paths(attempt: Path, identity: dict) -> list[tuple[Path, dict]]:
    contracts = identity["output_contract"]
    dataset = identity["metadata"]["dataset"]
    tool = identity["metadata"]["tool"]
    selected = [contracts[0]]
    if dataset == "scirnaseq3" and tool == "splitcode":
        selected = [contracts[1]]
    if dataset == "lr_splitseq_dual" and tool == "splitcode":
        # Contract zero is the dual-pass JSON report; use one extracted FASTQ
        # from each orientation (the other components have the same ID set).
        selected = [contracts[1], contracts[4]]
    resolved = []
    for contract in selected:
        path = Path(contract["path"].replace("{run_dir}", str(attempt)))
        resolved.append((path, contract))
    return resolved


def frozen_output_hashes(attempt: Path) -> dict[str, str]:
    hashes = {}
    for line in (attempt / "outputs.sha256").read_text().splitlines():
        digest, path = line.split(None, 1)
        hashes[str(Path(path.strip()))] = digest
    return hashes


def emitted_bitmap(attempt: Path, identity: dict, contract_exceptions: list[str]) -> tuple[int, dict]:
    metadata = identity["metadata"]
    dataset, tool = metadata["dataset"], metadata["tool"]
    parts: list[int] = []
    sources = []
    output_hashes = frozen_output_hashes(attempt)
    numeric_id_max = int(metadata["input_records"])
    prefix = None
    for index, (fastq, contract) in enumerate(materialized_paths(attempt, identity), 1):
        canonical = BITMAPS / f"{dataset}.{tool}.part{index}.bitmap"
        if not canonical.exists():
            command = [
                str(AUDITOR),
                str(fastq),
                str(contract["mate"]),
                str(contract["numeric_id_max"]),
                str(canonical),
            ]
            completed = subprocess.run(command, check=True, text=True, capture_output=True)
            audit = json.loads(completed.stdout)
        else:
            audit = None
        bitmap_meta, bitmap = parse_canonical(canonical)
        if bitmap_meta["numeric_id_max"] != numeric_id_max:
            raise RuntimeError(f"numeric maximum mismatch for {dataset}/{tool}")
        if len(bitmap) != (numeric_id_max + 7) // 8:
            raise RuntimeError(f"bitmap length mismatch for {dataset}/{tool}")
        if prefix is None:
            prefix = bitmap_meta["accession_prefix"]
        elif prefix != bitmap_meta["accession_prefix"]:
            raise RuntimeError(f"accession prefix mismatch for {dataset}/{tool}")
        parts.append(int.from_bytes(bitmap, "little"))
        sources.append(
            {
                "fastq": str(fastq),
                "fastq_sha256": output_hashes[str(fastq)],
                "canonical_bitmap": str(canonical),
                "canonical_bitmap_sha256": sha256(canonical),
                "audit": audit,
            }
        )
    combined = 0
    for part in parts:
        combined |= part
    return combined, {
        "attempt": str(attempt),
        "condition_id": metadata["condition_id"],
        "input_records": numeric_id_max,
        "accession_prefix": prefix,
        "identity_sha256": sha256(attempt / "identity.json"),
        "command_sha256": sha256(attempt / "command.json"),
        "repositories": identity.get("repositories", []),
        "executables": identity.get("executables", []),
        "configs": identity.get("configs", []),
        "contract_exception_reason": KNOWN_CONTRACT_EXCEPTIONS.get((dataset, tool)),
        "frozen_contract_exceptions": contract_exceptions,
        "sources": sources,
    }


def text_reference_bitmap(path: Path, numeric_id_max: int) -> tuple[int, int]:
    canonical = BITMAPS / f"reference.{path.parent.name}.{path.name}.bitmap.raw"
    if canonical.exists():
        raw = canonical.read_bytes()
    else:
        raw = bytearray((numeric_id_max + 7) // 8)
        with path.open("rb") as handle:
            for line_number, line in enumerate(handle, 1):
                token = line.strip().split(None, 1)[0]
                if not token:
                    continue
                try:
                    numeric_id = int(token.rsplit(b".", 1)[1])
                except (IndexError, ValueError) as error:
                    raise RuntimeError(f"bad reference ID at {path}:{line_number}") from error
                if not 1 <= numeric_id <= numeric_id_max:
                    raise RuntimeError(f"out-of-range reference ID {numeric_id} in {path}")
                bit = numeric_id - 1
                mask = 1 << (bit % 8)
                if raw[bit // 8] & mask:
                    raise RuntimeError(f"duplicate reference ID {numeric_id} in {path}")
                raw[bit // 8] |= mask
        canonical.write_bytes(raw)
        raw = bytes(raw)
    value = int.from_bytes(raw, "little")
    return value, popcount(value)


def complete_reference_bitmap(numeric_id_max: int) -> tuple[int, int]:
    return (1 << numeric_id_max) - 1, numeric_id_max


def ratio(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


def popcount(value: int) -> int:
    """Python 3.9-compatible population count for the campaign host."""
    return value.bit_count() if hasattr(int, "bit_count") else bin(value).count("1")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    BITMAPS.mkdir(parents=True, exist_ok=True)
    attempts = successful_attempts()
    missing = [(dataset, tool) for dataset in PRIMARY_DATASETS for tool in TOOLS if (dataset, tool) not in attempts]
    if missing:
        raise SystemExit(f"campaign incomplete; missing successful attempts: {missing}")

    emitted_sets: dict[tuple[str, str], int] = {}
    provenance = {}
    for key, (attempt, identity, contract_exceptions) in sorted(attempts.items()):
        if key[0] not in PRIMARY_DATASETS:
            continue
        emitted_sets[key], provenance[f"{key[0]}/{key[1]}"] = emitted_bitmap(
            attempt, identity, contract_exceptions
        )

    rows = []
    references = {}
    reference_sets = {}
    for dataset in PRIMARY_DATASETS:
        input_records = attempts[(dataset, TOOLS[0])][1]["metadata"]["input_records"]
        if dataset == "tenx_v2":
            reference, reference_count = complete_reference_bitmap(input_records)
            reference_info = {
                "definition": "all input pairs; every R1 is at least 26 nt",
                "count": reference_count,
            }
        else:
            reference_path = V_TOTAL_FILES[dataset]
            if not reference_path.exists():
                raise SystemExit(f"missing reference set: {reference_path}")
            reference, reference_count = text_reference_bitmap(reference_path, input_records)
            reference_info = {
                "path": str(reference_path),
                "sha256": sha256(reference_path),
                "count": reference_count,
            }
        references[dataset] = reference_info
        reference_sets[dataset] = reference
        for tool in TOOLS:
            emitted = emitted_sets[(dataset, tool)]
            emitted_count = popcount(emitted)
            intersection = popcount(emitted & reference)
            precision = ratio(intersection, emitted_count)
            recall = ratio(intersection, reference_count)
            f1 = ratio(2 * precision * recall, precision + recall)
            rows.append(
                {
                    "dataset": dataset,
                    "tool": tool,
                    "input_records": input_records,
                    "reference_records": reference_count,
                    "emitted_records": emitted_count,
                    "intersection_records": intersection,
                    "emitted_fraction": ratio(emitted_count, input_records),
                    "precision": precision,
                    "recall": recall,
                    "f1": f1,
                }
            )

    pairwise = []
    for dataset in PRIMARY_DATASETS:
        for left, right in combinations(TOOLS, 2):
            left_set = emitted_sets[(dataset, left)]
            right_set = emitted_sets[(dataset, right)]
            intersection = popcount(left_set & right_set)
            union = popcount(left_set | right_set)
            pairwise.append(
                {
                    "dataset": dataset,
                    "left": left,
                    "right": right,
                    "intersection_records": intersection,
                    "union_records": union,
                    "left_only_records": popcount(left_set & ~right_set),
                    "right_only_records": popcount(right_set & ~left_set),
                    "jaccard": ratio(intersection, union),
                }
            )

    venn = []
    for dataset in PRIMARY_DATASETS:
        seqproc = emitted_sets[(dataset, "seqproc")]
        matchbox = emitted_sets[(dataset, "matchbox")]
        splitcode = emitted_sets[(dataset, "splitcode")]
        union = seqproc | matchbox | splitcode
        categories = {
            "all_three": seqproc & matchbox & splitcode,
            "seqproc_splitcode": seqproc & splitcode & ~matchbox,
            "seqproc_matchbox": seqproc & matchbox & ~splitcode,
            "splitcode_matchbox": splitcode & matchbox & ~seqproc,
            "seqproc_only": seqproc & ~matchbox & ~splitcode,
            "splitcode_only": splitcode & ~seqproc & ~matchbox,
            "matchbox_only": matchbox & ~seqproc & ~splitcode,
        }
        row = {"dataset": dataset, "union_records": popcount(union)}
        for name, value in categories.items():
            count = popcount(value)
            row[f"{name}_records"] = count
            row[f"{name}_fraction_of_union"] = ratio(count, row["union_records"])
        splitcode_exclusive = splitcode & ~(seqproc | matchbox)
        splitcode_exclusive_valid = popcount(splitcode_exclusive & reference_sets[dataset])
        row["splitcode_exclusive_valid_records"] = splitcode_exclusive_valid
        row["splitcode_exclusive_invalid_fraction"] = 1.0 - ratio(
            splitcode_exclusive_valid, popcount(splitcode_exclusive)
        )
        venn.append(row)

    payload = {
        "schema_version": "1.0.0",
        "campaign_root": str(ROOT),
        "manifest": {
            "path": str(ROOT / "specs" / "publication-core-correctness-t32.yaml"),
            "sha256": sha256(ROOT / "specs" / "publication-core-correctness-t32.yaml"),
        },
        "schedule": {
            "path": str(ROOT / "specs" / "publication-core-correctness-t32.schedule.json"),
            "sha256": sha256(ROOT / "specs" / "publication-core-correctness-t32.schedule.json"),
        },
        "auditor": {"path": str(AUDITOR), "sha256": sha256(AUDITOR)},
        "summarizer": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
        "references": references,
        "metrics": rows,
        "pairwise": pairwise,
        "venn": venn,
        "provenance": provenance,
    }
    (OUT / "validity_metrics.json").write_text(json.dumps(payload, indent=2) + "\n")
    with (OUT / "validity_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (OUT / "pairwise_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pairwise[0]))
        writer.writeheader()
        writer.writerows(pairwise)
    with (OUT / "venn_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(venn[0]))
        writer.writeheader()
        writer.writerows(venn)

    print("dataset\ttool\temitted\temitted%\tintersection\tprecision%\trecall%\tF1")
    for row in rows:
        print(
            f"{row['dataset']}\t{row['tool']}\t{row['emitted_records']}\t"
            f"{100 * row['emitted_fraction']:.4f}\t{row['intersection_records']}\t"
            f"{100 * row['precision']:.4f}\t{100 * row['recall']:.4f}\t{row['f1']:.6f}"
        )


if __name__ == "__main__":
    main()
