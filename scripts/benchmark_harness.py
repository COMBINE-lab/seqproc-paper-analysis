#!/usr/bin/env python3
"""Content-addressed benchmark runner for development and release experiments.

The runner deliberately accepts argv arrays rather than shell command strings.
It records enough provenance to compare a baseline and a candidate commit while
the implementation is still changing; it does not require a final release tag.

Run specifications are JSON objects. See ``--write-example`` for the schema's
smallest useful example.
"""

from __future__ import annotations

import argparse
import base64
import concurrent.futures
import hashlib
import gzip
import heapq
import json
import os
import platform
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCHEMA_VERSION = "1.4.0"
DEFAULT_ENV_ALLOWLIST = (
    "PATH",
    "LD_LIBRARY_PATH",
    "LIBRARY_PATH",
    "LANG",
    "LC_ALL",
    "TMPDIR",
    "RUST_BACKTRACE",
    "RUST_LOG",
    "RAYON_NUM_THREADS",
    "ANTISEQ_CHUNK_SIZE",
)


class HarnessError(RuntimeError):
    """A specification, provenance, execution, or validation error."""


def canonical_json(value: Any) -> bytes:
    """Serialize *value* deterministically for IDs and checksums."""
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve(path: str | os.PathLike[str], cwd: Path) -> Path:
    resolved = Path(path)
    if not resolved.is_absolute():
        resolved = cwd / resolved
    return resolved.resolve()


def describe_input(entry: str | Mapping[str, Any], cwd: Path) -> dict[str, Any]:
    """Resolve and describe an input/config file.

    A mapping may provide a previously verified ``sha256``. This avoids hashing
    a multi-gigabyte FASTQ for every replicate. Set ``verify: true`` to hash and
    compare it during this invocation.
    """
    if isinstance(entry, str):
        item: dict[str, Any] = {"path": entry, "verify": True}
    elif isinstance(entry, Mapping):
        item = dict(entry)
    else:
        raise HarnessError(f"file entry must be a string or object, got {entry!r}")

    if "path" not in item:
        raise HarnessError(f"file entry has no path: {entry!r}")
    path = _resolve(str(item["path"]), cwd)
    if not path.is_file():
        raise HarnessError(f"required file does not exist: {path}")

    declared = item.get("sha256")
    verify = bool(item.get("verify", declared is None))
    observed = sha256_file(path) if verify or declared is None else str(declared)
    if declared is not None and verify and observed != declared:
        raise HarnessError(
            f"SHA-256 mismatch for {path}: expected {declared}, observed {observed}"
        )
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": observed,
        "digest_source": "verified" if verify else "declared",
    }


def _git(path: Path, *args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", "-C", str(path), *args],
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def describe_repository(
    entry: Mapping[str, Any], cwd: Path
) -> tuple[dict[str, Any], bytes, dict[str, bytes]]:
    name = str(entry.get("name", "repository"))
    if "path" not in entry:
        raise HarnessError(f"repository {name!r} has no path")
    path = _resolve(str(entry["path"]), cwd)
    try:
        commit = _git(path, "rev-parse", "HEAD").stdout.decode().strip()
        status = _git(path, "status", "--porcelain=v1").stdout.decode()
        diff = _git(path, "diff", "--binary", "HEAD").stdout
        untracked_output = _git(
            path, "ls-files", "--others", "--exclude-standard", "-z"
        ).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        raise HarnessError(f"cannot inspect git repository {path}: {error}") from error

    dirty = bool(status)
    expected_commit = entry.get("commit")
    if expected_commit is not None and commit != str(expected_commit):
        raise HarnessError(
            f"repository {name!r} is at {commit}, expected {expected_commit}"
        )
    if dirty and not bool(entry.get("allow_dirty", False)):
        raise HarnessError(
            f"repository {name!r} is dirty; commit it or set allow_dirty=true "
            "to archive its patch with this development run"
        )
    untracked: dict[str, bytes] = {}
    untracked_records: list[dict[str, Any]] = []
    for raw_relative in untracked_output.split(b"\0"):
        if not raw_relative:
            continue
        relative = raw_relative.decode("utf-8")
        untracked_path = (path / relative).resolve()
        try:
            untracked_path.relative_to(path)
        except ValueError as error:
            raise HarnessError(f"untracked path escapes repository: {relative}") from error
        if not untracked_path.is_file():
            continue
        content = untracked_path.read_bytes()
        untracked[relative] = content
        untracked_records.append(
            {
                "path": relative,
                "bytes": len(content),
                "sha256": sha256_bytes(content),
            }
        )

    state = {
        "name": name,
        "path": str(path),
        "commit": commit,
        "dirty": dirty,
        "status": status.splitlines(),
        "diff_sha256": sha256_bytes(diff),
        "untracked": untracked_records,
    }
    return state, diff, untracked


def resolve_executable(command0: str, cwd: Path, env: Mapping[str, str]) -> Path:
    candidate = Path(command0)
    if candidate.is_absolute() or "/" in command0:
        path = _resolve(command0, cwd)
    else:
        found = shutil.which(command0, path=env.get("PATH"))
        if found is None:
            raise HarnessError(f"executable not found: {command0}")
        path = Path(found).resolve()
    if not path.is_file():
        raise HarnessError(f"executable is not a file: {path}")
    return path


def describe_executable(
    entry: str | Mapping[str, Any], cwd: Path, env: Mapping[str, str]
) -> dict[str, Any]:
    if isinstance(entry, str):
        item: dict[str, Any] = {"path": entry, "verify": True}
    elif isinstance(entry, Mapping):
        item = dict(entry)
    else:
        raise HarnessError(f"executable entry must be a string or object, got {entry!r}")
    if "path" not in item:
        raise HarnessError(f"executable entry has no path: {entry!r}")
    path = resolve_executable(str(item["path"]), cwd, env)
    declared = item.get("sha256")
    verify = bool(item.get("verify", declared is None))
    observed = sha256_file(path) if verify or declared is None else str(declared)
    if declared is not None and verify and observed != declared:
        raise HarnessError(
            f"SHA-256 mismatch for executable {path}: expected {declared}, observed {observed}"
        )
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": observed,
        "digest_source": "verified" if verify else "declared",
    }


def build_environment(spec: Mapping[str, Any]) -> tuple[dict[str, str], dict[str, str]]:
    mode = spec.get("environment_mode", "allowlist")
    overrides = {str(k): str(v) for k, v in spec.get("environment", {}).items()}
    if mode == "inherit":
        execution = dict(os.environ)
    elif mode == "allowlist":
        allowlist = tuple(spec.get("environment_allowlist", DEFAULT_ENV_ALLOWLIST))
        execution = {key: os.environ[key] for key in allowlist if key in os.environ}
    else:
        raise HarnessError("environment_mode must be 'allowlist' or 'inherit'")
    execution.update(overrides)
    execution["LC_ALL"] = "C"

    # Record only the allowlisted/explicit environment, never unrelated secrets
    # inherited from an interactive shell or scheduler.
    recorded_keys = set(spec.get("environment_allowlist", DEFAULT_ENV_ALLOWLIST))
    recorded_keys.update(overrides)
    recorded_keys.add("LC_ALL")
    recorded = {key: execution[key] for key in sorted(recorded_keys) if key in execution}
    return execution, recorded


def host_snapshot() -> dict[str, Any]:
    snapshot: dict[str, Any] = {
        "hostname": platform.node(),
        "platform": platform.platform(),
        "python": sys.version,
        "cpu_count": os.cpu_count(),
        "process_affinity": sorted(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else None,
        "load_average": list(os.getloadavg()) if hasattr(os, "getloadavg") else None,
    }
    for name, path in {
        "cpu_governor": Path("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
        "cpu_boost": Path("/sys/devices/system/cpu/cpufreq/boost"),
        "smt_active": Path("/sys/devices/system/cpu/smt/active"),
    }.items():
        try:
            snapshot[name] = path.read_text().strip()
        except OSError:
            snapshot[name] = None
    try:
        memory_fields = {}
        for line in Path("/proc/meminfo").read_text().splitlines():
            key, value = line.split(":", 1)
            if key in {"MemTotal", "MemAvailable", "SwapTotal", "SwapFree"}:
                memory_fields[key] = value.strip()
        snapshot["memory"] = memory_fields
    except OSError:
        snapshot["memory"] = None
    try:
        result = subprocess.run(
            ["lscpu", "--json"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        snapshot["lscpu"] = json.loads(result.stdout)
    except (OSError, subprocess.CalledProcessError, json.JSONDecodeError):
        snapshot["lscpu"] = None
    try:
        result = subprocess.run(
            ["findmnt", "--json", "--target", "/scratch1"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        snapshot["scratch_mount"] = json.loads(result.stdout)
    except (OSError, subprocess.CalledProcessError, json.JSONDecodeError):
        snapshot["scratch_mount"] = None
    return snapshot


def parse_gnu_time(path: Path) -> dict[str, Any]:
    fields: dict[str, Any] = {}
    if not path.is_file():
        return fields
    conversions = {
        "User time (seconds)": ("user_cpu_seconds", float),
        "System time (seconds)": ("system_cpu_seconds", float),
        "Maximum resident set size (kbytes)": ("peak_rss_kib", int),
        "File system inputs": ("filesystem_inputs", int),
        "File system outputs": ("filesystem_outputs", int),
        "Exit status": ("time_exit_status", int),
    }
    for line in path.read_text(errors="replace").splitlines():
        stripped = line.strip()
        for label, (key, converter) in conversions.items():
            prefix = f"{label}:"
            if stripped.startswith(prefix):
                try:
                    fields[key] = converter(stripped[len(prefix) :].strip())
                except ValueError:
                    fields[key] = None
                break
    return fields


def _next_attempt(run_root: Path) -> tuple[int, Path]:
    run_root.mkdir(parents=True, exist_ok=True)
    existing = sorted(
        int(path.name.split("-", 1)[1])
        for path in run_root.glob("attempt-*")
        if path.is_dir() and path.name.split("-", 1)[1].isdigit()
    )
    attempt = (existing[-1] + 1) if existing else 1
    path = run_root / f"attempt-{attempt:04d}"
    path.mkdir()
    return attempt, path


def _substitute_run_dir(values: Iterable[str], run_dir: Path) -> list[str]:
    return [str(value).replace("{run_dir}", str(run_dir)) for value in values]


def _prepare_outputs(
    entries: Iterable[str | Mapping[str, Any]],
    cwd: Path,
    run_dir: Path,
    default_retain: bool,
) -> list[dict[str, Any]]:
    outputs: list[dict[str, Any]] = []
    for entry in entries:
        if isinstance(entry, str):
            item: dict[str, Any] = {"path": entry}
        elif isinstance(entry, Mapping):
            item = dict(entry)
        else:
            raise HarnessError(f"output entry must be a string or object, got {entry!r}")
        if "path" not in item:
            raise HarnessError(f"output entry has no path: {entry!r}")
        substituted = str(item["path"]).replace("{run_dir}", str(run_dir))
        output_format = item.get("format")
        if output_format not in (None, "fastq"):
            raise HarnessError(f"unsupported output format {output_format!r}")
        normalization = item.get("normalize")
        if normalization not in (
            None,
            "fastq_multiset",
            "fastq_id_multiset",
            "fastq_numeric_accession_set",
        ):
            raise HarnessError(f"unsupported output normalization {normalization!r}")
        if normalization is not None and output_format != "fastq":
            raise HarnessError("FASTQ normalization requires format: fastq")
        mate = int(item.get("mate", 0))
        if mate < 0:
            raise HarnessError("output mate cannot be negative")
        numeric_id_max = item.get("numeric_id_max")
        numeric_audit_executable = item.get("numeric_audit_executable")
        if normalization == "fastq_numeric_accession_set":
            if numeric_id_max is None or int(numeric_id_max) <= 0:
                raise HarnessError(
                    "fastq_numeric_accession_set requires a positive numeric_id_max"
                )
            numeric_id_max = int(numeric_id_max)
            if numeric_audit_executable is not None:
                numeric_audit_executable = _resolve(
                    str(numeric_audit_executable), cwd
                )
        elif numeric_id_max is not None:
            raise HarnessError("numeric_id_max requires fastq_numeric_accession_set")
        elif numeric_audit_executable is not None:
            raise HarnessError(
                "numeric_audit_executable requires fastq_numeric_accession_set"
            )
        min_bytes = int(item.get("min_bytes", 0))
        if min_bytes < 0:
            raise HarnessError("output min_bytes cannot be negative")
        outputs.append(
            {
                "path": _resolve(substituted, cwd),
                "min_bytes": min_bytes,
                "format": output_format,
                "normalize": normalization,
                "mate": mate,
                "numeric_id_max": numeric_id_max,
                "numeric_audit_executable": numeric_audit_executable,
                "retain": bool(item.get("retain", default_retain)),
            }
        )
    return outputs


def inspect_fastq(path: Path) -> dict[str, Any]:
    """Count and structurally validate a four-line FASTQ, optionally gzip-compressed."""
    records = 0
    for _ in _fastq_records(path):
        records += 1
    return {"format": "fastq", "records": records, "valid": True}


def _fastq_records_from_handle(handle, path: Path):
    record = 0
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline()
        separator = handle.readline()
        quality = handle.readline()
        record += 1
        if not sequence or not separator or not quality:
            raise HarnessError(f"truncated FASTQ record {record} in {path}")
        if not header.startswith(b"@") or not separator.startswith(b"+"):
            raise HarnessError(f"malformed FASTQ record {record} in {path}")
        if len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n")):
            raise HarnessError(
                f"sequence/quality length mismatch in FASTQ record {record} in {path}"
            )
        yield header, sequence, separator, quality


def _fastq_records(path: Path):
    opener = gzip.open if path.name.endswith((".gz", ".gzip")) else open
    with opener(path, "rb") as handle:
        yield from _fastq_records_from_handle(handle, path)


class _HashingReader:
    """Minimal read-only proxy that hashes the exact bytes consumed."""

    def __init__(self, handle, digest):
        self._handle = handle
        self._digest = digest

    def read(self, size: int = -1) -> bytes:
        value = self._handle.read(size)
        self._digest.update(value)
        return value

    def readline(self, size: int = -1) -> bytes:
        value = self._handle.readline(size)
        self._digest.update(value)
        return value

    def seek(self, *args):
        return self._handle.seek(*args)

    def tell(self) -> int:
        return self._handle.tell()

    def seekable(self) -> bool:
        return self._handle.seekable()

    def readable(self) -> bool:
        return True


def _normalized_fastq_line(lines: tuple[bytes, ...], mate: int) -> bytes:
    header = lines[0].rstrip(b"\r\n")
    normalized_id = header[1:].split(None, 1)[0]
    if normalized_id.endswith((b"/1", b"/2")):
        normalized_id = normalized_id[:-2]
    record = b"".join(lines)
    return (
        base64.b64encode(normalized_id)
        + b"\t"
        + str(mate).encode()
        + b"\t"
        + base64.b64encode(record)
        + b"\n"
    )


def _normalized_fastq_id_line(lines: tuple[bytes, ...], mate: int) -> bytes:
    header = lines[0].rstrip(b"\r\n")
    normalized_id = header[1:].split(None, 1)[0]
    if normalized_id.endswith((b"/1", b"/2")):
        normalized_id = normalized_id[:-2]
    return base64.b64encode(normalized_id) + b"\t" + str(mate).encode() + b"\n"


def _externally_sorted_fastq_digest(
    path: Path,
    mate: int,
    encoder,
    chunk_bytes: int,
    temp_dir: Path | None,
) -> str:
    if chunk_bytes <= 0:
        raise HarnessError("normalization chunk_bytes must be positive")
    with tempfile.TemporaryDirectory(
        prefix="seqproc-normalize-",
        dir=None if temp_dir is None else str(temp_dir),
    ) as temporary:
        root = Path(temporary)
        chunks: list[Path] = []
        buffered: list[bytes] = []
        buffered_bytes = 0

        def flush_chunk() -> None:
            nonlocal buffered_bytes
            if not buffered:
                return
            buffered.sort()
            chunk = root / f"chunk-{len(chunks):06d}"
            with chunk.open("wb") as handle:
                handle.writelines(buffered)
            chunks.append(chunk)
            buffered.clear()
            buffered_bytes = 0

        for lines in _fastq_records(path):
            encoded = encoder(lines, mate)
            buffered.append(encoded)
            buffered_bytes += len(encoded)
            if buffered_bytes >= chunk_bytes:
                flush_chunk()
        flush_chunk()

        digest = hashlib.sha256()
        handles = [chunk.open("rb") for chunk in chunks]
        try:
            for line in heapq.merge(*handles):
                digest.update(line)
        finally:
            for handle in handles:
                handle.close()
        return digest.hexdigest()


def inspect_and_normalize_fastq(
    path: Path,
    mate: int = 0,
    normalization: str | None = None,
    chunk_bytes: int = 64 * 1024 * 1024,
    temp_dir: Path | None = None,
    numeric_id_max: int | None = None,
    numeric_audit_executable: Path | None = None,
) -> dict[str, Any]:
    """Hash, validate, count, and optionally normalize a FASTQ in one input pass."""
    if normalization not in (
        None,
        "fastq_multiset",
        "fastq_id_multiset",
        "fastq_numeric_accession_set",
    ):
        raise HarnessError(f"unsupported output normalization {normalization!r}")
    if chunk_bytes <= 0:
        raise HarnessError("normalization chunk_bytes must be positive")
    if (
        normalization == "fastq_numeric_accession_set"
        and numeric_audit_executable is not None
        and not path.name.endswith((".gz", ".gzip"))
    ):
        return _inspect_with_numeric_auditor(
            path,
            mate,
            numeric_id_max,
            numeric_audit_executable,
            temp_dir,
        )

    encoder = {
        "fastq_multiset": _normalized_fastq_line,
        "fastq_id_multiset": _normalized_fastq_id_line,
    }.get(normalization)
    numeric_bitmap: bytearray | None = None
    numeric_prefix: bytes | None = None
    if normalization == "fastq_numeric_accession_set":
        if numeric_id_max is None or numeric_id_max <= 0:
            raise HarnessError(
                "fastq_numeric_accession_set requires a positive numeric_id_max"
            )
        numeric_bitmap = bytearray((numeric_id_max + 8) // 8)
    raw_digest = hashlib.sha256()
    records = 0

    with tempfile.TemporaryDirectory(
        prefix="seqproc-normalize-",
        dir=None if temp_dir is None else str(temp_dir),
    ) as temporary:
        root = Path(temporary)
        chunks: list[Path] = []
        buffered: list[bytes] = []
        buffered_bytes = 0

        def flush_chunk() -> None:
            nonlocal buffered_bytes
            if not buffered:
                return
            buffered.sort()
            chunk = root / f"chunk-{len(chunks):06d}"
            with chunk.open("wb") as handle:
                handle.writelines(buffered)
            chunks.append(chunk)
            buffered.clear()
            buffered_bytes = 0

        with path.open("rb") as raw_handle:
            hashing_handle = _HashingReader(raw_handle, raw_digest)
            if path.name.endswith((".gz", ".gzip")):
                stream = gzip.GzipFile(fileobj=hashing_handle, mode="rb")
            else:
                stream = hashing_handle
            try:
                for lines in _fastq_records_from_handle(stream, path):
                    records += 1
                    if encoder is not None:
                        encoded = encoder(lines, mate)
                        buffered.append(encoded)
                        buffered_bytes += len(encoded)
                        if buffered_bytes >= chunk_bytes:
                            flush_chunk()
                    elif numeric_bitmap is not None:
                        header_id = lines[0][1:].split(None, 1)[0]
                        if header_id.endswith((b"/1", b"/2")):
                            header_id = header_id[:-2]
                        prefix, separator, numeric_text = header_id.rpartition(b".")
                        if not separator or not prefix or not numeric_text.isdigit():
                            raise HarnessError(
                                f"FASTQ record {records} does not have a numeric accession ID: "
                                f"{header_id!r}"
                            )
                        if numeric_prefix is None:
                            numeric_prefix = prefix
                        elif prefix != numeric_prefix:
                            raise HarnessError(
                                f"FASTQ record {records} changes accession prefix from "
                                f"{numeric_prefix!r} to {prefix!r}"
                            )
                        numeric_id = int(numeric_text)
                        if numeric_id <= 0 or numeric_id > numeric_id_max:
                            raise HarnessError(
                                f"FASTQ record {records} numeric ID {numeric_id} is outside "
                                f"1..{numeric_id_max}"
                            )
                        byte_index, bit_index = divmod(numeric_id - 1, 8)
                        mask = 1 << bit_index
                        if numeric_bitmap[byte_index] & mask:
                            raise HarnessError(
                                f"duplicate numeric accession ID {header_id!r} in {path}"
                            )
                        numeric_bitmap[byte_index] |= mask
            finally:
                if stream is not hashing_handle:
                    stream.close()
        flush_chunk()

        result: dict[str, Any] = {
            "format": "fastq",
            "records": records,
            "valid": True,
            "sha256": raw_digest.hexdigest(),
        }
        if encoder is not None:
            normalized_digest = hashlib.sha256()
            handles = [chunk.open("rb") for chunk in chunks]
            try:
                for line in heapq.merge(*handles):
                    normalized_digest.update(line)
            finally:
                for handle in handles:
                    handle.close()
            result["normalized_sha256"] = normalized_digest.hexdigest()
        elif numeric_bitmap is not None:
            normalized_digest = hashlib.sha256()
            normalized_digest.update(b"fastq_numeric_accession_set_v1\0")
            normalized_digest.update(str(mate).encode("ascii") + b"\0")
            normalized_digest.update(str(numeric_id_max).encode("ascii") + b"\0")
            normalized_digest.update(numeric_prefix or b"")
            normalized_digest.update(b"\0")
            normalized_digest.update(numeric_bitmap)
            result["normalized_sha256"] = normalized_digest.hexdigest()
        return result


def _inspect_with_numeric_auditor(
    path: Path,
    mate: int,
    numeric_id_max: int | None,
    executable: Path,
    temp_dir: Path | None,
) -> dict[str, Any]:
    """Run the pinned compiled validator while hashing raw FASTQ bytes in parallel."""
    if numeric_id_max is None or numeric_id_max <= 0:
        raise HarnessError(
            "fastq_numeric_accession_set requires a positive numeric_id_max"
        )
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise HarnessError(f"numeric FASTQ auditor is not executable: {executable}")
    temporary = tempfile.NamedTemporaryFile(
        prefix="seqproc-numeric-set-",
        dir=None if temp_dir is None else str(temp_dir),
        delete=False,
    )
    canonical_path = Path(temporary.name)
    temporary.close()
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
            raw_digest_future = executor.submit(sha256_file, path)
            process = subprocess.run(
                [
                    str(executable),
                    str(path),
                    str(mate),
                    str(numeric_id_max),
                    str(canonical_path),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            raw_digest = raw_digest_future.result()
        if process.returncode != 0:
            detail = process.stderr.strip() or process.stdout.strip()
            raise HarnessError(
                f"numeric FASTQ audit failed for {path}: {detail or process.returncode}"
            )
        try:
            payload = json.loads(process.stdout)
            records = int(payload["records"])
        except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
            raise HarnessError(
                f"numeric FASTQ auditor returned invalid JSON for {path}: "
                f"{process.stdout!r}"
            ) from error
        return {
            "format": "fastq",
            "records": records,
            "valid": True,
            "sha256": raw_digest,
            "normalized_sha256": sha256_file(canonical_path),
            "validator": "fastq-numeric-audit-v1",
        }
    finally:
        canonical_path.unlink(missing_ok=True)


def normalized_fastq_multiset_sha256(
    path: Path,
    mate: int = 0,
    chunk_bytes: int = 64 * 1024 * 1024,
    temp_dir: Path | None = None,
) -> str:
    """Hash complete FASTQ records after a bounded-memory external sort."""
    return _externally_sorted_fastq_digest(
        path, mate, _normalized_fastq_line, chunk_bytes, temp_dir
    )


def normalized_fastq_id_multiset_sha256(
    path: Path,
    mate: int = 0,
    chunk_bytes: int = 64 * 1024 * 1024,
    temp_dir: Path | None = None,
) -> str:
    """Hash normalized read IDs after a bounded-memory external sort."""
    return _externally_sorted_fastq_digest(
        path, mate, _normalized_fastq_id_line, chunk_bytes, temp_dir
    )


def _terminate_process_group(process: subprocess.Popen[bytes], sig: int) -> None:
    try:
        os.killpg(process.pid, sig)
    except ProcessLookupError:
        pass


def prepare_spec(spec: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(spec.get("command"), list) or not spec["command"]:
        raise HarnessError("command must be a non-empty argv array")
    if any(not isinstance(token, str) for token in spec["command"]):
        raise HarnessError("every command token must be a string")
    cwd = Path(str(spec.get("cwd", "."))).resolve()
    if not cwd.is_dir():
        raise HarnessError(f"working directory does not exist: {cwd}")

    execution_env, recorded_env = build_environment(spec)
    repositories: list[dict[str, Any]] = []
    repository_snapshots: dict[str, dict[str, Any]] = {}
    for repo_entry in spec.get("repositories", []):
        state, diff, untracked = describe_repository(repo_entry, cwd)
        repositories.append(state)
        if state["dirty"]:
            repository_snapshots[state["name"]] = {
                "diff": diff,
                "untracked": untracked,
            }

    inputs = [describe_input(item, cwd) for item in spec.get("inputs", [])]
    configs = [describe_input(item, cwd) for item in spec.get("configs", [])]
    binary = describe_executable(spec["command"][0], cwd, execution_env)
    executables = [binary]
    for entry in spec.get("executables", []):
        described = describe_executable(entry, cwd, execution_env)
        if described["path"] not in {item["path"] for item in executables}:
            executables.append(described)
    identity = {
        "schema_version": SCHEMA_VERSION,
        "name": str(spec.get("name", "benchmark")),
        "command_template": list(spec["command"]),
        "cwd": str(cwd),
        "binary": binary,
        "executables": executables,
        "inputs": inputs,
        "configs": configs,
        "repositories": repositories,
        "environment": recorded_env,
        "output_contract": spec.get("outputs", []),
        "retain_outputs": bool(spec.get("retain_outputs", True)),
        "timeout_seconds": spec.get("timeout_seconds"),
        "metadata": spec.get("metadata", {}),
    }
    return {
        "identity": identity,
        "run_id": sha256_bytes(canonical_json(identity)),
        "cwd": cwd,
        "execution_env": execution_env,
        "repository_snapshots": repository_snapshots,
    }


def execute_spec(spec: Mapping[str, Any], output_root: Path) -> tuple[dict[str, Any], Path]:
    prepared = prepare_spec(spec)
    run_id = prepared["run_id"]
    run_root = output_root.resolve() / run_id
    attempt, run_dir = _next_attempt(run_root)
    command = _substitute_run_dir(spec["command"], run_dir)
    outputs = _prepare_outputs(
        spec.get("outputs", []),
        prepared["cwd"],
        run_dir,
        bool(spec.get("retain_outputs", True)),
    )
    stale_outputs = [str(item["path"]) for item in outputs if item["path"].exists()]
    if stale_outputs:
        raise HarnessError(
            "refusing to benchmark with pre-existing declared outputs: "
            + ", ".join(stale_outputs)
        )

    for name, snapshot in prepared["repository_snapshots"].items():
        safe_name = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in name)
        (run_dir / f"repository-{safe_name}.diff").write_bytes(snapshot["diff"])
        untracked_root = run_dir / f"repository-{safe_name}-untracked"
        for relative, content in snapshot["untracked"].items():
            destination = untracked_root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(content)

    command_record = {
        "argv": command,
        "cwd": str(prepared["cwd"]),
        "environment": prepared["identity"]["environment"],
    }
    (run_dir / "command.json").write_text(
        json.dumps(command_record, indent=2, sort_keys=True) + "\n"
    )
    (run_dir / "identity.json").write_text(
        json.dumps(prepared["identity"], indent=2, sort_keys=True) + "\n"
    )

    stdout_path = run_dir / "stdout.txt"
    stderr_path = run_dir / "stderr.txt"
    time_path = run_dir / "time.txt"
    time_binary = Path("/usr/bin/time")
    if not time_binary.is_file():
        raise HarnessError("GNU time is required at /usr/bin/time")
    timed_command = [str(time_binary), "-v", "-o", str(time_path), "--", *command]

    started_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    start = time.monotonic_ns()
    timeout_seconds = spec.get("timeout_seconds")
    if timeout_seconds is not None and float(timeout_seconds) <= 0:
        raise HarnessError("timeout_seconds must be positive")
    timed_out = False
    interrupted = False
    with stdout_path.open("wb") as stdout_handle, stderr_path.open("wb") as stderr_handle:
        process = subprocess.Popen(
            timed_command,
            cwd=prepared["cwd"],
            env=prepared["execution_env"],
            stdout=stdout_handle,
            stderr=stderr_handle,
            start_new_session=True,
        )
        try:
            exit_code = process.wait(
                timeout=None if timeout_seconds is None else float(timeout_seconds)
            )
        except subprocess.TimeoutExpired:
            timed_out = True
            _terminate_process_group(process, signal.SIGTERM)
            try:
                exit_code = process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                _terminate_process_group(process, signal.SIGKILL)
                exit_code = process.wait()
        except KeyboardInterrupt:
            interrupted = True
            _terminate_process_group(process, signal.SIGINT)
            try:
                exit_code = process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                _terminate_process_group(process, signal.SIGKILL)
                exit_code = process.wait()
    wall_seconds = (time.monotonic_ns() - start) / 1_000_000_000
    finished_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    output_records: list[dict[str, Any]] = []
    missing_outputs: list[str] = []
    invalid_outputs: list[str] = []
    indexed_output_counts: dict[int, dict[str, Any]] = {}
    fastq_jobs: list[tuple[int, dict[str, Any], dict[str, Any]]] = []
    for output_index, declaration in enumerate(outputs):
        path = declaration["path"]
        if path.is_file():
            size = path.stat().st_size
            record = {
                "path": str(path),
                "bytes": size,
                "retained": True,
            }
            output_records.append(record)
            if size < declaration["min_bytes"]:
                invalid_outputs.append(
                    f"{path}: {size} bytes is below minimum {declaration['min_bytes']}"
                )
            if declaration["format"] == "fastq":
                fastq_jobs.append((output_index, declaration, record))
            else:
                record["sha256"] = sha256_file(path)
        else:
            missing_outputs.append(str(path))

    futures: list[tuple[int, dict[str, Any], dict[str, Any], Any]] = []
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=max(1, min(len(fastq_jobs), 4))
    ) as executor:
        for output_index, declaration, record in fastq_jobs:
            future = executor.submit(
                inspect_and_normalize_fastq,
                declaration["path"],
                declaration["mate"],
                declaration["normalize"],
                64 * 1024 * 1024,
                run_dir,
                declaration["numeric_id_max"],
                declaration["numeric_audit_executable"],
            )
            futures.append((output_index, declaration, record, future))
        for output_index, declaration, record, future in futures:
            path = declaration["path"]
            try:
                counts = future.result()
                record["sha256"] = counts.pop("sha256")
                normalized_sha256 = counts.pop("normalized_sha256", None)
                validator = counts.pop("validator", None)
                if validator is not None:
                    record["validator"] = validator
                if normalized_sha256 is not None:
                    record["normalized_sha256"] = normalized_sha256
                    record["normalization"] = (
                        "fastq_multiset_v1"
                        if declaration["normalize"] == "fastq_multiset"
                        else (
                            "fastq_id_multiset_v1"
                            if declaration["normalize"] == "fastq_id_multiset"
                            else "fastq_numeric_accession_set_v1"
                        )
                    )
                    record["mate"] = declaration["mate"]
                    if declaration["numeric_id_max"] is not None:
                        record["numeric_id_max"] = declaration["numeric_id_max"]
                indexed_output_counts[output_index] = {"path": str(path), **counts}
            except (OSError, EOFError, HarnessError) as error:
                # Preserve the raw output identity even when structural
                # validation fails. The extra scan is confined to failed
                # runs; valid publication outputs remain single-pass.
                record["sha256"] = sha256_file(path)
                invalid_outputs.append(str(error))
                indexed_output_counts[output_index] = {
                    "path": str(path),
                    "format": "fastq",
                    "valid": False,
                    "error": str(error),
                }
    output_counts = [
        indexed_output_counts[index] for index in sorted(indexed_output_counts)
    ]

    timing = parse_gnu_time(time_path)
    termination_signal = -exit_code if exit_code < 0 else None
    preliminary_success = (
        exit_code == 0
        and not timed_out
        and not interrupted
        and not missing_outputs
        and not invalid_outputs
    )
    if preliminary_success:
        for declaration, record in zip(outputs, output_records):
            if not declaration["retain"]:
                try:
                    declaration["path"].unlink()
                    record["retained"] = False
                except OSError as error:
                    invalid_outputs.append(
                        f"cannot remove generated output {declaration['path']}: {error}"
                    )
    success = preliminary_success and not invalid_outputs
    result = {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_id,
        "attempt": attempt,
        "success": success,
        "exit_code": exit_code,
        "termination_signal": termination_signal,
        "timed_out": timed_out,
        "interrupted": interrupted,
        "started_utc": started_utc,
        "finished_utc": finished_utc,
        "wall_seconds": wall_seconds,
        **timing,
        "missing_outputs": missing_outputs,
        "invalid_outputs": invalid_outputs,
        "outputs": output_records,
        "output_counts": output_counts,
        "stdout_sha256": sha256_file(stdout_path),
        "stderr_sha256": sha256_file(stderr_path),
        "host": host_snapshot(),
        "identity": prepared["identity"],
    }
    (run_dir / "run.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    (run_dir / "outputs.sha256").write_text(
        "".join(f"{item['sha256']}  {item['path']}\n" for item in output_records)
    )
    (run_dir / "outputs.normalized.sha256").write_text(
        "".join(
            f"{item['normalized_sha256']}  {item['path']}\n"
            for item in output_records
            if "normalized_sha256" in item
        )
    )
    (run_dir / "output-counts.json").write_text(
        json.dumps(output_counts, indent=2, sort_keys=True) + "\n"
    )
    if interrupted:
        raise KeyboardInterrupt
    return result, run_dir


def example_spec() -> dict[str, Any]:
    return {
        "name": "seqproc-fixed-slice-baseline",
        "cwd": "/absolute/path/to/seqproc-paper-analysis",
        "command": [
            "/absolute/path/to/seqproc",
            "--geom",
            "/absolute/path/to/config.geom",
            "--file1",
            "/absolute/path/to/input.fastq",
            "--out1",
            "{run_dir}/output.fastq",
            "--threads",
            "8",
        ],
        "inputs": [
            {
                "path": "/absolute/path/to/input.fastq",
                "sha256": "verified-dataset-digest",
                "verify": False,
            }
        ],
        "configs": ["/absolute/path/to/config.geom"],
        "outputs": ["{run_dir}/output.fastq"],
        "repositories": [
            {"name": "seqproc", "path": "/absolute/path/to/seqproc"},
            {"name": "antisequence", "path": "/absolute/path/to/antisequence"},
        ],
        "environment_mode": "allowlist",
        "environment": {"ANTISEQ_CHUNK_SIZE": "512"},
        "metadata": {
            "dataset": "fixed-slice-synthetic",
            "threads": 8,
            "replicate": 1,
            "candidate": "baseline",
        },
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spec", nargs="?", type=Path, help="JSON run specification")
    parser.add_argument("--output-root", type=Path, default=Path("benchmark_runs"))
    parser.add_argument("--dry-run", action="store_true", help="validate and print run identity")
    parser.add_argument("--write-example", type=Path, help="write an example JSON specification")
    args = parser.parse_args(argv)

    if args.write_example:
        args.write_example.write_text(json.dumps(example_spec(), indent=2) + "\n")
        return 0
    if args.spec is None:
        parser.error("spec is required unless --write-example is used")

    try:
        spec = json.loads(args.spec.read_text())
        if args.dry_run:
            prepared = prepare_spec(spec)
            print(json.dumps({"run_id": prepared["run_id"], "identity": prepared["identity"]}, indent=2))
            return 0
        result, run_dir = execute_spec(spec, args.output_root)
    except (OSError, json.JSONDecodeError, HarnessError) as error:
        print(f"benchmark-harness: {error}", file=sys.stderr)
        return 2

    print(json.dumps({"run_id": result["run_id"], "attempt": result["attempt"], "success": result["success"], "run_dir": str(run_dir)}))
    if result["success"]:
        return 0
    if result["timed_out"]:
        return 124
    if result["termination_signal"]:
        return min(255, 128 + int(result["termination_signal"]))
    exit_code = int(result["exit_code"])
    return exit_code if 1 <= exit_code <= 125 else 1


if __name__ == "__main__":
    raise SystemExit(main())
