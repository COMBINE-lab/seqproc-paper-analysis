# Reproducible benchmark specifications

`scripts/benchmark_harness.py` executes one content-addressed run specification.
`scripts/run_frozen_schedule.py` is the publication-matrix coordinator: it
validates a frozen YAML/JSON manifest, deterministically randomizes conditions
within recorded blocks, protects the schedule with a SHA-256 sidecar, preserves
every failed attempt, and resumes without repeating successful runs.

The minimal scheduler manifest has this shape:

```yaml
schema_version: "1.0.0"
mode: development                 # or publication
study:
  name: example
  random_seed: 741211
artifacts:
  - path: /absolute/path/to/input.fastq
    sha256: 64_lowercase_hex_characters
execution:
  timeout_seconds: 3600
  sanitized_environment_allowlist: [PATH, LD_LIBRARY_PATH, LANG, TMPDIR]
runs:
  - id: dataset-t8-r1-seqproc
    block: dataset-t8-r1
    spec:
      name: seqproc
      cwd: /absolute/path/to/analysis
      command: [/absolute/path/to/seqproc, run, --geom, /path/config.geom,
                --file1, /path/input.fastq, --out1, "{run_dir}/out.fastq",
                --threads, "8"]
      inputs:
        - {path: /path/input.fastq, sha256: 64_lowercase_hex_characters, verify: false}
      configs:
        - {path: /path/config.geom, sha256: 64_lowercase_hex_characters, verify: false}
      outputs:
        - {path: "{run_dir}/out.fastq", format: fastq, min_bytes: 1}
      repositories:
        - {name: seqproc, path: /path/seqproc, commit: 40_lowercase_hex_characters}
      metadata: {dataset: example, tool: seqproc, threads: 8, replicate: 1}
```

All `artifacts` are hashed once before schedule generation or execution. Run
inputs/configs may therefore use `verify: false` to avoid hashing multi-GB files
for every replicate, but publication mode still requires their declared
SHA-256 values. Publication mode also requires full repository commits, clean
trees, and forbids `--max-runs`.

Generate and validate the immutable schedule:

```bash
python scripts/run_frozen_schedule.py \
  --manifest benchmark-manifest.yaml \
  --schedule schedule.json \
  --generate

python scripts/run_frozen_schedule.py \
  --manifest benchmark-manifest.yaml \
  --schedule schedule.json \
  --validate-only
```

Execute or resume it:

```bash
python scripts/run_frozen_schedule.py \
  --manifest benchmark-manifest.yaml \
  --schedule schedule.json \
  --output /benchmark/results
```

For a local seqproc-only integration tier, generate the manifest with
`scripts/make_seqproc_small_tier_manifest.py`. This is a development gate for
the harness and current binary; it is not a substitute for the full randomized
cross-tool publication matrix.

The analysis environment is hash-locked in `requirements.lock`. A clean check
can be run without installing into the repository:

```bash
uv run --with-requirements requirements.lock python -m pytest -q
```
