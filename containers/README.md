# split-pipe 1.4.0 reproduction environment

`split-pipe-1.4.0.podman.Dockerfile` is a rootless-Podman-compatible adaptation
of ParseNIP commit `f1ba1291aee95b4c6fa85bf0ce6678ffd078f6a1`, file
`dockerfiles/split-pipe.Dockerfile`. It downloads the same split-pipe 1.4.0
archive and Miniconda installer, verifies both SHA-256 digests, installs all
runtime dependencies in the Python 3.10 `spipe` environment, and records conda
and pip locks inside the image.

The split-pipe archive is not redistributed. Its bundled Parse Biosciences
Software License Agreement governs use; confirm that the intended use and host
comply before building or running the image.

## Build and verify

```bash
podman build \
  -f containers/split-pipe-1.4.0.podman.Dockerfile \
  -t seqproc-split-pipe:pbp-1.4.0 .

podman run --rm --entrypoint split-pipe \
  seqproc-split-pipe:pbp-1.4.0 --version
```

The calibration image built on 2026-08-19 had image ID
`a5eb92b407db654d072d020dde7616c4cf832face99015504ed080cb3ac474a8`
and digest
`sha256:bb55241a68c6d4fecf776b921c3fbee2540134ed0e8445574889dd0f32987dd8`.
Rebuilds should be recorded by immutable digest rather than assumed to have this
digest.

## Preprocessing-only reference

split-pipe 1.4.0 validates a split-pipe-formatted genome reference during
startup even when the only requested stage is `pre`. Build the supplied tiny
reference with split-pipe's own `mkref` command:

```bash
podman run --rm \
  -v "$PWD/containers/split-pipe-reference:/reference-source:ro" \
  -v "$PWD/work:/work" \
  seqproc-split-pipe:pbp-1.4.0 \
  --mode mkref \
  --parfile /reference-source/mkref.par \
  --gfasta dummy /reference-source/dummy-gene.fa \
  --output_dir /work/reference \
  --nthreads 1 --no_keep_going
```

The reference is only a startup requirement: `--mode pre --one_step` never
aligns reads, so its sequence cannot affect barcode acceptance.

## SRR6750041 vendor-set run

Mount gzipped R1 and R2 inputs as `/input`, the generated reference as
`/work/reference`, and the repository configuration as `/config`:

```bash
podman run --rm \
  -v "/path/to/gzipped-inputs:/input:ro" \
  -v "$PWD/work:/work" \
  -v "$PWD/configs/split-pipe:/config:ro" \
  seqproc-split-pipe:pbp-1.4.0 \
  --mode pre --one_step \
  --parfile /config/splitseq_pe_v1.par \
  --chemistry v1 --kit custom \
  --bc_round_set 1 v1 --bc_round_set 2 v1 --bc_round_set 3 v1 \
  --sample_bc_rounds 1 \
  --fq1 /input/SRR6750041_1.fastq.gz \
  --fq2 /input/SRR6750041_2.fastq.gz \
  --output_dir /work/output \
  --genome_dir /work/reference \
  --nthreads 32 --no_keep_going
```

The custom layout is required: the bundled chemistry-v1 template has a 22-nt
second linker, whereas SRR6750041 has the original 30-nt linker. The ENA headers
also use `N/1` and `N/2` in the second token; the configuration retains paired
identifier checking but disables split-pipe's incompatible assumption that the
first character of that token is the mate number.

Before the full run, this command was calibrated on the first 10,000,000 pairs.
It reproduced all 7,539,920 archived split-pipe IDs exactly: zero IDs occurred
in either side of the symmetric difference, and the decompressed sorted-ID
SHA-256 was
`b363cb3b804ca302851f2efdb0324fae062109291e961094b2b0cf8e00b66872`.

Use `scripts/splitpipe_full_concordance.py` to validate the resulting FASTQ,
construct a compact accepted-ID bitmap, compare it with all three final tool
bitmaps, and write JSON/CSV provenance artifacts.
