# Cross-tool semantic audit (2026-08-05)

This audit precedes the full timed matrix. It distinguishes exact semantic
equivalence from the best practical approximation supported by each tool.
Timing runs use `--match-mode one-best` for matchbox and retain output read-ID
digests so determinism can be checked within each tool independently of
cross-tool concordance.

## Findings that affect the benchmark specification

1. Matchbox's default `all` mode can emit the same input read more than once.
   On the 100,000-read protocol-shaped SPLiT-seq PE input, the legacy-style
   script emitted 230,224 paired output records. Explicit `one-best` mode
   emitted 77,709. Publication runs therefore specify `one-best`; output counts
   are records, not raw action invocations.
2. The legacy LR splitcode configuration contained only 20 whitelist entries
   per barcode round, versus 96 in the paired-end configuration. The
   publication LR comparison instead uses two linker tags and performs no
   barcode-whitelist predicate, matching the forward-orientation seqproc and
   matchbox task.
3. A splitcode linker distance of `3` is a practical substitution-distance
   index. On the adversarial LR input it accepts 60,000/100,000 reads, missing
   insertion/deletion cases. Setting `3:3:3` to permit substitutions and indels
   caused index construction to remain incomplete beyond the 30-second
   diagnostic window even for a 1,000-read cap. It is therefore not labeled
   semantically identical to the edit-distance implementations.
4. Seqproc and matchbox accepted the same 80,000 forward-oriented reads in the
   100,000-read LR adversarial input. Their normalized read-ID multiset digest
   was `4f47481127fcb22830818ba557fe1f4150b78b6144b82d6df9c57db5754f7830`.
   Splitcode accepted 60,000 and produced digest
   `f9f3208725209cee46552a23b30ff8938bf6c44a505686658ea2239ca4fab360`.
5. An approximate matchbox pattern ending exactly at the read boundary reaches
   `not yet implemented` at `src/core/matcher/read_matcher.rs:992` in matchbox
   0.3.2. The real sci-RNA-seq3 pattern safely retains its trailing wildcard.
   On the checked 1M real subset it emitted 897,666 paired reads; seqproc's
   prior validated count on the same subset was 892,555.
6. The checked seqproc BC2/BC3 whitelist contains the same 96 entries twice.
   Seqproc normalizes duplicate whitelist entries internally. Publication
   matchbox and splitcode configurations explicitly use one copy of each entry
   so duplicate input rows cannot create artificial ambiguity or work.

## Protocol-shaped audit inputs

The four inputs were generated deterministically by
`scripts/protocol_integration_benchmark.py` and are recorded in
`benchmark_specs/integration/current-small-tier-2026-08-05.yaml`. They include
exact, substitution, insertion, deletion, and (for LR-SPLiT-seq) reverse-
orientation cases. They are semantic/adversarial fixtures, not biological
truth sets and not performance inputs.

## Publication interpretation

The 10x v2 fixed-position task is directly equivalent. Other cross-tool rows
are labeled **best practical supported configuration**. Exact output identity
is required across threads and replicates within a tool/mode. Cross-tool
read-set equality is reported as a result, not assumed as a prerequisite for a
valid timing observation.
