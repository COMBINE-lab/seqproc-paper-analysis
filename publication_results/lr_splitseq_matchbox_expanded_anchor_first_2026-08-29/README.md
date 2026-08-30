# Matchbox LR-SPLiT-seq expanded-list anchor-first replacement

This archive replaces the obsolete LR-SPLiT-seq Matchbox parenthetical that
used the August 11 barcode-first expanded-list diagnostic.  The replacement
uses the same exact linkers, externally Hamming-distance-one-expanded barcode
lists, dual-orientation search, and `one-best` policy, but transcribes them with
the linker-first execution plan used by the final primary Matchbox
configuration.

The frozen Matchbox 0.3.2 binary has SHA-256
`992ebc5b0f447162bc5d999f51d490df3c38b7f0c1fb9b9f8e31aeda4fe6b2ec`.
The configuration is
`configs/matchbox/sensitivity_lr_splitseq_ham1_expanded_anchor_first.mb`
(SHA-256
`afacbbaa8eddd9fdefde25b8f4fcdb4811df772ccb43badfc453c870a2395df6`).

## Replacement values

The single frozen 32-thread `/dev/null` condition completed in
512.660708628 seconds and used 355,444 KiB (347.113 MiB) peak RSS.  It replaces
the former 4,308-second, 434.1-MiB materialized-output diagnostic.  The old and
new numbers do not form a strictly controlled speedup experiment because the
old condition materialized output, but output was only about 40 MiB and the
8.4-fold wall-time difference is dominated by the changed Matchbox search
plan.

The separate materialized 32-thread run completed in 525.14 seconds, used
358,388 KiB peak RSS, and emitted 406,132 records, all exactly 32 nt.  All
406,132 IDs occur in the 560,699-read conservative structural reference:

| Metric | Value |
|---|---:|
| Emitted input fraction | 7.045495% |
| Precision | 100.000000% |
| Recall | 72.433159% |
| F1 | 0.840130 |

## Old/new plan adjudication

The former barcode-first product emitted 406,257 reads.  The two products
share 406,101 read IDs (Jaccard 0.999539735): 156 occur only in the old product
and 31 only in the new product.  Every old-only and new-only read occurs in the
structural reference.  For all shared IDs, transformed sequence and quality
are identical.

On the 187 discordant reads, Matchbox `all` mode with the new plan emitted 174
distinct IDs: all 31 new-only IDs and 143 of the 156 old-only IDs.  Thus, most
old-only exclusions arise because Matchbox selects an equally scoring
anchor-first placement before the post-match barcode guard is evaluated.  The
remaining 13 old-only and 31 new-only calls reflect pattern-boundary and
candidate-selection differences between the two declarative transcriptions.
The publication retains deterministic `one-best`, consistent with the primary
LR configuration; `all` is not adopted because it can emit duplicate records
on multi-cassette reads in the complete dataset.

`old_new_reference_metrics.json` contains the machine-readable accuracy and
pairwise counts.  The frozen manifests, schedules, per-run reports, GNU time
records, output-length audit, and complete discordant-ID lists are stored in
this directory.  Raw FASTQ outputs remain in the campaign directory:

`/scratch1/seqproc-ecosystem/campaigns/journal-matchbox-lr-expanded-anchor-first-2026-08-29`
