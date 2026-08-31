# Matchbox fuzzy-linker boundary-guard diagnostic (2026-08-28)

## Question

Does exact barcode matching make approximate linker matching boundary-safe in
Matchbox 0.3.2, or are exact linkers required for the SPLiT-seq PE and
LR-SPLiT-seq comparisons?

This is a diagnostic configuration study, not a replacement performance
campaign.  All Matchbox runs used the pinned release binary from the journal
campaign (`matchbox 0.3.2`, repository commit
`76afc9941bd2bdc901a3cb5667f8785704c8c40d`, binary SHA-256
`992ebc5b0f447162bc5d999f51d490df3c38b7f0c1fb9b9f8e31aeda4fe6b2ec`),
`--match-mode one-best`, and the final protocol-specific structural references
from `structural_reference_2026-08-17`.

## Configurations

For each technology, four configurations were screened:

1. canonical barcode lists, approximate linkers, projected-component length
   guards;
2. canonical barcode lists, approximate linkers, projected-component and
   linker-span guards;
3. externally Hamming-one-expanded barcode lists, approximate linkers,
   projected-component length guards; and
4. externally Hamming-one-expanded barcode lists, approximate linkers,
   projected-component and linker-span guards.

The PE linker error rate is `0.10` for both 30-nt linkers, giving an edit budget
of three.  LR uses `0.10` for the 30-nt linker and `0.14` for the 22-nt linker,
giving the final edit-distance-three/edit-distance-three comparison.  A linker
span guard requires the matched input interval to have the nominal linker
length.  It admits substitutions and net-zero indel combinations but rejects a
match whose net indel changes the following boundary.

The eight scripts are under `configs/matchbox/` and have names beginning with
`diagnostic_splitseq_pe_*_fuzzy_linkers_` or
`diagnostic_lr_splitseq_*_fuzzy_linkers_`.

## Subset screen

PE used the first 10,000,000 read pairs; LR used the first 1,000,000 reads.  The
four conditions for each list type were run concurrently on disjoint 16-core
CPU sets.  Elapsed times from this screen are therefore useful only for
identifying the extreme expanded-list cost and are not publication-quality
performance measurements.

| Technology | Barcode list | Guard | Emitted | Precision | Recall | F1 | Output geometry |
|---|---|---:|---:|---:|---:|---:|---|
| SPLiT-seq PE | canonical | components | 6,257,913 | 0.997880 | 0.861742 | 0.924828 | all 34 nt |
| SPLiT-seq PE | canonical | components + linker spans | 5,835,387 | 0.999898 | 0.805183 | 0.892039 | all 34 nt |
| SPLiT-seq PE | Hamming-expanded | components | 7,210,925 | 0.995876 | 0.990982 | 0.993423 | all 34 nt |
| SPLiT-seq PE | Hamming-expanded | components + linker spans | 6,692,233 | 0.999889 | 0.923406 | 0.960127 | all 34 nt |
| LR-SPLiT-seq | canonical | components | 15,200 | 0.144408 | 0.022294 | 0.038625 | all 32 nt |
| LR-SPLiT-seq | canonical | components + linker spans | 1,599 | 1.000000 | 0.016241 | 0.031962 | all 32 nt |
| LR-SPLiT-seq | Hamming-expanded | components | 90,122 | 0.578627 | 0.529648 | 0.553055 | all 32 nt |
| LR-SPLiT-seq | Hamming-expanded | components + linker spans | 42,958 | 0.994343 | 0.433849 | 0.604113 | all 32 nt |

Every subset output contained unique read identifiers.  For LR, component-only
guards do not make fuzzy linkers safe.  Linker-span guards restore high
precision, but their sensitivity is worse than the already evaluated
expanded-barcode/exact-linker control (full-data F1 0.840).  No full LR rerun
was therefore warranted.

## Full SPLiT-seq PE validation

The two canonical-list PE configurations were promoted to one full-data,
32-thread materialized-output run.  Both emitted paired FASTQs with byte-
identical accepted-ID bitmaps, no duplicate identifiers, uniform 66-nt R1, and
uniform 34-nt projected R2.

| Configuration | Emitted | Precision | Recall | F1 | Materialized wall time | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|
| Existing exact-linker primary | 35,160,366 | 1.000000 | 0.612150 | 0.759421 | not rerun here | not rerun here |
| Fuzzy linkers + component guards | 50,266,141 | 0.997777 | 0.873200 | 0.931341 | 356.40 s | 11,079.5 MiB |
| Fuzzy linkers + component and linker-span guards | 46,835,210 | 0.999887 | 0.815319 | 0.898220 | 350.39 s | 10,324.6 MiB |

The component-guard and linker-span-guard rows intersect the conservative
reference in 50,154,400 and 46,829,916 pairs, respectively.  The complete
structural-reference and pairwise metrics are in `full_accuracy.json`.

Against the independent full split-pipe vendor accepted-read set, the
component-guard configuration has precision 0.985744, recall 0.830010, F1
0.901198, and Jaccard 0.820165.  The linker-span configuration has precision
0.999740, recall 0.784337, F1 0.879035, and Jaccard 0.784177.  This remains an
agreement comparison, not biological ground truth.

## Single `/dev/null` diagnostic timing

One additional full-data run of the preferred canonical-list component-guard
configuration used 32 pinned physical cores and directed both outputs to
`/dev/null`:

- wall time: 347.92 s;
- user CPU time: 5,059.35 s;
- system CPU time: 14.61 s; and
- peak RSS: 11,345,024 KiB (11,079.1 MiB).

This is a single diagnostic measurement, not a replacement for the randomized
three-replicate publication grid.  Relative to the existing exact-linker
Matchbox 32-thread result (278.445 s and 7,751.8 MiB), approximate linkers plus
component guards cost approximately 25% more wall time and 43% more memory,
while raising structural-reference F1 from 0.759 to 0.931 without external
whitelist expansion.

## Conclusion

Exact barcode identity alone is not a general guarantee of boundary safety in
Matchbox 0.3.2.  The answer is geometry-dependent:

- For anchored, fixed-layout SPLiT-seq PE, explicit projected-component length
  guards are sufficient in practice on the full dataset.  Requiring nominal
  linker spans is unnecessarily conservative.  The canonical-list
  component-guard configuration is the best primary Matchbox candidate tested.
- For unbounded, dual-orientation LR-SPLiT-seq, component guards remain unsafe.
  Linker-span guards restore precision but discard too many valid reads.  The
  existing exact-linker configuration remains the defensible primary LR
  configuration.
- Hamming-expanded lists can produce excellent PE sensitivity but remain an
  expensive, user-visible preprocessing workaround and are not needed to obtain
  the large PE improvement above.

Before replacing the paper's PE Matchbox row, run the preferred configuration
through the frozen randomized performance harness for three replicates at 1,
4, 16, and 32 threads.  Its full-data accuracy result is already sufficient to
replace the current exact-linker PE accuracy row if the authors adopt this more
qualitatively aligned configuration.
