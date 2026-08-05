# Publication benchmark prerequisites

The publication manifest uses an exact, bounded-memory auditor for uncompressed
FASTQ outputs whose read names are numeric ENA accessions. It validates every
FASTQ record, rejects duplicate/out-of-range/mixed-accession IDs, and writes a
canonical bitset whose SHA-256 is independent of output ordering. Raw FASTQ
SHA-256 is computed concurrently and retained separately.

Build the dependency-free auditor before freezing a manifest:

```bash
python3 scripts/build_fastq_numeric_audit.py
```

The build receipt records the compiler version, exact argv, source checksum,
and binary checksum. The manifest independently pins both source and binary.
The generic Python external-sort normalizers remain available for non-ENA IDs.
