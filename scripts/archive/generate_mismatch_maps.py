#!/usr/bin/env python3
"""
Generate expanded barcode TSV files with all mismatch variants.
For 8bp barcodes with ≤2 mismatches, each barcode generates ~277 entries:
- 1 exact match
- 24 single-mismatch variants (8 positions × 3 alternatives)
- 252 double-mismatch variants (C(8,2) × 9)

This allows using exact hash lookup instead of hamming distance matching.
"""

from pathlib import Path
from itertools import combinations

NUCLEOTIDES = ['A', 'C', 'G', 'T']

def generate_variants(barcode: str, max_mismatches: int = 2) -> set:
    """Generate all variants of a barcode within max_mismatches hamming distance."""
    variants = {barcode}  # Include exact match
    
    # Single mismatches
    if max_mismatches >= 1:
        for i in range(len(barcode)):
            for nuc in NUCLEOTIDES:
                if nuc != barcode[i]:
                    variant = barcode[:i] + nuc + barcode[i+1:]
                    variants.add(variant)
    
    # Double mismatches
    if max_mismatches >= 2:
        for i, j in combinations(range(len(barcode)), 2):
            for nuc1 in NUCLEOTIDES:
                if nuc1 != barcode[i]:
                    for nuc2 in NUCLEOTIDES:
                        if nuc2 != barcode[j]:
                            variant = barcode[:i] + nuc1 + barcode[i+1:j] + nuc2 + barcode[j+1:]
                            variants.add(variant)
    
    return variants


def load_original_tsv(path: Path) -> dict:
    """Load original TSV file: returns {sequence: id}"""
    barcodes = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    bc_id, sequence = parts[0], parts[1]
                    barcodes[sequence] = bc_id
    return barcodes


def generate_expanded_tsv(input_path: Path, output_path: Path, max_mismatches: int = 2):
    """Generate expanded TSV with all mismatch variants."""
    barcodes = load_original_tsv(input_path)
    
    # Generate all variants
    expanded = {}
    collisions = 0
    
    for sequence, bc_id in barcodes.items():
        variants = generate_variants(sequence, max_mismatches)
        for variant in variants:
            if variant in expanded:
                # Collision - mark as ambiguous
                if expanded[variant] != bc_id:
                    collisions += 1
                    expanded[variant] = "AMBIGUOUS"
            else:
                expanded[variant] = bc_id
    
    # Write expanded TSV (id\tsequence format for seqproc)
    with open(output_path, 'w') as f:
        for sequence, bc_id in sorted(expanded.items()):
            if bc_id != "AMBIGUOUS":
                f.write(f"{bc_id}\t{sequence}\n")
    
    return len(barcodes), len(expanded), collisions


def main():
    project_root = Path(__file__).parent.parent
    configs_dir = project_root / 'configs/seqproc'
    
    print("Generating expanded barcode maps with mismatch variants...")
    print("=" * 60)
    
    for bc_name in ['bc1', 'bc2', 'bc3']:
        input_path = configs_dir / f'splitseq_{bc_name}_map.tsv'
        output_path = configs_dir / f'splitseq_{bc_name}_expanded.tsv'
        
        if not input_path.exists():
            print(f"  {bc_name}: SKIP (input not found)")
            continue
        
        orig_count, expanded_count, collisions = generate_expanded_tsv(
            input_path, output_path, max_mismatches=2
        )
        
        print(f"  {bc_name}: {orig_count} barcodes → {expanded_count} variants ({collisions} collisions)")
        print(f"         Saved: {output_path}")
    
    print("=" * 60)
    print("Done! Use these expanded files with exact 'map' (no mismatch tolerance needed)")


if __name__ == '__main__':
    main()
