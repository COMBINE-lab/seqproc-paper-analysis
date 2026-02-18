#!/usr/bin/env python3
"""
End-to-end zebrafish pineal dataset analysis with seqproc.

This script:
1. Downloads the zebrafish pineal dataset (10x Chromium v2)
2. Runs seqproc to extract barcodes/UMIs
3. Outputs transformed FASTQs ready for salmon alevin

Dataset: SRR8315379, SRR8315380 (zebrafish pineal gland, 10x v2)
Reference: Alevin-fry paper (Nature Methods 2022)

Usage:
    python scripts/run_zebrafish_analysis.py [options]

Options:
    --data-dir      Directory to store downloaded data
    --output-dir    Directory for seqproc output
    --subset        Number of reads to subset (0 = all)
    --threads       Number of threads for seqproc
    --skip-download Skip downloading if files exist
"""

import argparse
import os
import subprocess
import sys
import time
import json
from pathlib import Path

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Run end-to-end zebrafish analysis with seqproc',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--data-dir', default='data/zebrafish',
                        help='Directory to store downloaded data')
    parser.add_argument('--output-dir', default='results/zebrafish_analysis',
                        help='Directory for analysis output')
    parser.add_argument('--subset', type=int, default=100000,
                        help='Number of reads to subset (0 = all, default 100K for testing)')
    parser.add_argument('--threads', type=int, default=4,
                        help='Number of threads for seqproc (multi-threading bug fixed)')
    parser.add_argument('--skip-download', action='store_true',
                        help='Skip downloading if files exist')
    parser.add_argument('--seqproc-bin', default=None,
                        help='Path to seqproc binary')
    return parser.parse_args()


def download_zebrafish_data(data_dir, skip_if_exists=False):
    """Download zebrafish pineal dataset from EBI.
    
    Returns paths to R1 and R2 files.
    """
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Dataset info from alevin-fry paper
    samples = [
        ('SRR8315379', 'ftp://ftp.ebi.ac.uk/vol1/fastq/SRR831/009/SRR8315379/SRR8315379_1.fastq.gz',
                      'ftp://ftp.ebi.ac.uk/vol1/fastq/SRR831/009/SRR8315379/SRR8315379_2.fastq.gz'),
        ('SRR8315380', 'ftp://ftp.ebi.ac.uk/vol1/fastq/SRR831/000/SRR8315380/SRR8315380_1.fastq.gz',
                      'ftp://ftp.ebi.ac.uk/vol1/fastq/SRR831/000/SRR8315380/SRR8315380_2.fastq.gz'),
    ]
    
    downloaded_files = {'r1': [], 'r2': []}
    
    for srr, url_r1, url_r2 in samples:
        r1_file = data_dir / f'{srr}_R1.fastq.gz'
        r2_file = data_dir / f'{srr}_R2.fastq.gz'
        
        for url, dest in [(url_r1, r1_file), (url_r2, r2_file)]:
            if skip_if_exists and dest.exists():
                print(f"  Skipping (exists): {dest.name}")
            else:
                print(f"  Downloading: {dest.name}...")
                try:
                    subprocess.run(['curl', '-L', '-o', str(dest), url], 
                                   check=True, capture_output=True)
                except subprocess.CalledProcessError as e:
                    print(f"    Warning: curl failed, trying wget...")
                    subprocess.run(['wget', '-O', str(dest), url], 
                                   check=True, capture_output=True)
        
        downloaded_files['r1'].append(str(r1_file))
        downloaded_files['r2'].append(str(r2_file))
    
    return downloaded_files


def create_10x_v2_geom(output_dir):
    """Create seqproc geometry file for 10x Chromium v2.
    
    10x v2 structure:
    - R1: 16bp barcode + 10bp UMI
    - R2: cDNA
    
    seqproc syntax:
    - b[N] = barcode of N bases
    - u[N] = UMI of N bases
    - r: = rest of read (biological)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    geom_content = """# 10x Chromium v2 geometry for seqproc
# R1: 16bp cell barcode + 10bp UMI
# R2: cDNA (biological sequence)

bc = b[16]
umi = u[10]
bio = r:

# Read 1 contains barcode and UMI
1{<bc><umi>}

# Read 2 contains biological sequence
2{<bio>}

# Output: barcode, UMI, and biological read
-> 1{<bc><umi>} 2{<bio>}
"""
    geom_file = output_dir / '10x_v2.geom'
    with open(geom_file, 'w') as f:
        f.write(geom_content)
    
    return str(geom_file)


def subset_fastq(input_file, output_file, num_reads):
    """Subset a FASTQ file to first N reads."""
    import gzip
    
    input_path = Path(input_file)
    output_path = Path(output_file)
    
    opener = gzip.open if str(input_file).endswith('.gz') else open
    
    with opener(input_file, 'rt') as fin, open(output_file, 'w') as fout:
        read_count = 0
        line_count = 0
        for line in fin:
            fout.write(line)
            line_count += 1
            if line_count % 4 == 0:
                read_count += 1
                if read_count >= num_reads:
                    break
    
    return read_count


def run_seqproc(r1_files, r2_files, geom_file, output_dir, threads, seqproc_bin=None):
    """Run seqproc on the zebrafish data.
    
    seqproc takes single input files, so we process each sample separately
    and concatenate results.
    
    Returns dict with timing and output info.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine seqproc binary
    if seqproc_bin:
        bin_path = seqproc_bin
    else:
        bin_path = os.environ.get('SEQPROC_BIN', '../seqproc/target/release/seqproc')
    
    bin_path = str(Path(bin_path).resolve())
    if not Path(bin_path).exists():
        raise FileNotFoundError(f"seqproc binary not found: {bin_path}")
    
    geom_path = str(Path(geom_file).resolve())
    
    total_runtime = 0
    total_reads = 0
    temp_outputs = []
    
    # Process each sample separately
    for i, (r1, r2) in enumerate(zip(r1_files, r2_files)):
        sample_out_r1 = output_dir / f'sample{i+1}_processed_R1.fq'
        sample_out_r2 = output_dir / f'sample{i+1}_processed_R2.fq'
        
        cmd = [
            bin_path,
            '--geom', geom_path,
            '--file1', str(Path(r1).resolve()),
            '--file2', str(Path(r2).resolve()),
            '--out1', str(sample_out_r1.resolve()),
            '--out2', str(sample_out_r2.resolve()),
            '--threads', str(threads)
        ]
        
        print(f"  Sample {i+1}: {Path(r1).name}")
        
        start_time = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True)
        end_time = time.time()
        
        total_runtime += (end_time - start_time)
        
        if result.returncode != 0:
            print(f"    Warning: seqproc exited with code {result.returncode}")
            if result.stderr:
                print(f"    Stderr: {result.stderr[:200]}")
        else:
            # Count reads
            if sample_out_r1.exists():
                with open(sample_out_r1, 'r') as f:
                    count = sum(1 for line in f if line.startswith('@'))
                total_reads += count
                print(f"    Processed: {count:,} reads in {end_time - start_time:.2f}s")
        
        temp_outputs.append((sample_out_r1, sample_out_r2))
    
    # Concatenate all outputs
    out_r1 = output_dir / 'zebrafish_processed_R1.fq'
    out_r2 = output_dir / 'zebrafish_processed_R2.fq'
    
    with open(out_r1, 'w') as fout:
        for tmp_r1, _ in temp_outputs:
            if tmp_r1.exists():
                with open(tmp_r1, 'r') as fin:
                    fout.write(fin.read())
    
    with open(out_r2, 'w') as fout:
        for _, tmp_r2 in temp_outputs:
            if tmp_r2.exists():
                with open(tmp_r2, 'r') as fin:
                    fout.write(fin.read())
    
    return {
        'runtime': total_runtime,
        'reads_out': total_reads,
        'output_r1': str(out_r1),
        'output_r2': str(out_r2),
        'returncode': 0 if total_reads > 0 else 1,
        'stderr': ''
    }


def main():
    args = parse_args()
    
    print("=" * 70)
    print("Zebrafish Pineal Dataset Analysis with seqproc")
    print("=" * 70)
    print(f"Data directory:   {args.data_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Subset size:      {args.subset if args.subset > 0 else 'ALL'}")
    print(f"Threads:          {args.threads}")
    print("=" * 70)
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Download data
    print("\n[1/4] Downloading zebrafish data...")
    downloaded = download_zebrafish_data(args.data_dir, skip_if_exists=args.skip_download)
    print(f"  Downloaded {len(downloaded['r1'])} sample pairs")
    
    # Step 2: Subset if requested
    if args.subset > 0:
        print(f"\n[2/4] Subsetting to {args.subset:,} reads...")
        subset_dir = output_dir / 'subset'
        subset_dir.mkdir(parents=True, exist_ok=True)
        
        r1_subset = []
        r2_subset = []
        for i, (r1, r2) in enumerate(zip(downloaded['r1'], downloaded['r2'])):
            r1_out = subset_dir / f'sample{i+1}_R1.fq'
            r2_out = subset_dir / f'sample{i+1}_R2.fq'
            
            count = subset_fastq(r1, r1_out, args.subset)
            subset_fastq(r2, r2_out, args.subset)
            print(f"  Sample {i+1}: {count:,} reads")
            
            r1_subset.append(str(r1_out))
            r2_subset.append(str(r2_out))
        
        input_r1 = r1_subset
        input_r2 = r2_subset
    else:
        print("\n[2/4] Using full dataset (no subset)")
        input_r1 = downloaded['r1']
        input_r2 = downloaded['r2']
    
    # Step 3: Create geometry file
    print("\n[3/4] Creating 10x v2 geometry file...")
    geom_file = create_10x_v2_geom(output_dir)
    print(f"  Created: {geom_file}")
    
    # Step 4: Run seqproc
    print("\n[4/4] Running seqproc...")
    try:
        results = run_seqproc(
            input_r1, input_r2, geom_file, output_dir,
            args.threads, args.seqproc_bin
        )
        
        print(f"\n  Runtime:    {results['runtime']:.2f} seconds")
        print(f"  Reads out:  {results['reads_out']:,}")
        print(f"  Output R1:  {results['output_r1']}")
        print(f"  Output R2:  {results['output_r2']}")
        
        if results['returncode'] != 0:
            print(f"\n  Warning: seqproc exited with code {results['returncode']}")
            if results['stderr']:
                print(f"  Stderr: {results['stderr'][:500]}")
        
        # Save results
        results_file = output_dir / 'analysis_results.json'
        with open(results_file, 'w') as f:
            json.dump({
                'dataset': 'zebrafish_pineal',
                'samples': ['SRR8315379', 'SRR8315380'],
                'subset': args.subset,
                'threads': args.threads,
                'runtime_seconds': results['runtime'],
                'reads_processed': results['reads_out'],
                'output_files': {
                    'r1': results['output_r1'],
                    'r2': results['output_r2']
                }
            }, f, indent=2)
        print(f"\n  Results saved to: {results_file}")
        
    except FileNotFoundError as e:
        print(f"\n  Error: {e}")
        print("  Please set SEQPROC_BIN environment variable or use --seqproc-bin")
        sys.exit(1)
    
    print("\n" + "=" * 70)
    print("Analysis Complete!")
    print("=" * 70)
    print("\nNext steps:")
    print("  1. The processed FASTQs can be used with salmon alevin")
    print("  2. Run: salmon alevin -l ISR -i <index> -1 <R1> -2 <R2> --chromiumV2 ...")
    print("  3. Then: alevin-fry generate-permit-list, collate, quant")
    print("=" * 70)


if __name__ == "__main__":
    main()
