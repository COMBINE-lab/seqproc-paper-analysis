import sys

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def process_fastq(input_path, output_path):
    with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
        while True:
            header = fin.readline()
            if not header: break
            seq = fin.readline().strip()
            plus = fin.readline()
            qual = fin.readline().strip()
            
            fout.write(header)
            fout.write(reverse_complement(seq) + '\n')
            fout.write(plus)
            fout.write(qual[::-1] + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rc_fastq.py <input.fq> <output.fq>")
        sys.exit(1)
    
    process_fastq(sys.argv[1], sys.argv[2])
