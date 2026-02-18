import sys

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def main():
    config_path = "/home/ubuntu/combine-lab/seqproc-paper-analysis-clean/configs/splitcode/splitseq_paper.config"
    output_path = "/home/ubuntu/combine-lab/seqproc-paper-analysis-clean/configs/splitcode/splitseq_singleend_rc.config"

    bc1_seqs = []
    bc2_seqs = []
    bc3_seqs = []
    linker1_seq = ""
    linker2_seq = ""

    with open(config_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('@') or line.startswith('ID'):
                continue
            
            parts = line.split('\t')
            name = parts[0]
            seq = parts[1]
            
            if name == 'linker1':
                linker1_seq = seq
            elif name == 'linker2':
                linker2_seq = seq
            elif name.startswith('bc1_'):
                bc1_seqs.append((name, seq))
            elif name.startswith('bc2_'):
                bc2_seqs.append((name, seq))
            elif name.startswith('bc3_'):
                bc3_seqs.append((name, seq))

    # RC sequences
    l1_rc = reverse_complement(linker1_seq)
    l2_rc = reverse_complement(linker2_seq)

    with open(output_path, 'w') as out:
        out.write("# SPLiT-seq single-end RC config\n")
        out.write("# RC Structure: bc1_rc -> linker2_rc -> bc2_rc -> linker1_rc -> bc3_rc -> umi_rc\n")
        out.write("ID\tTAG\tGROUP\tDISTANCE\tNEXT\n")

        # Linkers
        # Linker2 RC leads to bc2_rc
        out.write(f"linker2_rc\t{l2_rc}\tlinker2_rc\t3\t{{{{bc2_rc}}}}\n")
        # Linker1 RC leads to bc3_rc
        out.write(f"linker1_rc\t{l1_rc}\tlinker1_rc\t3\t{{{{bc3_rc}}}}\n")

        # BC1 RC -> Linker2 RC
        for name, seq in bc1_seqs:
            rc_seq = reverse_complement(seq)
            out.write(f"{name}_rc\t{rc_seq}\tbc1_rc\t2\t{{linker2_rc}}\n")

        # BC2 RC -> Linker1 RC
        for name, seq in bc2_seqs:
            rc_seq = reverse_complement(seq)
            out.write(f"{name}_rc\t{rc_seq}\tbc2_rc\t2\t{{linker1_rc}}\n")

        # BC3 RC -> End (UMI extraction handles the rest)
        for name, seq in bc3_seqs:
            rc_seq = reverse_complement(seq)
            out.write(f"{name}_rc\t{rc_seq}\tbc3_rc\t2\t-\n")

        out.write("\n")
        # Extract line: Explicitly define the chain
        # Note: We rely on finding bc1_rc first.
        out.write("@extract {{bc1_rc}}{linker2_rc}{{bc2_rc}}{linker1_rc}{{bc3_rc}}<umi_rc[10]>\n")

if __name__ == "__main__":
    main()
