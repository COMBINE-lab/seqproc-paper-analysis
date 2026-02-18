import gzip
import sys

def fix_fastq(input_file, output_file, read_num):
    print(f"Fixing {input_file} -> {output_file} (Read {read_num})...")
    with gzip.open(input_file, 'rt') as f_in, gzip.open(output_file, 'wt') as f_out:
        for i, line in enumerate(f_in):
            if i % 4 == 0:
                # Header line
                parts = line.strip().split()
                # Ensure read number is correct in the comment
                # @SRR... <read_num> length=...
                # Keep ID, change comment to match read_num
                if len(parts) > 1:
                    # check if parts[1] is a number
                    if parts[1].isdigit():
                        parts[1] = str(read_num)
                    # reconstruct
                    new_header = " ".join(parts)
                    f_out.write(new_header + "\n")
                else:
                    f_out.write(line)
            else:
                f_out.write(line)
    print("Done.")

if __name__ == "__main__":
    fix_fastq("SRR6750041_1M_R1.fastq.gz", "SRR6750041_1M_R1.fixed.fastq.gz", 1)
    fix_fastq("SRR6750041_1M_R2.fastq.gz", "SRR6750041_1M_R2.fixed.fastq.gz", 2)
