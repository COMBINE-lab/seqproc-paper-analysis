primer_rc = anchor_relative(hamming(f[AGATCGGAAGAGCGTCGTGTAG], 3))
bc_rc = b[16]
umi_rc = u[12]
rest_rc = r:
1{<rest_rc><umi_rc><bc_rc><primer_rc>}
-> 1{revcomp(<bc_rc>)revcomp(<umi_rc>)revcomp(<rest_rc>)}
