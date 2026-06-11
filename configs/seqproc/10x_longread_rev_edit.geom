# 10x Chromium long-read barcode extraction (Reverse Complement) using EDIT DISTANCE
# Structure: [_][barcode_rc:16][primer_rc:22][_]
# Primer: CTACACGACGCTCTTCCGATCT -> RC: AGATCGGAAGAGCGTCGTGTAG

#[search(relative)] #[edit(3)] primer_rc = f[AGATCGGAAGAGCGTCGTGTAG]
bc_rc = b[16]
umi_rc = u[12]
rest_rc = r:

1{<rest_rc><umi_rc><bc_rc><primer_rc>}
-> 1{revcomp(<bc_rc>)revcomp(<umi_rc>)revcomp(<rest_rc>)}
