#[search(relative)] #[edit(4)] primer = f[CTACACGACGCTCTTCCGATCT]
bc = b[16]
umi = u[12]
rest = r:
1{<primer><bc><umi><rest>}
-> 1{<bc><umi><rest>}
