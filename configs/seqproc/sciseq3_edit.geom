# definitions
#[edit(1)] anchor = f[CAGAGC]
# Preserve the protocol's natural 9- or 10-nt BC1 rather than padding it.
# This makes every tool project the same biological fields at 27/28 nt.
brc1 = b[9-10]
brc2 = b[10]
umi = u[8]

# read structure
1{
  <brc1><anchor><umi><brc2>
}2{r<read>:}
-> 1{<brc1><brc2><umi>}2{<read>}
