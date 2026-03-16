# definitions
#[hamming(1)] anchor = f[CAGAGC]
brc1 = norm(b[9-10])
brc2 = b[10]
umi = u[8]

# read structure
1{
  <brc1><anchor><umi><brc2>
}2{r<read>:}
-> 1{<brc1><brc2><umi>}2{<read>}
