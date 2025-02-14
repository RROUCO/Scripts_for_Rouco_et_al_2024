library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm9)

options(scipen = 999)

# load genome
genome <- getBSgenome("BSgenome.Mmusculus.UCSC.mm9", masked=FALSE)

# Source: http://rebase.neb.com/rebase/rebase.html
#5' --- GATC--- 3'
#3' ---CTAG --- 5'

# Search works only with motif whose reverse complement is exactly the same GATC (reverse complement = GATC)
pattern.DpnII <- "GATC"
seqInfo <- seqinfo(Mmusculus)

chroms <- seqlevels(seqInfo)

outfile <- "DpnII_starts_mm9.txt"
  
cat("" , file = outfile, sep = " ", fill = FALSE, labels = NULL, append = FALSE)

for(chrom in chroms) {
  print(chrom)
  res.DpnII <- matchPattern(pattern.DpnII, genome[[chrom]], max.mismatch=0)
  # start position of cutting sites
  pos <- start(res.DpnII)
  seqlen <- seqlengths(seqInfo)[chrom]
    
  cat(chrom, pos, seqlen, "\n" , file = outfile, sep = " ", fill = FALSE, labels = NULL, append = TRUE)
}
