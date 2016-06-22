install.packages("roxygen2", repos="http://cran.rstudio.com/")
### SRAdb
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite("GenomicAlignments")

library(roxygen2)
library(devtools)
library(GenomicRanges)
library(GenomicAlignments)
load_all()
library(dplyr)
sortAlignments = function(g) {
  g2 = transform(g, n=nchar(as.character(seqnames)))
  # h = select(g2[order(g2["n"], g2["seqnames"], g["start"], -g["end"]),], -c(9))
  rownames(h) <- seq_along(h[,1])
  g = GAlignments(seqnames=Rle(factor(h$seqnames)), pos=h$start, cigar=as.character(h$cigar), strand=Rle(factor(as.character(h$strand), levels=c('+', '-', '*'))))
  g
}


