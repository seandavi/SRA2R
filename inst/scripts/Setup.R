library(dplyr)
sortAlignments = function(g) {
  g2 = transform(g, n=nchar(as.character(seqnames)))
  h = select(g2[order(g2["n"], g2["seqnames"], g["start"], -g["end"]),], -c(9))
  rownames(h) <- seq_along(h)
  h
}

install.packages("roxygen2")
### SRAdb
source("https://bioconductor.org/biocLite.R")
biocLite("SRAdb")
biocLite("ShortRead")
biocLite("GenomicRanges")
biocLite("GenomicAlignments")

library(roxygen2)
library(devtools)
library(SRAdb)
library(ShortRead)
library(GenomicRanges)
library(GenomicAlignments)
load_all()
