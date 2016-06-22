getGAlignments = function(acc, seqnames = '*', start = 1, end = 0, ranges = NULL, sort = TRUE, dataframe = FALSE) {
  if (is.null(ranges)) {
    ranges = GRanges(seqnames = seqnames, IRanges(start = start, end = end))
  }
  if (class(ranges) == class(GRanges())){
    if (dataframe) {
      lapply(ranges, callGAlignments, acc=acc, sort=sort, dataframe=dataframe)
    } else {
      GAlignmentsList(lapply(ranges, callGAlignments, acc=acc, sort=sort, dataframe=dataframe))
    }
  } else if (class(ranges) == class(GRangesList())) {
    lapply(ranges, getGAlignments, acc=acc, sort=sort, dataframe=dataframe)
  }
}

callGAlignments = function(range, acc, sort, dataframe) {
  g = cpp_getGAlignments(acc, seqname = as.character(seqnames(range)), low_bound=start(range), up_bound = end(range))
  g = transform(g, n=nchar(as.character(seqnames)))
  if (sort) {
    g = select(g[order(g["n"], g["seqnames"], g["start"], -g["end"]),], -c(9))
    rownames(g) <- seq_along(g[,1])
  }
  if (dataframe) {
    g
  } else {
    GAlignments(seqnames=Rle(factor(g$seqnames)), pos=g$start, cigar=as.character(g$cigar), strand=Rle(factor(as.character(g$strand), levels=c('+', '-', '*'))))
  }
}