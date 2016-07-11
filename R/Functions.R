#'
#'
#' This returns the alignments of the specified accession represented as a Data.Frame or GAlignments object
#'
#' @param acc An accession or a path to an actual SRA file (with .sra suffix)
#' @param seqnames The reference names (Optional)
#' @param start The starting positions (Optional)
#' @param end The ending positions (Optional)
#' @param ranges GRanges object instead of seqnames, start, and end (Optional)
#' @param sort Boolean value representing whether or not alignments should be sorted
#' @param dataframe Boolean value representing whether return should be of type Data.Frame or not (GAlignments)
#' @return the aligned reads with GAlignments data
#' @export
#' @examples
#' getGAlignments("ERR1189688")
#' getGAlignments("ERR1189688", seqnames = 'X', start = 66943491, end = 6694392)
#' getGAlignments("ERR1189688", ranges = GRanges(seqnames = 'X', ranges = IRanges(start = 6694391, end = 6694392)))
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

#'
#'
#' Helper function to call cpp_getGAlignments and sort results, if applicable
#' 
#' @param range A range of locations to retrieve alignments from
#' @param acc An accession
#' @param sort Boolean value representing whether or not alignments should be sorted
#' @param dataframe Boolean value representing whether return should be of type Data.Frame or not (GAlignments)
#' @return GAlignments or DataFrame object representing reads that overlap with given range
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


#'
#'
#' Function to search SRA
#' 
#' @param search_terms String of terms to search
#' @param num Integer value of how many accession IDs to return
#' @param public Boolean value representing whether to return public or private access reads
#' @return a character vector of accessions that match given search terms
#' @export
searchSRA = function(search_terms, num = 20, public = TRUE) {
  library(rentrez)
  if (public) {
    id <- entrez_search(db="sra", term=paste(search_terms, " AND cluster_public[prop]", sep = ""), retmax = num)
  } else {
    id <- entrez_search(db="sra", term=paste(search_terms, " AND cluster_dbgap[prop]", sep = ""), retmax = num)
  }
  acclist <- entrez_fetch("sra", id=id$ids, rettype = "acclist")
  srr <- strsplit(acclist, '\n')[[1]]
  srr <- srr[grepl( "^<Acc", srr, perl=TRUE)]
  srr <- gsub( '<Acc>|</Acc>', '', srr, perl = T)
  return(srr)
}
