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

#' Search SRA using Eutils
#' 
#' @param term The text used for searching
#' @return A list of SRR ids
#' 
#' @export
#' @importFrom xml2 read_xml xml_children xml_text
#' @importFrom rentrez entrez_search entrez_fetch
#' @import magrittr
#' 
#' @examples 
#' ids = searchSRA('breast cancer[WORD] AND cluster_public[prop]')
#' head(ids)
searchSRA = function(term) {
  res = entrez_search(db="sra", term=term, use_history=TRUE)
  ids = entrez_fetch('sra', web_history=res$web_history, rettype='acclist') %>% 
    read_xml() %>% 
    xml_children() %>%
    xml_text()
  return(ids)
  }