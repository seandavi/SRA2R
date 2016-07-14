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