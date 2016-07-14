#' Get fastq equivalent from SRR record
#' 
#' This method will fetch reads from SRA
#' for the given SRR record. By default, all
#' reads are returned. 
#' 
#' @param object an SRR class object
#' 
#' @return a ShortReadQ object
#' @exportMethod getReads
setMethod('getReads',
          signature(object='SRR'),
          function(object,...) {
            tmpvals = getFastqReadsWithQuality(object)
            return(ShortReadQ(
              sread=DNAStringSet(tmpvals$reads),
              quality=BStringSet(tmpvals$qualities)
            ))
          }
)