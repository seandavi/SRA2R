##################################################################
#
# Contains all class definitions
#
# Follows similar pattern to ShortRead package
#
##################################################################

.sra2rValidity = function(object) TRUE

# VIRTUAL base class
setClass('.SRA2RBase')

.srrRegex = '^[DSE]RR\\d{6}$'

.srrValidity = function(object) {
  valid = grepl(.srrRegex,object,perl=TRUE)
  if(!valid) {
    message(paste('SRR should be a single string and match the regex',.srrRegex))
  }
  return(valid)
}

#' An S4 class representing an SRR accession
#' 
#' This class simply extends "character" and
#' checks that the accession matches a regex
#' for SRR, ERR, or DRR accessions.
#' 
setClass('SRR',
  contains=c('.SRA2RBase','character'),
  validity = .srrValidity
)

setGeneric('SRRGAlignments',function(object,regions,...) {
  standardGeneric('SRRGAlignments')
})

setGeneric('getReads',function(object,...) {
  standardGeneric('getReads')
})

