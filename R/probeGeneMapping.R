########################
## Benjamin Haibe-Kains, Zhaleh Safikhani
## All rights Reserved
## July 16, 2015
########################

#' Identify the best probes to use for each RNA gene, and reduce the data to
#' gene level
#' 
#' This function allows the summarization of RNA microarray data by the probe
#' most suitable for each gene. This function uses either the jetset package
#' to select the best probe to keep for each ENTREZ gene id, and
#' returns the dataset at a gene level. The different platforms are identified
#' using their GEO id.
#' 
#' @references 
#'  Li Q, Birkbak NJ, Gyorffy B, Szallasi Z and Eklund AC (2011). "Jetset:
#'  selecting the optimal microarray probe set to represent a gene." BMC
#'  Bioinformatics, 12, pp. 474.
#' 
#' @examples 
#' data(CGPsmall)
#' CGPsmall <- probeGeneMapping(CGPsmall)
#'
#' @param pSet [PharmacoSet] The PharmacoSet on which to preform the mapping platform
#' @export

probeGeneMapping <- function (pSet){

  pSet@molecularData$rna <- pSet@molecularData$rna[geneInfo(pSet)[,"BEST"] == TRUE, ]
  rownames(geneInfo(pSet)) <- paste('geneid', geneInfo(pSet)[,'GENEID'], sep='.')
  return (pSet)
}

