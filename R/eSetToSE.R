#` Eset to SE
#`
#' Converts all Expression Set objects within a PharmacoSet to SummarizedExperiments
#'
#' @param {S4} A PharmacoSet containing molecular data in an ExpressionSet
#'
#' @return {S4} A PharmacoSet containing molecular data in a SummarizedExperiment
#' 
#' @export
#' @importFrom parallel mclapply
#' @importFrom PharmacoGx PharmacoSetClass
#' 
eSetToSE <- function(pSet) {
  
  eSets <- pSet@molecularProfiles # Extract eSet data
  
  pSet@molecularProfiles <- 
    mclapply(eSets,
           FUN=function(eSet){
             SE <- SummarizedExperiment( # Build summarized experiment
               assays = Biobase::exprs(eSet),
               rowData = Biobase::fData(eSet),
               colData = Biobase::pData(eSet),
               metadata = list(eSet@experimentData, Biobase::annotation(eSet), 
                               Biobase::protocolData(eSet), eSet@.__classVersion__)
              )
              mDataType <- Biobase::annotation(eSet)
              pSet@molecularProfiles[[mDataType]] <- SE # Assign SE to pSet
              }
           )
  return(pSet)
}
