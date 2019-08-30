#' PSet molecularProfiles from ESets to SEs
#'
#' Converts all ExpressionSet objects within the molecularProfiles slot of a 
#'   PharmacoSet to SummarizedExperiments
#'
#' @param \code{S4} A PharmacoSet containing molecular data in ExpressionSets
#'
#' @return \code{S4} A PharmacoSet containing molecular data in a SummarizedExperiments
#' 
#' @importFrom parallel mclapply
#' @importFrom PharmacoGx PharmacoSetClass
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom Biobase exprs fData pData annotation protocolData
#' @improtFrom S4Vectors SimpleList
#' 
#' @export
convertPsetMolecularProfilesToSE <- function(pSet) {
  
  eSets <- pSet@molecularProfiles # Extract eSet data
  
  pSet@molecularProfiles <- 
    mclapply(eSets,
           FUN=function(eSet){
             
             # Change rownames form probes to EnsembleGeneId
             if (Biobase::annotation(eSet) == "rna") {
               rownames(eSet) <- fData(eSet)$EnsemblGeneId
             }
             
             # Build summarized experiment from eSet
             SE <- SummarizedExperiment(
               # Convert to SimpleList class to meet SE requirements
               assays = S4Vectors::SimpleList(as.list(eSet@assayData)),
               # Switch rearrange columns so that IDs are first, probes second
               rowData = S4Vectors::DataFrame(Biobase::fData(eSet),
                                              rownames=rownames(eSet) 
                                              ),
               colData = S4Vectors::DataFrame(Biobase::pData(eSet),
                                              rownames=colnames(eSet)
                                              ),
               metadata = list("experimentData" = eSet@experimentData, 
                               "annotation" = Biobase::annotation(eSet), 
                               "protocolData" = Biobase::protocolData(eSet), 
                               ".__classVersion__" = eSet@.__classVersion__)
               )
               # Assign SE to pSet
               mDataType <- Biobase::annotation(eSet)
               pSet@molecularProfiles[[mDataType]] <- SE
          })
  pSet
}