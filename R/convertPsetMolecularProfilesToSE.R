#' PSet molecularProfiles from ESets to SEs
#'
#' Converts all ExpressionSet objects within the molecularProfiles slot of a 
#'   PharmacoSet to SummarizedExperiments
#'
#' @param pSet \code{S4} A PharmacoSet containing molecular data in ExpressionSets
#'
#' @return \code{S4} A PharmacoSet containing molecular data in a SummarizedExperiments
#' 
#' @importFrom parallel mclapply
#' @importFrom SummarizedExperiment SummarizedExperiment Assays
#' @importFrom Biobase exprs fData pData annotation protocolData
#' @importFrom S4Vectors SimpleList DataFrame
#' 
#' @export
convertPsetMolecularProfilesToSE <- function(pSet) {
  
  eSets <- pSet@molecularProfiles # Extract eSet data
  
  pSet@molecularProfiles <- 
    #parallel::mc
    lapply(eSets,
           FUN=function(eSet){
             
             # Change rownames from probes to EnsemblGeneId for rna data type
             if (grepl("^rna$", Biobase::annotation(eSet))) {
               rownames(eSet) <- Biobase::fData(eSet)$EnsemblGeneId
             }
             
             # Build summarized experiment from eSet
             SE <- SummarizedExperiment::SummarizedExperiment(
               # Convert to SimpleList class to meet SE requirements
               assays = S4Vectors::SimpleList(eSet@assayData$exprs
                                              ),
               # Switch rearrange columns so that IDs are first, probes second
               rowData = S4Vectors::DataFrame(Biobase::fData(eSet),
                                              rownames=rownames(Biobase::fData(eSet)) 
                                              ),
               colData = S4Vectors::DataFrame(Biobase::pData(eSet),
                                              rownames=rownames(Biobase::pData(eSet))
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
  setNames(pSet@molecularProfiles, names(eSets))
  pSet
}