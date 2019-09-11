##' Validate PSet molecularProfiles Conversion
##' 
##' Checks that all the information contained in an ExpressionSet molecularProfile 
##'   was successfully tranferred to the SummarizedExperiment molecularProfile
##'   
##' @param [S4] a PSet containing molecularProfiles as SummarizedExperiments
##' @param [S4] a PSet containing molecularProfiles as ExpressionSets
##' 
##' @return [string] Any slots which are not the same
##' 
##' @importFrom 
##' 
##' @export
##' 
#validatePsetMolecularProfilesToSEConversion <- function()