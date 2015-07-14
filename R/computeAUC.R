########################
## Zhaleh Safikhani
## All rights Reserved
## July 8, 2015
## Function to compute AUC for given concentration and viability vectors
########################


#' Return AUC (Area Under the drug response curve) for an experiment of a pSet by taking 
#' its concentration and viability as input.
#' 
#' @param concentration [vector] A concentration range that the AUC should be computed for that range.
#' Concentration range by default considered as not logarithmic scaled.
#' @param viability [vector] Viablities correspondant to the concentration range passed as first parameter.
#' The range of viablity values by definition should be between 0 and 100. But the viabalities greater than
#' 100 and lower than 0 are also accepted.
#' @param trunc [binary] A flag that identify if the viabality values should be truncated to be in the
#' range of (0,100)
#' @export
#' @import caTools

computeAUC <-
  function(concentration, viability, trunc = TRUE) {
    
    if(length(concentration) < 2){
      return(NA)
    }
    
    if(trunc) {viability = pmin(as.numeric(viability), 100); viability = pmax(as.numeric(viability), 0)}
    trapezoid.integral <- caTools::trapz(log10(as.numeric(concentration)) ,as.numeric(viability))
    AUC <- round(1- (trapezoid.integral/trapz(log10(as.numeric(concentration)), rep(100, length(viability)))), digits=2)
    return (AUC)
  }

