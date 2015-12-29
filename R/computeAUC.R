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
#' @param area.type Should the area be computed using the actual data ("Actual"), or a fitted curve ("Fitted")
#' @param verbose [boolean] If 'TRUE' the function will retrun warnings and other infomrative messages.
#' @return Numeric AUC value
#' @export
#' @import caTools

computeAUC <- function(concentration, viability, trunc=TRUE, area.type=c("Fitted","Actual"), verbose=TRUE) {
  concentration <- as.numeric(concentration[!is.na(concentration)])
  viability <- as.numeric(viability[!is.na(viability)])
  ii <- which(concentration == 0)
  if(length(ii) > 0) {
    concentration <- concentration[-ii]
    viability <- viability[-ii]
  }
  
  if(missing(area.type)){
    area.type <- "Fitted"
  }
  if(length(concentration) < 2){
    return(NA)
  }
  if(area.type == "Actual"){
    if(trunc) {viability = pmin(as.numeric(viability), 100); viability = pmax(as.numeric(viability), 0)}
    trapezoid.integral <- caTools::trapz(log10(as.numeric(concentration)) ,as.numeric(viability))
    AUC <- round(1- (trapezoid.integral/caTools::trapz(log10(as.numeric(concentration)), rep(100, length(viability)))), digits=2)
  }else{
    AUC <- .computeAUCUnderFittedCurve(concentration, viability, trunc)
  }
  return (AUC)
}
  
