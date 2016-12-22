#' Return Slope (normalized slope of the drug response curve) for an experiment of a pSet by taking 
#' its concentration and viability as input.
#' 
#' @examples
#' dose <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8") 
#' viability <- c("108.67","111","102.16","100.27","90","87","74","57")
#' computeSlope(dose, viability)
#' 
#' @param concentration [vector] A concentration range that the AUC should be computed for that range.
#' Concentration range by default considered as not logarithmic scaled. Converted to numeric by function if necessary.
#' @param viability [vector] Viablities corresponding to the concentration range passed as first parameter.
#' The range of viablity values by definition should be between 0 and 100. But the viabalities greater than
#' 100 and lower than 0 are also accepted.
#' @param trunc [binary] A flag that identify if the viabality values should be truncated to be in the
#' range of (0,100)
#' @param verbose [boolean] If 'TRUE' the function will retrun warnings and other infomrative messages.
#' @return Returns the normalized linear slope of the drug response curve
#' @export

computeSlope <- function(concentration, viability, trunc=TRUE, verbose=TRUE) {
  concentration <- as.numeric(concentration[!is.na(concentration)])
  viability <- as.numeric(viability[!is.na(viability)])
  ii <- which(concentration == 0)
  if(length(ii) > 0) {
    concentration <- concentration[-ii]
    viability <- viability[-ii]
  }
  ##convert to nanomolar with the assumption that always concentrations are in micro molar
  concentration <- concentration 
  concentration <- log10(concentration) + 6
  if(trunc) {
    viability <- pmin(viability, 100)
    viability <- pmax(viability, 0)
  }
  
  most.sensitive <- NULL
  for(dose in concentration)
  {
    most.sensitive <- rbind(most.sensitive, cbind(dose,0))
  }
  
  slope.prime <- .optimizeRegression(x = most.sensitive[,1], y = most.sensitive[,2])
  slope <- .optimizeRegression(x = concentration, y = viability)
  slope <- round(slope/abs(slope.prime),digits=2)
  return(-slope) 
}
