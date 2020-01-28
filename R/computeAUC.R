#' Computes the AUC for a Drug Dose Viability Curve
#' 
#' Returns the AUC (Area Under the drug response Curve) given concentration and viability as input, normalized by the concentration
#' range of the experiment. The area returned is the response (1-Viablility) area, i.e. area under the curve when the response curve 
#' is plotted on a log10 concentration scale, with high AUC implying high sensitivity to the drug. The function can calculate both 
#' the area under a fitted Hill Curve to the data, and a trapz numeric integral of the actual data provided. Alternatively, the parameters
#' of a Hill Slope returned by logLogisticRegression can be passed in if they already known. 
#' 
#' @examples
#' dose <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8") 
#' viability <- c("108.67","111","102.16","100.27","90","87","74","57")
#' computeAUC(dose, viability)
#' 
#' 
#' @param concentration [vector] is a vector of drug concentrations.
#' @param viability [vector] is a vector whose entries are the viability values observed in the presence of the
#' drug concentrations whose logarithms are in the corresponding entries of conc, where viability 0
#' indicates that all cells died, and viability 1 indicates that the drug had no effect on the cells. 
#' @param Hill_fit [list or vector] In the order: c("Hill Slope", "E_inf", "EC50"), the parameters of a Hill Slope 
#' as returned by logLogisticRegression. If conc_as_log is set then the function assumes logEC50 is passed in, and if
#' viability_as_pct flag is set, it assumes E_inf is passed in as a percent. Otherwise, E_inf is assumed to be a decimal, 
#' and EC50 as a concentration. 
#' @param conc_as_log [logical], if true, assumes that log10-concentration data has been given rather than concentration data.
#' @param viability_as_pct [logical], if false, assumes that viability is given as a decimal rather
#' than a percentage, and returns AUC as a decimal. Otherwise, viability is interpreted as percent, and AUC is returned 0-100.
#' @param trunc [logical], if true, causes viability data to be truncated to lie between 0 and 1 before
#' curve-fitting is performed.
#' @param area.type Should the area be computed using the actual data ("Actual"), or a fitted curve ("Fitted")
#' @param verbose [logical], if true, causes warnings thrown by the function to be printed.
#' @return Numeric AUC value
#' @export
#' @import caTools

computeAUC <- function (concentration,
   viability,
   Hill_fit,
   conc_as_log = FALSE,
   viability_as_pct = TRUE,
   trunc = TRUE,
   area.type = c("Fitted", "Actual"),
   verbose = TRUE
   #, ...
   ) {

  if(missing(concentration)){

    stop("The concentration values to integrate over must always be provided.")

  }
if (missing(area.type)) {
    area.type <- "Fitted"
} else {
    area.type <- match.arg(area.type)
}
if (area.type == "Fitted" && missing(Hill_fit)) {

    Hill_fit <- logLogisticRegression(concentration,
      viability,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose)
    cleanData <- sanitizeInput(conc=concentration, 
      Hill_fit=Hill_fit,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose)
    pars <- cleanData[["Hill_fit"]]
    concentration <- cleanData[["log_conc"]]
} else if (area.type == "Fitted" && !missing(Hill_fit)){

  cleanData <- sanitizeInput(conc = concentration, 
    viability = viability,
    Hill_fit = Hill_fit,
    conc_as_log = conc_as_log,
    viability_as_pct = viability_as_pct,
    trunc = trunc,
    verbose = verbose)
  pars <- cleanData[["Hill_fit"]]
  concentration <- cleanData[["log_conc"]]
} else if (area.type == "Actual" && !missing(viability)){
  cleanData <- sanitizeInput(conc = concentration,
     viability = viability,
     conc_as_log = conc_as_log,
     viability_as_pct = viability_as_pct,
     trunc = trunc,
     verbose = verbose)
  concentration <- cleanData[["log_conc"]]
  viability <- cleanData[["viability"]]
} else if (area.type == "Actual" && missing(viability)){

  stop("To calculate the actual area using a trapezoid integral, the raw viability values are needed!")
}

if (length(concentration) < 2) {
  return(NA)
}

a <- min(concentration)
b <- max(concentration)
if (area.type == "Actual") {
  trapezoid.integral <- caTools::trapz(concentration, viability)
  AUC <- 1 - trapezoid.integral / (b - a)
}
else {
    if(pars[2]==1){
      AUC <- 0
    }else if(pars[1]==0){
      AUC <- (1-pars[2])/2
    } else {
      AUC <- as.numeric((1 - pars[2]) / (pars[1] * (b - a)) * log10((1 + (10 ^ (b - pars[3])) ^ pars[1]) / (1 + (10 ^ (a - pars[3])) ^ pars[1])))
    }
}

if(viability_as_pct){

  AUC <- AUC*100

}

return(AUC)
}