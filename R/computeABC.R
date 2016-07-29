########################
## Mark Freeman
## All rights Reserved
## August 18, 2015
## Function to calculate area between dose-response curves over their common concentration range
########################

#' Fits dose-response curves to data given by the user
#' and returns the ABC of the fitted curves.
#' 
#' @examples
#' dose1 <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8") 
#' viability1 <- c("108.67","111","102.16","100.27","90","87","74","57")
#' dose2 <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8") 
#' viability2 <- c("100.94","112.5","86","104.16","75","68","48","29")
#' computeABC(dose1, dose2, viability1, viability2)
#' 
#' @param conc1 [vector] is a vector of drug concentrations.
#' @param conc2 [vector] is a vector of drug concentrations.
#' @param viability1 [vector] is a vector whose entries are the viability values observed in the presence of the
#' drug concentrations whose logarithms are in the corresponding entries of conc1, expressed as percentages
#' of viability in the absence of any drug.
#' @param viability2 [vector] is a vector whose entries are the viability values observed in the presence of the
#' drug concentrations whose logarithms are in the corresponding entries of conc2, expressed as percentages
#' of viability in the absence of any drug.
#' @param Hill_fit1 [list or vector] In the order: c("Hill Slope", "E_inf", "EC50"), the parameters of a Hill Slope 
#' as returned by logLogisticRegression. If conc_as_log is set then the function assumes logEC50 is passed in, and if
#' viability_as_pct flag is set, it assumes E_inf is passed in as a percent. Otherwise, E_inf is assumed to be a decimal, 
#' and EC50 as a concentration. 
#' @param Hill_fit2 [list or vector] In the order: c("Hill Slope", "E_inf", "EC50"), the parameters of a Hill Slope 
#' as returned by logLogisticRegression. If conc_as_log is set then the function assumes logEC50 is passed in, and if
#' viability_as_pct flag is set, it assumes E_inf is passed in as a percent. Otherwise, E_inf is assumed to be a decimal, 
#' and EC50 as a concentration. 
#' @param conc_as_log [logical], if true, assumes that log10-concentration data has been given rather than concentration data.
#' @param viability_as_pct [logical], if false, assumes that viability is given as a decimal rather
#' than a percentage, and returns ABC as a decimal. Otherwise, viability is interpreted as percent, and AUC is returned 0-100.
#' @param verbose [logical], if true, causes warnings thrown by the function to be printed.
#' @param trunc [logical], if true, causes viability data to be truncated to lie between 0 and 1 before
#' curve-fitting is performed.
#' @return The numeric area of the absolute difference between the two hill slopes
#' @export 
computeABC <- function(conc1, conc2, viability1, viability2,
                        Hill_fit1,
                        Hill_fit2,
                        conc_as_log = FALSE,
                        viability_as_pct = TRUE,
                        trunc = TRUE,
                        verbose=TRUE) {

if(missing(conc1) | missing(conc2)){

    stop("Both Concentration vectors the drugs were tested on must always be provided.")

}
if (missing(Hill_fit1) | missing(Hill_fit2)) {

    Hill_fit1 <- logLogisticRegression(conc1,
      viability1,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose)
    cleanData <- sanitizeInput(conc=conc1, 
      Hill_fit=Hill_fit1,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose)
    pars1 <- cleanData[["Hill_fit"]]
    log_conc1 <- cleanData[["log_conc"]]
    Hill_fit2 <- logLogisticRegression(conc2,
      viability2,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose)
    cleanData <- sanitizeInput(conc=conc2, 
      Hill_fit=Hill_fit2,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose)
    pars2 <- cleanData[["Hill_fit"]]
    log_conc2 <- cleanData[["log_conc"]]

} else {

  cleanData <- sanitizeInput(conc = conc1, 
    viability = viability1,
    Hill_fit = Hill_fit1,
    conc_as_log = conc_as_log,
    viability_as_pct = viability_as_pct,
    trunc = trunc,
    verbose = verbose)
  pars1 <- cleanData[["Hill_fit"]]
  log_conc1 <- cleanData[["log_conc"]]
  cleanData <- sanitizeInput(conc = conc2, 
    viability = viability2,
    Hill_fit = Hill_fit2,
    conc_as_log = conc_as_log,
    viability_as_pct = viability_as_pct,
    trunc = trunc,
    verbose = verbose)
  pars2 <- cleanData[["Hill_fit"]]
  log_conc2 <- cleanData[["log_conc"]]
}
  
  #FIT CURVE AND CALCULATE IC50
  if (max(log_conc1) < min(log_conc2) | max(log_conc2) < min(log_conc1)) {
    return(NA)
  } else {
    extrema <- sort(c(min(log_conc1), max(log_conc1), min(log_conc2), max(log_conc2)))
    support <- .GetSupportVec(c(extrema[2], extrema[3]))
    ABC <- as.numeric(caTools::trapz(support, abs(.Hill(support, pars1) - .Hill(support, pars2))) / (extrema[3] - extrema[2]))
    if(viability_as_pct){
      ABC <- ABC*100
    }
    return(ABC)
  }
}