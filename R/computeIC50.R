#'Fits dose-response curves to data given by the user
#'and returns the IC50 of the fitted curve.
#'
#'@param conc [vector] is a vector of drug concentrations.
#'@param viability [vector] is a vector whose entries are the viability values observed in the presence of the
#'drug concentrations whose logarithms are in the corresponding entries of the log_conc, expressed as percentages
#'of viability in the absence of any drug.
#'@param trunc [logical], if true, causes viability data to be truncated to lie between 0 and 1 before
#'curve-fitting is performed.
#'@param verbose [logical] should diagnostic messages be printed? (default=FALSE)
#'@return An estimate of the IC50 for the concentrations and viabilities provided
#'@export

computeIC50 <- function(conc, viability, trunc = TRUE, verbose=FALSE) {
  
  conc <- as.numeric(conc[!is.na(conc)])
  viability <- as.numeric(viability[!is.na(viability)])
  ii <- which(conc == 0)
  if(length(ii) > 0) {
    conc <- conc[-ii]
    viability <- viability[-ii]
  }
  
  #CHECK THAT FUNCTION INPUTS ARE APPROPRIATE
  if (prod(is.finite(conc)) != 1) {
    print(conc)
    stop("Concentration vector contains elements which are not real numbers.")
  }
  
  if (prod(is.finite(viability)) != 1) {
    print(viability)
    stop("Viability vector contains elements which are not real numbers.")
  }
  
  if (is.logical(trunc) == FALSE) {
    print(trunc)
    stop("'trunc' is not a logical.")
  }
  
  if (length(conc) != length(viability)) {
    print(conc)
    print(viability)
    stop("Concentration vector is not of same length as viability vector.")
  }
  
  if (min(conc) < 0) {
    stop("Concentration vector contains negative data.")
  }
  
  if (min(viability) < 0 & verbose) {
    warning("Warning: Negative viability data.")
  }
  
  if (max(viability) > 100 & verbose) {
    warning("Warning: Viability data exceeds negative control.")
  }
  
  #CONVERT DOSE-RESPONSE DATA TO APPROPRIATE INTERNAL REPRESENTATION
  log_conc <- log10(conc)
  viability <- viability / 100
  
  if (trunc == TRUE) {
    viability[which(viability < 0)] <- 0
    viability[which(viability > 1)] <- 1
  }
  
  #FIT CURVE AND CALCULATE IC50
  pars <- unlist(logLogisticRegression(log_conc,
                                       viability,
                                       conc_as_log = TRUE,
                                       viability_as_pct = FALSE,
                                       trunc = trunc))
  return(as.numeric(ifelse(pars[2] < 1 / 2, 10 ^ pars[3] * (1 - 2 * pars[2]) ^ (-1 / pars[1]), Inf)))
}
