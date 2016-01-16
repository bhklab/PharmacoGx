########################
## Mark Freeman
## All rights Reserved
## August 18, 2015
## Function to calculate area between dose-response curves over their common concentration range
########################

#'  Fits dose-response curves to data given by the user
#'  and returns the ABC of the fitted curves.
#'  
#'  @param conc1 [vector] is a vector of drug concentrations.
#'  
#'  @param conc2 [vector] is a vector of drug concentrations.
#'  
#'  @param viability1 [vector] is a vector whose entries are the viability values observed in the presence of the
#'  drug concentrations whose logarithms are in the corresponding entries of conc1, expressed as percentages
#'  of viability in the absence of any drug.
#'  
#'  @param viability2 [vector] is a vector whose entries are the viability values observed in the presence of the
#'  drug concentrations whose logarithms are in the corresponding entries of conc2, expressed as percentages
#'  of viability in the absence of any drug.
#'  
#'  @param trunc [logical], if true, causes viability data to be truncated to lie between 0 and 1 before
#'  curve-fitting is performed.

computeABC <- function(conc1, conc2, viability1, viability2, trunc = TRUE) {
  
  conc1 <- na.omit(as.numeric(conc1))
  conc2 <- na.omit(as.numeric(conc2))
  viability1 <- na.omit(as.numeric(viability1))
  viability2 <- na.omit(as.numeric(viability2))
  
  #CHECK THAT FUNCTION INPUTS ARE APPROPRIATE
  if (prod(is.finite(conc1)) != 1) {
    print(conc1)
    stop("Concentration vector conc1 contains elements which are not real numbers.")
  }
  
  if (prod(is.finite(conc2)) != 1) {
    print(conc2)
    stop("Concentration vector conc2 contains elements which are not real numbers.")
  }
  
  if (prod(is.finite(viability1)) != 1) {
    print(viability1)
    stop("Viability vector viability1 contains elements which are not real numbers.")
  }
  
  if (prod(is.finite(viability2)) != 1) {
    print(viability2)
    stop("Viability vector viability2 contains elements which are not real numbers.")
  }
  
  if (is.logical(trunc) == FALSE) {
    print(trunc)
    stop("'trunc' is not a logical.")
  }
  
  if (length(conc1) != length(viability1)) {
    print(conc1)
    print(viability1)
    stop("Concentration vector conc1 is not of same length as viability vector viability1.")
  }
  
  if (length(conc2) != length(viability2)) {
    print(conc2)
    print(viability2)
    stop("Concentration vector conc2 is not of same length as viability vector viability2.")
  }
  
  if (min(conc1) < 0) {
    stop("Concentration vector conc1 contains negative data.")
  }
  
  if (min(conc2) < 0) {
    stop("Concentration vector conc2 contains negative data.")
  }
  
  if (min(viability1) < 0) {
    warning("Warning: Negative viability data in viability vector viability1.")
  }
  
  if (min(viability2) < 0) {
    warning("Warning: Negative viability data in viability vector viability2.")
  }
  
  if (max(viability1) > 100) {
    warning("Warning: Viability data in vector viability1 exceeds negative control.")
  }
  
  if (max(viability2) > 100) {
    warning("Warning: Viability data in vector viability2 exceeds negative control.")
  }
  
  
  #CONVERT DOSE-RESPONSE DATA TO APPROPRIATE INTERNAL REPRESENTATION
  log_conc1 <- log10(conc1)
  log_conc2 <- log10(conc2)
  viability1 <- viability1 / 100
  viability2 <- viability2 / 100
  
  if (trunc == TRUE) {
    viability1[which(viability1 < 0)] <- 0
    viability2[which(viability2 < 0)] <- 0
    viability1[which(viability1 > 1)] <- 1
    viability2[which(viability2 > 1)] <- 1
  }
  
  #FIT CURVE AND CALCULATE IC50
  if (max(log_conc1) < min(log_conc2) | max(log_conc2) < min(log_conc1)) {
    return(NA)
  } else {
    pars1 <- unlist(logLogisticRegression(log_conc1,
                                         viability1,
                                         conc_as_log = TRUE,
                                         viability_as_pct = FALSE,
                                         trunc = trunc))
    pars2 <- unlist(logLogisticRegression(log_conc2,
                                         viability2,
                                         conc_as_log = TRUE,
                                         viability_as_pct = FALSE,
                                         trunc = trunc))
    extrema <- sort(c(min(log_conc1), max(log_conc1), min(log_conc2), max(log_conc2)))
    support <- .GetSupportVec(c(extrema[2], extrema[3]))
    return(caTools::trapz(support, abs(.Hill(support, pars1) - .Hill(support, pars2))) / (extrema[3] - extrema[2]))
  }
}