#' Fits dose-response curves to data given by the user
#' and returns the Amax of the fitted curve.
#' Amax: 100 - viability at maximum concentarion (in fitted curve)
#'
#' @examples
#' dose <- c(0.0025,0.008,0.025,0.08,0.25,0.8,2.53,8)
#' viability <- c(108.67,111,102.16,100.27,90,87,74,57)
#' computeAmax(dose, viability)
#'
#' @param concentration `numeric` is a vector of drug concentrations.
#'
#' @param viability `numeric` is a vector whose entries are the viability values observed in the presence of the
#' drug concentrations whose logarithms are in the corresponding entries of the log_conc, expressed as percentages
#' of viability in the absence of any drug.
#'
#' @param trunc `logical`, if true, causes viability data to be truncated to lie between 0 and 1 before
#' curve-fitting is performed.
#' @param verbose `logical` should warnings be printed
#' @return The numerical Amax
#' @export
computeAmax <- function(
  concentration,
  viability,
  trunc = TRUE,
  verbose = FALSE
) {
  concentration <- as.numeric(concentration[!is.na(concentration)])
  viability <- as.numeric(viability[!is.na(viability)])
  ii <- which(concentration == 0)
  if (length(ii) > 0) {
    concentration <- concentration[-ii]
    viability <- viability[-ii]
  }
  if (length(concentration) < 2) {
    if (verbose) {
      warning("Insufficient non-zero concentrations for curve fitting")
    }
    x <- NA_real_
    names(x) <- "Amax"
    return(x)
  }

  #CHECK THAT FUNCTION INPUTS ARE APPROPRIATE
  if (!all(is.finite(concentration))) {
    stop(
      "Concentration vector contains non-finite values: ",
      toString(concentration[!is.finite(concentration)])
    )
  }

  if (!all(is.finite(viability))) {
    stop(
      "Viability vector contains non-finite values: ",
      toString(viability[!is.finite(viability)])
    )
  }

  if (!is.logical(trunc)) {
    stop("'trunc' must be a logical value.")
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be a logical value.")
  }

  if (length(concentration) != length(viability)) {
    stop(
      "Concentration vector is not the same length as the viability vector.",
      " Lengths: ",
      length(concentration),
      " vs ",
      length(viability)
    )
  }

  if (min(concentration) < 0) {
    stop("Concentration vector contains negative data.")
  }

  if (min(viability) < 0 && verbose) {
    warning("Negative viability data detected.")
  }

  if (max(viability) > 100 && verbose) {
    warning("Viability data exceeds negative control.")
  }

  #CONVERT DOSE-RESPONSE DATA TO APPROPRIATE INTERNAL REPRESENTATION
  log_conc <- log10(concentration)
  viability <- viability / 100

  if (trunc == TRUE) {
    viability[which(viability < 0)] <- 0
    viability[which(viability > 1)] <- 1
  }

  #FIT CURVE AND CALCULATE IC50
  pars <- unlist(logLogisticRegression(
    log_conc,
    viability,
    conc_as_log = TRUE,
    viability_as_pct = FALSE,
    trunc = trunc
  ))
  x <- 100 - .Hill(max(log_conc), pars) * 100
  names(x) <- "Amax"
  return(x)
}
