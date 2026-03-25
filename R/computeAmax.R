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
#' of viability in the absence of any drug when `viability_as_pct = TRUE`, or
#' as decimals when `viability_as_pct = FALSE`.
#'
#' @param viability_as_pct `logical(1)` whether the viability values are given
#' as percentages or decimals.
#' @param trunc `logical`, if true, causes viability data to be truncated to lie between 0 and 1 before
#' curve-fitting is performed.
#' @param verbose `logical` should warnings be printed
#' @return The numerical Amax expressed as a percentage.
#' @export
computeAmax <- function(
  concentration,
  viability,
  viability_as_pct = TRUE,
  trunc = TRUE,
  verbose = FALSE
) {
  if (
    !is.logical(viability_as_pct) ||
      length(viability_as_pct) != 1L ||
      is.na(viability_as_pct)
  ) {
    stop("'viability_as_pct' must be a logical value.")
  }

  if (!is.logical(trunc) || length(trunc) != 1L || is.na(trunc)) {
    stop("'trunc' must be a logical value.")
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be a logical value.")
  }

  clean_data <- sanitizeInput(
    conc = concentration,
    viability = viability,
    conc_as_log = FALSE,
    viability_as_pct = viability_as_pct,
    trunc = trunc,
    verbose = verbose
  )

  log_conc <- clean_data[["log_conc"]]
  viability_clean <- clean_data[["viability"]]

  if (length(log_conc) < 2L) {
    if (verbose) {
      warning("Insufficient non-zero concentrations for curve fitting")
    }
    x <- NA_real_
    names(x) <- "Amax"
    return(x)
  }

  #FIT CURVE AND CALCULATE IC50
  pars <- unlist(logLogisticRegression(
    log_conc,
    viability_clean,
    conc_as_log = TRUE,
    viability_as_pct = FALSE,
    trunc = FALSE,
    verbose = FALSE
  ))
  internal <- .normalizeHillPars(
    hill_fit = pars,
    conc_as_log = TRUE,
    viability_as_pct = FALSE
  )
  x <- 100 -
    .pgx_hill_curve(
      max(log_conc),
      unname(internal[c("HS", "E0", "E_inf", "log10EC50")])
    ) *
      100
  names(x) <- "Amax"
  return(x)
}
