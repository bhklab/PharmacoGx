#' Compute the linear slope of a drug response curve
#'
#' @examples
#' dose <- c(0.0025,0.008,0.025,0.08,0.25,0.8,2.53,8)
#' viability <- c(108.67,111,102.16,100.27,90,87,74,57)
#' computeSlope(dose, viability)
#'
#' @param concentration `numeric` vector of drug concentrations. Concentrations
#'   are assumed to be provided on the raw scale and are transformed to
#'   `log10(concentration)` internally before fitting.
#' @param viability `numeric` vector of percent viability values aligned to
#'   `concentration`. Values outside the range `[0, 100]` are accepted and can
#'   be truncated with `trunc`.
#' @param trunc `logical(1)` If `TRUE`, clip viability values to lie in
#'   `[0, 100]` before fitting the regression.
#' @param verbose `logical(1)` If `TRUE`, emit warnings and informative messages
#'   produced while sanitizing the inputs.
#' @return Returns the negative ordinary least-squares slope of normalized
#'   viability (0-1 scale) regressed on `log10(concentration)`, rounded to two
#'   decimals.
#'
#' @export
computeSlope <- function(
  concentration,
  viability,
  trunc = TRUE,
  verbose = TRUE
) {
  clean_data <- sanitizeInput(
    conc = concentration,
    viability = viability,
    conc_as_log = FALSE,
    viability_as_pct = TRUE,
    trunc = trunc,
    verbose = verbose
  )

  log_conc <- clean_data[["log_conc"]]
  viability_clean <- clean_data[["viability"]]

  if (length(log_conc) < 2L || length(unique(log_conc)) < 2L) {
    return(NA_real_)
  }

  slope <- stats::coef(stats::lm(viability_clean ~ log_conc))[["log_conc"]]
  if (!is.finite(slope)) {
    return(NA_real_)
  }

  round(-unname(slope), digits = 2)
}
