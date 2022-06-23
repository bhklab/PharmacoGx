#' @describeIn computeICn Returns the IC50 of a Drug Dose response curve
#'
#' @examples
#' dose <- c(0.0025,0.008,0.025,0.08,0.25,0.8,2.53,8)
#' viability <- c(108.67,111,102.16,100.27,90,87,74,57)
#' computeICn(dose, viability)
#'
#' @return `numeric(1)` The ICn of the Hill curve over the specified dose
#' range.
#'
#' @export
computeIC50 <- function(concentration,
                       viability,
                       Hill_fit,
                       conc_as_log = FALSE,
                       viability_as_pct = TRUE,
                       verbose = TRUE,
                       trunc = TRUE) {

  return(computeICn(concentration = concentration,
                    viability = viability,
                    Hill_fit = Hill_fit,
                    n = ifelse(viability_as_pct, 50, .5),
                    conc_as_log = conc_as_log,
                    viability_as_pct = viability_as_pct,
                    verbose=TRUE,
                    trunc=TRUE))
}
