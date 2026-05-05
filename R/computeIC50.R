#' @describeIn computeICn Returns the relative IC50 of a drug dose-response
#'   curve. Use `computeAC50()` for an absolute 50% viability threshold.
#'
#' @return `numeric(1)` The IC50 of the Hill curve over the specified dose
#'   range.
#'
#' @export
computeIC50 <- function(
  concentration,
  viability,
  Hill_fit,
  reference = c("relative", "absolute"),
  conc_as_log = FALSE,
  viability_as_pct = TRUE,
  verbose = TRUE,
  trunc = TRUE
) {
  reference <- match.arg(reference)

  computeICn(
    concentration = concentration,
    viability = viability,
    Hill_fit = Hill_fit,
    n = ifelse(viability_as_pct, 50, .5),
    reference = reference,
    conc_as_log = conc_as_log,
    viability_as_pct = viability_as_pct,
    verbose = verbose,
    trunc = trunc
  )
}

#' @describeIn computeICn Returns the AC50 of a drug dose-response curve using
#'   an absolute 50% viability threshold. This is the appropriate helper for
#'   comparisons against published AC50 values.
#' @return `numeric(1)` The AC50 of the Hill curve over the specified dose
#'   range.
#' @export
computeAC50 <- function(
  concentration,
  viability,
  Hill_fit,
  conc_as_log = FALSE,
  viability_as_pct = TRUE,
  verbose = TRUE,
  trunc = TRUE
) {
  computeACn(
    concentration = concentration,
    viability = viability,
    Hill_fit = Hill_fit,
    n = ifelse(viability_as_pct, 50, .5),
    conc_as_log = conc_as_log,
    viability_as_pct = viability_as_pct,
    verbose = verbose,
    trunc = trunc
  )
}
