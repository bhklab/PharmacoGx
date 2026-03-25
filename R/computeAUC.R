#' Computes the normalized response area for a drug dose viability curve
#'
#' Returns the normalized response area over the experiment's concentration
#' range. Historically this function is named `computeAUC()`, but the quantity
#' returned is the area under the response curve `(1 - viability)` on the log10
#' concentration scale, rather than the classical area under the viability curve
#' itself. Larger values therefore imply greater drug sensitivity. The function
#' can calculate this quantity using either a fitted Hill curve or a trapezoidal
#' integral of the observed data. Alternatively, parameters returned by
#' logLogisticRegression can be supplied directly.
#'
#' @examples
#' dose <- c(0.0025,0.008,0.025,0.08,0.25,0.8,2.53,8)
#' viability <- c(108.67,111,102.16,100.27,90,87,74,57)
#' computeAUC(dose, viability)
#'
#'
#' @param concentration `numeric` vector of drug concentrations.
#' @param viability `numeric` vector of observed viabilities aligned to
#'   `concentration`. Viability can be supplied as percentages or proportions,
#'   depending on `viability_as_pct`.
#' @param Hill_fit `list` or `vector` of Hill-curve parameters as returned by
#'   logLogisticRegression. When `conc_as_log = TRUE`, EC50 values are assumed to
#'   already be on the log10 scale. When `viability_as_pct = TRUE`, response
#'   parameters such as `E0` and `E_inf` are assumed to be expressed as
#'   percentages; otherwise they are assumed to be proportions.
#' @param conc_as_log `logical`, if `TRUE`, assumes that log10-concentration data
#'   has been given rather than concentration data.
#' @param viability_as_pct `logical`, if `FALSE`, assumes that viability is
#'   given as a proportion rather than a percentage, and returns the normalized
#'   response area on the 0-1 scale. Otherwise, viability is interpreted as a
#'   percentage and the result is returned on the 0-100 scale.
#' @param trunc `logical`, if `TRUE`, clips viability data to lie between 0 and
#'   1 after scale normalization and before either the actual or fitted area is
#'   computed.
#' @param area.type Character string indicating whether to compute the normalized
#'   response area using the observed data (`"Actual"`) or a fitted curve
#'   (`"Fitted"`).
#' @param verbose `logical`, if true, causes warnings thrown by the function to be printed.
#' @return Numeric normalized response area value.
#'
#' @export
#' @import caTools
#' @importFrom stats integrate
computeAUC <- function(
  concentration,
  viability,
  Hill_fit,
  conc_as_log = FALSE,
  viability_as_pct = TRUE,
  trunc = TRUE,
  area.type = c("Fitted", "Actual"),
  verbose = TRUE
  #, ...
) {
  if (missing(concentration)) {
    stop("The concentration values to integrate over must always be provided.")
  }
  if (missing(area.type)) {
    area.type <- "Fitted"
  } else {
    area.type <- match.arg(area.type)
  }
  if (area.type == "Fitted" && missing(Hill_fit)) {
    if (missing(viability)) {
      stop(
        "To fit a curve (area.type='Fitted'), supply raw viability or provide Hill_fit."
      )
    }
    Hill_fit <- logLogisticRegression(
      concentration,
      viability,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose
    )
    cleanData <- sanitizeInput(
      conc = concentration,
      Hill_fit = Hill_fit,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose
    )
    pars <- cleanData[["Hill_fit"]]
    concentration <- cleanData[["log_conc"]]
  } else if (area.type == "Fitted" && !missing(Hill_fit)) {
    cleanData <- sanitizeInput(
      conc = concentration,
      viability = viability,
      Hill_fit = Hill_fit,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose
    )
    pars <- cleanData[["Hill_fit"]]
    concentration <- cleanData[["log_conc"]]
  } else if (area.type == "Actual" && !missing(viability)) {
    cleanData <- sanitizeInput(
      conc = concentration,
      viability = viability,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose
    )
    concentration <- cleanData[["log_conc"]]
    viability <- cleanData[["viability"]]
  } else if (area.type == "Actual" && missing(viability)) {
    stop(
      "To calculate the actual area using a trapezoid integral, the raw viability values are needed!"
    )
  }

  if (length(concentration) < 2) {
    return(NA)
  }

  a <- min(concentration)
  b <- max(concentration)
  if (area.type == "Actual") {
    trapezoid.integral <- caTools::trapz(concentration, viability)
    AUC <- 1 - trapezoid.integral / (b - a)
  } else {
    if (b == a) {
      return(NA_real_)
    }

    if (
      length(pars) >= 8 &&
        all(
          c(
            "HS1",
            "E0",
            "E_inf1",
            "HS2",
            "E_inf2",
            "log10EC50_1",
            "log10EC50_2",
            "Frac"
          ) %in%
            names(pars)
        )
    ) {
      curve_fun <- .pgx_biphasic_curve
      param_vec <- unname(pars[c(
        "HS1",
        "E0",
        "E_inf1",
        "HS2",
        "E_inf2",
        "log10EC50_1",
        "log10EC50_2",
        "Frac"
      )])
    } else {
      curve_fun <- .pgx_hill_curve
      param_vec <- unname(pars[c("HS", "E0", "E_inf", "log10EC50")])
    }

    integral <- try(
      stats::integrate(
        f = function(x) curve_fun(x, param_vec),
        lower = a,
        upper = b
      ),
      silent = TRUE
    )

    if (inherits(integral, "try-error")) {
      AUC <- NA_real_
    } else {
      AUC <- 1 - integral$value / (b - a)
    }
  }

  if (viability_as_pct) {
    AUC <- AUC * 100
  }

  return(AUC)
}
