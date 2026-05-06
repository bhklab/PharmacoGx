# Internal helpers for relative and absolute inhibition calculations.
.compute_inhibition_fraction <- function(n, viability_as_pct) {
  if (viability_as_pct) {
    n / 100
  } else {
    n
  }
}

.compute_target_viability <- function(
  pars,
  n,
  viability_as_pct,
  reference = c("relative", "absolute")
) {
  reference <- match.arg(reference)
  frac <- .compute_inhibition_fraction(
    n = n,
    viability_as_pct = viability_as_pct
  )

  if (reference == "absolute") {
    return(1 - frac)
  }

  e0 <- pars[["E0"]]
  einf <- pars[["E_inf"]]
  e0 - frac * (e0 - einf)
}

.compute_threshold_dose <- function(pars, target, conc_as_log) {
  e0 <- pars[["E0"]]
  einf <- pars[["E_inf"]]
  hs <- pars[["HS"]]
  log_ec50 <- pars[["log10EC50"]]

  top <- max(e0, einf)
  bottom <- min(e0, einf)
  tol <- sqrt(.Machine$double.eps)

  if (target >= (top - tol)) {
    return(ifelse(conc_as_log, -Inf, 0))
  }
  if (target <= (bottom + tol)) {
    return(Inf)
  }
  if (hs <= 0) {
    return(NA_real_)
  }

  ratio <- (e0 - einf) / (target - einf) - 1
  if (!is.finite(ratio) || ratio <= 0) {
    return(NA_real_)
  }

  log_icn <- log_ec50 + (1 / hs) * log10(ratio)
  if (conc_as_log) {
    log_icn
  } else {
    10^log_icn
  }
}

.is_biphasic_fit <- function(pars) {
  required_names <- c(
    "HS1",
    "E0",
    "E_inf1",
    "HS2",
    "E_inf2",
    "log10EC50_1",
    "log10EC50_2",
    "Frac"
  )
  length(pars) >= length(required_names) &&
    all(required_names %in% names(pars))
}

#' @describeIn computeICn Returns the ACn of a drug dose response curve using an
#'   absolute viability threshold.
#' @export
computeACn <- function(
  concentration,
  viability,
  Hill_fit,
  n,
  conc_as_log = FALSE,
  viability_as_pct = TRUE,
  verbose = TRUE,
  trunc = TRUE
) {
  .computeICn_impl(
    concentration = concentration,
    viability = viability,
    Hill_fit = Hill_fit,
    n = n,
    reference = "absolute",
    conc_as_log = conc_as_log,
    viability_as_pct = viability_as_pct,
    verbose = verbose,
    trunc = trunc
  )
}

.computeICn_impl <- function(
  concentration,
  viability,
  Hill_fit,
  n,
  reference,
  conc_as_log = FALSE,
  viability_as_pct = TRUE,
  verbose = TRUE,
  trunc = TRUE,
  warn_absolute = FALSE
) {
  if (missing(Hill_fit) && !missing(concentration) && !missing(viability)) {
    Hill_fit <- logLogisticRegression(
      conc = concentration,
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
  } else if (!missing(Hill_fit)) {
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
  } else {
    stop(
      "Insufficient information to calculate ICn. Please enter concentration and viability or Hill parameters."
    )
  }

  if (.is_biphasic_fit(pars)) {
    stop(
      "computeICn() and computeIC50() currently support Hill fits only. ",
      "For biphasic fits, fit a Hill model for IC metrics or use an AUC-based ",
      "summary.",
      call. = FALSE
    )
  }

  if (reference == "absolute" && warn_absolute) {
    warning(
      paste(
        "`reference = \"absolute\"` is deprecated for computeICn().",
        "Use computeACn() instead."
      ),
      call. = FALSE
    )
  }

  target <- .compute_target_viability(
    pars = pars,
    n = n,
    viability_as_pct = viability_as_pct,
    reference = reference
  )
  .compute_threshold_dose(
    pars = pars,
    target = target,
    conc_as_log = conc_as_log
  )
}

#' Computes inhibition concentrations for a drug dose viability curve
#'
#' Returns the concentration corresponding to the requested percent inhibition.
#' By default, inhibition is measured relative to the fitted dynamic range of the
#' curve, so `computeIC50()` returns the half-maximal inhibitory concentration.
#' Absolute-threshold helpers `computeACn()` and `computeAC50()` are provided for
#' the historical "absolute viability" behaviour. Use `computeAC50()`, not the
#' default `computeIC50()`, when comparing against published AC50 values that
#' represent an absolute 50% viability threshold.
#'
#' @examples
#' dose <- c(0.0025, 0.008, 0.025, 0.08, 0.25, 0.8, 2.53, 8)
#' viability <- c(108.67, 111, 102.16, 100.27, 90, 87, 74, 57)
#'
#' # Relative half-maximal inhibition concentration.
#' computeIC50(dose, viability)
#'
#' # Absolute 50% viability threshold, suitable for published AC50 comparisons.
#' computeAC50(dose, viability)
#'
#' computeICn(dose, viability, n = 10)
#'
#' @param concentration `numeric` is a vector of drug concentrations.
#' @param viability `numeric` is a vector whose entries are the viability values observed in the presence of the
#' drug concentrations whose logarithms are in the corresponding entries of conc, where viability 0
#' indicates that all cells died, and viability 1 indicates that the drug had no effect on the cells.
#' @param Hill_fit `list` or `vector` In the order: c("Hill Slope", "E_inf", "EC50"), the parameters of a Hill Slope
#' as returned by logLogisticRegression. If conc_as_log is set then the function assumes logEC50 is passed in, and if
#' viability_as_pct flag is set, it assumes E_inf is passed in as a percent. Otherwise, E_inf is assumed to be a decimal,
#' and EC50 as a concentration.
#' @param n `numeric` The inhibition level to compute. If `viability_as_pct =
#'   TRUE` it is treated as a percent inhibition; otherwise it is assumed to be
#'   a decimal fraction.
#' @param reference `character(1)` Whether inhibition should be measured
#'   relative to the fitted dynamic range (`"relative"`, default) or against an
#'   absolute viability threshold (`"absolute"`). The absolute-threshold mode is
#'   deprecated in favour of `computeACn()`.
#' @param conc_as_log `logical`, if true, assumes that log10-concentration data has been given rather than concentration data,
#' and that log10(ICn) should be returned instead of ICn.
#' @param viability_as_pct `logical`, if false, assumes that viability is given as a decimal rather
#' than a percentage, and that E_inf passed in as decimal.
#' @param trunc `logical`, if true, causes viability data to be truncated to lie between 0 and 1 before
#' curve-fitting is performed.
#' @param verbose `logical`, if true, causes warnings thrown by the function to be printed.
#' @return a numeric value for the concentration of the requested inhibition
#'   level.
#' @export
computeICn <- function(
  concentration,
  viability,
  Hill_fit,
  n,
  reference = c("relative", "absolute"),
  conc_as_log = FALSE,
  viability_as_pct = TRUE,
  verbose = TRUE,
  trunc = TRUE
) {
  reference <- match.arg(reference)

  .computeICn_impl(
    concentration = concentration,
    viability = viability,
    Hill_fit = Hill_fit,
    n = n,
    reference = reference,
    conc_as_log = conc_as_log,
    viability_as_pct = viability_as_pct,
    verbose = verbose,
    trunc = trunc,
    warn_absolute = reference == "absolute"
  )
}
