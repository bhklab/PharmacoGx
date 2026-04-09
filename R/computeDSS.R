##TODO:: Add function documentation
computeDSS <- function(
  concentration,
  viability,
  Hill_fit,
  t_param = 10,
  conc_as_log = FALSE,
  viability_as_pct = TRUE,
  trunc = TRUE,
  verbose = TRUE,
  dss_type = 3,
  censor = FALSE
  #, ...
) {
  if (missing(concentration)) {
    stop("The concentration values to integrate over must always be provided.")
  }
  if (missing(Hill_fit)) {
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
  } else {
    # sanitizeInput normalizes concentration, returning log10(conc) when
    # conc_as_log = FALSE and leaving it unchanged otherwise
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
  }

  if ("log10EC50" %in% names(pars)) {
    log_ec50 <- as.numeric(pars[["log10EC50"]])
  } else if ("EC50" %in% names(pars)) {
    ec50 <- as.numeric(pars[["EC50"]])
    if (!is.finite(ec50) || length(ec50) != 1 || ec50 <= 0) {
      stop("Hill fit EC50 must be a single positive finite value.")
    }
    log_ec50 <- log10(ec50)
  } else {
    stop("Hill fit parameters must include 'log10EC50' or 'EC50'.")
  }

  if (!is.finite(log_ec50) || length(log_ec50) != 1) {
    stop("Failed to derive a valid log10 EC50 from Hill fit parameters.")
  }

  if (log_ec50 > max(concentration)) {
    return(0)
  }

  if (!viability_as_pct && t_param <= 1) {
    warning(
      paste(
        "Fractional `t_param` values are deprecated when",
        "`viability_as_pct = FALSE`.",
        "Interpreting the input as a fractional inhibition threshold."
      ),
      call. = FALSE
    )
    t_param_pct <- t_param * 100
  } else {
    t_param_pct <- t_param
  }
  t_param_fraction <- t_param_pct / 100

  hill_external <- list(
    HS = pars[["HS"]],
    E0 = pars[["E0"]],
    E_inf = pars[["E_inf"]],
    log10EC50 = log_ec50
  )

  x2 <- max(concentration)
  x1 <- computeACn(
    concentration = concentration,
    Hill_fit = hill_external,
    n = t_param_fraction,
    conc_as_log = TRUE,
    viability_as_pct = FALSE
  )
  if (!is.finite(x1)) {
    return(0)
  }

  x1 <- min(max(x1, min(concentration)), x2)

  e_inf_pct <- hill_external$E_inf * 100
  if (censor) {
    if (e_inf_pct > 50) {
      return(NA)
    } else if (all(concentration < log_ec50)) {
      return(0)
    }
  }

  auc_fraction <- computeAUC(
    concentration = c(x1, x2),
    Hill_fit = hill_external,
    conc_as_log = TRUE,
    viability_as_pct = FALSE,
    verbose = verbose,
    trunc = trunc
  )
  AUC <- auc_fraction * 100

  DSS <- (AUC * (x2 - x1) - t_param_pct * (x2 - x1)) /
    ((100 - t_param_pct) * (max(concentration) - min(concentration)))
  if (dss_type == 1) {
    return(DSS)
  }
  delta <- 100 - e_inf_pct
  if (!is.finite(delta)) {
    if (verbose) {
      warning("Asymptotic viability is non-finite; returning NA for DSS.")
    }
    return(NA_real_)
  }
  if (delta <= 0) {
    if (verbose) {
      warning(
        "E_inf exceeds or equals 100%; DSS type 2/3 undefined, returning NA."
      )
    }
    return(NA_real_)
  }
  DSS <- DSS / log(delta)
  if (dss_type == 2) {
    return(DSS)
  }
  DSS <- DSS * (x2 - x1) / (max(concentration) - min(concentration))
  if (dss_type == 3) {
    return(DSS)
  } else {
    stop("Invalid DSS type entered.")
  }
}
