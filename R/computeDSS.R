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

  ec50 <- if ("EC50" %in% names(pars)) pars[["EC50"]] else NA_real_
  log_ec50 <- if ("log10EC50" %in% names(pars)) {
    pars[["log10EC50"]]
  } else {
    NA_real_
  }

  if (!is.na(ec50)) {
    log_ec50 <- if (conc_as_log) ec50 else log10(ec50)
  }

  if (is.na(log_ec50)) {
    stop("Unable to determine EC50 from Hill fit parameters.")
  }

  if (log_ec50 > max(concentration)) {
    return(0)
  }

  t_param_pct <- if (viability_as_pct) t_param else t_param * 100
  t_param_fraction <- t_param_pct / 100

  hill_external <- list(
    HS = pars[["HS"]],
    E0 = pars[["E0"]],
    E_inf = pars[["E_inf"]],
    log10EC50 = log_ec50
  )

  x2 <- max(concentration)
  x1 <- computeICn(
    concentration = concentration,
    Hill_fit = hill_external,
    n = t_param_fraction,
    conc_as_log = TRUE,
    viability_as_pct = FALSE
  )
  if (!is.finite(x1)) {
    return(0)
  }

  x1 <- max(x1, min(concentration))

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
  DSS <- DSS / log(100 - e_inf_pct)
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
