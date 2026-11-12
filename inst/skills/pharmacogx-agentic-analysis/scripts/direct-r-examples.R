# Bounded development API examples: bundled/synthetic data only; no downloads.
# Run after selecting the verified package environment. Returns an in-memory list.
pgx_direct_r_examples <- function() {
  stopifnot(requireNamespace("PharmacoGx", quietly = TRUE))
  data_env <- new.env()
  utils::data("CCLEsmall", package = "PharmacoGx", envir = data_env)
  pset <- data_env$CCLEsmall
  methods::validObject(pset)
  drug <- PharmacoGx::treatmentNames(pset)[1]
  feature <- rownames(PharmacoGx::featureInfo(pset, "rna"))[1]
  response <- PharmacoGx::summarizeSensitivityProfiles(
    pset,
    drugs = drug,
    sensitivity.measure = "aac_recomputed",
    summary.stat = "mean",
    verbose = FALSE
  )
  stopifnot(length(response) > 0L, any(is.finite(response)))
  available_cores <- parallel::detectCores()
  if (length(available_cores) != 1L || is.na(available_cores)) {
    stop(paste(
      "CPU detection is unavailable in this execution environment.",
      "The tested drugSensitivitySig implementation requires a finite core count",
      "even with nthread = 1. Run in an approved R runtime with CPU detection;",
      "do not interpret this environment failure as a data or scientific result."
    ))
  }
  signature <- PharmacoGx::drugSensitivitySig(
    pset,
    mDataType = "rna",
    drugs = drug,
    features = feature,
    cells = PharmacoGx::sampleNames(pset),
    sensitivity.measure = "aac_recomputed",
    molecular.summary.stat = "mean",
    sensitivity.summary.stat = "mean",
    standardize = "SD",
    modeling.method = "anova",
    inference.method = "analytic",
    returnValues = c("estimate", "pvalue", "fdr"),
    nthread = 1,
    parallel.on = "drug",
    verbose = FALSE
  )
  stopifnot("estimate" %in% dimnames(signature)[[3]])
  dose <- c(0.1, 1, 10)
  viability <- c(0.9, 0.7, 0.4)
  auc <- PharmacoGx::computeAUC(
    dose,
    viability,
    area.type = "Actual",
    viability_as_pct = FALSE
  )
  aac <- PharmacoGx::computeAAC(
    dose,
    viability,
    area.type = "Actual",
    viability_as_pct = FALSE
  )
  stopifnot(isTRUE(all.equal(auc + aac, 1, tolerance = 1e-8)))
  fit_dose <- 10^seq(-2, 2, length.out = 9)
  fit_viability <- 0.2 + 0.8 / (1 + fit_dose)
  fit <- PharmacoGx::logLogisticRegression(
    fit_dose,
    fit_viability,
    fit_type = "hill",
    curve_direction = "decreasing",
    conc_as_log = FALSE,
    viability_as_pct = FALSE
  )
  prediction <- fit$E_inf +
    (fit$E0 - fit$E_inf) /
      (1 + (fit_dose / fit$EC50)^fit$HS)
  stopifnot(max(abs(prediction - fit_viability)) < 0.05)
  relative_ic50 <- PharmacoGx::computeIC50(
    fit_dose,
    Hill_fit = fit,
    reference = "relative",
    viability_as_pct = FALSE
  )
  absolute_ac50 <- PharmacoGx::computeAC50(
    fit_dose,
    Hill_fit = fit,
    viability_as_pct = FALSE
  )
  stopifnot(
    is.finite(relative_ic50),
    is.finite(absolute_ac50),
    absolute_ac50 > relative_ic50
  )
  mono_a <- c(0.8, 0.6)
  mono_b <- c(0.7, 0.5)
  bliss <- PharmacoGx::computeBliss(mono_a, mono_b)
  hsa <- PharmacoGx::computeHSA(mono_a, mono_b)
  stopifnot(
    isTRUE(all.equal(bliss, c(0.56, 0.30))),
    isTRUE(all.equal(hsa, c(0.7, 0.5)))
  )
  list(
    label = "development API checks; not biological validation",
    response = response,
    signature = signature,
    auc = auc,
    aac = aac,
    fit = fit,
    prediction = prediction,
    relative_ic50 = relative_ic50,
    absolute_ac50 = absolute_ac50,
    bliss_reference = bliss,
    hsa_reference = hsa,
    session = utils::sessionInfo()
  )
}
