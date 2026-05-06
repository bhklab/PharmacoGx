#' Fit dose-response curves for screening data
#'
#' Fits dose-response curves for each cell/drug experiment contained in
#' `scrn_data` using the updated four-parameter Hill or biphasic model.
#'
#' @param scrn_data Data frame containing at least the columns `cell_id`,
#'   `drug_id`, `conc` (numeric concentrations) and `viability` (fractional
#'   viability).
#' @param drug_subset Optional character or numeric vector of drug identifiers to
#'   retain.
#' @param cell_subset Optional character vector of cell identifiers to retain.
#' @param output_type Character string indicating which result to return. One of
#'   `"curves"`, `"metrics"`, or `"all"`.
#' @param main_fit_func Character string selecting the primary model family.
#'   One of `"hill"` or `"biphasic"`.
#'
#' @return Depending on `output_type`, either a `data.table` with fitted curve
#'   coordinates, a `data.table` of summary metrics, or a list containing both.
#'   Returns `NULL` when no experiments remain after filtering or minimum
#'   measurement checks.
#'
#' @examples
#' scrn <- data.frame(
#'   cell_id = rep(c("cellA", "cellB"), each = 3),
#'   drug_id = rep("drugX", 6),
#'   conc = rep(c(0.1, 1, 10), 2),
#'   viability = c(0.95, 0.75, 0.5, 0.9, 0.7, 0.45)
#' )
#' curveFittingPGX(scrn, output_type = "metrics", main_fit_func = "hill")
#'
#' @export
#' @importFrom data.table as.data.table data.table rbindlist setorder
#' @importFrom CoreGx .getSupportVec
curveFittingPGX <- function(
  scrn_data,
  drug_subset = NULL,
  cell_subset = NULL,
  output_type = c("curves", "metrics", "all"),
  main_fit_func = c("hill", "biphasic")
) {
  output_type <- match.arg(output_type)
  main_fit_func <- match.arg(tolower(main_fit_func), c("hill", "biphasic"))

  required_cols <- c("cell_id", "drug_id", "conc", "viability")
  missing_cols <- setdiff(required_cols, names(scrn_data))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dt <- data.table::as.data.table(scrn_data)[, ..required_cols]
  dt[, experiment := paste(cell_id, drug_id, sep = "_")]
  data.table::setorder(dt, experiment, conc)
  dt <- unique(dt)

  if (!is.null(drug_subset)) {
    dt <- dt[drug_id %in% drug_subset]
  }
  if (!is.null(cell_subset)) {
    dt <- dt[cell_id %in% cell_subset]
  }
  if (!nrow(dt)) {
    return(NULL)
  }

  keep_experiments <- dt[, .N, by = experiment][N >= 3, experiment]
  dt <- dt[experiment %in% keep_experiments]
  if (!nrow(dt)) {
    return(NULL)
  }

  experiment_list <- split(dt, dt$experiment)

  fit_results <- lapply(experiment_list, function(exp_dt) {
    conc <- exp_dt$conc
    viability_pct <- pmin(exp_dt$viability * 100, 100)

    fit <- try(
      logLogisticRegression(
        conc = conc,
        viability = viability_pct,
        conc_as_log = FALSE,
        viability_as_pct = TRUE,
        trunc = TRUE,
        fit_type = main_fit_func
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      return(NULL)
    }

    rsq <- attr(fit, "Rsquare")
    auc <- tryCatch(
      computeAUC(
        concentration = conc,
        Hill_fit = fit,
        conc_as_log = FALSE,
        viability_as_pct = TRUE,
        trunc = TRUE
      ),
      warning = function(w) NA_real_
    )
    aac <- tryCatch(
      computeAAC(
        concentration = conc,
        Hill_fit = fit,
        conc_as_log = FALSE,
        viability_as_pct = TRUE,
        trunc = TRUE
      ),
      warning = function(w) NA_real_
    )

    support <- CoreGx::.getSupportVec(log10(conc))

    if (main_fit_func == "hill") {
      ic50 <- tryCatch(
        computeIC50(
          Hill_fit = fit,
          conc_as_log = FALSE,
          viability_as_pct = TRUE,
          trunc = TRUE
        ),
        warning = function(w) NA_real_
      )
      fit_internal <- .normalizeHillPars(
        hill_fit = fit,
        conc_as_log = FALSE,
        viability_as_pct = TRUE
      )
      curve_fun <- .pgx_hill_curve
      curve_pars <- unname(fit_internal[c("HS", "E0", "E_inf", "log10EC50")])
      metrics <- data.table::data.table(
        experiment = exp_dt$experiment[1],
        min_conc = min(conc),
        max_conc = max(conc),
        HS = fit$HS,
        E0 = fit$E0,
        E_inf = fit$E_inf,
        EC50 = fit$EC50,
        emax = curve_fun(log10(max(conc)), curve_pars) * 100,
        Rsquare = rsq,
        auc = auc,
        aac = aac,
        ic50 = ic50
      )
    } else {
      hill_fit <- try(
        logLogisticRegression(
          conc = conc,
          viability = viability_pct,
          conc_as_log = FALSE,
          viability_as_pct = TRUE,
          trunc = TRUE,
          fit_type = "hill"
        ),
        silent = TRUE
      )
      ic50 <- if (inherits(hill_fit, "try-error")) {
        NA_real_
      } else {
        tryCatch(
          computeIC50(
            Hill_fit = hill_fit,
            conc_as_log = FALSE,
            viability_as_pct = TRUE,
            trunc = TRUE
          ),
          warning = function(w) NA_real_
        )
      }
      fit_internal <- .normalizeBiphasicPars(
        hill_fit = fit,
        conc_as_log = FALSE,
        viability_as_pct = TRUE
      )
      curve_fun <- .pgx_biphasic_curve
      curve_pars <- unname(fit_internal[c(
        "HS1",
        "E0",
        "E_inf1",
        "HS2",
        "E_inf2",
        "log10EC50_1",
        "log10EC50_2",
        "Frac"
      )])
      metrics <- data.table::data.table(
        experiment = exp_dt$experiment[1],
        min_conc = min(conc),
        max_conc = max(conc),
        HS1 = fit$HS1,
        E0 = fit$E0,
        E_inf1 = fit$E_inf1,
        HS2 = fit$HS2,
        E_inf2 = fit$E_inf2,
        EC50_1 = fit$EC50_1,
        EC50_2 = fit$EC50_2,
        Frac = fit$Frac,
        emax = curve_fun(log10(max(conc)), curve_pars) * 100,
        Rsquare = rsq,
        auc = auc,
        aac = aac,
        ic50 = ic50
      )
    }

    fitted_vals <- curve_fun(support, curve_pars) * 100

    curves <- data.table::data.table(
      experiment = exp_dt$experiment[1],
      log10_conc = support,
      conc = 10^support,
      fitted_viability = fitted_vals,
      drug_id = exp_dt$drug_id[1],
      cell_id = exp_dt$cell_id[1]
    )

    observed <- data.table::data.table(
      experiment = exp_dt$experiment[1],
      conc = conc,
      observed_viability = viability_pct
    )
    curves <- merge(
      curves,
      observed,
      by = c("experiment", "conc"),
      all.x = TRUE
    )

    list(metrics = metrics, curves = curves)
  })

  fit_results <- Filter(Negate(is.null), fit_results)
  if (!length(fit_results)) {
    return(NULL)
  }

  metrics_dt <- data.table::rbindlist(
    lapply(fit_results, `[[`, "metrics"),
    use.names = TRUE,
    fill = TRUE
  )
  curves_dt <- data.table::rbindlist(
    lapply(fit_results, `[[`, "curves"),
    use.names = TRUE,
    fill = TRUE
  )

  if (output_type == "metrics") {
    return(metrics_dt)
  }
  if (output_type == "curves") {
    return(curves_dt)
  }

  list(metrics = metrics_dt, curves = curves_dt)
}
