#' @importFrom BiocParallel bplapply
.calculateSensitivitiesStar <- function(
  pSets = list(),
  exps = NULL,
  cap = NA,
  na.rm = TRUE,
  area.type = c("Fitted", "Actual"),
  nthread = 1
) {
  if (missing(area.type)) {
    area.type <- "Fitted"
  }
  if (is.null(exps)) {
    stop("expriments is empty!")
  }
  for (study in names(pSets)) {
    sensitivityProfiles(pSets[[study]])$auc_recomputed_star <- NA
  }
  if (!is.na(cap)) {
    trunc <- TRUE
  } else {
    trunc <- FALSE
  }

  for (i in seq_len(nrow(exps))) {
    ranges <- list()
    for (study in names(pSets)) {
      ranges[[study]] <- as.numeric(sensitivityRaw(pSets[[study]])[
        exps[i, study],
        ,
        "Dose"
      ])
    }
    ranges <- .getCommonConcentrationRange(ranges)
    names(ranges) <- names(pSets)
    for (study in names(pSets)) {
      myx <- as.numeric(
        sensitivityRaw(pSets[[study]])[exps[i, study], , "Dose"]
      ) %in%
        ranges[[study]]
      sensitivityRaw(pSets[[study]])[exps[i, study], !myx, ] <- NA
    }
  }

  op <- options()
  options(mc.cores = nthread)
  on.exit(options(op))

  for (study in names(pSets)) {
    auc_recomputed_star <- unlist(
      bplapply(
        rownames(sensitivityRaw(pSets[[study]])),
        FUN = function(experiment, exps, study, dataset, area.type) {
          if (!experiment %in% exps[, study]) {
            return(NA_real_)
          }
          return(
            computeAUC(
              concentration = as.numeric(dataset[experiment, , 1]),
              viability = as.numeric(dataset[experiment, , 2]),
              trunc = trunc,
              conc_as_log = FALSE,
              viability_as_pct = TRUE,
              area.type = area.type
            ) /
              100
          )
        },
        exps = exps,
        study = study,
        dataset = sensitivityRaw(pSets[[study]]),
        area.type = area.type
      )
    )
    sensitivityProfiles(pSets[[study]])$auc_recomputed_star <-
      auc_recomputed_star
  }
  return(pSets)
}

## This function computes AUC for the whole raw sensitivity data of a pset
.calculateFromRaw <- function(
  raw.sensitivity,
  cap = NA,
  nthread = 1,
  family = c("normal", "Cauchy"),
  scale = 0.07,
  n = 1
) {
  family <- match.arg(family)

  AUC <- vector(length = dim(raw.sensitivity)[1])
  names(AUC) <- dimnames(raw.sensitivity)[[1]]

  IC50 <- vector(length = dim(raw.sensitivity)[1])
  names(IC50) <- dimnames(raw.sensitivity)[[1]]

  trunc <- !is.na(cap)

  if (nthread == 1) {
    pars <- lapply(
      names(AUC),
      FUN = function(exp, raw.sensitivity, family, scale, n) {
        if (
          length(grep("///", raw.sensitivity[exp, , "Dose"])) > 0 ||
            all(is.na(raw.sensitivity[exp, , "Dose"]))
        ) {
          NA
        } else {
          logLogisticRegression(
            raw.sensitivity[exp, , "Dose"],
            raw.sensitivity[exp, , "Viability"],
            trunc = trunc,
            conc_as_log = FALSE,
            viability_as_pct = TRUE,
            family = family,
            scale = scale,
            median_n = n
          )
        }
      },
      raw.sensitivity = raw.sensitivity,
      family = family,
      scale = scale,
      n = n
    )
    names(pars) <- dimnames(raw.sensitivity)[[1]]
    AUC <- unlist(lapply(
      names(pars),
      FUN = function(exp, raw.sensitivity, pars) {
        if (any(is.na(pars[[exp]]))) {
          NA
        } else {
          computeAUC(
            concentration = raw.sensitivity[exp, , "Dose"],
            Hill_fit = pars[[exp]],
            trunc = trunc,
            conc_as_log = FALSE,
            viability_as_pct = TRUE
          )
        }
      },
      raw.sensitivity = raw.sensitivity,
      pars = pars
    ))
    IC50 <- unlist(lapply(
      names(pars),
      function(exp, pars) {
        if (any(is.na(pars[[exp]]))) {
          NA
        } else {
          computeIC50(
            Hill_fit = pars[[exp]],
            trunc = trunc,
            conc_as_log = FALSE,
            viability_as_pct = TRUE
          )
        }
      },
      pars = pars
    ))
  } else {
    pars <- parallel::mclapply(
      names(AUC),
      FUN = function(exp, raw.sensitivity, family, scale, n, trunc) {
        if (
          length(grep("///", raw.sensitivity[exp, , "Dose"])) > 0 ||
            all(is.na(raw.sensitivity[exp, , "Dose"]))
        ) {
          NA
        } else {
          logLogisticRegression(
            raw.sensitivity[exp, , "Dose"],
            raw.sensitivity[exp, , "Viability"],
            trunc = trunc,
            conc_as_log = FALSE,
            viability_as_pct = TRUE,
            family = family,
            scale = scale,
            median_n = n
          )
        }
      },
      raw.sensitivity = raw.sensitivity,
      family = family,
      scale = scale,
      n = n,
      trunc = trunc,
      mc.cores = nthread
    )
    names(pars) <- dimnames(raw.sensitivity)[[1]]
    AUC <- unlist(parallel::mclapply(
      names(pars),
      FUN = function(exp, raw.sensitivity, pars, trunc) {
        if (any(is.na(pars[[exp]]))) {
          NA
        } else {
          computeAUC(
            concentration = raw.sensitivity[exp, , "Dose"],
            Hill_fit = pars[[exp]],
            trunc = trunc,
            conc_as_log = FALSE,
            viability_as_pct = TRUE
          )
        }
      },
      raw.sensitivity = raw.sensitivity,
      pars = pars,
      trunc = trunc,
      mc.cores = nthread
    ))
    IC50 <- unlist(parallel::mclapply(
      names(pars),
      FUN = function(exp, pars, trunc) {
        if (any(is.na(pars[[exp]]))) {
          NA
        } else {
          computeIC50(
            Hill_fit = pars[[exp]],
            trunc = trunc,
            conc_as_log = FALSE,
            viability_as_pct = TRUE
          )
        }
      },
      pars = pars,
      trunc = trunc,
      mc.cores = nthread
    ))
  }
  names(AUC) <- dimnames(raw.sensitivity)[[1]]
  names(IC50) <- dimnames(raw.sensitivity)[[1]]

  return(list("AUC" = AUC, "IC50" = IC50, "pars" = pars))
}


## This function computes intersected concentration range between a list of
## concentration ranges
.getCommonConcentrationRange <- function(doses) {
  min.dose <- 0
  max.dose <- 10^100
  for (i in seq_len(length(doses))) {
    di <- sort(as.numeric(doses[[i]]))
    min.dose <- max(min.dose, min(di, na.rm = TRUE), na.rm = TRUE)
    max.dose <- min(max.dose, max(di, na.rm = TRUE), na.rm = TRUE)
    doses[[i]] <- di
  }
  common.ranges <- vector("list", length(doses))
  if (!is.finite(min.dose) || !is.finite(max.dose) || min.dose > max.dose) {
    return(common.ranges)
  }
  for (i in seq_len(length(doses))) {
    di <- doses[[i]]
    common.ranges[[i]] <- di[
      seq(
        which.min(abs(di - min.dose)),
        max(
          which(
            abs(di - max.dose) == min(abs(di - max.dose), na.rm = TRUE)
          )
        )
      )
    ]
  }
  return(common.ranges)
}

## predict viability from concentration data and curve parameters
.Hill <- function(x, pars) {
  internal <- .normalizeHillPars(
    hill_fit = pars,
    conc_as_log = TRUE,
    viability_as_pct = FALSE
  )

  .pgx_hill_curve(
    x,
    unname(internal[c("HS", "E0", "E_inf", "log10EC50")])
  )
}

## calculate residual of fit
## FIXME:: Why is this different from CoreGx?
#' @importFrom CoreGx .dmedncauchys .edmedncauchys
.residual <- function(
  x,
  y,
  n,
  pars,
  scale = 0.07,
  family = c("normal", "Cauchy"),
  trunc = FALSE
) {
  family <- match.arg(family)
  internal <- .normalizeHillPars(
    hill_fit = pars,
    conc_as_log = TRUE,
    viability_as_pct = FALSE
  )
  sum(.pgx_curve_residual(
    x = x,
    y = y,
    n = n,
    pars = unname(internal[c("HS", "E0", "E_inf", "log10EC50")]),
    f = .pgx_hill_curve,
    scale = scale,
    family = family,
    trunc = trunc,
    delta = 1
  ))
}

##FIXME:: Why is this different from CoreGx?
.meshEval <- function(
  log_conc,
  viability,
  lower_bounds = NULL,
  upper_bounds = NULL,
  density = NULL,
  scale = 0.07,
  n = 1,
  family = c("normal", "Cauchy"),
  trunc = FALSE
) {
  family <- match.arg(family)
  bounds <- .pgx_prepare_bounds(
    log_conc = log_conc,
    density = density,
    step = NULL,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    fit_type = "hill"
  )
  guess <- .pgx_initial_guess(
    log_conc = log_conc,
    viability = viability,
    lower_bounds = bounds$lower,
    upper_bounds = bounds$upper,
    fit_type = "hill"
  )
  .pgx_mesh_eval(
    x = log_conc,
    y = viability,
    f = .pgx_hill_curve,
    guess = unname(guess),
    lower_bounds = bounds$lower,
    upper_bounds = bounds$upper,
    density = bounds$density,
    n = n,
    scale = scale,
    family = family,
    trunc = trunc
  )
}

## FIXME:: Documentation?
#  Fits dose-response curves to data given by the user
#  and returns the AUC of the fitted curve, normalized to the length of the concentration range.
#
#  @param concentration `numeric` is a vector of drug concentrations.
#
#  @param viability `numeric` is a vector whose entries are the viability values observed in the presence of the
#  drug concentrations whose logarithms are in the corresponding entries of the log_conc, expressed as percentages
#  of viability in the absence of any drug.
#
#  @param trunc `logical`, if true, causes viability data to be truncated to lie between 0 and 1 before
#  curve-fitting is performed.
#' @importFrom CoreGx .getSupportVec
#' @export
#' @keywords internal
.computeAUCUnderFittedCurve <- function(
  concentration,
  viability,
  trunc = TRUE,
  verbose = FALSE
) {
  log_conc <- concentration
  #FIT CURVE AND CALCULATE IC50
  pars <- unlist(logLogisticRegression(
    log_conc,
    viability,
    conc_as_log = TRUE,
    viability_as_pct = FALSE,
    trunc = trunc
  ))
  internal <- .normalizeHillPars(
    hill_fit = pars,
    conc_as_log = TRUE,
    viability_as_pct = FALSE
  )
  x <- .getSupportVec(log_conc)
  y_hat <- .pgx_hill_curve(
    x,
    unname(internal[c("HS", "E0", "E_inf", "log10EC50")])
  )
  1 -
    caTools::trapz(x, y_hat) /
      (log_conc[length(log_conc)] - log_conc[1])
}

#This function is being used in computeSlope
.optimizeRegression <- function(x, y, x0 = -3, y0 = 100) {
  beta1 <- (sum(x * y) - y0 * sum(x)) / (sum(x * x) - x0 * sum(x))
  return(beta1)
}

updateMaxConc <- function(pSet) {
  sensitivityInfo(pSet)$max.conc <- apply(
    sensitivityRaw(pSet)[,, "Dose"],
    1,
    max,
    na.rm = TRUE
  )
  return(pSet)
}
