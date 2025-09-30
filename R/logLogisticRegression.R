#' Fit dose-response curves with flexible Hill or biphasic models
#'
#' Fits dose-response data using either a four-parameter Hill curve or an
#' eight-parameter biphasic curve. Compared to the historical implementation,
#' the updated fitting routine estimates the top asymptote (`E0`) in addition to
#' the bottom asymptote (`E_inf`), improves the optimisation strategy through a
#' robust residual, and supports biphasic response profiles.
#'
#' @param conc Numeric vector of drug concentrations.
#' @param viability Numeric vector of viabilities aligned to `conc`.
#' @param density Optional numeric vector controlling the coarse mesh search. If
#'   omitted, sensible defaults are derived from the data and the selected
#'   `fit_type`.
#' @param step Optional numeric vector defining the initial step sizes for the
#'   pattern search used when optimisation stalls. Defaults are based on
#'   `density` when not supplied.
#' @param precision Positive numeric tolerance controlling termination of the
#'   pattern search.
#' @param lower_bounds,upper_bounds Optional numeric vectors specifying bounds
#'   for the optimisation. Missing entries are substituted with data-driven
#'   defaults.
#' @param scale Positive numeric value used by the median-based likelihood when
#'   `family = "Cauchy"`.
#' @param family Character string selecting the residual family, either
#'   `"normal"` (default) or `"Cauchy"`.
#' @param median_n Integer number of replicate measurements represented by each
#'   viability value.
#' @param conc_as_log Logical; when `TRUE`, `conc` is assumed to already be on a
#'   log10 scale and `EC50` is returned on that scale.
#' @param viability_as_pct Logical; when `TRUE`, viabilities are treated as
#'   percentages and the returned `E0`/`E_inf` are expressed in percent.
#' @param trunc Logical flag indicating whether viability values should be
#'   truncated to the unit interval prior to fitting.
#' @param verbose Logical flag enabling diagnostic warnings.
#' @param fit_type Character string selecting the curve family to fit. One of
#'   `"hill"` (default) or `"biphasic"`.
#'
#' @return A named list of fitted parameters with an `Rsquare` attribute. For
#'   Hill fits: `HS`, `E0`, `E_inf`, and `EC50`. For biphasic fits: `HS1`, `E0`,
#'   `E_inf1`, `HS2`, `E_inf2`, `EC50_1`, `EC50_2`, and `Frac`.
#'
#' @examples
#' dose <- c(0.0025, 0.008, 0.025, 0.08, 0.25, 0.8, 2.53, 8)
#' viability <- c(108.67, 111, 102.16, 100.27, 90, 87, 74, 57)
#' logLogisticRegression(dose, viability)
#'
#' @export
#' @importFrom stats approx integrate optim var
#' @importFrom CoreGx .dmedncauchys .edmedncauchys
logLogisticRegression <- function(
  conc,
  viability,
  density = NULL,
  step = NULL,
  precision = 1e-4,
  lower_bounds = NULL,
  upper_bounds = NULL,
  scale = 0.07,
  family = c("normal", "Cauchy"),
  median_n = 1,
  conc_as_log = FALSE,
  viability_as_pct = TRUE,
  trunc = TRUE,
  verbose = FALSE,
  fit_type = c("hill", "biphasic")
) {
  fit_type <- match.arg(tolower(fit_type), c("hill", "biphasic"))
  family <- match.arg(tolower(family), c("normal", "cauchy"))
  family <- if (family == "cauchy") "Cauchy" else "normal"

  if (!is.numeric(median_n) || median_n != round(median_n)) {
    stop("`median_n` must be an integer.")
  }
  if (median_n < 1) {
    stop("`median_n` must be greater than or equal to 1.")
  }
  if (precision <= 0) {
    stop("`precision` must be strictly positive.")
  }
  if (scale <= 0) {
    stop("`scale` must be strictly positive.")
  }

  clean_data <- sanitizeInput(
    conc = conc,
    viability = viability,
    conc_as_log = conc_as_log,
    viability_as_pct = viability_as_pct,
    trunc = trunc,
    verbose = verbose
  )

  log_conc <- clean_data[["log_conc"]]
  viability_clean <- clean_data[["viability"]]

  if (length(log_conc) < 3L) {
    stop("At least three unique concentrations are required to fit a curve.")
  }

  bounds <- .pgx_prepare_bounds(
    log_conc = log_conc,
    density = density,
    step = step,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    fit_type = fit_type
  )

  guess <- .pgx_initial_guess(
    log_conc = log_conc,
    viability = viability_clean,
    lower_bounds = bounds$lower,
    upper_bounds = bounds$upper,
    fit_type = fit_type
  )

  curve_fun <- if (fit_type == "hill") .pgx_hill_curve else .pgx_biphasic_curve

  fitted <- .pgx_fit_curve(
    x = log_conc,
    y = viability_clean,
    f = curve_fun,
    density = bounds$density,
    step = bounds$step,
    precision = precision,
    lower_bounds = bounds$lower,
    upper_bounds = bounds$upper,
    scale = scale,
    family = family,
    median_n = median_n,
    trunc = trunc,
    gritty_guess = guess
  )

  params <- fitted$pars

  if (fit_type == "hill") {
    result <- list(
      HS = params[["HS"]],
      E0 = if (viability_as_pct) params[["E0"]] * 100 else params[["E0"]],
      E_inf = if (viability_as_pct) {
        params[["E_inf"]] * 100
      } else {
        params[["E_inf"]]
      },
      EC50 = if (conc_as_log) {
        params[["log10EC50"]]
      } else {
        10^params[["log10EC50"]]
      }
    )
  } else {
    result <- list(
      HS1 = params[["HS1"]],
      E0 = if (viability_as_pct) params[["E0"]] * 100 else params[["E0"]],
      E_inf1 = if (viability_as_pct) {
        params[["E_inf1"]] * 100
      } else {
        params[["E_inf1"]]
      },
      HS2 = params[["HS2"]],
      E_inf2 = if (viability_as_pct) {
        params[["E_inf2"]] * 100
      } else {
        params[["E_inf2"]]
      },
      EC50_1 = if (conc_as_log) {
        params[["log10EC50_1"]]
      } else {
        10^params[["log10EC50_1"]]
      },
      EC50_2 = if (conc_as_log) {
        params[["log10EC50_2"]]
      } else {
        10^params[["log10EC50_2"]]
      },
      Frac = params[["Frac"]]
    )
  }

  if (trunc) {
    upper_bound <- if (viability_as_pct) 100 else 1
    lower_bound <- 0
    if (fit_type == "hill") {
      result$E0 <- pmin(pmax(result$E0, lower_bound), upper_bound)
      result$E_inf <- pmin(pmax(result$E_inf, lower_bound), upper_bound)
    } else {
      result$E0 <- pmin(pmax(result$E0, lower_bound), upper_bound)
      result$E_inf1 <- pmin(pmax(result$E_inf1, lower_bound), upper_bound)
      result$E_inf2 <- pmin(pmax(result$E_inf2, lower_bound), upper_bound)
    }
  }

  attr(result, "Rsquare") <- fitted$Rsquare
  result
}

.pgx_prepare_bounds <- function(
  log_conc,
  density,
  step,
  lower_bounds,
  upper_bounds,
  fit_type
) {
  if (fit_type == "hill") {
    default_lower <- c(0.2, 0, 0, min(log_conc) - 2)
    default_upper <- c(4.5, 1.2, 1.2, max(log_conc) + 2)
    default_density <- c(2, 10, 10, 5)
  } else {
    default_lower <- c(
      0.1,
      0,
      0,
      min(log_conc) - 3,
      0.1,
      0,
      min(log_conc) - 3,
      0
    )
    default_upper <- c(
      5,
      1.5,
      1.5,
      max(log_conc) + 3,
      5,
      1.5,
      max(log_conc) + 3,
      1
    )
    default_density <- c(2, 10, 10, 5, 10, 5, 5, 5)
  }

  if (is.null(density)) {
    density <- default_density
  } else {
    density <- as.numeric(density)
    if (fit_type == "hill" && length(density) == 3) {
      density <- c(density[1], default_density[2], density[2], density[3])
    }
    if (length(density) != length(default_density)) {
      stop(
        "`density` must match the parameter count for the selected fit type."
      )
    }
  }

  if (is.null(step)) {
    step <- 0.5 / density
  } else {
    step <- as.numeric(step)
    if (fit_type == "hill" && length(step) == 3) {
      step <- c(step[1], 0.5 / density[2], step[2], step[3])
    }
    if (length(step) != length(default_density)) {
      stop("`step` must match the parameter count for the selected fit type.")
    }
  }

  if (is.null(lower_bounds)) {
    lower_bounds <- default_lower
  } else {
    lower_bounds <- as.numeric(lower_bounds)
    if (fit_type == "hill" && length(lower_bounds) == 3) {
      lower_bounds <- c(lower_bounds[1], 0, lower_bounds[2], lower_bounds[3])
    }
    if (length(lower_bounds) != length(default_lower)) {
      stop(
        "`lower_bounds` must match the parameter count for the selected fit type."
      )
    }
  }

  if (is.null(upper_bounds)) {
    upper_bounds <- default_upper
  } else {
    upper_bounds <- as.numeric(upper_bounds)
    if (fit_type == "hill" && length(upper_bounds) == 3) {
      upper_bounds <- c(upper_bounds[1], 1.2, upper_bounds[2], upper_bounds[3])
    }
    if (length(upper_bounds) != length(default_upper)) {
      stop(
        "`upper_bounds` must match the parameter count for the selected fit type."
      )
    }
  }

  if (any(upper_bounds - lower_bounds <= 0)) {
    stop("Upper bounds must exceed lower bounds for all parameters.")
  }

  if (any(density <= 0) || any(step <= 0)) {
    stop("`density` and `step` must contain positive values.")
  }

  list(
    lower = lower_bounds,
    upper = upper_bounds,
    density = density,
    step = step
  )
}


.pgx_initial_guess <- function(
  log_conc,
  viability,
  lower_bounds,
  upper_bounds,
  fit_type
) {
  finite_mask <- is.finite(log_conc) & is.finite(viability)
  log_conc_valid <- log_conc[finite_mask]
  viability_valid <- viability[finite_mask]
  if (!length(log_conc_valid)) {
    log_conc_valid <- log_conc
    viability_valid <- viability
  }

  viab_span <- max(viability_valid) - min(viability_valid)
  if (!is.finite(viab_span) || viab_span < .Machine$double.eps) {
    viab_span <- 1
  }
  viab_norm <- (viability_valid - min(viability_valid)) / viab_span

  ord <- order(viab_norm, log_conc_valid)
  viab_norm <- viab_norm[ord]
  log_conc_valid <- log_conc_valid[ord]

  unique_mask <- !duplicated(viab_norm)
  viab_norm_unique <- viab_norm[unique_mask]
  log_conc_unique <- log_conc_valid[unique_mask]

  ec50_guess <- stats::median(log_conc_valid)
  q25 <- q75 <- NA_real_
  if (length(viab_norm_unique) >= 2L) {
    ec50_guess_interp <- stats::approx(
      x = viab_norm_unique,
      y = log_conc_unique,
      xout = 0.5,
      rule = 2
    )$y
    if (is.finite(ec50_guess_interp)) {
      ec50_guess <- ec50_guess_interp
    }

    q25 <- stats::approx(
      x = viab_norm_unique,
      y = log_conc_unique,
      xout = 0.25,
      rule = 2
    )$y
    q75 <- stats::approx(
      x = viab_norm_unique,
      y = log_conc_unique,
      xout = 0.75,
      rule = 2
    )$y
  }
  if (is.finite(q25) && is.finite(q75) && abs(q75 - q25) > 1e-3) {
    hs_guess <- log(81) / (q75 - q25)
  } else {
    hs_guess <- 1
  }

  e0_guess <- max(viability_valid)
  einf_guess <- min(viability_valid)

  if (fit_type == "hill") {
    guess <- c(
      pmin(pmax(hs_guess, lower_bounds[1]), upper_bounds[1]),
      pmin(pmax(e0_guess, lower_bounds[2]), upper_bounds[2]),
      pmin(pmax(einf_guess, lower_bounds[3]), upper_bounds[3]),
      pmin(pmax(ec50_guess, lower_bounds[4]), upper_bounds[4])
    )
    names(guess) <- c("HS", "E0", "E_inf", "log10EC50")
  } else {
    guess <- c(
      pmin(pmax(hs_guess, lower_bounds[1]), upper_bounds[1]),
      pmin(pmax(e0_guess, lower_bounds[2]), upper_bounds[2]),
      pmin(pmax(einf_guess, lower_bounds[3]), upper_bounds[3]),
      pmin(pmax(1, lower_bounds[4]), upper_bounds[4]),
      pmin(
        pmax((min(viability_valid) + einf_guess) / 2, lower_bounds[5]),
        upper_bounds[5]
      ),
      pmin(
        pmax(stats::median(log_conc_valid) - 0.5, lower_bounds[6]),
        upper_bounds[6]
      ),
      pmin(
        pmax(stats::median(log_conc_valid) + 0.5, lower_bounds[7]),
        upper_bounds[7]
      ),
      pmin(pmax(0.5, lower_bounds[8]), upper_bounds[8])
    )
    names(guess) <- c(
      "HS1",
      "E0",
      "E_inf1",
      "HS2",
      "E_inf2",
      "log10EC50_1",
      "log10EC50_2",
      "Frac"
    )
  }
  guess
}

.pgx_fit_curve <- function(
  x,
  y,
  f,
  density,
  step,
  precision,
  lower_bounds,
  upper_bounds,
  scale,
  family,
  median_n,
  trunc,
  gritty_guess,
  span = 1,
  delta = 1
) {
  names(gritty_guess) <- names(gritty_guess)
  objective <- function(pars) {
    residual <- .pgx_curve_residual(
      x = x,
      y = y,
      n = median_n,
      pars = pars,
      f = f,
      scale = scale,
      family = family,
      trunc = trunc,
      delta = delta
    )
    if (identical(f, .pgx_biphasic_curve)) {
      weights <- rep(1, length(y))
      edge_idx <- c(1, 2, length(y) - 1, length(y))
      edge_idx <- edge_idx[edge_idx >= 1 & edge_idx <= length(y)]
      weights[edge_idx] <- 10
      sum(weights * residual)
    } else {
      sum(residual)
    }
  }

  opt <- try(
    stats::optim(
      par = gritty_guess,
      fn = objective,
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(
        factr = 1e-08,
        ndeps = rep(1e-4, length(gritty_guess)),
        trace = 0
      ),
      method = "L-BFGS-B"
    ),
    silent = TRUE
  )

  if (inherits(opt, "try-error")) {
    failed <- TRUE
    guess <- gritty_guess
  } else {
    failed <- FALSE
    guess <- opt$par
  }

  guess_residual <- sum(.pgx_curve_residual(
    x = x,
    y = y,
    n = median_n,
    pars = guess,
    f = f,
    scale = scale,
    family = family,
    trunc = trunc,
    delta = delta
  ))
  gritty_residual <- sum(.pgx_curve_residual(
    x = x,
    y = y,
    n = median_n,
    pars = gritty_guess,
    f = f,
    scale = scale,
    family = family,
    trunc = trunc,
    delta = delta
  ))

  if (failed || any(!is.finite(guess)) || guess_residual >= gritty_residual) {
    guess <- .pgx_mesh_eval(
      x = x,
      y = y,
      f = f,
      guess = gritty_guess,
      lower_bounds = lower_bounds,
      upper_bounds = upper_bounds,
      density = density,
      n = median_n,
      scale = scale,
      family = family,
      trunc = trunc
    )
    guess_residual <- sum(.pgx_curve_residual(
      x = x,
      y = y,
      n = median_n,
      pars = guess,
      f = f,
      scale = scale,
      family = family,
      trunc = trunc,
      delta = delta
    ))
    guess <- .pgx_pattern_search(
      x = x,
      y = y,
      f = f,
      guess = guess,
      n = median_n,
      guess_residual = guess_residual,
      lower_bounds = lower_bounds,
      upper_bounds = upper_bounds,
      span = span,
      precision = precision,
      step = step,
      scale = scale,
      family = family,
      trunc = trunc
    )
  }

  fitted_vals <- f(x, guess)
  rsq <- if (stats::var(y) == 0) {
    NA_real_
  } else {
    1 - (stats::var(y - fitted_vals) / stats::var(y))
  }

  names(guess) <- names(gritty_guess)
  attr(guess, "Rsquare") <- rsq
  list(pars = guess, Rsquare = rsq)
}

.pgx_curve_residual <- function(
  x,
  y,
  n,
  pars,
  f,
  scale,
  family,
  trunc,
  delta = 1
) {
  diffs <- f(x, pars) - y
  if (family == "Cauchy") {
    if (!trunc) {
      -log(CoreGx::.dmedncauchys(diffs, n, scale))
    } else {
      down_truncated <- abs(y) >= 1
      up_truncated <- abs(y) <= 0
      c(
        -log(CoreGx::.dmedncauchys(
          diffs[!(down_truncated | up_truncated)],
          n,
          scale
        )),
        -log(CoreGx::.edmedncauchys(
          -diffs[up_truncated | down_truncated],
          n,
          scale
        ))
      )
    }
  } else {
    huber_loss <- ifelse(
      abs(diffs) <= delta,
      0.5 * diffs^2,
      delta * (abs(diffs) - 0.5 * delta)
    )
    if (!trunc) {
      huber_loss
    } else {
      down_truncated <- abs(y) >= 1
      up_truncated <- abs(y) <= 0
      c(
        huber_loss[!(down_truncated | up_truncated)],
        delta * (abs(-diffs[up_truncated | down_truncated]) - 0.5 * delta)
      )
    }
  }
}

.pgx_mesh_eval <- function(
  x,
  y,
  f,
  guess,
  lower_bounds,
  upper_bounds,
  density,
  n,
  scale,
  family,
  trunc
) {
  guess_residual <- sum(.pgx_curve_residual(
    x = x,
    y = y,
    n = n,
    pars = guess,
    f = f,
    scale = scale,
    family = family,
    trunc = trunc
  ))

  periods <- rep(1, length(guess))
  if (length(guess) > 1) {
    for (idx in 2:length(guess)) {
      periods[idx] <- periods[idx - 1] *
        (density[idx - 1] * (upper_bounds[idx - 1] - lower_bounds[idx - 1]) + 1)
    }
  }

  current <- lower_bounds
  total_points <- prod((upper_bounds - lower_bounds) * density + 1)

  for (point in seq_len(total_points)) {
    test_residual <- sum(.pgx_curve_residual(
      x = x,
      y = y,
      n = n,
      pars = current,
      f = f,
      scale = scale,
      family = family,
      trunc = trunc
    ))
    if (is.finite(test_residual) && test_residual < guess_residual) {
      guess <- current
      guess_residual <- test_residual
    }
    for (idx in seq_along(guess)) {
      if (point %% periods[idx] == 0) {
        current[idx] <- current[idx] + 1 / density[idx]
        if (current[idx] > upper_bounds[idx]) {
          current[idx] <- lower_bounds[idx]
        }
      }
    }
  }
  guess
}

.pgx_pattern_search <- function(
  x,
  y,
  f,
  guess,
  n,
  guess_residual,
  lower_bounds,
  upper_bounds,
  span,
  precision,
  step,
  scale,
  family,
  trunc
) {
  neighbours <- matrix(nrow = 2 * length(guess), ncol = length(guess))
  neighbour_residuals <- rep(NA_real_, nrow(neighbours))

  while (span > precision) {
    for (idx in seq_len(nrow(neighbours))) {
      neighbours[idx, ] <- guess
      dimension <- ceiling(idx / 2)
      if (idx %% 2 == 1) {
        neighbours[idx, dimension] <- min(
          guess[dimension] + span * step[dimension],
          upper_bounds[dimension]
        )
      } else {
        neighbours[idx, dimension] <- max(
          guess[dimension] - span * step[dimension],
          lower_bounds[dimension]
        )
      }
      neighbour_residuals[idx] <- sum(.pgx_curve_residual(
        x = x,
        y = y,
        f = f,
        pars = neighbours[idx, ],
        n = n,
        scale = scale,
        family = family,
        trunc = trunc
      ))
    }

    min_residual <- min(neighbour_residuals, na.rm = TRUE)
    if (min_residual < guess_residual) {
      guess <- neighbours[which.min(neighbour_residuals), ]
      guess_residual <- min_residual
    } else {
      span <- span / 2
    }
  }
  guess
}

.pgx_hill_curve <- function(x, pars) {
  hs <- pars[1]
  e0 <- pars[2]
  einf <- pars[3]
  log_ec50 <- pars[4]
  einf + (e0 - einf) / (1 + (10^x / 10^log_ec50)^hs)
}

.pgx_biphasic_curve <- function(x, pars) {
  hs1 <- pars[1]
  e0 <- pars[2]
  einf1 <- pars[3]
  hs2 <- pars[4]
  einf2 <- pars[5]
  log_ec50_1 <- pars[6]
  log_ec50_2 <- pars[7]
  frac <- pars[8]

  span <- max(einf1, einf2) - e0
  section1 <- span * frac / (1 + 10^((log_ec50_1 - x) * hs1))
  section2 <- span * (1 - frac) / (1 + 10^((log_ec50_2 - x) * hs2))
  e0 + section1 + section2
}
