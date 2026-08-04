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
#' @param curve_direction Character string controlling asymptote ordering.
#'   `"unconstrained"` preserves the existing behavior. `"decreasing"`
#'   requires fitted viability to be non-increasing with dose by constraining
#'   every fitted high-dose asymptote to be no greater than `E0`.
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
  fit_type = c("hill", "biphasic"),
  curve_direction = c("unconstrained", "decreasing")
) {
  fit_type <- match.arg(tolower(fit_type), c("hill", "biphasic"))
  curve_direction <- match.arg(
    tolower(curve_direction),
    c("unconstrained", "decreasing")
  )
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

  starting_guesses <- .pgx_starting_guesses(
    log_conc = log_conc,
    viability = viability_clean,
    initial_guess = guess,
    lower_bounds = bounds$lower,
    upper_bounds = bounds$upper,
    fit_type = fit_type
  )

  parameterization <- .pgx_parameterization(
    lower_bounds = bounds$lower,
    upper_bounds = bounds$upper,
    fit_type = fit_type,
    curve_direction = curve_direction
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
    starting_guesses = starting_guesses
  )

  if (curve_direction == "decreasing") {
    asymptote_indices <- if (fit_type == "hill") 3L else c(3L, 5L)
    unconstrained_is_decreasing <- all(
      fitted$pars[asymptote_indices] <= fitted$pars[2L]
    )

    if (!unconstrained_is_decreasing) {
      starting_guesses <- t(apply(
        starting_guesses,
        1,
        parameterization$to_optimizer
      ))
      colnames(starting_guesses) <- names(guess)

      fitted <- .pgx_fit_curve(
        x = log_conc,
        y = viability_clean,
        f = curve_fun,
        density = bounds$density,
        step = bounds$step,
        precision = precision,
        lower_bounds = parameterization$lower,
        upper_bounds = parameterization$upper,
        scale = scale,
        family = family,
        median_n = median_n,
        trunc = trunc,
        starting_guesses = starting_guesses,
        to_native = parameterization$to_native
      )
    }
  }

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
      0.1,
      0,
      min(log_conc) - 3,
      min(log_conc) - 3,
      0
    )
    default_upper <- c(
      5,
      1.5,
      1.5,
      5,
      1.5,
      max(log_conc) + 3,
      max(log_conc) + 3,
      1
    )
    default_density <- c(2, 10, 10, 5, 10, 5, 5, 5)
  }

  default_step <- 0.5 / default_density

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
    invalid_density <- !is.finite(density)
    if (any(invalid_density)) {
      density[invalid_density] <- default_density[invalid_density]
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
    invalid_step <- !is.finite(step)
    if (any(invalid_step)) {
      step[invalid_step] <- default_step[invalid_step]
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
    invalid_lower <- !is.finite(lower_bounds)
    if (any(invalid_lower)) {
      lower_bounds[invalid_lower] <- default_lower[invalid_lower]
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
    invalid_upper <- !is.finite(upper_bounds)
    if (any(invalid_upper)) {
      upper_bounds[invalid_upper] <- default_upper[invalid_upper]
    }
  }

  if (!all(is.finite(upper_bounds)) || !all(is.finite(lower_bounds))) {
    stop("`lower_bounds` and `upper_bounds` must contain finite values.")
  }

  if (any((upper_bounds - lower_bounds) <= 0)) {
    stop("Upper bounds must exceed lower bounds for all parameters.")
  }

  if (!all(is.finite(density)) || !all(is.finite(step))) {
    stop("`density` and `step` must contain finite positive values.")
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

  if (!length(log_conc_valid) || !length(viability_valid)) {
    stop(
      "Unable to compute initial guesses: all concentration or viability values are non-finite. ",
      "Please verify input data."
    )
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

.pgx_starting_guesses <- function(
  log_conc,
  viability,
  initial_guess,
  lower_bounds,
  upper_bounds,
  fit_type
) {
  dose_quantiles <- as.numeric(stats::quantile(
    log_conc,
    probs = c(0.25, 0.5, 0.75),
    names = FALSE
  ))

  if (fit_type == "hill") {
    canonical_guesses <- cbind(
      HS = 1,
      E0 = max(viability),
      E_inf = min(viability),
      log10EC50 = dose_quantiles
    )
  } else {
    ec50_pairs <- rbind(
      c(dose_quantiles[1], dose_quantiles[2]),
      c(dose_quantiles[1], dose_quantiles[3]),
      c(dose_quantiles[2], dose_quantiles[3])
    )
    canonical_guesses <- cbind(
      HS1 = 1,
      E0 = max(viability),
      E_inf1 = min(viability),
      HS2 = 1,
      E_inf2 = min(viability),
      log10EC50_1 = ec50_pairs[, 1],
      log10EC50_2 = ec50_pairs[, 2],
      Frac = 0.5
    )
  }

  guesses <- rbind(initial_guess, canonical_guesses)
  guesses <- sweep(guesses, 2, lower_bounds, pmax)
  guesses <- sweep(guesses, 2, upper_bounds, pmin)
  guesses <- unique(guesses)
  colnames(guesses) <- names(initial_guess)
  guesses
}

.pgx_parameterization <- function(
  lower_bounds,
  upper_bounds,
  fit_type,
  curve_direction
) {
  identity_transform <- function(parameters) {
    parameters
  }
  if (curve_direction == "unconstrained") {
    return(list(
      lower = lower_bounds,
      upper = upper_bounds,
      to_optimizer = identity_transform,
      to_native = identity_transform
    ))
  }

  e0_index <- 2L
  asymptote_indices <- if (fit_type == "hill") 3L else c(3L, 5L)
  required_e0 <- max(lower_bounds[asymptote_indices])
  if (upper_bounds[e0_index] < required_e0) {
    stop(
      "Decreasing curves require the E0 upper bound to be at least every high-dose asymptote lower bound."
    )
  }

  optimizer_lower <- lower_bounds
  optimizer_upper <- upper_bounds
  optimizer_lower[e0_index] <- max(
    optimizer_lower[e0_index],
    required_e0
  )
  optimizer_lower[asymptote_indices] <- 0
  optimizer_upper[asymptote_indices] <- 1

  to_native <- function(parameters) {
    native <- parameters
    e0 <- native[e0_index]
    for (index in asymptote_indices) {
      asymptote_upper <- min(upper_bounds[index], e0)
      native[index] <- lower_bounds[index] +
        parameters[index] * (asymptote_upper - lower_bounds[index])
    }
    names(native) <- names(parameters)
    native
  }

  to_optimizer <- function(parameters) {
    optimizer <- parameters
    optimizer[e0_index] <- pmin(
      pmax(optimizer[e0_index], optimizer_lower[e0_index]),
      optimizer_upper[e0_index]
    )
    e0 <- optimizer[e0_index]
    for (index in asymptote_indices) {
      asymptote_upper <- min(upper_bounds[index], e0)
      asymptote_range <- asymptote_upper - lower_bounds[index]
      if (asymptote_range <= .Machine$double.eps) {
        optimizer[index] <- 0
      } else {
        native_asymptote <- pmin(
          pmax(parameters[index], lower_bounds[index]),
          asymptote_upper
        )
        optimizer[index] <-
          (native_asymptote - lower_bounds[index]) / asymptote_range
      }
    }
    names(optimizer) <- names(parameters)
    optimizer
  }

  list(
    lower = optimizer_lower,
    upper = optimizer_upper,
    to_optimizer = to_optimizer,
    to_native = to_native
  )
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
  starting_guesses,
  to_native = identity,
  span = 1,
  delta = 1
) {
  starting_guesses <- as.matrix(starting_guesses)
  objective <- function(pars) {
    native_pars <- to_native(pars)
    residual <- .pgx_curve_residual(
      x = x,
      y = y,
      n = median_n,
      pars = native_pars,
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

  raw_residuals <- apply(starting_guesses, 1, objective)
  canonical_indices <- seq_len(nrow(starting_guesses))[-1]
  canonical_order <- canonical_indices[order(
    raw_residuals[canonical_indices],
    canonical_indices
  )]
  selected_indices <- unique(c(1L, utils::head(canonical_order, 2L)))
  optimization_scales <- abs(raw_residuals[selected_indices])
  optimization_scales <- optimization_scales[
    is.finite(optimization_scales) & optimization_scales > .Machine$double.eps
  ]
  optimization_scale <- if (length(optimization_scales)) {
    min(optimization_scales)
  } else {
    1
  }

  optimized <- lapply(selected_indices, function(index) {
    result <- try(
      stats::optim(
        par = starting_guesses[index, ],
        fn = objective,
        lower = lower_bounds,
        upper = upper_bounds,
        control = list(
          factr = 1e-08,
          fnscale = optimization_scale,
          maxit = 500L,
          ndeps = rep(1e-4, ncol(starting_guesses)),
          trace = 0
        ),
        method = "L-BFGS-B"
      ),
      silent = TRUE
    )
    if (inherits(result, "try-error") || any(!is.finite(result$par))) {
      return(NULL)
    }
    result$residual <- objective(result$par)
    if (!is.finite(result$residual)) {
      return(NULL)
    }
    result$start_index <- index
    result
  })
  optimized <- Filter(Negate(is.null), optimized)

  raw_best <- min(raw_residuals)
  raw_tolerance <- sqrt(.Machine$double.eps) * max(1, abs(raw_best))
  raw_best_index <- which(raw_residuals <= raw_best + raw_tolerance)[1]
  raw_guess <- starting_guesses[raw_best_index, ]

  if (length(optimized)) {
    optimized_residuals <- vapply(
      optimized,
      function(result) result$residual,
      numeric(1)
    )
    optimized_best <- min(optimized_residuals)
    optimized_tolerance <- sqrt(.Machine$double.eps) *
      max(1, abs(optimized_best))
    optimized_best_index <- which(
      optimized_residuals <= optimized_best + optimized_tolerance
    )[1]
    guess <- optimized[[optimized_best_index]]$par
    guess_residual <- optimized_residuals[optimized_best_index]
  } else {
    guess <- raw_guess
    guess_residual <- raw_best
  }

  use_fallback <- !length(optimized) ||
    !is.finite(guess_residual) ||
    guess_residual >= raw_best - raw_tolerance

  if (use_fallback) {
    density_vec <- rep_len(density, length(lower_bounds))
    grid_counts <- ceiling(pmax(
      1,
      (upper_bounds - lower_bounds) * density_vec + 1
    ))
    grid_size <- prod(grid_counts)
    max_grid_points <- getOption("PharmacoGx.max_mesh_points", 1e5)

    if (grid_size <= max_grid_points) {
      guess <- .pgx_mesh_eval(
        x = x,
        y = y,
        f = f,
        guess = raw_guess,
        lower_bounds = lower_bounds,
        upper_bounds = upper_bounds,
        density = density,
        n = median_n,
        scale = scale,
        family = family,
        trunc = trunc,
        to_native = to_native
      )
      guess_residual <- objective(guess)
    } else {
      warning(sprintf(
        "Skipping full mesh evaluation: grid size (%s) exceeds max_grid_points (%s). Adjust density or option 'PharmacoGx.max_mesh_points' to enable mesh search.",
        format(grid_size, scientific = TRUE),
        format(max_grid_points, scientific = TRUE)
      ))
      guess <- raw_guess
      guess_residual <- raw_best
    }

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
      trunc = trunc,
      to_native = to_native
    )
  }

  guess <- to_native(guess)
  fitted_vals <- f(x, guess)
  rsq <- if (stats::var(y) == 0) {
    NA_real_
  } else {
    1 - (stats::var(y - fitted_vals) / stats::var(y))
  }

  names(guess) <- colnames(starting_guesses)
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
      truncated_idx <- down_truncated | up_truncated
      residual <- numeric(length(y))
      if (any(!truncated_idx)) {
        residual[!truncated_idx] <- -log(CoreGx::.dmedncauchys(
          diffs[!truncated_idx],
          n,
          scale
        ))
      }
      if (any(truncated_idx)) {
        residual[truncated_idx] <- -log(CoreGx::.edmedncauchys(
          -diffs[truncated_idx],
          n,
          scale
        ))
      }
      residual
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
      truncated_idx <- down_truncated | up_truncated
      residual <- numeric(length(y))
      if (any(!truncated_idx)) {
        residual[!truncated_idx] <- huber_loss[!truncated_idx]
      }
      if (any(truncated_idx)) {
        residual[truncated_idx] <- huber_loss[truncated_idx]
      }
      residual
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
  trunc,
  to_native = identity
) {
  guess_residual <- sum(.pgx_curve_residual(
    x = x,
    y = y,
    n = n,
    pars = to_native(guess),
    f = f,
    scale = scale,
    family = family,
    trunc = trunc
  ))

  density_vec <- rep_len(density, length(guess))
  ranges <- upper_bounds - lower_bounds
  grid_counts <- pmax(1L, as.integer(round(density_vec * ranges)) + 1L)
  step_sizes <- numeric(length(guess))
  multi_step <- grid_counts > 1L
  step_sizes[multi_step] <- ranges[multi_step] / (grid_counts[multi_step] - 1L)

  total_points <- prod(grid_counts)
  indices <- integer(length(guess))
  current <- lower_bounds + indices * step_sizes

  for (point in seq_len(total_points)) {
    test_residual <- sum(.pgx_curve_residual(
      x = x,
      y = y,
      n = n,
      pars = to_native(current),
      f = f,
      scale = scale,
      family = family,
      trunc = trunc
    ))
    if (is.finite(test_residual) && test_residual < guess_residual) {
      guess <- current
      guess_residual <- test_residual
    }
    if (point < total_points) {
      for (idx in seq_along(indices)) {
        indices[idx] <- indices[idx] + 1L
        if (indices[idx] < grid_counts[idx]) {
          current[idx] <- lower_bounds[idx] + indices[idx] * step_sizes[idx]
          break
        }
        indices[idx] <- 0L
        current[idx] <- lower_bounds[idx]
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
  trunc,
  to_native = identity
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
        pars = to_native(neighbours[idx, ]),
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

#' Evaluate a Hill curve at log10 concentrations
#'
#' Computes viability predictions for a four-parameter Hill model given log10
#' concentration inputs and parameter vector.
#'
#' @param x `numeric` vector of log10-transformed concentrations.
#' @param pars `numeric` vector containing `HS`, `E0`, `E_inf`, and `log10EC50`
#'   in that order.
#'
#' @return `numeric` vector of predicted viabilities on the 0–1 scale.
#'
#' @keywords internal
.pgx_hill_curve <- function(x, pars) {
  hs <- pars[1]
  e0 <- pars[2]
  einf <- pars[3]
  log_ec50 <- pars[4]
  einf + (e0 - einf) / (1 + (10^x / 10^log_ec50)^hs)
}

#' Evaluate a biphasic Hill curve at log10 concentrations
#'
#' Produces viability predictions for a two-component (biphasic) Hill model
#' using log10 dose inputs and the combined parameter vector.
#'
#' @param x `numeric` vector of log10-transformed concentrations.
#' @param pars `numeric` vector containing biphasic parameters in the order
#'   `HS1`, `E0`, `E_inf1`, `HS2`, `E_inf2`, `log10EC50_1`, `log10EC50_2`,
#'   `Frac`.
#'
#' @return `numeric` vector of predicted viabilities on the 0–1 scale.
#'
#' @keywords internal
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
