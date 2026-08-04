library(PharmacoGx)

context("Testing LogLogisticRegression.")

##TO-DO::Supress print to console from this test file

test_that("Errors are checked.", {
  expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60))) #should complain
  expect_warning({
    res <- logLogisticRegression(
      c(1, 2, 3),
      c(70, 60, 50),
      viability_as_pct = FALSE
    )
    expect_named(res, c("HS", "E0", "E_inf", "EC50"))
  })
  expect_error(logLogisticRegression(
    c(-1, 2, 3),
    c(70, 60, 50),
    conc_as_log = FALSE
  )) #should complain

  expect_error(logLogisticRegression(c(1, 2, 3), c(70, 60, 50), median_n = 0)) #should complain
  expect_error(logLogisticRegression(
    c(1, 2, 3),
    c(50, 60, 70),
    median_n = 3 / 2
  )) #should complain
  expect_error(logLogisticRegression(
    c(1, 2, 3),
    c(50, 60, 70),
    density = c(1, 1)
  )) #should complain
  expect_error(logLogisticRegression(
    c(1, 2, 3),
    c(50, 60, 70),
    density = c(1, 1, -1)
  )) #should complain
  expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), precision = 0)) #should complain
  expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), scale = 0)) #should complain
  expect_error(logLogisticRegression(
    c(1, 2, 3),
    c(50, 60, 70),
    lower_bounds = c(0, 0, 0),
    upper_bounds = c(1, 1, -1)
  )) #should complain
  expect_error(logLogisticRegression(
    c(1, 2, 3),
    c(50, 60, 70),
    family = "The Addams Family"
  )) #should complain
  expect_error(logLogisticRegression(c(1, 2), c(70, 60)))
  expect_error(logLogisticRegression(
    c(1, 2, 3),
    c(70, 60, 50),
    curve_direction = "sideways"
  ))
})

test_that("Hill and biphasic fits recover known parameters", {
  hill_fun <- PharmacoGx:::.pgx_hill_curve
  log_conc <- seq(-3, 3, length.out = 25)
  hill_pars <- c(HS = 1.4, E0 = 1, E_inf = 0.2, log10EC50 = -0.2)
  conc <- 10^log_conc
  viability <- hill_fun(log_conc, hill_pars)

  hill_fit <- logLogisticRegression(
    conc = conc,
    viability = viability,
    conc_as_log = FALSE,
    viability_as_pct = FALSE
  )

  expect_named(hill_fit, c("HS", "E0", "E_inf", "EC50"))
  expect_equal(hill_fit$HS, hill_pars[["HS"]], tolerance = 0.1)
  expect_equal(hill_fit$E0, hill_pars[["E0"]], tolerance = 0.05)
  expect_equal(hill_fit$E_inf, hill_pars[["E_inf"]], tolerance = 0.05)
  expect_equal(log10(hill_fit$EC50), hill_pars[["log10EC50"]], tolerance = 0.1)

  viability_pct <- viability * 100
  hill_fit_pct <- logLogisticRegression(
    conc = conc,
    viability = viability_pct,
    conc_as_log = FALSE,
    viability_as_pct = TRUE
  )
  expect_equal(hill_fit_pct$E0, hill_pars[["E0"]] * 100, tolerance = 5)
  expect_equal(hill_fit_pct$E_inf, hill_pars[["E_inf"]] * 100, tolerance = 5)

  hill_fit_cauchy <- logLogisticRegression(
    conc = conc,
    viability = viability,
    conc_as_log = FALSE,
    viability_as_pct = FALSE,
    family = "Cauchy"
  )
  expect_named(hill_fit_cauchy, c("HS", "E0", "E_inf", "EC50"))

  # check truncation bounds values above 100%
  inflated <- viability_pct
  inflated[1] <- 150
  fit_trunc <- logLogisticRegression(
    conc = conc,
    viability = inflated,
    conc_as_log = FALSE,
    viability_as_pct = TRUE,
    trunc = TRUE
  )
  expect_lte(fit_trunc$E0, 100 + 1e-6)

  biphasic_fun <- PharmacoGx:::.pgx_biphasic_curve
  biphasic_pars <- c(
    HS1 = 1.3,
    E0 = 1,
    E_inf1 = 0.3,
    HS2 = 0.8,
    E_inf2 = 0.1,
    log10EC50_1 = -0.5,
    log10EC50_2 = 0.5,
    Frac = 0.6
  )
  biphasic_viab <- biphasic_fun(log_conc, biphasic_pars)
  biphasic_fit <- logLogisticRegression(
    conc = conc,
    viability = biphasic_viab,
    conc_as_log = FALSE,
    viability_as_pct = FALSE,
    fit_type = "biphasic"
  )
  expect_named(
    biphasic_fit,
    c("HS1", "E0", "E_inf1", "HS2", "E_inf2", "EC50_1", "EC50_2", "Frac")
  )
  expect_gt(attr(biphasic_fit, "Rsquare"), 0.99)
  expect_equal(
    unname(sort(log10(c(biphasic_fit$EC50_1, biphasic_fit$EC50_2)))),
    unname(sort(biphasic_pars[c("log10EC50_1", "log10EC50_2")])),
    tolerance = 0.35
  )
})

test_that("biphasic fitting respects parameter-specific bounds", {
  log_conc <- seq(-2, 2, length.out = 25)
  true_pars <- c(
    HS1 = 1.2,
    E0 = 1,
    E_inf1 = 0.35,
    HS2 = 0.7,
    E_inf2 = 0.15,
    log10EC50_1 = -0.6,
    log10EC50_2 = 0.7,
    Frac = 0.55
  )
  viability <- PharmacoGx:::.pgx_biphasic_curve(log_conc, true_pars)
  lower <- c(0.4, 0.8, 0, 0.3, 0, -1, 0.2, 0.25)
  upper <- c(2, 1.1, 0.6, 1.5, 0.6, -0.2, 1.2, 0.8)

  fit <- logLogisticRegression(
    conc = 10^log_conc,
    viability = viability,
    viability_as_pct = FALSE,
    fit_type = "biphasic",
    lower_bounds = lower,
    upper_bounds = upper
  )

  fitted_native <- c(
    fit$HS1,
    fit$E0,
    fit$E_inf1,
    fit$HS2,
    fit$E_inf2,
    log10(fit$EC50_1),
    log10(fit$EC50_2),
    fit$Frac
  )
  expect_true(all(fitted_native >= lower - 1e-8))
  expect_true(all(fitted_native <= upper + 1e-8))
  expect_gt(attr(fit, "Rsquare"), 0.99)
})

test_that("multi-start fitting avoids divergent CTRPv2 minima", {
  fixtures <- list(
    list(
      conc = c(
        0.002,
        0.0041,
        0.0081,
        0.016,
        0.032,
        0.065,
        0.13,
        0.26,
        0.52,
        1,
        2.1,
        4.2,
        8.3,
        17,
        33,
        66
      ),
      viability = c(
        126.2,
        125.8,
        137.1,
        105.5,
        113.8,
        125.3,
        125.6,
        109.7,
        93.52,
        73.86,
        118.3,
        111.1,
        94.28,
        74.31,
        9.685,
        2.315
      ),
      minimum_r_squared = 0.94
    ),
    list(
      conc = c(
        0.009,
        0.018,
        0.036,
        0.072,
        0.14,
        0.29,
        0.58,
        1.2,
        2.3,
        4.6,
        9.2,
        18,
        37,
        74,
        150,
        300
      ),
      viability = c(
        108.2,
        114,
        110.5,
        113.2,
        107.6,
        103.6,
        97.74,
        114.6,
        119,
        112.9,
        114.3,
        126.1,
        64.41,
        106.3,
        75.89,
        3.613
      ),
      minimum_r_squared = 0.85
    ),
    list(
      conc = c(
        0.0046,
        0.0091,
        0.018,
        0.036,
        0.073,
        0.15,
        0.29,
        0.58,
        1.2,
        2.3,
        4.7,
        9.3,
        19,
        37,
        75,
        150
      ),
      viability = c(
        108.2,
        102.9,
        103,
        101.4,
        103.2,
        102.5,
        101.9,
        102.4,
        80.48,
        94.81,
        61.18,
        98.23,
        72.42,
        83.14,
        77.86,
        2.857
      ),
      minimum_r_squared = 0.74
    ),
    list(
      conc = c(
        0.001,
        0.002,
        0.0041,
        0.0081,
        0.016,
        0.032,
        0.065,
        0.13,
        0.26,
        0.52,
        1,
        2.1,
        4.2,
        8.3,
        17,
        33
      ),
      viability = c(
        111.9,
        115,
        117.2,
        112.8,
        110.1,
        42.16,
        116,
        114.4,
        113.1,
        113.9,
        104.9,
        99.83,
        81.74,
        95.05,
        57.8,
        3.198
      ),
      minimum_r_squared = 0.70
    ),
    list(
      conc = c(
        0.001,
        0.002,
        0.0041,
        0.0081,
        0.016,
        0.032,
        0.065,
        0.13,
        0.26,
        0.52,
        1,
        2.1,
        4.2,
        8.3,
        17,
        33
      ),
      viability = c(
        105.8,
        108,
        102.3,
        101.7,
        103.2,
        106.4,
        102.6,
        104.7,
        104.4,
        104.5,
        101.9,
        101.9,
        105.3,
        108.1,
        99.97,
        100.7
      ),
      minimum_r_squared = 0.40
    ),
    list(
      conc = c(
        0.002,
        0.0041,
        0.0081,
        0.016,
        0.032,
        0.065,
        0.13,
        0.26,
        0.52,
        1,
        2.1,
        4.2,
        8.3,
        17,
        33,
        66
      ),
      viability = c(
        101.1,
        94.53,
        97.12,
        97.52,
        93.17,
        92.48,
        98.43,
        98.93,
        100.5,
        97.02,
        94.49,
        96.09,
        95.6,
        106.3,
        95.04,
        96.56
      ),
      minimum_r_squared = 0.12
    )
  )

  for (fixture in fixtures) {
    fit <- logLogisticRegression(
      fixture$conc,
      fixture$viability,
      upper_bounds = c(4.5, NA, NA, NA)
    )
    perturbed_fit <- logLogisticRegression(
      fixture$conc,
      fixture$viability +
        rep(c(-1, 1), length.out = length(fixture$viability)) *
          1e-12,
      upper_bounds = c(4.5, NA, NA, NA)
    )

    expect_gte(attr(fit, "Rsquare"), fixture$minimum_r_squared)
    expect_gte(attr(perturbed_fit, "Rsquare"), fixture$minimum_r_squared)
    expect_equal(
      unlist(fit),
      unlist(perturbed_fit),
      tolerance = 1e-4
    )
  }
})

test_that("unconstrained remains the default curve direction", {
  conc <- 10^seq(-2, 2, length.out = 16)
  viability <- seq(95, 15, length.out = 16)

  default_fit <- logLogisticRegression(conc, viability)
  explicit_fit <- logLogisticRegression(
    conc,
    viability,
    curve_direction = "unconstrained"
  )

  expect_equal(default_fit, explicit_fit)
  expect_equal(attr(default_fit, "Rsquare"), attr(explicit_fit, "Rsquare"))
})

test_that("decreasing Hill fits are non-increasing across data modes", {
  log_conc <- seq(-3, 3, length.out = 25)
  native_pars <- c(HS = 1.4, E0 = 0.95, E_inf = 0.15, log10EC50 = -0.2)
  fractional_viability <- PharmacoGx:::.pgx_hill_curve(log_conc, native_pars)

  for (family in c("normal", "Cauchy")) {
    for (viability_as_pct in c(FALSE, TRUE)) {
      viability <- if (viability_as_pct) {
        fractional_viability * 100
      } else {
        fractional_viability
      }
      fit <- logLogisticRegression(
        conc = 10^log_conc,
        viability = viability,
        viability_as_pct = viability_as_pct,
        family = family,
        curve_direction = "decreasing"
      )
      scale_factor <- if (viability_as_pct) 100 else 1
      fitted_native <- c(
        fit$HS,
        fit$E0 / scale_factor,
        fit$E_inf / scale_factor,
        log10(fit$EC50)
      )
      dense_predictions <- PharmacoGx:::.pgx_hill_curve(
        seq(-4, 4, length.out = 501),
        fitted_native
      )

      expect_named(fit, c("HS", "E0", "E_inf", "EC50"))
      expect_lte(fit$E_inf, fit$E0 + 1e-8)
      expect_true(all(diff(dense_predictions) <= 1e-10))
      expect_gt(attr(fit, "Rsquare"), 0.99)
    }
  }
})

test_that("decreasing biphasic fits constrain both asymptotes", {
  log_conc <- seq(-3, 3, length.out = 31)
  native_pars <- c(
    HS1 = 1.3,
    E0 = 1,
    E_inf1 = 0.25,
    HS2 = 0.8,
    E_inf2 = 0.1,
    log10EC50_1 = -0.7,
    log10EC50_2 = 0.6,
    Frac = 0.6
  )
  viability <- PharmacoGx:::.pgx_biphasic_curve(log_conc, native_pars)

  for (family in c("normal", "Cauchy")) {
    for (viability_as_pct in c(FALSE, TRUE)) {
      observed_viability <- if (viability_as_pct) viability * 100 else viability
      fit <- logLogisticRegression(
        conc = 10^log_conc,
        viability = observed_viability,
        viability_as_pct = viability_as_pct,
        family = family,
        fit_type = "biphasic",
        curve_direction = "decreasing"
      )
      scale_factor <- if (viability_as_pct) 100 else 1
      fitted_native <- c(
        fit$HS1,
        fit$E0 / scale_factor,
        fit$E_inf1 / scale_factor,
        fit$HS2,
        fit$E_inf2 / scale_factor,
        log10(fit$EC50_1),
        log10(fit$EC50_2),
        fit$Frac
      )
      dense_predictions <- PharmacoGx:::.pgx_biphasic_curve(
        seq(-4, 4, length.out = 501),
        fitted_native
      )

      expect_named(
        fit,
        c("HS1", "E0", "E_inf1", "HS2", "E_inf2", "EC50_1", "EC50_2", "Frac")
      )
      expect_lte(fit$E_inf1, fit$E0 + 1e-8)
      expect_lte(fit$E_inf2, fit$E0 + 1e-8)
      expect_true(all(diff(dense_predictions) <= 1e-10))
      expect_gt(attr(fit, "Rsquare"), 0.99)
    }
  }
})

test_that("decreasing fits preserve feasible native custom bounds", {
  log_conc <- seq(-2, 2, length.out = 25)
  viability <- PharmacoGx:::.pgx_hill_curve(
    log_conc,
    c(HS = 1.1, E0 = 0.9, E_inf = 0.25, log10EC50 = 0.1)
  )
  lower <- c(0.5, 0.7, 0.1, -1)
  upper <- c(2, 1, 0.5, 1)
  fit <- logLogisticRegression(
    conc = 10^log_conc,
    viability = viability,
    viability_as_pct = FALSE,
    lower_bounds = lower,
    upper_bounds = upper,
    curve_direction = "decreasing"
  )
  fitted_native <- c(fit$HS, fit$E0, fit$E_inf, log10(fit$EC50))

  expect_true(all(fitted_native >= lower - 1e-8))
  expect_true(all(fitted_native <= upper + 1e-8))
  expect_lte(fit$E_inf, fit$E0 + 1e-8)

  biphasic_pars <- c(
    HS1 = 1.2,
    E0 = 0.95,
    E_inf1 = 0.3,
    HS2 = 0.7,
    E_inf2 = 0.15,
    log10EC50_1 = -0.5,
    log10EC50_2 = 0.6,
    Frac = 0.55
  )
  biphasic_viability <- PharmacoGx:::.pgx_biphasic_curve(
    log_conc,
    biphasic_pars
  )
  biphasic_lower <- c(0.4, 0.8, 0.1, 0.3, 0.05, -1, 0.1, 0.25)
  biphasic_upper <- c(2, 1.1, 0.6, 1.5, 0.5, 0, 1.2, 0.8)
  biphasic_fit <- logLogisticRegression(
    conc = 10^log_conc,
    viability = biphasic_viability,
    viability_as_pct = FALSE,
    fit_type = "biphasic",
    lower_bounds = biphasic_lower,
    upper_bounds = biphasic_upper,
    curve_direction = "decreasing"
  )
  biphasic_fitted_native <- c(
    biphasic_fit$HS1,
    biphasic_fit$E0,
    biphasic_fit$E_inf1,
    biphasic_fit$HS2,
    biphasic_fit$E_inf2,
    log10(biphasic_fit$EC50_1),
    log10(biphasic_fit$EC50_2),
    biphasic_fit$Frac
  )
  expect_true(all(biphasic_fitted_native >= biphasic_lower - 1e-8))
  expect_true(all(biphasic_fitted_native <= biphasic_upper + 1e-8))
  expect_lte(biphasic_fit$E_inf1, biphasic_fit$E0 + 1e-8)
  expect_lte(biphasic_fit$E_inf2, biphasic_fit$E0 + 1e-8)

  expect_error(
    logLogisticRegression(
      conc = 10^log_conc,
      viability = viability,
      viability_as_pct = FALSE,
      lower_bounds = c(0.5, 0.1, 0.6, -1),
      upper_bounds = c(2, 0.5, 1, 1),
      curve_direction = "decreasing"
    ),
    "E0 upper bound"
  )
  expect_error(
    logLogisticRegression(
      conc = 10^log_conc,
      viability = viability,
      viability_as_pct = FALSE,
      fit_type = "biphasic",
      lower_bounds = c(0.2, 0.1, 0.6, 0.2, 0.7, -2, -2, 0),
      upper_bounds = c(2, 0.5, 1, 2, 1, 2, 2, 1),
      curve_direction = "decreasing"
    ),
    "E0 upper bound"
  )
})

test_that("decreasing mode handles flat and increasing profiles", {
  conc <- 10^seq(-2, 2, length.out = 21)
  flat_fit <- logLogisticRegression(
    conc,
    rep(50, length(conc)),
    curve_direction = "decreasing"
  )
  expect_lte(flat_fit$E_inf, flat_fit$E0 + 1e-8)
  expect_true(is.na(attr(flat_fit, "Rsquare")))

  increasing_viability <- seq(15, 90, length.out = length(conc))
  unconstrained_fit <- logLogisticRegression(conc, increasing_viability)
  decreasing_fit <- logLogisticRegression(
    conc,
    increasing_viability,
    curve_direction = "decreasing"
  )

  expect_gt(unconstrained_fit$E_inf, unconstrained_fit$E0)
  expect_lte(decreasing_fit$E_inf, decreasing_fit$E0 + 1e-8)
})
