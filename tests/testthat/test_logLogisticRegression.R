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
})
