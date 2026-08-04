library(PharmacoGx)

context("Testing computeDSS.")

test_that("computeDSS is consistent across percent and fractional viabilities", {
  concentration <- 10^(seq(-2, 2))
  fit_pct <- c(HS = 1, E0 = 100, E_inf = 20, EC50 = 1)
  fit_fraction <- c(HS = 1, E0 = 1, E_inf = 0.2, EC50 = 1)

  dss_pct <- computeDSS(
    concentration = concentration,
    Hill_fit = fit_pct,
    viability_as_pct = TRUE
  )
  dss_fraction <- computeDSS(
    concentration = concentration,
    Hill_fit = fit_fraction,
    viability_as_pct = FALSE
  )

  expect_equal(dss_fraction, dss_pct, tolerance = 1e-8)
  expect_gt(dss_fraction, 0)
})

test_that("computeDSS supports deprecated fractional t_param inputs with a warning", {
  concentration <- 10^(seq(-2, 2))
  fit_fraction <- c(HS = 1, E0 = 1, E_inf = 0.2, EC50 = 1)

  expect_warning(
    deprecated_threshold <- computeDSS(
      concentration = concentration,
      Hill_fit = fit_fraction,
      viability_as_pct = FALSE,
      t_param = 0.1
    ),
    "deprecated"
  )

  expect_equal(
    deprecated_threshold,
    computeDSS(
      concentration = concentration,
      Hill_fit = fit_fraction,
      viability_as_pct = FALSE,
      t_param = 10
    ),
    tolerance = 1e-8
  )
})
