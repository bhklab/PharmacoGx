library(PharmacoGx)

context("Checking computeAmax.")

test_that("computeAmax treats percent and decimal viabilities consistently", {
  dose <- c(0.001, 0.004, 0.016, 0.064, 0.256, 1.024, 4.096, 16.384)
  viability_pct <- c(99, 95, 80, 50, 25, 16, 15, 15)
  viability_dec <- c(0.99, 0.95, 0.80, 0.50, 0.25, 0.16, 0.15, 0.15)

  amax_pct <- computeAmax(
    concentration = dose,
    viability = viability_pct,
    viability_as_pct = TRUE,
    trunc = FALSE
  )
  amax_dec <- computeAmax(
    concentration = dose,
    viability = viability_dec,
    viability_as_pct = FALSE,
    trunc = FALSE
  )

  expect_equal(amax_dec, amax_pct, tolerance = 1e-4)
})

test_that("computeAmax validates viability_as_pct", {
  expect_error(
    computeAmax(
      concentration = c(1, 2, 3),
      viability = c(50, 40, 30),
      viability_as_pct = NA
    ),
    "'viability_as_pct' must be a logical value."
  )
})
