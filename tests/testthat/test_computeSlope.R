library(PharmacoGx)

context("Checking computeSlope.")

test_that("computeSlope matches ordinary regression for reported examples", {
  dose <- c(0.0025, 0.008, 0.025, 0.08, 0.25, 0.8, 2.53, 8)

  expect_equal(
    computeSlope(
      concentration = dose,
      viability = c(120, 111, 104, 100, 80, 56, 30, 20),
      trunc = TRUE,
      verbose = FALSE
    ),
    0.25
  )
  expect_equal(
    computeSlope(
      concentration = dose,
      viability = c(120, 111, 104, 100, 80, 56, 30, 20),
      trunc = FALSE,
      verbose = FALSE
    ),
    0.30
  )

  expect_equal(
    computeSlope(
      concentration = dose,
      viability = c(108.67, 111, 102.16, 100.27, 70, 56, 30, 10),
      trunc = TRUE,
      verbose = FALSE
    ),
    0.27
  )
  expect_equal(
    computeSlope(
      concentration = dose,
      viability = c(108.67, 111, 102.16, 100.27, 70, 56, 30, 10),
      trunc = FALSE,
      verbose = FALSE
    ),
    0.30
  )

  expect_equal(
    computeSlope(
      concentration = dose,
      viability = c(108.67, 111, 102.16, 100.27, 90, 87, 74, 57),
      trunc = TRUE,
      verbose = FALSE
    ),
    0.11
  )
  expect_equal(
    computeSlope(
      concentration = dose,
      viability = c(108.67, 111, 102.16, 100.27, 90, 87, 74, 57),
      trunc = FALSE,
      verbose = FALSE
    ),
    0.14
  )
})
