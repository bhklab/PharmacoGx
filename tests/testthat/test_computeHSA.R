library(PharmacoGx)

context("Testing computeHSA.")

test_that("computeHSA handles missing values according to na.rm", {
  expect_equal(
    computeHSA(c(0.75, NA, 0.4), c(0.65, 0.5, NA)),
    c(0.65, NA, NA)
  )

  expect_equal(
    computeHSA(c(0.75, NA, 0.4), c(0.65, 0.5, NA), na.rm = TRUE),
    c(0.65, 0.5, 0.4)
  )
})
