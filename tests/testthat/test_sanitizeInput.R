library(PharmacoGx)

context("Checking the sanitization of input to curve fitting and sensitivity summary funcitons")

test_that("Function sanitizeInput handles no input correctly.", {
  expect_error(sanitizeInput(), "Both 'Hill_fit' and 'viability'")
})