library(PharmacoGx)

context("Checking the sanitization of input to curve fitting and sensitivity summary funcitons")

test_that("Function sanitizeInput handles no input correctly.", {
  expect_error(sanitizeInput(), "Both 'Hill_fit' and 'viability'")
})

test_that("Function yells at user sufficiently.", {
  
  expect_error(PharmacoGx:::sanitizeInput(conc=c(-3,-2,-1,0,1,2,3),viability=c(100,90,80,70,60,50,40),conc_as_log = FALSE,viability_as_pct = TRUE, verbose=TRUE))
  expect_warning(PharmacoGx:::sanitizeInput(conc=c(-3,-2,-1,0,1,2,3),viability=c(100,90,80,70,60,50,40),conc_as_log = TRUE,viability_as_pct = FALSE, verbose=TRUE))
  
})

test_that("Function returns correct values.", {
  
  expect_equal(PharmacoGx:::sanitizeInput(conc=c(-3,-2,-1,0,1,2,3),
                                          viability=c(100,90,80,70,60,50,40),
                                          conc_as_log = TRUE,
                                          viability_as_pct = TRUE,
                                          verbose=TRUE), list(log_conc=c(-3,-2,-1,0,1,2,3),viability=c(100,90,80,70,60,50,40)/100))
  
})