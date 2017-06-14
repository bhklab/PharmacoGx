library(PharmacoGx)

context("Checking the sanitization of input to curve fitting and sensitivity summary funcitons")

test_that("Function sanitizeInput handles no input correctly.", {
  expect_error(sanitizeInput(), "Both 'Hill_fit' and 'viability'")
})

test_that("Function sanitizeInput yells at user sufficiently.", {
  expect_error(PharmacoGx:::sanitizeInput(conc=c(-3,-2,-1,0,1,2,3),viability=c(100,90,80,70,60,50,40),conc_as_log = FALSE,viability_as_pct = TRUE, verbose=TRUE))
  expect_error(PharmacoGx:::sanitizeInput(conc=c(-3,-2,-1),viability=c(100,90,80,70,60,50,40),conc_as_log = FALSE,viability_as_pct = TRUE, verbose=TRUE))
  expect_error(PharmacoGx:::sanitizeInput(conc=c(-3,-2,-1,0,1,2,3),viability=c(100,90,80,70),conc_as_log = FALSE,viability_as_pct = TRUE, verbose=TRUE))
  expect_warning(PharmacoGx:::sanitizeInput(conc=c(-3,-2,-1,0,1,2,3),viability=c(100,90,80,70,60,50,40),conc_as_log = TRUE,viability_as_pct = FALSE, verbose=TRUE), "'viability_as_pct' flag may be set incorrectly")
})

test_that("Function sanitizeInput returns correct values.", {
  expect_equal(PharmacoGx:::sanitizeInput(conc=c(-3,-2,-1,0,1,2,3),
                                          viability=c(100,90,80,70,60,50,40),
                                          conc_as_log = TRUE,
                                          viability_as_pct = TRUE,
                                          verbose=TRUE), list(log_conc=c(-3,-2,-1,0,1,2,3),viability=c(100,90,80,70,60,50,40)/100))
  expect_equal(PharmacoGx:::sanitizeInput(conc=c(0,1,2,3),
                                          viability=c(100,90,80,70),
                                          conc_as_log = FALSE,
                                          viability_as_pct = TRUE,
                                          verbose=TRUE), list(log_conc=log10(c(1,2,3)),viability=c(90,80,70)/100))
  expect_equal(PharmacoGx:::sanitizeInput(conc=c(0,1,2,3),
                                          viability=c(100,90,80,70)/100,
                                          conc_as_log = FALSE,
                                          viability_as_pct = FALSE,
                                          verbose=TRUE), list(log_conc=log10(c(1,2,3)),viability=c(90,80,70)/100))
  expect_equal(PharmacoGx:::sanitizeInput(conc=c(-3,-2,-1,0,1,2,3),
                                          viability=c(110,90,80,70,60,50,40),
                                          conc_as_log = TRUE,
                                          viability_as_pct = TRUE,
                                          verbose=FALSE), list(log_conc=c(-3,-2,-1,0,1,2,3),viability=c(100,90,80,70,60,50,40)/100))
  expect_equal(PharmacoGx:::sanitizeInput(conc=c(-3,-2,-1,0,1,2,3),
                                          viability=c(110,90,80,70,60,50,40)/100,
                                          conc_as_log = TRUE,
                                          viability_as_pct = FALSE,
                                          verbose=FALSE), list(log_conc=c(-3,-2,-1,0,1,2,3),viability=c(100,90,80,70,60,50,40)/100))
  expect_equal(PharmacoGx:::sanitizeInput(conc=c(-3,-2,-1,0,1,2,3),
                                          viability=c(110,90,80,70,60,50,40),
                                          conc_as_log = TRUE,
                                          viability_as_pct = TRUE, trunc = FALSE,
                                          verbose=FALSE), list(log_conc=c(-3,-2,-1,0,1,2,3),viability=c(110,90,80,70,60,50,40)/100))
})