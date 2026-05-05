library(PharmacoGx)

context("Checking computeIC50/ICn.")

test_that("Function complains when given insensible input", {
  expect_error(
    computeIC50(
      concentration = c(1, 2, 3),
      viability = c(50, 60, 70),
      Hill_fit = c(1, 0, 0.1)
    ),
    "Please pass in only one"
  )
  # expect_silent(computeIC50(concentration = c(1, 2, 3),
  #     #   viability1 = c(50, 60, 70),
  #   Hill_fit2 = c(0.5, 0.2, 1)))

  expect_error(
    computeIC50(concentration = c(1, 2, 3, 5), viability = c(50, 60, 70)),
    "is not of same length"
  ) #should complain
  expect_error(
    computeIC50(
      concentration = c(-1, 2, 3),
      viability = c(50, 60, 70),
      conc_as_log = FALSE
    ),
    "'conc_as_log' flag may be set incorrectly"
  ) #should complain
  ##TO-DO:: Add wanring strings to expect_warning call
  expect_error(
    computeIC50(
      concentration = c(NA, "cat", 3),
      viability = c(50, 60, 70),
      conc_as_log = FALSE
    ),
    "At least three unique concentrations"
  ) #should complain
  expect_error(
    computeIC50(concentration = c(1, 2, Inf), viability = c(50, 60, 70)),
    "non-finite values"
  ) #should complain
  expect_warning(
    computeIC50(
      concentration = c(1, 2, 3),
      viability = c(.50, .60, .70),
      viability_as_pct = TRUE
    ),
    "as_pct"
  ) #should complain
  expect_error(computeIC50()) #should complain
})

test_that("Functions return right values", {
  expect_equal(
    computeIC50(
      concentration = seq(-3, 3),
      Hill_fit = c(1, 0, 0),
      conc_as_log = TRUE,
      viability_as_pct = FALSE
    ),
    0
  )
  expect_equal(
    computeIC50(
      concentration = seq(1, 3),
      Hill_fit = c(1, 0, 0),
      conc_as_log = TRUE,
      viability_as_pct = FALSE
    ),
    0
  )
  expect_equal(
    computeIC50(
      concentration = seq(1, 3),
      Hill_fit = c(1, .9, 0),
      conc_as_log = TRUE,
      viability_as_pct = FALSE
    ),
    0
  )
  expect_equal(
    computeIC50(
      concentration = seq(1, 3),
      Hill_fit = c(1, .5, 0),
      conc_as_log = TRUE,
      viability_as_pct = FALSE
    ),
    0
  )
  expect_equal(
    computeAC50(
      concentration = seq(1, 3),
      Hill_fit = c(1, .6, .55),
      conc_as_log = TRUE,
      viability_as_pct = FALSE
    ),
    Inf
  )
  expect_equal(
    .Hill(
      computeICn(
        concentration = seq(1, 3),
        Hill_fit = c(1, 0, 0),
        n = .7,
        conc_as_log = TRUE,
        viability_as_pct = FALSE
      ),
      c(1, 0, 0)
    ),
    .3
  )
  expect_equal(
    computeICn(
      concentration = seq(1, 3),
      Hill_fit = c(1, 0, 0),
      n = 0,
      conc_as_log = TRUE,
      viability_as_pct = FALSE
    ),
    -Inf
  )
})

test_that("Relative IC50 follows fitted half-max while AC50 keeps absolute target", {
  fit <- c(HS = 1, E0 = 1, E_inf = 0.2, EC50 = 1)

  expect_equal(
    computeIC50(
      concentration = 10^(seq(-2, 2)),
      Hill_fit = fit,
      conc_as_log = FALSE,
      viability_as_pct = FALSE
    ),
    1
  )

  expect_equal(
    computeAC50(
      concentration = 10^(seq(-2, 2)),
      Hill_fit = fit,
      conc_as_log = FALSE,
      viability_as_pct = FALSE
    ),
    5 / 3
  )

  expect_warning(
    computeIC50(
      concentration = 10^(seq(-2, 2)),
      Hill_fit = fit,
      reference = "absolute",
      conc_as_log = FALSE,
      viability_as_pct = FALSE
    ),
    "deprecated"
  )
})

test_that("ICn helpers reject biphasic fits with a clear error", {
  biphasic_fit <- c(
    HS1 = 1.3,
    E0 = 1,
    E_inf1 = 0.3,
    HS2 = 0.8,
    E_inf2 = 0.1,
    EC50_1 = 0.3,
    EC50_2 = 3,
    Frac = 0.6
  )

  expect_error(
    computeIC50(
      concentration = 10^seq(-2, 2),
      Hill_fit = biphasic_fit,
      viability_as_pct = FALSE
    ),
    "currently support Hill fits only"
  )
  expect_error(
    computeICn(
      concentration = 10^seq(-2, 2),
      Hill_fit = biphasic_fit,
      n = 0.5,
      viability_as_pct = FALSE
    ),
    "currently support Hill fits only"
  )
})
