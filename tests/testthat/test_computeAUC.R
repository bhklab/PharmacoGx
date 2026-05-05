library(PharmacoGx)

context("Testing computeAUC/computeAAC.")

test_that("Actual and fitted AUC/AAC values are complementary", {
  concentration <- c(0.1, 1, 10)
  viability <- c(0.9, 0.7, 0.4)
  fit <- c(HS = 1, E0 = 1, E_inf = 0.2, EC50 = 1)

  expect_equal(
    computeAUC(
      concentration = concentration,
      viability = viability,
      area.type = "Actual",
      viability_as_pct = FALSE
    ) +
      computeAAC(
        concentration = concentration,
        viability = viability,
        area.type = "Actual",
        viability_as_pct = FALSE
      ),
    1,
    tolerance = 1e-8
  )

  expect_equal(
    computeAUC(
      concentration = concentration,
      Hill_fit = fit,
      viability_as_pct = FALSE
    ) +
      computeAAC(
        concentration = concentration,
        Hill_fit = fit,
        viability_as_pct = FALSE
      ),
    1,
    tolerance = 1e-8
  )
})

test_that("curveFittingPGX exposes both AUC and AAC metrics", {
  concentration <- c(0.1, 0.3, 1, 3, 10)
  viability <- PharmacoGx:::.pgx_hill_curve(
    log10(concentration),
    c(1, 1, 0.2, 0)
  )
  scrn <- data.frame(
    cell_id = rep("cellA", length(concentration)),
    drug_id = rep("drugX", length(concentration)),
    conc = concentration,
    viability = viability
  )

  metrics <- curveFittingPGX(
    scrn,
    output_type = "metrics",
    main_fit_func = "hill"
  )

  expect_true(all(c("auc", "aac", "EC50", "ic50") %in% colnames(metrics)))
  expect_equal(metrics$auc + metrics$aac, 100, tolerance = 1e-6)
  expect_equal(metrics$ic50, metrics$EC50, tolerance = 0.1)
})

test_that("raw fitted AUC and AAC forward fit_type", {
  concentration <- 10^seq(-3, 3, length.out = 25)
  log_conc <- log10(concentration)
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
  viability <- PharmacoGx:::.pgx_biphasic_curve(log_conc, biphasic_pars)
  biphasic_fit <- logLogisticRegression(
    conc = concentration,
    viability = viability,
    viability_as_pct = FALSE,
    fit_type = "biphasic"
  )

  expect_equal(
    computeAUC(
      concentration = concentration,
      viability = viability,
      viability_as_pct = FALSE,
      fit_type = "biphasic"
    ),
    computeAUC(
      concentration = concentration,
      Hill_fit = biphasic_fit,
      viability_as_pct = FALSE
    ),
    tolerance = 1e-8
  )
  expect_equal(
    computeAAC(
      concentration = concentration,
      viability = viability,
      viability_as_pct = FALSE,
      fit_type = "biphasic"
    ),
    computeAAC(
      concentration = concentration,
      Hill_fit = biphasic_fit,
      viability_as_pct = FALSE
    ),
    tolerance = 1e-8
  )
  expect_equal(
    computeAUC(
      concentration = concentration,
      viability = viability,
      viability_as_pct = FALSE
    ),
    computeAUC(
      concentration = concentration,
      viability = viability,
      viability_as_pct = FALSE,
      fit_type = "hill"
    )
  )
})
