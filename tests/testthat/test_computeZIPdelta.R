library(PharmacoGx)

context("Testing ZIP delta helpers.")

test_that("estimateProjParams returns the same projected parameters with or without Rsqr", {
  dose_to <- c(0.1, 0.3, 1, 3, 10)
  dose_add <- 1
  HS_add <- 1.2
  EC50_add <- 0.7
  E_inf_add <- 0.3
  E_ninf_proj <- .Hill(log10(dose_add), c(HS_add, E_inf_add, log10(EC50_add)))

  combo_viability <- hillCurve(
    dose = log10(dose_to),
    HS = 0.9,
    EC50 = log10(1.5),
    E_inf = 0.15,
    E_ninf = E_ninf_proj
  )

  with_rsq <- estimateProjParams(
    dose_to = dose_to,
    combo_viability = combo_viability,
    dose_add = dose_add,
    EC50_add = EC50_add,
    HS_add = HS_add,
    E_inf_add = E_inf_add,
    show_Rsqr = TRUE
  )
  without_rsq <- estimateProjParams(
    dose_to = dose_to,
    combo_viability = combo_viability,
    dose_add = dose_add,
    EC50_add = EC50_add,
    HS_add = HS_add,
    E_inf_add = E_inf_add,
    show_Rsqr = FALSE
  )

  expect_equal(without_rsq$HS_proj, with_rsq$HS_proj, tolerance = 1e-6)
  expect_equal(without_rsq$EC50_proj, with_rsq$EC50_proj, tolerance = 1e-6)
  expect_equal(without_rsq$E_inf_proj, with_rsq$E_inf_proj, tolerance = 1e-6)
  expect_equal(without_rsq$E_ninf_proj, with_rsq$E_ninf_proj, tolerance = 1e-6)
})

test_that(".computeZIPdelta matches between show_Rsqr branches", {
  doses1 <- rep(c(0.1, 1, 10), each = 3)
  doses2 <- rep(c(0.1, 1, 10), times = 3)
  HS_1 <- rep(1, 9)
  HS_2 <- rep(1.2, 9)
  EC50_1 <- rep(0.5, 9)
  EC50_2 <- rep(0.7, 9)
  E_inf_1 <- rep(0.2, 9)
  E_inf_2 <- rep(0.3, 9)
  zip <- computeZIP(
    treatment1dose = doses1,
    HS_1 = HS_1,
    EC50_1 = EC50_1,
    E_inf_1 = E_inf_1,
    treatment2dose = doses2,
    HS_2 = HS_2,
    EC50_2 = EC50_2,
    E_inf_2 = E_inf_2
  )

  with_rsq <- .computeZIPdelta(
    treatment1id = rep("d1", 9),
    treatment2id = rep("d2", 9),
    treatment1dose = doses1,
    treatment2dose = doses2,
    sampleid = rep("s1", 9),
    HS_1 = HS_1,
    HS_2 = HS_2,
    EC50_1 = EC50_1,
    EC50_2 = EC50_2,
    E_inf_1 = E_inf_1,
    E_inf_2 = E_inf_2,
    combo_viability = zip,
    ZIP = zip,
    show_Rsqr = TRUE
  )
  without_rsq <- .computeZIPdelta(
    treatment1id = rep("d1", 9),
    treatment2id = rep("d2", 9),
    treatment1dose = doses1,
    treatment2dose = doses2,
    sampleid = rep("s1", 9),
    HS_1 = HS_1,
    HS_2 = HS_2,
    EC50_1 = EC50_1,
    EC50_2 = EC50_2,
    E_inf_1 = E_inf_1,
    E_inf_2 = E_inf_2,
    combo_viability = zip,
    ZIP = zip,
    show_Rsqr = FALSE
  )

  expect_equal(without_rsq$delta_score, with_rsq$delta_score, tolerance = 1e-6)
  expect_equal(without_rsq$treatment1dose, with_rsq$treatment1dose)
  expect_equal(without_rsq$treatment2dose, with_rsq$treatment2dose)
})
