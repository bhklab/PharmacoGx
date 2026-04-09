library(PharmacoGx)

context("drugSensitivitySig defaults")

test_that("drugSensitivitySig works with explicit optional arguments", {
  data(CCLEsmall)
  res <- drugSensitivitySig(
    CCLEsmall,
    mDataType = "rna",
    features = rownames(featureInfo(CCLEsmall, "rna"))[1],
    sensitivity.measure = "aac_recomputed",
    molecular.summary.stat = "mean",
    sensitivity.summary.stat = "mean",
    returnValues = "estimate",
    modeling.method = "anova",
    inference.method = "analytic",
    parallel.on = "drug",
    nthread = 1,
    verbose = FALSE
  )
  expect_s4_class(res, "PharmacoSig")
  expect_equal(dim(res), c(1, length(treatmentNames(CCLEsmall)), 8))
})

test_that("drugSensitivitySig warns on legacy AUC sensitivity aliases", {
  data(CCLEsmall)
  expect_warning(
    drugSensitivitySig(
      CCLEsmall,
      mDataType = "rna",
      features = rownames(featureInfo(CCLEsmall, "rna"))[1],
      sensitivity.measure = "auc_recomputed",
      returnValues = "estimate",
      parallel.on = "drug",
      nthread = 1,
      verbose = TRUE
    ),
    "deprecated"
  )
})

test_that("drugSensitivitySig works with default optional arguments", {
  data(CCLEsmall)
  res <- drugSensitivitySig(
    CCLEsmall,
    mDataType = "rna",
    features = rownames(featureInfo(CCLEsmall, "rna"))[1],
    parallel.on = "drug",
    nthread = 1,
    verbose = FALSE
  )
  expect_s4_class(res, "PharmacoSig")
  expect_equal(dim(res), c(1, length(treatmentNames(CCLEsmall)), 8))
})
