library(PharmacoGx)

context("drugSensitivitySig defaults")

test_that("drugSensitivitySig works with defaulted optional arguments", {
  data(CCLEsmall)
  res <- drugSensitivitySig(
    CCLEsmall,
    mDataType = "rna",
    features = rownames(featureInfo(CCLEsmall, "rna"))[1],
    sensitivity.measure = "auc_recomputed",
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
