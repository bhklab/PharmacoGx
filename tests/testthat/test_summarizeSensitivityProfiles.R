library(PharmacoGx)

context("Checking summarizeSensitivityProfiles function.")

data("GDSCsmall")

## FIXME:: No S4 method for summarizeSensitivityProfiles with class 'missing'
#test_that("Summarize Sensitivity Profiles fails gracefully.", {
#  expect_error(summarizeSensitivityProfiles(), "argument \"pSet\" is missing")
#})


test_that("Summarize Sensitivity Profiles function outputs data with right dimensions and dimnames, class", {
  testSummary <- summarizeSensitivityProfiles(GDSCsmall)
  expect_equal(colnames(testSummary), sampleNames(GDSCsmall))
  expect_equal(rownames(testSummary), treatmentNames(GDSCsmall))
  expect_equivalent(is(testSummary, "matrix"), TRUE)
})

test_that("summarizeSensitivityProfiles produces correct values.",{

  GDSCsmall2 <- subsetTo(GDSCsmall, drugs="AZD6482")
  testCells <- sensitivityProfiles(GDSCsmall2)[order(sensitivityInfo(GDSCsmall2)$sampleid),"auc_recomputed", drop=FALSE]

  testSummary <- summarizeSensitivityProfiles(GDSCsmall2, summary.stat = "median", fill.missing=FALSE)
  testSummary <- testSummary[order(names(testSummary))]
  names(testSummary)<- NULL
  expect_equivalent(testSummary, mapply(function(x,y) {median(c(x,y))}, testCells[seq(1,18,2),], testCells[seq(1,18,2)+1,]))

  testSummary <- summarizeSensitivityProfiles(GDSCsmall2, summary.stat = "mean", fill.missing=FALSE)
  testSummary <- testSummary[order(names(testSummary))]
  names(testSummary)<- NULL
  expect_equivalent(testSummary, mapply(function(x,y) {mean(c(x,y))}, testCells[seq(1,18,2),], testCells[seq(1,18,2)+1,]))

  testSummary <- summarizeSensitivityProfiles(GDSCsmall2, summary.stat = "first", fill.missing=FALSE)
  testSummary <- testSummary[order(names(testSummary))]
  names(testSummary)<- NULL
  expect_equivalent(testSummary, mapply(function(x,y) {x}, testCells[seq(1,18,2),], testCells[seq(1,18,2)+1,]))

  testSummary <- summarizeSensitivityProfiles(GDSCsmall2, summary.stat = "last", fill.missing=FALSE)
  testSummary <- testSummary[order(names(testSummary))]
  names(testSummary)<- NULL
  expect_equivalent(testSummary, mapply(function(x,y) {y}, testCells[seq(1,18,2),], testCells[seq(1,18,2)+1,]))

})


test_that("Summarize Sensitivity Profiles parameters work as expected", {
  expect_silent(summarizeSensitivityProfiles(GDSCsmall, verbose = FALSE))
  expect_equal(ncol(summarizeSensitivityProfiles(GDSCsmall, fill.missing = FALSE)), length(unique(sensitivityInfo(GDSCsmall)$sampleid)))
  expect_equal(ncol(summarizeSensitivityProfiles(GDSCsmall, fill.missing = TRUE)), length(sampleNames(GDSCsmall)))
})
