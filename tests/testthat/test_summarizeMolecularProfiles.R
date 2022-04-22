library(PharmacoGx)
library(SummarizedExperiment)
library(S4Vectors)

context("Checking summarizeMolecularProfiles function.")

data("GDSCsmall")

test_that("Summarize Molecular Profiles fails gracefully.",{
  ## FIXME:: No method defition for summarizeMolecularProfiles with class 'missing'
  #expect_error(summarizeMolecularProfiles(), "argument \"pSet\" is missing")
  expect_error(summarizeMolecularProfiles(GDSCsmall), "argument \"mDataType\" is missing")
  expect_error(summarizeMolecularProfiles(GDSCsmall, "rnaseq"), "Invalid mDataType")
})

test_that("Summarize Molecular Profiles function outputs data with right dimensions and dimnames, class", {
  testSummary <- summarizeMolecularProfiles(GDSCsmall, "rna")
  expect_equal(colnames(testSummary), sampleNames(GDSCsmall))
  expect_equivalent(is(testSummary, "SummarizedExperiment"), TRUE)
  expect_length(rownames(testSummary), 300)
})

test_that("Summarize Molecular Profiles correctly summarizes replicates", {
  myx <- "647-V" == colData(molecularProfilesSlot(GDSCsmall)$rna)$sampleid
  testCells <- SummarizedExperiment::assay(molecularProfilesSlot(GDSCsmall)$rna, 1)[,myx]
  testSummary <- summarizeMolecularProfiles(GDSCsmall, "rna", summary.stat = "median")
  expect_equal(SummarizedExperiment::assay(testSummary, 1)[,"647-V"], apply(testCells, 1, median))
  testSummary <- summarizeMolecularProfiles(GDSCsmall, "rna", summary.stat = "mean")
  expect_equal(SummarizedExperiment::assay(testSummary, 1)[,"647-V"], apply(testCells, 1, mean))
  testSummary <- summarizeMolecularProfiles(GDSCsmall, "rna", summary.stat = "first")
  expect_equal(SummarizedExperiment::assay(testSummary, 1)[,"647-V"], testCells[,1])
  testSummary <- summarizeMolecularProfiles(GDSCsmall, "rna", summary.stat = "last")
  expect_equal(SummarizedExperiment::assay(testSummary, 1)[,"647-V"], testCells[,-1])

  GDSCsmall2 <- subsetTo(GDSCsmall, cells = c("22RV1", "23132-87"))
  colData(molecularProfilesSlot(GDSCsmall2)$mutation)$sampleid <- "22RV1"
  testCells <- SummarizedExperiment::assay(molecularProfilesSlot(GDSCsmall2)$mutation, 1)

  testSummary <- summarizeMolecularProfiles(GDSCsmall2, "mutation", summary.stat = "or")
  expect_equal(sum(as.numeric(SummarizedExperiment::assay(testSummary, 1)), na.rm=TRUE), 2)

  testSummary <- summarizeMolecularProfiles(GDSCsmall2, "mutation", summary.stat = "and")
  expect_equal(sum(as.numeric(SummarizedExperiment::assay(testSummary, 1)), na.rm=TRUE), 0)

})


test_that("Summarize Molecular Profiles parameters work as expected", {
  expect_equal(summarizeMolecularProfiles(GDSCsmall, "rna", summarize=FALSE), molecularProfilesSlot(GDSCsmall)$rna)
  expect_silent(summarizeMolecularProfiles(GDSCsmall, "rna", verbose = FALSE))
  GDSCsmall2 <- GDSCsmall
  molecularProfilesSlot(GDSCsmall2)$rna <- molecularProfilesSlot(GDSCsmall2)$rna[,1]
  expect_equivalent(ncol(summarizeMolecularProfiles(GDSCsmall2, "rna", fill.missing = FALSE)), 1)
  expect_equivalent(ncol(summarizeMolecularProfiles(GDSCsmall2, "rna", fill.missing = TRUE)), length(sampleNames(GDSCsmall2)))
})
