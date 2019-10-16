library(PharmacoGx)

context("Checking summarizeMolecularProfiles function.")

data("GDSCsmall")

test_that("Summarize Molecular Profiles fails gracefully.",{
  expect_error(summarizeMolecularProfiles(), "argument \"pSet\" is missing")
  expect_error(summarizeMolecularProfiles(GDSCsmall), "argument \"mDataType\" is missing")
  expect_error(summarizeMolecularProfiles(GDSCsmall, "rnaseq"), "Invalid mDataType")
})

test_that("Summarize Molecular Profiles function outputs data with right dimensions and dimnames, class", {
  testSummary <- summarizeMolecularProfiles(GDSCsmall, "rna")
  expect_equal(colnames(testSummary), cellNames(GDSCsmall))
  expect_equivalent(class(testSummary), "ExpressionSet")
  expect_length(rownames(testSummary), 300)
})

test_that("Summarize Molecular Profiles correctly summarizes replicates", {
  myx <- "647-V" == Biobase::pData(GDSCsmall@molecularProfiles$rna)$cellid
  testCells <- Biobase::exprs(GDSCsmall@molecularProfiles$rna)[,myx]
  testSummary <- summarizeMolecularProfiles(GDSCsmall, "rna", summary.stat = "median")
  expect_equal(Biobase::exprs(testSummary)[,"647-V"], apply(testCells, 1, median))
  testSummary <- summarizeMolecularProfiles(GDSCsmall, "rna", summary.stat = "mean")
  expect_equal(Biobase::exprs(testSummary)[,"647-V"], apply(testCells, 1, mean))
  testSummary <- summarizeMolecularProfiles(GDSCsmall, "rna", summary.stat = "first")
  expect_equal(Biobase::exprs(testSummary)[,"647-V"], testCells[,1])
  testSummary <- summarizeMolecularProfiles(GDSCsmall, "rna", summary.stat = "last")
  expect_equal(Biobase::exprs(testSummary)[,"647-V"], testCells[,-1])
  
  GDSCsmall2 <- subsetTo(GDSCsmall, cells = c("22RV1", "23132-87"))
  Biobase::pData(GDSCsmall2@molecularProfiles$mutation)$cellid <- "22RV1"
  testCells <- Biobase::exprs(GDSCsmall2@molecularProfiles$mutation)
  
  testSummary <- summarizeMolecularProfiles(GDSCsmall2, "mutation", summary.stat = "or")
  expect_equal(sum(as.numeric(exprs(testSummary)), na.rm=TRUE), 2)
  
  testSummary <- summarizeMolecularProfiles(GDSCsmall2, "mutation", summary.stat = "and")
  expect_equal(sum(as.numeric(exprs(testSummary)), na.rm=TRUE), 0)
  
})


test_that("Summarize Molecular Profiles parameters work as expected", {
  expect_equal(summarizeMolecularProfiles(GDSCsmall, "rna", summarize=FALSE), GDSCsmall@molecularProfiles$rna)
  expect_silent(summarizeMolecularProfiles(GDSCsmall, "rna", verbose = FALSE))
  GDSCsmall2 <- GDSCsmall
  GDSCsmall2@molecularProfiles$rna <- GDSCsmall2@molecularProfiles$rna[,1]
  expect_equivalent(ncol(summarizeMolecularProfiles(GDSCsmall2, "rna", fill.missing = FALSE)), 1)
  expect_equivalent(ncol(summarizeMolecularProfiles(GDSCsmall2, "rna", fill.missing = TRUE)), length(cellNames(GDSCsmall2)))
})