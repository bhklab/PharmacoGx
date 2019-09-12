#library(PharmacoGx)
#
#context("Checking subsetTo.")
#
#test_that("Intersection result did not change since last time", {
# data(CCLEsmall)
# CCLEsmaller <- subsetTo(CCLEsmall, drugs=drugNames(CCLEsmall), cells=cellNames(CCLEsmall))
# expect_equal(CCLEsmaller@annotation, CCLEsmall@annotation)
# expect_equal(attributes(CCLEsmaller@molecularProfiles$rna), attributes(CCLEsmall$rna@assays))
#})