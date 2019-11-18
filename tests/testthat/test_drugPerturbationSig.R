library(PharmacoGx)
require(parallel)
context("Checking drugPerturbationSig.")

test_that("Perturbation result did not change since last time", {
	data(CMAPsmall)
	drug.perturbation <- drugPerturbationSig(CMAPsmall, mDataType="rna", nthread=1)	
	expect_equal_to_reference(drug.perturbation@.Data, "drug.perturbationSmall.rds")
})