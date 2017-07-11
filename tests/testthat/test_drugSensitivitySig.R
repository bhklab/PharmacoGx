library(PharmacoGx)
require(parallel)
context("Checking drugSensitivitySig.")

test_that("Sensitivity result did not change since last time", {
	data(GDSCsmall)
	drug.sensitivity <- drugSensitivitySig(GDSCsmall, mDataType="rna", nthread=1, features = fNames(GDSCsmall, "rna")[1:50])	
	expect_equal_to_reference(drug.sensitivity@.Data, "drug.sensitivityGDSCSmall.rds")
	drug.sensitivity <- drugSensitivitySig(GDSCsmall, mDataType="rna", nthread=1, features = fNames(GDSCsmall, "rna")[1:50], sensitivity.cutoff = 0.2, sensitivity.measure="auc_recomputed")
	expect_equal_to_reference(drug.sensitivity@.Data, "drug.sensitivity.discreteGDSCSmall.rds", tolerance = 0.02)
	drug.sensitivity <- drugSensitivitySig(GDSCsmall, mDataType="rna", nthread=1, drugs=drugNames(GDSCsmall)[1:2], features = fNames(GDSCsmall, "rna")[1:10], sensitivity.measure=c("auc_recomputed","auc_published"))
	expect_equal_to_reference(drug.sensitivity@.Data, "drug.sensitivity.MANOVAGDSCSmall.rds")

})