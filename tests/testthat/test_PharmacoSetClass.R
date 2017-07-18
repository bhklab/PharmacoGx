library(PharmacoGx)
require(parallel)
context("Checking PharmacoSet Class Methods.")


test_that("cellInfo result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(cellInfo(GDSCsmall), "cellInfo.GDSCsmall.rds")
})

test_that("drugInfo result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(drugInfo(GDSCsmall), "drugInfo.GDSCsmall.rds")
})

test_that("phenoInfo result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(phenoInfo(GDSCsmall, "rna"), "phenoInfo.GDSCsmall.rds")
})

test_that("molecularProfiles result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(molecularProfiles(GDSCsmall, "rna"), "molecularProfiles.GDSCsmall.rds")
})

test_that("featureInfo result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(featureInfo(GDSCsmall, "rna"), "featureInfo.GDSCsmall.rds")
})

test_that("sensitivityInfo result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(sensitivityInfo(GDSCsmall), "sensitivityInfo.GDSCsmall.rds")
})

test_that("sensitivityProfiles result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(sensitivityProfiles(GDSCsmall), "sensitivityProfiles.GDSCsmall.rds")
})


test_that("sensitivityMeasures result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(sensitivityMeasures(GDSCsmall), "sensitivityMeasures.GDSCsmall.rds")
})


test_that("drugNames result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(drugNames(GDSCsmall), "drugNames.GDSCsmall.rds")
})

test_that("cellNames result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(cellNames(GDSCsmall), "cellNames.GDSCsmall.rds")
})


test_that("fNames result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(fNames(GDSCsmall, "rna"), "fNames.GDSCsmall.rds")
})



test_that("pSetName result did not change since last time", {
	data(GDSCsmall)
	expect_equal(pSetName(GDSCsmall), "GDSC")
})


