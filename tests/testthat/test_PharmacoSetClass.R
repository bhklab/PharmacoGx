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

test_that("updateCellId works without duplicates", {
	data(GDSCsmall)
	newNames <- c("Test","Test2",cellNames(GDSCsmall)[3:length(cellNames(GDSCsmall))])


	cellNames(GDSCsmall) <- newNames

	expect_true(all(unique(sensitivityInfo(GDSCsmall)$cellid) %in% newNames))
	expect_true(all(unique(sensitivityInfo(GDSCsmall)$cellid) %in% newNames))
	expect_equal(sort(unique(rownames(cellInfo(GDSCsmall)))), sort(newNames))
	expect_equal(sort(rownames(sensNumber(GDSCsmall))), sort(newNames))

})


test_that("updateCellId works with duplicates", {
	data(GDSCsmall)
	newNames <- c("Test","Test",cellNames(GDSCsmall)[3:length(cellNames(GDSCsmall))])


	cellNames(GDSCsmall) <- newNames

	expect_true(all(unique(sensitivityInfo(GDSCsmall)$cellid) %in% newNames))
	expect_equal(sort(unique(rownames(cellInfo(GDSCsmall)))), sort(unique(newNames)))
	expect_equal(sort(rownames(sensNumber(GDSCsmall))), sort(unique(newNames)))

})



test_that("updateDrugId works without duplicates", {
	data(GDSCsmall)
	newNames <- c("Test","Test2",drugNames(GDSCsmall)[3:length(drugNames(GDSCsmall))])


	drugNames(GDSCsmall) <- newNames

	expect_true(all(unique(sensitivityInfo(GDSCsmall)$drugid) %in% newNames))
	expect_equal(sort(unique(rownames(drugInfo(GDSCsmall)))), sort(newNames))
	expect_equal(sort(colnames(sensNumber(GDSCsmall))), sort(newNames))

})


test_that("updateDrugId works without duplicates", {
    data(GDSCsmall)
    newNames <- c("Test","Test",drugNames(GDSCsmall)[3:length(drugNames(GDSCsmall))])


    drugNames(GDSCsmall) <- newNames

    expect_true(all(unique(sensitivityInfo(GDSCsmall)$drugid) %in% newNames))
    expect_equal(sort(unique(rownames(drugInfo(GDSCsmall)))), sort(unique(newNames)))
    expect_equal(sort(colnames(sensNumber(GDSCsmall))), sort(unique(newNames)))

})























