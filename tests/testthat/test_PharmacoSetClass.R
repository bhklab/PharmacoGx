library(PharmacoGx)
require(parallel)
context("Checking PharmacoSet Class Methods.")


test_that("cellInfo result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(sampleInfo(GDSCsmall), "cellInfo.GDSCsmall.rds")
})

test_that("drugInfo result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(treatmentInfo(GDSCsmall), "drugInfo.GDSCsmall.rds")
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
	expect_equal_to_reference(treatmentNames(GDSCsmall), "drugNames.GDSCsmall.rds")
})

test_that("cellNames result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(sampleNames(GDSCsmall), "cellNames.GDSCsmall.rds")
})


test_that("fNames result did not change since last time", {
	data(GDSCsmall)
	expect_equal_to_reference(fNames(GDSCsmall, "rna"), "fNames.GDSCsmall.rds")
})


test_that("name result did not change since last time", {
	data(GDSCsmall)
	expect_equal(name(GDSCsmall), "GDSC")
})

test_that("updateSampleId works without duplicates", {
	data(GDSCsmall)
	newNames <- c("Test","Test2",sampleNames(GDSCsmall)[3:length(sampleNames(GDSCsmall))])


	sampleNames(GDSCsmall) <- newNames

	expect_true(all(unique(sensitivityInfo(GDSCsmall)$sampleid) %in% newNames))
	expect_true(all(unique(sensitivityInfo(GDSCsmall)$sampleid) %in% newNames))
	expect_equal(sort(unique(rownames(sampleInfo(GDSCsmall)))), sort(newNames))
	expect_equal(sort(rownames(sensNumber(GDSCsmall))), sort(newNames))

})


test_that("updateSampleId works with duplicates", {
	data(GDSCsmall)
	newNames <- c("Test","Test",sampleNames(GDSCsmall)[3:length(sampleNames(GDSCsmall))])


	expect_warning(sampleNames(GDSCsmall) <- newNames, "Duplicated ids passed to updateSampleId. Merging old ids into the same identifier")

	expect_true(all(unique(sensitivityInfo(GDSCsmall)$sampleid) %in% newNames))
	expect_equal(sort(unique(rownames(sampleInfo(GDSCsmall)))), sort(unique(newNames)))
	expect_equal(sort(rownames(sensNumber(GDSCsmall))), sort(unique(newNames)))

})



test_that("updateTreatmentId works without duplicates", {
	data(GDSCsmall)
	newNames <- c("Test","Test2",treatmentNames(GDSCsmall)[3:length(treatmentNames(GDSCsmall))])

	treatmentNames(GDSCsmall) <- newNames

	expect_true(all(unique(sensitivityInfo(GDSCsmall)$treatmentid) %in% newNames))
	expect_equal(sort(unique(rownames(treatmentInfo(GDSCsmall)))), sort(newNames))
	expect_equal(sort(colnames(sensNumber(GDSCsmall))), sort(newNames))

})


test_that("updateTreatmentId works with duplicates", {
    data(GDSCsmall)
    newNames <- c("Test","Test",treatmentNames(GDSCsmall)[3:length(treatmentNames(GDSCsmall))])

    expect_warning(treatmentNames(GDSCsmall) <- newNames,
		"Duplicated ids passed to updateTreatmentId. Merging old ids into the same identifier")

    expect_true(all(unique(sensitivityInfo(GDSCsmall)$treatmentid) %in% newNames))
    expect_equal(sort(unique(rownames(treatmentInfo(GDSCsmall)))), sort(unique(newNames)))
    expect_equal(sort(colnames(sensNumber(GDSCsmall))), sort(unique(newNames)))
})
