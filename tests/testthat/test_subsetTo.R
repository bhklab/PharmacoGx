library(PharmacoGx)

context("Checking subsetTo.")

test_that("Intersection result did not change since last time", {
	data(CMAPsmall)
	CMAPsmaller <- subsetTo(CMAPsmall, drugs=drugNames(CMAPsmall)[1:5], cells=cellNames(CMAPsmall)[1:5])
	expect_equal_to_reference(CMAPsmaller, "subsetTo.CMAPsmaller.rds")
})

