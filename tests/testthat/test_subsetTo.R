library(PharmacoGx)

context("Checking subsetTo.")


test_that("Intersection result did not change since last time", {
	data(CMAPsmall)
	CMAPsmaller <- subsetTo(CMAPsmall, drugs=drugNames(CMAPsmall), cells=cellNames(CMAPsmall))
	expect_equal(CMAPsmaller, CMAPsmall)
})

