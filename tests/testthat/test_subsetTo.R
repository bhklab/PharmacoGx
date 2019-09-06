library(PharmacoGx)

context("Checking subsetTo.")


test_that("Intersection result did not change since last time", {
	data(CCLE)
	CCLEsmaller <- subsetTo(CCLE, drugs=drugNames(CCLE), cells=cellNames(CCLE))
	expect_equal(CCLEsmaller, CCLE)
})