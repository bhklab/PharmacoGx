library(PharmacoGx)

context("Checking intersectPSet.")

test_that("Intersection result did not change since last time",{

	data(GDSCsmall)
	data(CCLEsmall)
	common <- intersectPSet(list('GDSC'=GDSCsmall, 'CCLE'=CCLEsmall), intersectOn = c("drugs", "cell.lines", "concentrations"))
	expect_equal_to_reference(common, "intersectedSmallData.rds")

})