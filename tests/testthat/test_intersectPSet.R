library(PharmacoGx)
require(parallel)
context("Checking intersectPSet.")

test_that("Intersection result did not change since last time", {
	data(GDSCsmall)
	data(CCLEsmall)
	common <- intersectPSet(list('GDSC'=GDSCsmall, 'CCLE'=CCLEsmall), intersectOn = c("drugs", "cell.lines","concentrations"))	
	expect_equal_to_reference(common, "intersectedSmallData.rds")
	expect_equal(sum(!is.na(common$CCLE@sensitivity$profiles$auc_recomputed_star)),sum(!is.na(common$CCLE@sensitivity$profiles$auc_recomputed_star)))
	expect_equal(drugNames(common$CCLE), drugNames(common$GDSC))
	expect_equal(cellNames(common$CCLE), cellNames(common$GDSC))
})

test_that("Strict Intersection result did not change since last time", {
	data(GDSCsmall)
	data(CCLEsmall)
	common <- intersectPSet(list('GDSC'=GDSCsmall, 'CCLE'=CCLEsmall), intersectOn = c("drugs", "cell.lines","concentrations"), strictIntersect=TRUE)	
	expect_equal_to_reference(common, "intersectedSmallDataStrict.rds")
	expect_equal(sum(!is.na(common$CCLE@sensitivity$profiles$auc_recomputed_star)),sum(!is.na(common$CCLE@sensitivity$profiles$auc_recomputed_star)))
	expect_equal(drugNames(common$CCLE), drugNames(common$GDSC))
	expect_equal(cellNames(common$CCLE), cellNames(common$GDSC))
	expect_equal(rownames(common$GDSC@sensitivity$profiles), rownames(common$CCLE@sensitivity$profiles))
})