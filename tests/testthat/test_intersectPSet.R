library(PharmacoGx)
require(parallel)
context("Checking intersectPSet.")

test_that("Intersection result did not change since last time", {

	data(GDSCsmall)
	data(CCLEsmall)
	# cell <- intersect(cellNames(GDSCsmall), cellNames(CCLEsmall))[[1]]
	# CCLEsmall2 <- subsetTo(CCLEsmall, cells=cell)
	# GDSCsmall2 <- subsetTo(GDSCsmall, cells=cell)
	# common <- intersectPSet(list('GDSC'=GDSCsmall2, 'CCLE'=CCLEsmall2), intersectOn = c("drugs", "cell.lines", "concentrations"))
	# common <- intersectPSet(list('GDSC'=GDSCsmall, 'CCLE'=CCLEsmall), intersectOn = c("drugs", "cell.lines","concentrations"))	
	expect_equal_to_reference(intersectPSet(list('GDSC'=GDSCsmall, 'CCLE'=CCLEsmall), intersectOn = c("drugs", "cell.lines","concentrations")), "intersectedSmallData.rds")
	# expect_equal(sum(!is.na(common$CCLE@sensitivity$profiles$auc_recomputed_star)),sum(!is.na(common$CCLE@sensitivity$profiles$auc_recomputed_star)))
	# expect_equal(drugNames(common$CCLE), drugNames(common$GDSC))
	# expect_equal(cellNames(common$CCLE), cellNames(common$GDSC))

})