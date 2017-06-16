library(PharmacoGx)

context("Checking unionList function.")

test_that("union List works as union with arbitrary number of arguments",{
	expect_equal(unionList(1,2,3,4,2,2,1), c(1,2,3,4))
	expect_equal(unionList(1,2), c(1,2))
	expect_equal(unionList(1), c(1))
	expect_equal(unionList(), list())
})

test_that("unionList unlists things and unions properly", {

	expect_equal(unionList(list(1,2,3,4,3)), c(1,2,3,4))
	expect_equal(unionList(list(1,2,3), list(4,3)), c(1,2,3,4))
})