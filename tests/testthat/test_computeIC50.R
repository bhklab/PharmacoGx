library(PharmacoGx)

context("Checking computeIC50/ICn.")

test_that("Function complains when given insensible input",{
	expect_error(computeIC50(concentration = c(1, 2, 3),
		viability = c(50, 60, 70),
		Hill_fit = c(1, 0, 0.1)), "Please pass in only one")
	# expect_silent(computeIC50(concentration = c(1, 2, 3),
	# 		# 	viability1 = c(50, 60, 70),
	# 	Hill_fit2 = c(0.5, 0.2, 1)))

	expect_error(computeIC50(concentration = c(1, 2, 3, 5), viability = c(50, 60, 70)), "is not of same length") #should complain
	expect_error(computeIC50(concentration = c(-1, 2, 3),viability = c(50, 60, 70),conc_as_log = FALSE),"'x_as_log' flag may be set incorrectly") #should complain
	##TO-DO:: Add wanring strings to expect_warning call
	expect_error(computeIC50(concentration = c(NA, "cat", 3), viability = c(50, 60, 70), conc_as_log = FALSE), "real numbers") #should complain
	expect_error(computeIC50(concentration = c(1, 2, Inf), viability = c(50, 60, 70)), "real numbers, NA-values, and/or -Inf") #should complain
	expect_warning(computeIC50(concentration = c(1, 2, 3),
		viability = c(.50, .60, .70),
		viability_as_pct = TRUE), "as_pct") #should complain
	expect_error(computeIC50()) #should complain
})

test_that("Functions return right values",{
	expect_equal(computeIC50(concentration = seq(-3,3), Hill_fit=c(1,0,0), conc_as_log=TRUE, viability_as_pct=FALSE), 0)
	expect_equal(computeIC50(concentration = seq(1,3), Hill_fit=c(1,0,0), conc_as_log=TRUE, viability_as_pct=FALSE), 0)
	expect_equal(computeIC50(concentration = seq(1,3), Hill_fit=c(1,.9,0), conc_as_log=TRUE, viability_as_pct=FALSE), NA_real_)
	expect_equal(computeIC50(concentration = seq(1,3), Hill_fit=c(1,.5,0), conc_as_log=TRUE, viability_as_pct=FALSE), Inf)
	expect_equal(.Hill(computeICn(concentration = seq(1,3), Hill_fit=c(1,0,0), n=.7, conc_as_log=TRUE, viability_as_pct=FALSE), c(1,0,0)), .3)
	expect_equal(computeICn(concentration = seq(1,3), Hill_fit=c(1,0,0), n=0, conc_as_log=TRUE, viability_as_pct=FALSE), -Inf)
})
