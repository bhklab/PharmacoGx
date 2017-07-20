library(PharmacoGx)

context("Checking computeABC.")

test_that("Function complains when given insensible input",{
	
	expect_error(computeABC(conc1 = c(1, 2, 3),
		conc2 = c(1, 2, 3),
		viability1 = c(50, 60, 70),
		viability2 = c(40, 90, 10),
		Hill_fit1 = c(1, 0, 0.1),
		Hill_fit2 = c(0.5, 0.2, 1)), "Please pass in only one")
	# expect_silent(computeABC(conc1 = c(1, 2, 3),
	# 	conc2 = c(1, 2, 3),
	# 	viability1 = c(50, 60, 70),
	# 	Hill_fit2 = c(0.5, 0.2, 1)))

	expect_error(computeABC(conc1 = c(1, 2),
		conc2 = c(1, 2, 3),
		viability1 = c(50, 60, 70),
		viability2 = c(40, 90, 10)), "is not of same length") #should complain
	expect_error(computeABC(conc1 = c(-1, 2, 3),
		conc2 = c(1, -2, 3),
		viability1 = c(50, 60, 70),
		viability2 = c(40, 90, 10),
		conc_as_log = FALSE)) #should complain
	expect_error(computeABC(conc1 = c(NA, "cat", 3),
		conc2 = c(1, -2, 3),
		viability1 = c(50, 60, 70),
		viability2 = c(40, 90, 10),
		conc_as_log = FALSE)) #should complain
	expect_error(computeABC(conc1 = c(1, 2, 3),
		conc2 = c(1, -2, 3),
		viability1 = c(50, 60, 70),
		viability2 = c(40, 90, 10),
		verbose = NA)) #should complain
	expect_error(computeABC(conc1 = c(1, 2, 3),
		conc2 = c(1, -2, 3),
		viability1 = c(50, 60, 70))) #should complain
	expect_error(computeABC(conc1 = c(1, 2, Inf),
		conc2 = c(1, -2, 3),
		viability1 = c(50, 60, 70),
		viability2 = c(40, 90, 10))) #should complain
	expect_error(computeABC(conc1 = c(1, 2, 3),
		conc2 = c(1, -2, 3),
		viability1 = c(.50, .60, .70),
		viability2 = c(.40, .90, .10),
		viability_as_pct = TRUE)) #should complain
	expect_warning(computeABC(conc1 = c(1, 2, 3),
		conc2 = c(1, 2, 3),
		viability1 = c(.50, .60, .70),
		viability2 = c(.40, .90, .10),
		viability_as_pct = TRUE)) #should complain
	expect_error(computeABC()) #should complain
})

test_that("Function values make sense",{
	expect_equal(computeABC(conc1=c(-1,0,1), conc2=c(-1,0,1), Hill_fit1=c(0,1,0), Hill_fit2=c(1,0,0), conc_as_log=TRUE, viability_as_pct=FALSE), 0.5)
	expect_equal(computeABC(conc1=c(-1,0,1), conc2=c(-1,0,1), Hill_fit1=c(0,1,0), Hill_fit2=c(0,1,0), conc_as_log=TRUE, viability_as_pct=FALSE), 0)
})

