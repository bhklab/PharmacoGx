library(PharmacoGx)

context("Testing LogLogisticRegression.")


test_that("Errors are checked.",{
	
    expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60))) #should complain
    expect_warning(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), viability_as_pct = FALSE)) #should complain
    expect_error(logLogisticRegression(c(-1, 2, 3), c(50, 60, 70), conc_as_log = FALSE)) #should complain
    expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), median_n = 0)) #should complain
    expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), median_n = 3/2)) #should complain
    expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), density = c(1, 1))) #should complain
    expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), density = c(1, 1, -1))) #should complain
    expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), precision = 0)) #should complain
    expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), scale = 0)) #should complain
    expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), lower_bounds = c(0, 0, 0), upper_bounds = c(1, 1, -1))) #should complain
    expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), family = "The Addams Family")) #should complain
    
})

test_that("Values returned as expected (previous runs of function).",{
	
	expect_equivalent(
		logLogisticRegression(seq(-10,10,0.1), .Hill(seq(-10,10,0.1), c(1,0,0))
			, conc_as_log=TRUE, viability_as_pct = FALSE), list(1,0,0))
	expect_equivalent(
		logLogisticRegression(seq(-10,10,0.1), .Hill(seq(-10,10,0.1), c(1,0,0))
			, conc_as_log=TRUE, viability_as_pct = FALSE, family="Cauchy"), list(1,0,0))
	expect_equal(logLogisticRegression(c(1, 2, 3), c(50, 60, 70)), list(HS = 1.22388251239748, E_inf = 59.9999882231278, EC50 = 1e-06)) #should run with no objections
    expect_equal(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), family="Cauchy"), list(HS = 1.08959673823682, E_inf = 60.0004126580051, EC50 = 3.39074306428225e-06)) #should run with no objections
	expect_equal(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), trunc=FALSE), list(HS = 1.22388251239748, E_inf = 59.9999882231278, EC50 = 1e-06)) #should run with no objections
    expect_equal(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), family="Cauchy", trunc=FALSE), list(HS = 1.08959673823682, E_inf = 60.0004126580051, EC50 = 3.39074306428225e-06)) #should run with no objections
})