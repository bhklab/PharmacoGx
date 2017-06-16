library(PharmacoGx)

context("Testing LogLogisticRegression.")


test_that("Stuff works",{
	
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
    expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60, 70), family = “The Addams Family”)) #should complain
    logLogisticRegression(c(1, 2, 3), c(50, 60, 70)) #should run with no objections
})