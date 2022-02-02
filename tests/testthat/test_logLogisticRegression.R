library(PharmacoGx)

context("Testing LogLogisticRegression.")

##TO-DO::Supress print to console from this test file

test_that("Errors are checked.",{
	
    expect_error(logLogisticRegression(c(1, 2, 3), c(50, 60))) #should complain
    expect_warning(logLogisticRegression(c(1, 2, 3), c(70, 60, 50), viability_as_pct = FALSE)) #should complain
    expect_error(logLogisticRegression(c(-1, 2, 3), c(70, 60, 50), conc_as_log = FALSE)) #should complain

    expect_error(logLogisticRegression(c(1, 2, 3), c(70, 60, 50), median_n = 0)) #should complain
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

	expect_equal(logLogisticRegression(c(1, 2, 3), c(70, 60, 50)), 
        structure(list(HS = 0.766235008089789, E_inf = 0, EC50 = 3.14114781624406), Rsquare = 0.983872720888105), 
        tolerance=1e-3) #should run with no objections

    expect_equal(logLogisticRegression(c(1, 2, 3), c(70, 60, 50), family="Cauchy"), 
        structure(list(HS = 0.766099642757164, E_inf = 0, EC50 = 3.13883130512133), Rsquare = 0.983863295704847), 
        tolerance=1e-3) #should run with no objections

	expect_equal(logLogisticRegression(c(1, 2, 3), c(70, 60, 50), trunc=FALSE), 
        structure(list(HS = 0.766234869170442, E_inf = 0, EC50 = 3.14114817741207), Rsquare = 0.983872721139933), 
        tolerance=1e-3) #should run with no objections

    expect_equal(logLogisticRegression(c(1, 2, 3), c(70, 60, 50), family="Cauchy", trunc=FALSE), 
        structure(list(HS = 0.766099642757164, E_inf = 0, EC50 = 3.13883130512133), Rsquare = 0.983863295704847),
        tolerance=1e-3) #should run with no objections

    ## These next few tests make sure trunc is doing something sensible
    expect_equal(logLogisticRegression(c(0.1, 1, 2, 3), c(110, 70, 60, 50), family="Cauchy", trunc=FALSE), 
        structure(list(HS = 1.83941027802297, E_inf = 46.2252841409534, EC50 = 0.929240163785174), Rsquare = 0.950154102951421),
        tolerance=1e-3) #should run with no objections
    
    expect_equal(logLogisticRegression(c(0.1, 1, 2, 3), c(110, 70, 60, 50), family="Cauchy", trunc=TRUE), 
        structure(list(HS = 2.06741101065827, E_inf = 48.0764684303728, EC50 = 0.900808050726654), Rsquare = 0.986745954925997),
        tolerance=1e-3) #should run with no objections

    expect_equivalent(logLogisticRegression(c(0.1, 1, 2, 3), c(100, 70, 60, 50), trunc=TRUE), 
        logLogisticRegression(c(0.1, 1, 2, 3), c(500, 70, 60, 50), trunc=TRUE))

})
