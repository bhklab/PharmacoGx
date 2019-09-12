##' Validate PSet molecularProfiles Conversion
##' 
##' Checks that all the information contained in an ExpressionSet molecularProfile 
##'   was successfully tranferred to the SummarizedExperiment molecularProfile
##'   
##' @param pSet [S4] a PSet containing molecularProfiles as SummarizedExperiments
##' @param pSet [S4] a PSet containing molecularProfiles as ExpressionSets
##' 
##' @return [string] Any slots which are not the same
##' 
##' @importFrom testthat expect_equal context
##' 
##' @export
##' 
#
#validatePsetMolecularProfilesToSEConversion <- function(pSet_old, pSet_new) 
#{
#  test_that::context("Checking pSet names are unchanged.") {
#    expect_equal(pSet_old@names, pSet_new@annotation)
#    
#  test_that::context("Checking cellInfo is unchanged.") {
#    testthat::expect_equal()
#  }
#  
#  test_that::context("Checking cellInfo is unchanged.") {
#    testthat::expect_equal()
#  }
#  
#  test_that::context("Checking cellInfo is unchanged.") {
#    testthat::expect_equal()
#  }
#  
#  test_that::context("Checking cellInfo is unchanged.") {
#    testthat::expect_equal()
#  }
#  
#  test_that::context("Checking cellInfo is unchanged.") {
#    testthat::expect_equal()
#  }
#  
#  test_that::context("Checking cellInfo is unchanged.") {
#    testthat::expect_equal()
#  }
#  
#  test_that::context("Checking cellInfo is unchanged.") {
#    testthat::expect_equal()
#  }
#}#