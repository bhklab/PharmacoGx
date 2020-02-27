##' Validate PSet molecularProfiles Conversion
##' 
##' Checks that all the information contained in an ExpressionSet molecularProfile 
##'   was successfully tranferred to the SummarizedExperiment molecularProfile
##'   
##' @param pSet [S4] a PSet containing molecularProfiles as SummarizedExperiments
##' @param pSet [S4] a PSet containing molecularProfiles as ExpressionSets
##' 
##' @return [message] Any slots which are not the same
##' 
##' @importFrom testthat expect_equal test_that
##' 
##' @export
##'
validatePsetMolecularProfilesToSEConversion <- function(pSet_old, pSet_new) {
  
  # Testing that pSets are in correct order
  print("Checking is pSet structures are correct")
  
  testthat::expect_true(
    all(vapply(pSet_old@molecularProfiles, function(x) { is(x, "ExpressionSet") }, FUN.VALUE = logical(1))),
    info = "Old pSet doesn't contain ExpressionSet objects, maybe argument order is wrong?"
  )
  
  testthat::expect_true(
    all(vapply(pSet_new@molecularProfiles, function(x) { is(x, "SummarizedExperiment") }, FUN.VALUE = logical(1))),
    info = "New pSet doesn't contain SummarizedExperiment objects, maybe argument order is wrong?"
  )
  
  # Comparing molecularProfiles slot data  
  print("Checking molecularProfiles slots hold equivalent data.")
  
    for (i in seq_len(length(pSet_old@molecularProfiles))) {
      testthat::expect_true(
        all(
          exprs(pSet_old@molecularProfiles[[i]]) == 
          assays(pSet_new@molecularProfiles[[i]])[[1]],
          na.rm = TRUE
        ),
        info = "The assay data is not equivalent"
      )
    }
    ## TODO:: Rewrite this as an apply statement
    for (i in seq_len(length(pSet_old@molecularProfiles))) { # Have to compare like this due to NAs in data
      # Checking phenoData
      testthat::expect_true(
      if (nrow(pData(pSet_old@molecularProfiles[[i]])) > 0) {
        all(
          as(pSet_old@molecularProfiles[[i]]@phenoData, "data.frame") == 
            as.data.frame(pSet_new@molecularProfiles[[i]]@colData[
              seq_len(length(pSet_new@molecularProfiles[[i]]@colData) -1)]),
          na.rm = TRUE)
      } else { TRUE },
          info = "The phenoData is not equivalent",
        )
      # Checking featureData
      testthat::expect_true(
        if (nrow(fData(pSet_old@molecularProfiles[[i]])) > 0) {
          all(
            as(pSet_old@molecularProfiles[[i]]@featureData, "data.frame") == 
              as.data.frame(pSet_new@molecularProfiles[[i]]@elementMetadata[
                seq_len(length(pSet_new@molecularProfiles[[i]]@elementMetadata) -1)]),
            na.rm = TRUE)
        } else { TRUE },
        info = "The featureData is not equivalent",
      )
      # Checking protocolData
      testthat::expect_true(
        all(
          as(pSet_old@molecularProfiles[[i]]@protocolData, "data.frame") ==
            as(pSet_new@molecularProfiles[[i]]@metadata$protocolData, "data.frame"),
          na.rm = TRUE),
        info = "The protocolData is not equivalent"
      )
    }
      
    testthat::expect_equal(
      lapply(pSet_old@molecularProfiles, function(x) { x@annotation }), 
      lapply(pSet_new@molecularProfiles, function(x) { x@metadata$annotation }),
      info = "The annotation is not equivalent"
    )
      
    testthat::expect_equal(
      lapply(pSet_old@molecularProfiles, function(x) { x@experimentData }), 
      lapply(pSet_new@molecularProfiles, function(x) { x@metadata$experimentData }),
      info = "The experimentData is not equivalent"
    )
    
    testthat::expect_equal(
      lapply(pSet_old@molecularProfiles, function(x) { x@.__classVersion__ }), 
      lapply(pSet_new@molecularProfiles, function(x) { x@metadata$.__classVersion__}),
      info = "The .__classVersion__ is not equivalent"
    )
    
  # Comparing remainder of pSet slots; should not be affect by conversion
  print("Comparing remainder of pSet slots")
    
  testthat::test_that("Checking pSet@annotation slot is unchanged.", {
    testthat::expect_equal(pSet_old@annotation, pSet_new@annotation)
  })
  
  testthat::test_that("Checking pSet@cell slot is unchanged.", {
    testthat::expect_equal(pSet_old@cell, pSet_new@cell)
  })
  
  testthat::test_that("Checking pSet@drug slot is unchanged.", {
    testthat::expect_equal(pSet_old@drug, pSet_new@drug)
  })
  
  testthat::test_that("Checking pSet@sensitivity slot is unchanged.", {
    testthat::expect_equal(pSet_old@sensitivity, pSet_new@sensitivity)
  })
  
  testthat::test_that("Checking pSet@datasetType slot is unchanged.", {
    testthat::expect_equal(pSet_old@datasetType, pSet_new@datasetType)
  })
  
  testthat::test_that("Checking pSet@perturbation slot is unchanged.", {
    testthat::expect_equal(pSet_old@perturbation, pSet_new@perturbation)
  })
  
  testthat::test_that("Checking pSet@curation slot is unchanged.", {
    testthat::expect_equal(pSet_old@curation, pSet_new@curation)
  })
}