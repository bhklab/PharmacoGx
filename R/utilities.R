#' PSet molecularProfiles from ESets to SEs
#'
#' Converts all ExpressionSet objects within the molecularProfiles slot of a 
#'   PharmacoSet to SummarizedExperiments
#'
#' @param pSet \code{S4} A PharmacoSet containing molecular data in ExpressionSets
#'
#' @return \code{S4} A PharmacoSet containing molecular data in a SummarizedExperiments
#' 
#' @importFrom parallel mclapply
#' @importFrom SummarizedExperiment assay assays assayNames
#' @importClassesFrom SummarizedExperiment SummarizedExperiment Assays
#' @importFrom Biobase exprs fData pData annotation protocolData
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom stats setNames
#' 
#' @export
#' @keywords internal
.convertPsetMolecularProfilesToSE <- function(pSet) {

  if (!is.null(pSet@annotation$version) && pSet@annotation$version >= 2 ){
    return(pSet)
  }

  eSets <- pSet@molecularProfiles # Extract eSet data
  
  pSet@molecularProfiles <-
    lapply(eSets,
           function(eSet){
             
             # Change rownames from probes to EnsemblGeneId for rna data type
             if (grepl("^rna$", Biobase::annotation(eSet))) {
               rownames(eSet) <- Biobase::fData(eSet)$EnsemblGeneId
             }
             
             # Build summarized experiment from eSet
             SE <- SummarizedExperiment::SummarizedExperiment(
               ## TODO:: Do we want to pass an environment for better memory efficiency?
               assays=SimpleList(as.list(Biobase::assayData(eSet))
               ),
               # Switch rearrange columns so that IDs are first, probes second
               rowData=S4Vectors::DataFrame(Biobase::fData(eSet),
                                            rownames=rownames(Biobase::fData(eSet)) 
               ),
               colData=S4Vectors::DataFrame(Biobase::pData(eSet),
                                            rownames=rownames(Biobase::pData(eSet))
               ),
               metadata=list("experimentData" = eSet@experimentData, 
                             "annotation" = Biobase::annotation(eSet), 
                             "protocolData" = Biobase::protocolData(eSet)
               )
             )
             ## TODO:: Determine if this can be done in the SE constructor?
             # Extract names from expression set
             SummarizedExperiment::assayNames(SE) <- Biobase::assayDataElementNames(eSet)
             # Assign SE to pSet
             mDataType <- Biobase::annotation(eSet)
             pSet@molecularProfiles[[mDataType]] <- SE
           })
  setNames(pSet@molecularProfiles, names(eSets))
  pSet@annotation$version <- 2
  pSet
}

##' Validate PSet molecularProfiles Conversion
##' 
##' Checks that all the information contained in an ExpressionSet molecularProfile 
##'   was successfully tranferred to the SummarizedExperiment molecularProfile
##'   
##' @param pSet_new \code{S4} a PSet containing molecularProfiles as SummarizedExperiments
##' @param pSet_old \code{S4} a PSet containing molecularProfiles as ExpressionSets
##' 
##' @return \code{message} Any slots which are not the same
##' 
##' @importFrom testthat expect_equal test_that
##' @import SummarizedExperiment
##' @import Biobase
##' 
##' @export
##' @keywords internal
.validatePsetMolecularProfilesToSEConversion <- function(pSet_old, pSet_new) {
  
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
    for (j in seq_along(assays(pSet_new@molecularProfiles[[i]]))) {
      testthat::expect_true(
        all(
          as.list(assayData(pSet_old@molecularProfiles[[i]]))[[j]] == 
            assay(pSet_new@molecularProfiles[[i]], j),
          na.rm = TRUE
        ),
        info = "The assay data is not equivalent"
      )
    }
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
  
  ##TODO:: Removed .__classVersion__ from SE as it is a property specific to eSet
  # testthat::expect_equal(
  #   lapply(pSet_old@molecularProfiles, function(x) { x@.__classVersion__ }), 
  #   lapply(pSet_new@molecularProfiles, function(x) { x@metadata$.__classVersion__}),
  #   info = "The .__classVersion__ is not equivalent"
  # )
  
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
  message("Tests pass!")
}

##TODO:: Determine why CCLEsmall is 3x larger in memory after conversion?
#' Utility function to resave all datasets after modifying convertPSetMolecularProfiles
#' 
#' Converts all example dastasets specificed as an argument from 
#'   molecularProfiles as ExpressionSet to molecularProfiles as 
#'   SummarizedExperiment and saves them in the data folder
#' 
#' @param datasets \code{character} A list of the example datasets to update
#' 
#' @return \code{none} Works by side effects alone to resave all example 
#'   datasets in a package to have SummarizedExperiments for molecularProfiles
#' 
#' @export
#' @keywords internal
.resaveAllExampleDatasets <- function(datasets) {
  for (dataset in datasets) {
    dataDir <- paste0(grep('data', list.dirs(), value=TRUE))
    load(paste0(dataDir, '/', dataset, '_old.rda'))
    assign(dataset, .convertPsetMolecularProfilesToSE(get(dataset)))
    save(list=dataset, file=paste0(dataDir, '/', dataset, '.rda'), compress='xz')
  }
}
