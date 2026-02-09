library(PharmacoGx)
library(testthat)
data(CCLEsmall)

# --
context("Testing PharmacoSet subset methods...")

test_that('subsetByTreatment works...', {
  expect_true({
    treatments <- treatmentNames(CCLEsmall)[1:5]
    suppressMessages({
      CCLE_sub <- subsetByTreatment(CCLEsmall, treatments)
    })
    all(treatmentNames(CCLE_sub) %in% treatments)
  })
})

test_that('subsetBySample works...', {
  expect_true({
    samples <- sampleNames(CCLEsmall)[1:5]
    suppressMessages({
      CCLE_sub <- subsetBySample(CCLEsmall, samples)
    })
    all(sampleNames(CCLE_sub) %in% samples)
  })
})

test_that('subsetByFeature works...', {
  expect_true({
    features <- head(rownames(featureInfo(CCLEsmall, 'rna')), 5)
    suppressMessages({
      CCLE_sub <- subsetByFeature(
        CCLEsmall,
        features = features,
        mDataTypes = 'rna'
      )
    })
    assay_data <- SummarizedExperiment::assay(molecularProfiles(
      CCLE_sub,
      mDataType = 'rna'
    ))
    identical(sort(rownames(featureInfo(CCLE_sub, 'rna'))), sort(features)) &&
      nrow(assay_data) == length(features)
  })
})
