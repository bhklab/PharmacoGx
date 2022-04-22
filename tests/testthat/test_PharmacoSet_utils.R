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

})
