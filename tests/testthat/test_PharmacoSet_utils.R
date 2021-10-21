library(PharmacoGx)
library(testthat)
data(CCLEsmall)

# -- 
context("Testing PharmacoSet subset methods...")

test_that('subsetByTreatment works...', {
    expect_true({
        treatments <- drugNames(CCLEsmall)[1:5]
        suppressMessages({
            CCLE_sub <- subsetByTreatment(CCLEsmall, treatments)
        })
        all(drugNames(CCLE_sub) %in% treatments)
    })
})

test_that('subsetBySample works...', {
    expect_true({
        samples <- cellNames(CCLEsmall)[1:5]
        suppressMessages({
            CCLE_sub <- subsetBySample(CCLEsmall, samples)
        })
        all(cellNames(CCLE_sub) %in% samples)
    })
})

test_that('subsetByFeature works...', {
    
})
