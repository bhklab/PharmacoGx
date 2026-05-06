library(PharmacoGx)

context("[` method validation for PharmacoSet")

test_that("numeric indices subset correctly", {
  data(CCLEsmall)
  ps <- CCLEsmall[1, 1]
  expect_equal(sampleNames(ps), sampleNames(CCLEsmall)[1])
  expect_equal(treatmentNames(ps), treatmentNames(CCLEsmall)[1])
})

test_that("single missing dimension subsets all entries", {
  data(CCLEsmall)
  cell_subset <- CCLEsmall[1, ]
  expect_equal(sampleNames(cell_subset), sampleNames(CCLEsmall)[1])
  expect_equal(treatmentNames(cell_subset), treatmentNames(CCLEsmall))

  drug_subset <- CCLEsmall[, 1]
  expect_equal(sampleNames(drug_subset), sampleNames(CCLEsmall))
  expect_equal(treatmentNames(drug_subset), treatmentNames(CCLEsmall)[1])
})

test_that("integer(0) selections return empty subsets", {
  data(CCLEsmall)
  empty_cells <- CCLEsmall[integer(0), 1]
  expect_equal(length(sampleNames(empty_cells)), 0)
  expect_equal(treatmentNames(empty_cells), treatmentNames(CCLEsmall)[1])

  empty_drugs <- CCLEsmall[1, integer(0)]
  expect_equal(sampleNames(empty_drugs), sampleNames(CCLEsmall)[1])
  expect_equal(length(treatmentNames(empty_drugs)), 0)

  empty_both <- CCLEsmall[integer(0), integer(0)]
  expect_equal(length(sampleNames(empty_both)), 0)
  expect_equal(length(treatmentNames(empty_both)), 0)
})

test_that("invalid numeric indices raise informative errors", {
  data(CCLEsmall)
  expect_error(CCLEsmall[-1, 1], "positive")
  expect_error(CCLEsmall[0, 1], "positive")
  expect_error(CCLEsmall[NA_integer_, 1], "must not contain NA")
  expect_error(
    CCLEsmall[length(sampleNames(CCLEsmall)) + 1, 1],
    "exceed"
  )
  expect_error(CCLEsmall[1, NA_integer_], "must not contain NA")
})
