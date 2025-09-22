library(PharmacoGx)
library(checkmate)
context("Testing Download Functionality")

test_that("AvailablePset Works as Expected", {
  canonical_df <- availablePSets(canonical = TRUE)
  print(canonical_df)
  checkmate::expect_data_frame(
    canonical_df,
    types = c(
      "character",
      "character",
      "character",
      "character",
      "list",
      "character",
      "character"
    ),
    ncols = 8,
    min.rows = 1,
  )

  required_cols = c(
    "Dataset Name",
    "Date Created",
    "PSet Name",
    "DOI",
    "Download"
  )

  # the names of the df should include the required columns
  expect_true(all(required_cols %in% colnames(canonical_df)))

  non_cano_df <- availablePSets(canonical = FALSE)

  # non canonical should have more rows than canonical
  expect_true(nrow(non_cano_df) > nrow(canonical_df))
})
