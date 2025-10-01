library(PharmacoGx)

context("summarizeMolecularProfiles missing ucell handling")

matching_logic <- function(pp, dd, ucell) {
  ucell_idx <- match(ucell, pp[, "sampleid"])
  missing_ucell <- is.na(ucell_idx) & !is.na(ucell)
  if (any(missing_ucell)) {
    alt_idx <- match(ucell[missing_ucell], rownames(pp))
    ucell_idx[missing_ucell] <- alt_idx
  }
  dd_idx <- match(pp[ucell_idx, "sampleid"], colnames(dd))
  if (any(missing_ucell)) {
    alt_dd <- match(ucell[missing_ucell], colnames(dd))
    dd_idx[missing_ucell] <- alt_dd
  }
  list(ucell_idx = ucell_idx, dd_idx = dd_idx)
}

test_that("missing cells absent everywhere remain NA", {
  pp <- data.frame(sampleid = c("A", "B"))
  rownames(pp) <- c("rowA", "rowB")
  dd <- matrix(1:4, nrow = 2, dimnames = list(NULL, c("C", "D")))
  res <- matching_logic(pp, dd, c("Z"))
  expect_true(all(is.na(res$ucell_idx)))
  expect_true(all(is.na(res$dd_idx)))
})

test_that("duplicate rownames fallback picks first occurrence", {
  pp <- matrix(c(NA, NA, "C"), ncol = 1, dimnames = list(c("dup", "dup", "other"), "sampleid"))
  dd <- matrix(1:9, nrow = 3, dimnames = list(NULL, c("X", "Y", "C")))
  res <- matching_logic(pp, dd, c("dup"))
  expect_equal(res$ucell_idx, 1L)
  expect_true(is.na(res$dd_idx))
})

test_that("fallback to dd column names populates dd_idx", {
  pp <- data.frame(sampleid = c("A", NA))
  rownames(pp) <- c("rowA", "rowB")
  dd <- matrix(1:4, nrow = 2, dimnames = list(NULL, c("A", "only_dd")))
  res <- matching_logic(pp, dd, c("only_dd"))
  expect_true(is.na(res$ucell_idx))
  expect_equal(res$dd_idx, 2L)
})
