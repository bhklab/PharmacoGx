local_mcp_environment <- function() {
  skip_if_not_installed("ellmer")
  env <- new.env(parent = globalenv())
  source(
    system.file("mcp", "PharmacoGx-tools.R", package = "PharmacoGx"),
    local = env
  )
  env
}

test_that("legacy expression assays map to biological samples and aggregate", {
  env <- local_mcp_environment()
  mat <- matrix(
    c(1, NA, 3, Inf, 20, NA, 8, 4),
    nrow = 2,
    dimnames = list(c("gene1", "gene2"), paste0("assay", 1:4))
  )
  profile <- SummarizedExperiment::SummarizedExperiment(
    assays = list(exprs = mat),
    colData = S4Vectors::DataFrame(cellid = c("A", "A", "A", "B"))
  )
  pset <- structure(
    list(),
    cell = data.frame(cellid = c("A", "B")),
    molecularProfiles = list(rna = profile)
  )
  result <- env$pgx_summarize_molecular_profiles(
    pset,
    "rna",
    c("gene2", "gene1"),
    c("B", "A", "absent"),
    "median"
  )
  expect_identical(dimnames(result), list(c("gene2", "gene1"), c("B", "A")))
  expect_equal(unname(result), matrix(c(4, 8, NA, 3), nrow = 2))
  mean_result <- env$pgx_summarize_molecular_profiles(
    pset,
    "rna",
    "gene1",
    "A",
    "mean"
  )
  expect_equal(unname(mean_result[1, 1]), 8)
  first <- env$pgx_summarize_profile_by_sample(profile, "gene1", "A", "first")
  last <- env$pgx_summarize_profile_by_sample(profile, "gene1", "A", "last")
  expect_equal(unname(first[1, 1]), 1)
  expect_equal(unname(last[1, 1]), 20)
  expect_error(env$pgx_summarize_profile_by_sample(
    profile,
    "gene1",
    "A",
    "sum"
  ))
  expect_equal(
    dim(env$pgx_summarize_profile_by_sample(
      profile,
      "absent",
      "A",
      "median"
    )),
    c(0L, 1L)
  )
  expect_equal(
    dim(env$pgx_summarize_profile_by_sample(
      profile,
      "gene1",
      "absent",
      "median"
    )),
    c(1L, 0L)
  )
})

test_that("sampleid and already sample-labelled matrices remain supported", {
  env <- local_mcp_environment()
  mat <- matrix(c(2, 6), nrow = 1, dimnames = list("g", c("A", "B")))
  expect_equal(
    env$pgx_summarize_profile_by_sample(mat, "g", c("B", "A"), "mean"),
    mat[, c("B", "A"), drop = FALSE]
  )
  profile <- SummarizedExperiment::SummarizedExperiment(
    assays = list(exprs = mat),
    colData = S4Vectors::DataFrame(sampleid = c("X", "Y"))
  )
  result <- env$pgx_summarize_profile_by_sample(
    profile,
    "g",
    c("Y", "X"),
    "mean"
  )
  expect_equal(unname(result), matrix(c(6, 2), nrow = 1))
  expect_identical(colnames(result), c("Y", "X"))
  SummarizedExperiment::colData(profile)$sampleid[1] <- NA_character_
  expect_error(
    env$pgx_summarize_profile_by_sample(profile, "g", "Y", "mean"),
    "missing biological sample identifiers"
  )
})

test_that("legacy binary aggregation preserves missingness and rejects nonbinary values", {
  env <- local_mcp_environment()
  mat <- matrix(
    c("0", NA, "1", NA),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("a1", "a2"))
  )
  profile <- SummarizedExperiment::SummarizedExperiment(
    assays = list(exprs = mat),
    colData = S4Vectors::DataFrame(cellid = c("A", "A"))
  )
  expect_equal(
    unname(env$pgx_summarize_profile_by_sample(
      profile,
      c("g1", "g2"),
      "A",
      "or"
    )),
    matrix(c(1, NA), ncol = 1)
  )
  expect_equal(
    unname(env$pgx_summarize_profile_by_sample(
      profile,
      "g1",
      "A",
      "and"
    )),
    matrix(0, ncol = 1)
  )
  SummarizedExperiment::assay(profile)[1, 1] <- "2"
  expect_error(
    env$pgx_summarize_profile_by_sample(profile, "g1", "A", "or"),
    "binary 0/1"
  )
})

test_that("MCP correlations use finite pairs and explicit asymptotic inference", {
  env <- local_mcp_environment()
  x <- c(1, 2, 3, 4, 5, 6, Inf, NA)
  y <- c(1, 2, 4, 3, 6, 5, 7, 8)
  samples <- paste0("s", seq_along(x))
  mat <- matrix(x, nrow = 1, dimnames = list("gene", samples))
  result <- env$pgx_feature_response_association(
    mat,
    setNames(y, samples),
    "spearman",
    "BH"
  )
  reference <- cor.test(x[1:6], y[1:6], method = "spearman", exact = FALSE)
  expect_equal(result$n, 6L)
  expect_equal(result$effect_size, unname(reference$estimate))
  expect_equal(result$p_value, reference$p.value)
  expect_identical(result$p_value_method, "asymptotic (exact=FALSE)")
  constant <- env$pgx_feature_response_association(
    mat,
    setNames(rep(1, 8), samples),
    "spearman",
    "BH"
  )
  expect_true(is.na(constant$p_value))
  expect_true(is.na(constant$effect_size))
})
