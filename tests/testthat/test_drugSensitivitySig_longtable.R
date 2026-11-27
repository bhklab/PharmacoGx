make_signature_response_fixture <- function(layout) {
  data(CCLEsmall, envir = environment())
  pset <- CCLEsmall
  drug <- treatmentNames(pset)[1L]
  cells <- head(intersect(sampleNames(pset), phenoInfo(pset, "rna")$sampleid), 30L)
  features <- head(rownames(featureInfo(pset, "rna")), 3L)
  values <- seq_along(cells) / (length(cells) + 1L)
  samples <- data.frame(sampleid = cells)
  scores <- values[c(seq(2L, length(cells), 2L), seq(1L, length(cells), 2L))]
  if (identical(layout, "compact")) {
    response <- CoreGx::TreatmentResponseExperiment(
      rowData = data.frame(treatmentid = drug, dose = c(-7, -6)),
      rowIDs = c("treatmentid", "dose"),
      colData = samples,
      colIDs = "sampleid",
      assays = list(
        profiles = data.frame(
          treatmentid = drug,
          sampleid = cells,
          aac_recomputed = values
        )
      ),
      assayIDs = list(profiles = c("treatmentid", "sampleid"))
    )
  } else {
    # Transitional LongTable: the monotherapy arm uses treatment1id and the
    # stored response retains a historical metric name, next to a combo arm.
    response <- CoreGx::LongTable(
      rowData = data.frame(treatment1id = drug, treatment2id = c("", "partner")),
      rowIDs = c("treatment1id", "treatment2id"),
      colData = samples,
      colIDs = "sampleid",
      assays = list(
        profiles = data.frame(
          treatment1id = drug,
          treatment2id = rep(c("", "partner"), each = length(cells)),
          sampleid = rep(cells, 2L),
          auc_recomputed = c(values, 1 - values),
          SCORE = c(scores, 1 - scores)
        )
      ),
      assayIDs = list(profiles = c("treatment1id", "treatment2id", "sampleid"))
    )
  }
  # Only replace response storage in this fixture. The public setter does not
  # convert between list and LongTable layouts.
  methods::slot(pset, "treatmentResponse") <- response
  stopifnot(methods::validObject(pset))
  expected <- function(values) {
    matrix(values, nrow = 1L, dimnames = list(drug, cells))
  }
  list(
    args = list(
      object = pset, mDataType = "rna", drugs = drug, cells = cells,
      features = features, nthread = 1L, verbose = FALSE
    ),
    expected = expected(values),
    expected_score = expected(scores)
  )
}

test_that("drug signatures read stored metrics from compact LongTable assays", {
  # The compact fixture omits the dose row key. Only CoreGx >= 2.17.1
  # can construct it; the devel CI job checks with that version.
  skip_if_not_installed("CoreGx", minimum_version = "2.17.1")
  fixture <- make_signature_response_fixture("compact")
  response <- treatmentResponse(fixture$args$object)
  expect_true("dose" %in% CoreGx::rowIDs(response))
  expect_false("dose" %in% CoreGx::assayCols(response, "profiles"))
  direct <- do.call(drugSensitivitySig, fixture$args)
  supplied <- do.call(
    drugSensitivitySig, c(fixture$args, list(sProfiles = fixture$expected))
  )
  expect_equal(as(direct, "array"), as(supplied, "array"))
  expect_true(all(is.finite(as(direct, "array")[, , "estimate"])))
  expect_error(
    do.call(
      drugSensitivitySig, c(fixture$args, list(sensitivity.measure = "dose"))
    ),
    "Invalid sensitivity measure.*aac_recomputed"
  )
})

test_that("drug signatures keep transitional measures and monotherapies distinct", {
  fixture <- make_signature_response_fixture("transitional")
  direct_alias <- do.call(drugSensitivitySig, fixture$args)
  direct_score <- do.call(
    drugSensitivitySig, c(fixture$args, list(sensitivity.measure = "SCORE"))
  )
  legacy <- do.call(
    drugSensitivitySig, c(fixture$args, list(sProfiles = fixture$expected))
  )
  explicit_score <- do.call(
    drugSensitivitySig,
    c(fixture$args, list(sProfiles = fixture$expected_score))
  )
  expect_equal(as(direct_alias, "array"), as(legacy, "array"))
  expect_equal(as(direct_score, "array"), as(explicit_score, "array"))
  expect_false(isTRUE(all.equal(
    as(direct_score, "array"), as(direct_alias, "array")
  )))
  expect_true(all(is.finite(as(direct_score, "array")[, , "estimate"])))
})
