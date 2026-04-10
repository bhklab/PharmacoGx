library(PharmacoGx)

context("Checking summarizeSensitivityProfiles LongTable support.")

build_longtable_monotherapy_pset <- function() {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      expr = matrix(
        1,
        nrow = 1,
        ncol = 2,
        dimnames = list("gene1", c("c1", "c2"))
      )
    ),
    colData = S4Vectors::DataFrame(
      sampleid = c("cell1", "cell2"),
      batchid = c("b1", "b1"),
      row.names = c("c1", "c2")
    ),
    rowData = S4Vectors::DataFrame(row.names = "gene1")
  )
  S4Vectors::metadata(se)$annotation <- "rna"

  tre <- CoreGx::TreatmentResponseExperiment(
    rowData = data.frame(
      treatment1id = c("drugA", "drugB"),
      treatment2id = c("", "drugX"),
      stringsAsFactors = FALSE
    ),
    rowIDs = c("treatment1id", "treatment2id"),
    colData = data.frame(
      sampleid = c("cell1", "cell2"),
      stringsAsFactors = FALSE
    ),
    colIDs = "sampleid",
    assays = list(
      profiles = data.frame(
        treatment1id = c("drugA", "drugA", "drugB", "drugB"),
        treatment2id = c("", "", "drugX", "drugX"),
        sampleid = c("cell1", "cell2", "cell1", "cell2"),
        SCORE = c(1, 2, 3, 4),
        stringsAsFactors = FALSE
      )
    ),
    assayIDs = list(
      profiles = c("treatment1id", "treatment2id", "sampleid")
    )
  )

  PharmacoSet2(
    name = "toyLongTable",
    treatment = data.frame(row.names = "drugA"),
    sample = data.frame(
      tissueid = c("t1", "t2"),
      row.names = c("cell1", "cell2")
    ),
    molecularProfiles = MultiAssayExperiment::MultiAssayExperiment(
      experiments = list(rna = se)
    ),
    treatmentResponse = tre,
    curation = list(
      sample = data.frame(),
      treatment = data.frame(),
      tissue = data.frame()
    )
  )
}

test_that("summarizeSensitivityProfiles infers monotherapy treatment keys", {
  pset <- build_longtable_monotherapy_pset()

  res <- summarizeSensitivityProfiles(
    pset,
    sensitivity.measure = "SCORE",
    verbose = FALSE
  )

  expect_true(is.matrix(res))
  expect_identical(rownames(res), "drugA")
  expect_identical(colnames(res), c("cell1", "cell2"))
  expect_equal(unname(res["drugA", ]), c(1, 2))
})

test_that("summarizeSensitivityProfiles explains explicit measure requirements on LongTable datasets", {
  pset <- build_longtable_monotherapy_pset()

  expect_error(
    summarizeSensitivityProfiles(pset, verbose = FALSE),
    regexp = paste(
      "explicit sensitivity.measure",
      "SCORE",
      "sProfiles",
      sep = ".*"
    )
  )
})
