library(PharmacoGx)

context("drugSensitivitySig defaults")

test_that("drugSensitivitySig works with explicit optional arguments", {
  data(CCLEsmall)
  res <- drugSensitivitySig(
    CCLEsmall,
    mDataType = "rna",
    features = rownames(featureInfo(CCLEsmall, "rna"))[1],
    sensitivity.measure = "aac_recomputed",
    molecular.summary.stat = "mean",
    sensitivity.summary.stat = "mean",
    returnValues = "estimate",
    modeling.method = "anova",
    inference.method = "analytic",
    parallel.on = "drug",
    nthread = 1,
    verbose = FALSE
  )
  expect_s4_class(res, "PharmacoSig")
  expect_equal(dim(res), c(1, length(treatmentNames(CCLEsmall)), 8))
})

test_that("drugSensitivitySig warns on legacy AUC sensitivity aliases", {
  data(CCLEsmall)
  expect_warning(
    drugSensitivitySig(
      CCLEsmall,
      mDataType = "rna",
      features = rownames(featureInfo(CCLEsmall, "rna"))[1],
      sensitivity.measure = "auc_recomputed",
      returnValues = "estimate",
      parallel.on = "drug",
      nthread = 1,
      verbose = TRUE
    ),
    "deprecated"
  )
})

test_that("drugSensitivitySig works with default optional arguments", {
  data(CCLEsmall)
  res <- drugSensitivitySig(
    CCLEsmall,
    mDataType = "rna",
    features = rownames(featureInfo(CCLEsmall, "rna"))[1],
    parallel.on = "drug",
    nthread = 1,
    verbose = FALSE
  )
  expect_s4_class(res, "PharmacoSig")
  expect_equal(dim(res), c(1, length(treatmentNames(CCLEsmall)), 8))
})

test_that("drugSensitivitySig ignores invalid sensitivity.measure when sProfiles are supplied", {
  data(CCLEsmall)
  drugs <- treatmentNames(CCLEsmall)[1:2]
  cells <- sampleNames(CCLEsmall)[1:5]
  sProfiles <- summarizeSensitivityProfiles(
    CCLEsmall,
    sensitivity.measure = "aac_recomputed",
    drugs = drugs,
    cell.lines = cells,
    verbose = FALSE
  )

  res <- drugSensitivitySig(
    CCLEsmall,
    mDataType = "rna",
    drugs = drugs,
    cells = cells,
    features = rownames(featureInfo(CCLEsmall, "rna"))[1:2],
    sensitivity.measure = "definitely_not_a_measure",
    sProfiles = sProfiles,
    parallel.on = "drug",
    nthread = 1,
    verbose = FALSE
  )

  expect_s4_class(res, "PharmacoSig")
  expect_equal(dim(res), c(2, length(drugs), 8))
})

test_that("drugSensitivitySig supports bundled non-RNA modalities with supplied sProfiles", {
  data(CCLEsmall)
  drugs <- treatmentNames(CCLEsmall)[1:2]
  cells <- sampleNames(CCLEsmall)[1:10]
  sProfiles <- summarizeSensitivityProfiles(
    CCLEsmall,
    sensitivity.measure = "aac_recomputed",
    drugs = drugs,
    cell.lines = cells,
    verbose = FALSE
  )

  rnaseq_res <- suppressWarnings(drugSensitivitySig(
    CCLEsmall,
    mDataType = "rnaseq",
    drugs = drugs,
    cells = cells,
    features = rownames(featureInfo(CCLEsmall, "rnaseq"))[1:2],
    sProfiles = sProfiles,
    modeling.method = "pearson",
    parallel.on = "drug",
    nthread = 1,
    verbose = FALSE
  ))
  cnv_res <- suppressWarnings(drugSensitivitySig(
    CCLEsmall,
    mDataType = "cnv",
    drugs = drugs,
    cells = cells,
    features = rownames(featureInfo(CCLEsmall, "cnv"))[1:2],
    sProfiles = sProfiles,
    modeling.method = "pearson",
    parallel.on = "drug",
    nthread = 1,
    verbose = FALSE
  ))
  mutation_res <- suppressWarnings(drugSensitivitySig(
    CCLEsmall,
    mDataType = "mutation",
    drugs = drugs,
    cells = cells,
    features = rownames(featureInfo(CCLEsmall, "mutation"))[1:2],
    sProfiles = sProfiles,
    molecular.summary.stat = "or",
    modeling.method = "pearson",
    parallel.on = "drug",
    nthread = 1,
    verbose = FALSE
  ))

  expect_s4_class(rnaseq_res, "PharmacoSig")
  expect_s4_class(cnv_res, "PharmacoSig")
  expect_s4_class(mutation_res, "PharmacoSig")
  expect_equal(dim(rnaseq_res), c(2, length(drugs), 8))
  expect_equal(dim(cnv_res), c(2, length(drugs), 8))
  expect_equal(dim(mutation_res), c(2, length(drugs), 8))
})

test_that("drugSensitivitySig accepts continuous custom annotations and mirna", {
  data(CCLEsmall)
  drugs <- treatmentNames(CCLEsmall)[1:2]
  cells <- sampleNames(CCLEsmall)[1:10]
  sProfiles <- summarizeSensitivityProfiles(
    CCLEsmall,
    sensitivity.measure = "aac_recomputed",
    drugs = drugs,
    cell.lines = cells,
    verbose = FALSE
  )

  mirna_pset <- CCLEsmall
  S4Vectors::metadata(molecularProfilesSlot(mirna_pset)[[
    "rna"
  ]])$annotation <- "mirna"
  mirna_res <- suppressWarnings(drugSensitivitySig(
    mirna_pset,
    mDataType = "rna",
    drugs = drugs,
    cells = cells,
    features = rownames(featureInfo(mirna_pset, "rna"))[1:2],
    sProfiles = sProfiles,
    modeling.method = "pearson",
    parallel.on = "drug",
    nthread = 1,
    verbose = FALSE
  ))

  custom_pset <- CCLEsmall
  S4Vectors::metadata(molecularProfilesSlot(custom_pset)[[
    "rna"
  ]])$annotation <- "rnaseq.comp"
  custom_res <- suppressWarnings(drugSensitivitySig(
    custom_pset,
    mDataType = "rna",
    drugs = drugs,
    cells = cells,
    features = rownames(featureInfo(custom_pset, "rna"))[1:2],
    sProfiles = sProfiles,
    modeling.method = "pearson",
    parallel.on = "drug",
    nthread = 1,
    verbose = FALSE
  ))

  expect_s4_class(mirna_res, "PharmacoSig")
  expect_s4_class(custom_res, "PharmacoSig")
  expect_equal(dim(mirna_res), c(2, length(drugs), 8))
  expect_equal(dim(custom_res), c(2, length(drugs), 8))
})

test_that("drugSensitivitySig supports lm alias and spearman for continuous inputs", {
  data(CCLEsmall)
  drugs <- treatmentNames(CCLEsmall)[1:2]
  cells <- sampleNames(CCLEsmall)[1:10]
  features <- rownames(featureInfo(CCLEsmall, "rna"))[1:2]

  anova_res <- drugSensitivitySig(
    CCLEsmall,
    mDataType = "rna",
    drugs = drugs,
    cells = cells,
    features = features,
    modeling.method = "anova",
    parallel.on = "drug",
    nthread = 1,
    verbose = FALSE
  )
  lm_res <- drugSensitivitySig(
    CCLEsmall,
    mDataType = "rna",
    drugs = drugs,
    cells = cells,
    features = features,
    modeling.method = "lm",
    parallel.on = "drug",
    nthread = 1,
    verbose = FALSE
  )
  spearman_res <- suppressWarnings(drugSensitivitySig(
    CCLEsmall,
    mDataType = "rna",
    drugs = drugs,
    cells = cells,
    features = features,
    modeling.method = "spearman",
    parallel.on = "drug",
    nthread = 1,
    verbose = FALSE
  ))

  expect_equal(lm_res@.Data, anova_res@.Data)
  expect_s4_class(spearman_res, "PharmacoSig")
  expect_equal(dim(spearman_res), c(length(features), length(drugs), 8))
})
