make_custom_profile <- function(
  pset,
  source = "rna",
  annotation = "rna.custom"
) {
  profile <- molecularProfilesSlot(pset)[[source]]
  values <- SummarizedExperiment::assay(profile)
  SummarizedExperiment::assays(profile) <- list(
    decoy = if (is.numeric(values)) -values else values,
    chosen = values
  )
  S4Vectors::metadata(profile)$annotation <- annotation
  molecularProfilesSlot(pset) <- list("lab.profile+[v1]" = profile)
  pset
}

test_that("custom keys and inner assays give the original molecular summary", {
  data(GDSCsmall)
  custom <- make_custom_profile(GDSCsmall)
  baseline <- summarizeMolecularProfiles(GDSCsmall, "rna", verbose = FALSE)
  selected <- summarizeMolecularProfiles(
    custom,
    "lab.profile+[v1]",
    assay = "chosen",
    verbose = FALSE
  )
  expect_equal(
    SummarizedExperiment::assay(selected),
    SummarizedExperiment::assay(baseline)
  )
  expect_identical(SummarizedExperiment::assayNames(selected), "chosen")
  expect_equal(
    S4Vectors::metadata(selected),
    S4Vectors::metadata(molecularProfilesSlot(custom)[[1]])
  )
  expect_equal(
    SummarizedExperiment::assay(summarizeMolecularProfiles(
      custom,
      "lab.profile+[v1]",
      assay = 2L,
      verbose = FALSE
    )),
    SummarizedExperiment::assay(selected)
  )
  expect_equal(
    SummarizedExperiment::assay(summarizeMolecularProfiles(
      custom,
      "lab.profile+[v1]",
      verbose = FALSE
    )),
    -SummarizedExperiment::assay(selected)
  )
  expect_equal(
    summarizeMolecularProfiles(custom, "lab.profile+[v1]", summarize = FALSE),
    molecularProfilesSlot(custom)[[1]]
  )
})

test_that("assay selectors reject missing, ambiguous and invalid selections", {
  data(GDSCsmall)
  custom <- make_custom_profile(GDSCsmall)
  for (selector in list(
    "absent",
    NA_character_,
    character(),
    c("chosen", "decoy"),
    0,
    3,
    1.5,
    Inf,
    TRUE
  )) {
    expect_error(
      summarizeMolecularProfiles(custom, "lab.profile+[v1]", assay = selector),
      "Invalid assay"
    )
  }
  expect_error(summarizeMolecularProfiles(custom, "rna"), "Invalid mDataType")
  expect_error(drugPerturbationSig(custom, "absent"), "Invalid mDataType")
  profile <- molecularProfilesSlot(custom)[[1]]
  SummarizedExperiment::assayNames(profile) <- c("same", "same")
  molecularProfilesSlot(custom)[[1]] <- profile
  expect_error(
    summarizeMolecularProfiles(custom, "lab.profile+[v1]", assay = "same"),
    "Invalid assay"
  )
  names(SummarizedExperiment::assays(profile)) <- NULL
  molecularProfilesSlot(custom)[[1]] <- profile
  expect_null(SummarizedExperiment::assayNames(summarizeMolecularProfiles(
    custom,
    "lab.profile+[v1]",
    assay = 2,
    verbose = FALSE
  )))
})

test_that("sensitivity modeling forwards inner assay selection and records it", {
  data(CCLEsmall)
  custom <- make_custom_profile(CCLEsmall)
  args <- list(
    object = CCLEsmall,
    mDataType = "rna",
    features = rownames(featureInfo(CCLEsmall, "rna"))[1:2],
    drugs = treatmentNames(CCLEsmall)[1],
    verbose = FALSE
  )
  baseline <- do.call(drugSensitivitySig, args)
  args$object <- custom
  args$mDataType <- "lab.profile+[v1]"
  args$assay <- "chosen"
  selected <- do.call(drugSensitivitySig, args)
  expect_equal(selected@.Data, baseline@.Data)
  expect_identical(selected@Arguments$assay, "chosen")
  args$assay <- "absent"
  expect_error(do.call(drugSensitivitySig, args), "Invalid assay")
  args$assay <- "chosen"
  args$mProfiles <- SummarizedExperiment::assay(summarizeMolecularProfiles(
    CCLEsmall,
    "rna",
    features = args$features,
    verbose = FALSE
  ))
  expect_equal(do.call(drugSensitivitySig, args)@.Data, baseline@.Data)
})

test_that("perturbation modeling accepts custom RNA keys and dotted annotations", {
  data(CMAPsmall)
  custom <- make_custom_profile(CMAPsmall)
  features <- rownames(featureInfo(CMAPsmall, "rna"))[1:3]
  drug <- treatmentNames(CMAPsmall)[1]
  baseline <- drugPerturbationSig(
    CMAPsmall,
    "rna",
    drugs = drug,
    features = features
  )
  selected <- drugPerturbationSig(
    custom,
    "lab.profile+[v1]",
    drugs = drug,
    features = features,
    assay = "chosen"
  )
  expect_equal(selected@.Data, baseline@.Data)
  expect_identical(selected@Arguments$assay, "chosen")
  expect_error(
    drugPerturbationSig(custom, "lab.profile+[v1]", assay = "absent"),
    "Invalid assay"
  )
  nonstandard <- make_custom_profile(CMAPsmall, annotation = "custom abundance")
  expect_error(drugPerturbationSig(nonstandard, "lab.profile+[v1]"), "Only rna")
  expect_equal(
    drugPerturbationSig(
      nonstandard,
      "lab.profile+[v1]",
      drugs = drug,
      features = features,
      assay = "chosen",
      molecular.type = "rna"
    )@.Data,
    baseline@.Data
  )
  for (annotation in c(
    "rnaseq.custom",
    "mutation.custom",
    "fusion.custom",
    "cnv.custom"
  )) {
    altered <- make_custom_profile(CMAPsmall, annotation = annotation)
    expect_error(drugPerturbationSig(altered, "lab.profile+[v1]"), "Only rna")
    expect_error(
      drugPerturbationSig(altered, "lab.profile+[v1]", molecular.type = "rna"),
      "conflicts"
    )
  }
})

test_that("nonstandard metadata needs an explicit biological type for fitting", {
  data(CCLEsmall)
  for (annotation in list(
    NULL,
    NA_character_,
    "custom abundance",
    c("rna", "mutation")
  )) {
    custom <- make_custom_profile(CCLEsmall, annotation = annotation)
    args <- list(
      object = custom,
      mDataType = "lab.profile+[v1]",
      features = rownames(featureInfo(CCLEsmall, "rna"))[1],
      drugs = treatmentNames(CCLEsmall)[1],
      assay = "chosen",
      verbose = FALSE
    )
    expect_error(do.call(drugSensitivitySig, args), "supply molecular.type")
    args$molecular.type <- "rna"
    expect_s4_class(do.call(drugSensitivitySig, args), "PharmacoSig")
  }
  custom <- make_custom_profile(CCLEsmall)
  expect_error(
    summarizeMolecularProfiles(
      custom,
      "lab.profile+[v1]",
      molecular.type = "mutation"
    ),
    "conflicts"
  )
  expect_error(
    summarizeMolecularProfiles(
      custom,
      "lab.profile+[v1]",
      molecular.type = "anything"
    ),
    "Invalid molecular.type"
  )
})

test_that("dotted binary annotations preserve logical aggregation and encoding", {
  data(GDSCsmall)
  custom <- make_custom_profile(
    GDSCsmall,
    source = "mutation",
    annotation = "mutation.custom"
  )
  baseline <- summarizeMolecularProfiles(
    GDSCsmall,
    "mutation",
    summary.stat = "or",
    verbose = FALSE
  )
  selected <- summarizeMolecularProfiles(
    custom,
    "lab.profile+[v1]",
    summary.stat = "or",
    assay = "chosen",
    verbose = FALSE
  )
  expect_equal(
    SummarizedExperiment::assay(selected),
    SummarizedExperiment::assay(baseline)
  )
  expect_error(
    summarizeMolecularProfiles(
      custom,
      "lab.profile+[v1]",
      summary.stat = "mean"
    ),
    "choose among: and, or"
  )
})

test_that("CCLE annotation aliases retain their scientific data types", {
  expect_identical(
    PharmacoGx:::.normalizeMolecularAnnotation("isoforms.vendor"),
    "isoform"
  )
  expect_identical(
    PharmacoGx:::.normalizeMolecularAnnotation("mut"),
    "mutation"
  )
  data(GDSCsmall)
  custom <- make_custom_profile(GDSCsmall, annotation = "mut")
  profile <- molecularProfilesSlot(custom)[[1]][1, 1:3, drop = FALSE]
  values <- matrix(c(0, 1, NA_real_), nrow = 1, dimnames = dimnames(profile))
  SummarizedExperiment::assays(profile) <- list(binary = values)
  SummarizedExperiment::colData(profile)$sampleid <- rep(
    sampleNames(custom)[1],
    3
  )
  molecularProfilesSlot(custom)[[1]] <- profile
  for (stat in c("or", "and")) {
    summarized <- summarizeMolecularProfiles(
      custom,
      "lab.profile+[v1]",
      assay = "binary",
      summary.stat = stat,
      fill.missing = FALSE,
      verbose = FALSE
    )
    expect_identical(
      unname(SummarizedExperiment::assay(summarized)[1, 1]),
      if (stat == "or") "1" else "0"
    )
  }
  expect_error(
    summarizeMolecularProfiles(
      custom,
      "lab.profile+[v1]",
      summary.stat = "mean"
    ),
    "choose among: and, or"
  )
  expect_error(
    drugSensitivitySig(
      custom,
      "lab.profile+[v1]",
      molecular.summary.stat = "or",
      modeling.method = "spearman",
      verbose = FALSE
    ),
    "continuous molecular"
  )
  SummarizedExperiment::assay(profile)[1, 1] <- 2
  molecularProfilesSlot(custom)[[1]] <- profile
  expect_error(
    summarizeMolecularProfiles(custom, "lab.profile+[v1]", summary.stat = "or"),
    "only binary 0/1"
  )
})

test_that("inner assay selection is respected during single-feature binarization", {
  data(GDSCsmall)
  custom <- make_custom_profile(GDSCsmall)
  feature <- rownames(featureInfo(custom, "lab.profile+[v1]"))[1]
  for (direction in c("less", "greater")) {
    baseline <- summarizeMolecularProfiles(
      GDSCsmall,
      "rna",
      features = feature,
      binarize.threshold = 0,
      binarize.direction = direction,
      verbose = FALSE
    )
    selected <- summarizeMolecularProfiles(
      custom,
      "lab.profile+[v1]",
      assay = "chosen",
      features = feature,
      binarize.threshold = 0,
      binarize.direction = direction,
      verbose = FALSE
    )
    expect_equal(
      SummarizedExperiment::assay(selected),
      SummarizedExperiment::assay(baseline)
    )
  }
})

test_that("intersection uses biological groups with independent custom keys", {
  data(GDSCsmall)
  for (annotations in list(
    c("rna.custom", "rnaseq.custom"),
    c("mutation.custom", "mut.custom")
  )) {
    psets <- lapply(annotations, function(annot) {
      make_custom_profile(GDSCsmall, annotation = annot)
    })
    cells <- unique(phenoInfo(psets[[1]], "lab.profile+[v1]")$sampleid)[1:3]
    for (i in seq_along(psets)) {
      profile <- molecularProfilesSlot(psets[[i]])[[1]]
      keep <- SummarizedExperiment::colData(profile)$sampleid %in%
        cells[i:(i + 1)]
      molecularProfilesSlot(psets[[i]]) <- stats::setNames(
        list(profile[, keep]),
        paste0("profile", i, "+[v1]")
      )
    }
    for (strict in c(FALSE, TRUE)) {
      selected <- intersectPSet(
        psets,
        intersectOn = "cell.lines",
        strictIntersect = strict,
        verbose = FALSE
      )
      for (pset in selected) {
        expect_identical(
          unique(
            SummarizedExperiment::colData(molecularProfilesSlot(pset)[[
              1
            ]])$sampleid
          ),
          cells[2]
        )
        expect_identical(
          SummarizedExperiment::assayNames(molecularProfilesSlot(pset)[[1]]),
          c("decoy", "chosen")
        )
      }
    }
  }
})

test_that("molecular cell selectors prefer exact custom profile keys", {
  data(GDSCsmall)
  custom <- make_custom_profile(GDSCsmall)
  cells <- unique(phenoInfo(custom, "lab.profile+[v1]")$sampleid)[1:2]
  selected <- subsetTo(
    custom,
    molecular.data.cells = list("lab.profile+[v1]" = cells[1], rna = cells[2])
  )
  expect_identical(
    unique(phenoInfo(selected, "lab.profile+[v1]")$sampleid),
    cells[1]
  )
  legacy <- subsetTo(custom, molecular.data.cells = list(rna = cells[2]))
  expect_identical(
    unique(phenoInfo(legacy, "lab.profile+[v1]")$sampleid),
    cells[2]
  )
  expect_identical(
    SummarizedExperiment::assayNames(molecularProfilesSlot(selected)[[1]]),
    c("decoy", "chosen")
  )
})
