#' Creates a signature representing the association between gene expression (or
#' other molecular profile) and drug dose response, for use in drug sensitivity
#' analysis.
#'
#' Given a Pharmacoset of the sensitivity experiment type, and a list of drugs,
#' the function will compute a signature for the effect gene expression on the
#' molecular profile of a cell. The function returns the estimated coefficient,
#' the t-stat, the p-value and the false discovery rate associated with that
#' coefficient, in a 3 dimensional array, with genes in the first direction,
#' drugs in the second, and the selected return values in the third.
#'
#' @examples
#' data(GDSCsmall)
#' drug.sensitivity <- drugSensitivitySig(GDSCsmall,
#'   mDataType = "rna",
#'   nthread = 1, features = fNames(GDSCsmall, "rna")[1]
#' )
#' print(drug.sensitivity)
#'
#' @param object `PharmacoSet` a PharmacoSet of the perturbation experiment type
#' @param mDataType `character` which one of the molecular data types to use
#'   in the analysis, out of dna, rna, rnaseq, snp, cnv
#' @param drugs `character` a vector of drug names for which to compute the
#'   signatures. Should match the names used in the PharmacoSet.
#' @param features `character` a vector of features for which to compute the
#'   signatures. Should match the names used in correspondant molecular data in PharmacoSet.
#' @param cells `character` allows choosing exactly which cell lines to include for the signature fitting.
#'   Should be a subset of sampleNames(pSet)
#' @param tissues `character` a vector of which tissue types to include in the signature fitting.
#'   Should be a subset of sampleInfo(pSet)$tissueid
#' @param nthread `numeric` if multiple cores are available, how many cores
#'   should the computation be parallelized over?
#' @param returnValues `character` Which of estimate, t-stat, p-value and fdr
#'   should the function return for each gene drug pair?
#' @param sensitivity.measure `character` which measure of the drug dose
#'   sensitivity should the function use for its computations? Use the
#'   sensitivityMeasures function to find out what measures are available for each PSet.
#' @param molecular.summary.stat `character` What summary statistic should be used to
#'   summarize duplicates for cell line molecular profile measurements?
#' @param sensitivity.summary.stat `character` What summary statistic should be used to
#'   summarize duplicates for cell line sensitivity measurements?
#' @param sensitivity.cutoff `numeric` Allows the user to binarize the sensitivity data using this threshold.
#' @param standardize `character` One of "SD", "rescale", or "none", for the form of standardization of
#'   the data to use. If "SD", the the data is scaled so that SD = 1. If rescale, then the data is scaled so that the 95%
#'   interquantile range lies in \[0,1\]. If none no rescaling is done.
#' @param molecular.cutoff Allows the user to binarize the sensitivity data using this threshold.
#' @param molecular.cutoff.direction `character` One of "less" or "greater", allows to set direction of binarization.
#' @param verbose `logical` 'TRUE' if the warnings and other informative message shoud be displayed
#' @param parallel.on One of "gene" or "drug", chooses which level to parallelize computation (by gene, or by drug).
#' @param modeling.method One of "anova" or "pearson". If "anova", nested linear models (including and excluding the molecular feature) adjusted for
#'   are fit after the data is standardized, and ANOVA is used to estimate significance. If "pearson", partial correlation adjusted for tissue of origin are
#'   fit to the data, and a Pearson t-test (or permutation) test are used. Note that the difference is in whether standardization is done across the whole
#'   dataset (anova) or within each tissue (pearson), as well as the test applied.
#' @param inference.method Should "analytic" or "resampling" (permutation testing + bootstrap) inference be used to estimate significance.
#'   For permutation testing, QUICK-STOP is used to adaptively stop permutations. Resampling is currently only implemented for "pearson" modelling method.
#' @param ... additional arguments not currently fully supported by the function
#'
#' @return `array` a 3D array with genes in the first dimension, drugs in the
#'   second, and return values in the third.
#'
#' @importMethodsFrom CoreGx drugSensitivitySig
#' @export
setMethod(
  "drugSensitivitySig",
  signature(object = "PharmacoSet"),
  function(object, mDataType, drugs, features, cells, tissues, sensitivity.measure = "auc_recomputed",
           molecular.summary.stat = c("mean", "median", "first", "last", "or", "and"),
           sensitivity.summary.stat = c("mean", "median", "first", "last"),
           returnValues = c("estimate", "pvalue", "fdr"),
           sensitivity.cutoff, standardize = c("SD", "rescale", "none"), molecular.cutoff = NA,
           molecular.cutoff.direction = c("less", "greater"),
           nthread = 1, parallel.on = c("drug", "gene"), modeling.method = c("anova", "pearson"),
           inference.method = c("analytic", "resampling"), verbose = TRUE, ...) {
    .drugSensitivitySigPharmacoSet(
      object, mDataType, drugs, features, cells, tissues, sensitivity.measure,
      molecular.summary.stat, sensitivity.summary.stat, returnValues,
      sensitivity.cutoff, standardize, molecular.cutoff, molecular.cutoff.direction,
      nthread, parallel.on, modeling.method, inference.method, verbose, ...
    )
  }
)

#' @import parallel
#' @importFrom SummarizedExperiment assayNames assay
#' @keywords internal
.drugSensitivitySigPharmacoSet <- function(object,
                                           mDataType,
                                           drugs,
                                           features,
                                           cells,
                                           tissues,
                                           sensitivity.measure = "auc_recomputed",
                                           molecular.summary.stat = c("mean", "median", "first", "last", "or", "and"),
                                           sensitivity.summary.stat = c("mean", "median", "first", "last"),
                                           returnValues = c("estimate", "pvalue", "fdr"),
                                           sensitivity.cutoff, standardize = c("SD", "rescale", "none"),
                                           molecular.cutoff = NA,
                                           molecular.cutoff.direction = c("less", "greater"),
                                           nthread = 1,
                                           parallel.on = c("drug", "gene"),
                                           modeling.method = c("anova", "pearson"),
                                           inference.method = c("analytic", "resampling"),
                                           verbose = TRUE,
                                           ...) {

  ### This function needs to: Get a table of AUC values per cell line / drug
  ### Be able to recompute those values on the fly from raw data if needed to change concentration
  ### Be able to choose different summary methods on fly if needed (need to add annotation to table to tell what summary method previously used)
  ### Be able to extract genomic data
  ### Run rankGeneDrugSens in parallel at the drug level
  ### Return matrix as we had before

  # sensitivity.measure <- match.arg(sensitivity.measure)
  molecular.summary.stat <- match.arg(molecular.summary.stat)
  sensitivity.summary.stat <- match.arg(sensitivity.summary.stat)
  standardize <- match.arg(standardize)
  molecular.cutoff.direction <- match.arg(molecular.cutoff.direction)
  parallel.on <- match.arg(parallel.on)
  dots <- list(...)
  ndots <- length(dots)
  modeling.method <- match.arg(modeling.method)
  inference.method <- match.arg(inference.method)



  if (is.null(dots[["sProfiles"]]) & !all(sensitivity.measure %in% colnames(sensitivityProfiles(object)))) {
    stop(sprintf("Invalid sensitivity measure for %s, choose among: %s", annotation(object)$name, paste(colnames(sensitivityProfiles(object)), collapse = ", ")))
  }

  if (!(mDataType %in% names(molecularProfilesSlot(object)))) {
    stop(sprintf("Invalid mDataType for %s, choose among: %s", annotation(object)$name, paste(names(molecularProfilesSlot(object)), collapse = ", ")))
  }
  switch(S4Vectors::metadata(molecularProfilesSlot(object)[[mDataType]])$annotation,
    "mutation" = {
      if (!is.element(molecular.summary.stat, c("or", "and"))) {
        stop("Molecular summary statistic for mutation must be either 'or' or 'and'")
      }
    },
    "fusion" = {
      if (!is.element(molecular.summary.stat, c("or", "and"))) {
        stop("Molecular summary statistic for fusion must be either 'or' or 'and'")
      }
    },
    "rna" = {
      if (!is.element(molecular.summary.stat, c("mean", "median", "first", "last"))) {
        stop("Molecular summary statistic for rna must be either 'mean', 'median', 'first' or 'last'")
      }
    },
    "cnv" = {
      if (!is.element(molecular.summary.stat, c("mean", "median", "first", "last"))) {
        stop("Molecular summary statistic for cnv must be either 'mean', 'median', 'first' or 'last'")
      }
    },
    "rnaseq" = {
      if (!is.element(molecular.summary.stat, c("mean", "median", "first", "last"))) {
        stop("Molecular summary statistic for rna must be either 'mean', 'median', 'first' or 'last'")
      }
    },
    "isoform" = {
      if (!is.element(molecular.summary.stat, c("mean", "median", "first", "last"))) {
        stop("Molecular summary statistic for rna must be either 'mean', 'median', 'first' or 'last'")
      }
    },
    stop(sprintf("No summary statistic for %s has been implemented yet", S4Vectors::metadata(molecularProfilesSlot(object)[[mDataType]])$annotation))
  )

  if (!is.element(sensitivity.summary.stat, c("mean", "median", "first", "last"))) {
    stop("Sensitivity summary statistic for sensitivity must be either 'mean', 'median', 'first' or 'last'")
  }

  if (missing(sensitivity.cutoff)) {
    sensitivity.cutoff <- NA
  }
  if (missing(drugs)) {
    if(is.null(dots[["sProfiles"]])){
      drugn <- drugs <- treatmentNames(object)
    } else {
      drugn <- drugs <- rownames(dots[["sProfiles"]])
    }
  } else {
    drugn <- drugs
  }

  if (missing(cells)) {
    celln <- cells <- sampleNames(object)
  } else {
    celln <- cells
  }

  availcore <- parallel::detectCores()

  if (nthread > availcore) {
    nthread <- availcore
  }

  if (parallel.on == "drug") {
    nthread_drug <- nthread
    nthread_gene <- 1
  } else {
    nthread_gene <- nthread
    nthread_drug <- 1
  }

  if (missing(features)) {
    features <- rownames(featureInfo(object, mDataType))
  } else {
    fix <- is.element(features, rownames(featureInfo(object, mDataType)))
    if (verbose && !all(fix)) {
      warning(sprintf("%i/%i features can be found", sum(fix), length(features)))
    }
    features <- features[fix]
  }

  # if(missing(modeling.method)){
  #   modeling.method <- "anova"
  # }
  #
  # if(missing(inference.method)){
  #   inference.method <- "analytic"
  # }

  if (is.null(dots[["sProfiles"]])) {
    drugpheno.all <- lapply(sensitivity.measure, function(sensitivity.measure) {
      return(t(summarizeSensitivityProfiles(object,
        sensitivity.measure = sensitivity.measure,
        summary.stat = sensitivity.summary.stat,
        verbose = verbose
      )))
    })
  } else {
    sProfiles <- dots[["sProfiles"]]
    drugpheno.all <- list(t(sProfiles))
  }

  dix <- is.element(drugn, do.call(colnames, drugpheno.all))
  if (verbose && !all(dix)) {
    warning(sprintf("%i/%i drugs can be found", sum(dix), length(drugn)))
  }
  if (!any(dix)) {
    stop("None of the drugs were found in the dataset")
  }
  drugn <- drugn[dix]

  cix <- is.element(celln, do.call(rownames, drugpheno.all))
  if (verbose && !all(cix)) {
    warning(sprintf("%i/%i cells can be found", sum(cix), length(celln)))
  }
  if (!any(cix)) {
    stop("None of the cells were found in the dataset")
  }
  celln <- celln[cix]

  if (!missing(tissues)) {
    celln <- celln[sampleInfo(object)[celln, "tissueid"] %in% tissues]
  } else {
    tissues <- unique(sampleInfo(object)[celln, "tissueid"])
  }

  molecularProfilesSlot(object)[[mDataType]] <- summarizeMolecularProfiles(
    object = object,
    mDataType = mDataType,
    summary.stat = molecular.summary.stat,
    binarize.threshold = molecular.cutoff,
    binarize.direction = molecular.cutoff.direction,
    verbose = verbose
  )[features, ]

  if (!is.null(dots[["mProfiles"]])) {
    mProfiles <- dots[["mProfiles"]]
    SummarizedExperiment::assay(molecularProfilesSlot(object)[[mDataType]]) <- mProfiles[features, colnames(molecularProfilesSlot(object)[[mDataType]]), drop = FALSE]
  }

  drugpheno.all <- lapply(drugpheno.all, function(x) {
    x[intersect(phenoInfo(object, mDataType)[, "sampleid"], celln), , drop = FALSE]
  })

  molcellx <- phenoInfo(object, mDataType)[, "sampleid"] %in% celln

  type <- as.factor(sampleInfo(object)[phenoInfo(object, mDataType)[molcellx, "sampleid"], "tissueid"])

  if ("batchid" %in% colnames(phenoInfo(object, mDataType))) {
    batch <- phenoInfo(object, mDataType)[molcellx, "batchid"]
  } else {
    batch <- rep(NA, times = nrow(phenoInfo(object, mDataType)))
  }
  batch[!is.na(batch) & batch == "NA"] <- NA
  batch <- as.factor(batch)
  names(batch) <- phenoInfo(object, mDataType)[molcellx, "sampleid"]
  batch <- batch[rownames(drugpheno.all[[1]])]
  if (verbose) {
    message("Computing drug sensitivity signatures...")
  }

  ### Calculate approximate number of perms needed



  if (is.null(dots[["req_alpha"]])) {
    req_alpha <- 0.05 / (nrow(molecularProfilesSlot(object)[[mDataType]])) ## bonferonni correction
  } else {
    req_alpha <- dots[["req_alpha"]]
  }



  # splitix <- parallel::splitIndices(nx = length(drugn), ncl = nthread_drug)
  # splitix <- splitix[vapply(splitix, length, FUN.VALUE=numeric(1)) > 0]
  mcres <- parallel::mclapply(seq_along(drugn), function(x, drugn, expr, drugpheno, type, batch, standardize, nthread, modeling.method, inference.method, req_alpha) {
    res <- NULL
    for (i in drugn[x]) {
      ## using a linear model (x ~ concentration + cell + batch)
      dd <- lapply(drugpheno, function(x) x[, i])
      dd <- do.call(cbind, dd)
      colnames(dd) <- seq_len(ncol(dd))
      if (!is.na(sensitivity.cutoff)) {
        dd <- factor(ifelse(dd > sensitivity.cutoff, 1, 0), levels = c(0, 1))
      }
      rr <- rankGeneDrugSensitivity(data = expr, drugpheno = dd, type = type, batch = batch, single.type = FALSE, standardize = standardize, nthread = nthread, verbose = verbose, modeling.method = modeling.method, inference.method = inference.method, req_alpha)
      res <- c(res, list(rr$all))
    }
    names(res) <- drugn[x]
    return(res)
  },
  drugn = drugn, expr = t(molecularProfiles(object, mDataType)[features, molcellx, drop = FALSE]),
  drugpheno = drugpheno.all, type = type, batch = batch, nthread = nthread_gene, standardize = standardize,
  modeling.method = modeling.method, inference.method = inference.method,
  req_alpha = req_alpha, mc.cores = nthread_drug, mc.preschedule = FALSE
  )

  res <- do.call(c, mcres)
  res <- res[!vapply(res, is.null, FUN.VALUE = logical(1))]
  drug.sensitivity <- array(NA,
    dim = c(
      nrow(featureInfo(object, mDataType)[features, , drop = FALSE]),
      length(res), ncol(res[[1]])
    ),
    dimnames = list(rownames(featureInfo(object, mDataType)[features, , drop = FALSE]), names(res), colnames(res[[1]]))
  )
  for (j in seq_len(ncol(res[[1]]))) {
    ttt <- vapply(res, function(x, j, k) {
      xx <- array(NA, dim = length(k), dimnames = list(k))
      xx[rownames(x)] <- x[, j, drop = FALSE]
      return(xx)
    },
    j = j,
    k = rownames(featureInfo(object, mDataType)[features, , drop = FALSE]),
    FUN.VALUE = numeric(dim(drug.sensitivity)[1])
    )
    drug.sensitivity[rownames(featureInfo(object, mDataType)[features, , drop = FALSE]), names(res), j] <- ttt
  }

  drug.sensitivity <- PharmacoSig(drug.sensitivity,
    PSetName = name(object),
    Call = as.character(match.call()),
    SigType = "Sensitivity",
    Arguments = list(
      "mDataType" = mDataType,
      "drugs" = drugs,
      "features" = features,
      "cells" = cells,
      "tissues" = tissues,
      "sensitivity.measure" = sensitivity.measure,
      "molecular.summary.stat" = molecular.summary.stat,
      "sensitivity.summary.stat" = sensitivity.summary.stat,
      "returnValues" = returnValues,
      "sensitivity.cutoff" = sensitivity.cutoff,
      "standardize" = standardize,
      "molecular.cutoff" = molecular.cutoff,
      "molecular.cutoff.direction" = molecular.cutoff.direction,
      "nthread" = nthread,
      "verbose" = verbose
    )
  )

  return(drug.sensitivity)
}
