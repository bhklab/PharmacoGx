#' @importFrom stats complete.cases
#' @importFrom stats p.adjust

#################################################
## Rank genes based on drug effect in the Connectivity Map
##
## inputs:
##      - data: gene expression data matrix
##            - drugpheno: sensititivity values fo thr drug of interest
##            - type: cell or tissue type for each experiment
##            - duration: experiment duration in hours
##      - batch: experiment batches
##            - single.type: Should the statitsics be computed for each cell/tissue type separately?
##      - nthread: number of parallel threads (bound to the maximum number of cores available)
##
## outputs:
## list of datafraes with the statistics for each gene, for each type
##
## Notes:    duration is not taken into account as only 4 perturbations lasted 12h, the other 6096 lasted 6h
#################################################

rankGeneDrugSensitivity <- function(
  data,
  drugpheno,
  type,
  batch,
  single.type = FALSE,
  standardize = "SD",
  nthread = 1,
  verbose = FALSE,
  modeling.method = c("anova", "pearson", "lm", "spearman"),
  inference.method = c("analytic", "resampling"),
  req_alpha = 0.05
) {
  if (nthread != 1) {
    availcore <- parallel::detectCores()
    if ((missing(nthread) || nthread < 1 || nthread > availcore) && verbose) {
      warning(
        "nthread undefined, negative or larger than available cores. Resetting to maximum number of cores."
      )
      nthread <- availcore
    }
  }

  modeling.method <- match.arg(modeling.method)
  if (identical(modeling.method, "lm")) {
    modeling.method <- "anova"
  }
  inference.method <- match.arg(inference.method)

  if (modeling.method == "anova" && inference.method == "resampling") {
    stop("Resampling based inference for anova model is not yet implemented.")
  }

  if (is.null(dim(drugpheno))) {
    drugpheno <- data.frame(drugpheno)
  } else if (!is(drugpheno, "data.frame")) {
    drugpheno <- as.data.frame(drugpheno)
  }

  if (missing(type) || all(is.na(type))) {
    type <- array("other", dim = nrow(data), dimnames = list(rownames(data)))
  }
  if (missing(batch) || all(is.na(batch))) {
    batch <- array(1, dim = nrow(data), dimnames = list(rownames(data)))
  }
  if (any(c(nrow(drugpheno), length(type), length(batch)) != nrow(data))) {
    stop(
      "length of drugpheno, type, duration, and batch should be equal to the number of rows of data!"
    )
  }
  rownames(drugpheno) <- names(type) <- names(batch) <- rownames(data)
  if (modeling.method == "spearman") {
    if (any(unlist(lapply(drugpheno, is.factor)))) {
      stop(
        "Spearman modeling requires continuous sensitivity inputs. Use 'pearson', 'anova', or 'lm' for discrete sensitivity data."
      )
    }
    drugpheno <- data.frame(
      lapply(drugpheno, function(x) {
        rank(as.numeric(x), na.last = "keep", ties.method = "average")
      }),
      check.names = FALSE
    )
    rownames(drugpheno) <- rownames(data)
  }

  fit_feature <- function(feature_idx, data, type, batch, drugpheno) {
    feature_data <- data[, feature_idx]

    if (modeling.method == "anova") {
      return(geneDrugSensitivity(
        feature_data,
        type = type,
        batch = batch,
        drugpheno = drugpheno,
        verbose = verbose,
        standardize = standardize
      ))
    }

    if (modeling.method == "spearman") {
      if (!is.numeric(feature_data)) {
        stop(
          "Spearman modeling is only implemented for continuous molecular features. Use 'pearson', 'anova', or 'lm' for discrete features."
        )
      }

      return(geneDrugSensitivityPCorr(
        rank(
          as.numeric(feature_data),
          na.last = "keep",
          ties.method = "average"
        ),
        type = type,
        batch = batch,
        drugpheno = drugpheno,
        verbose = verbose,
        test = inference.method,
        req_alpha = req_alpha
      ))
    }

    if (!is.character(feature_data)) {
      return(geneDrugSensitivityPCorr(
        feature_data,
        type = type,
        batch = batch,
        drugpheno = drugpheno,
        verbose = verbose,
        test = inference.method,
        req_alpha = req_alpha
      ))
    }

    geneDrugSensitivityPBCorr(
      as.factor(feature_data),
      type = type,
      batch = batch,
      drugpheno = drugpheno,
      verbose = verbose,
      test = inference.method,
      req_alpha = req_alpha
    )
  }

  res <- NULL
  utype <- sort(unique(as.character(type)))
  ltype <- list("all" = utype)
  if (single.type) {
    ltype <- c(ltype, as.list(utype))
    names(ltype)[-1] <- utype
  }
  res <- NULL
  ccix <- complete.cases(data, type, batch, drugpheno)
  nn <- sum(ccix)

  if (modeling.method == "anova") {
    if (!any(unlist(lapply(drugpheno, is.factor)))) {
      if (ncol(drugpheno) > 1) {
        ##### FIX NAMES!!! This is important
        nc <- lapply(seq_len(ncol(drugpheno)), function(i) {
          est <- paste("estimate", i, sep = ".")
          se <- paste("se", i, sep = ".")
          tstat <- paste("tstat", i, sep = ".")

          nc <- c(est, se, tstat)
          return(nc)
        })
        nc <- c(nc, n = nn, "fstat" = NA, "pvalue" = NA, "fdr")
      } else {
        nc <- c("estimate", "se", "n", "tstat", "fstat", "pvalue", "df", "fdr")
      }
    } else {
      nc <- c("estimate", "se", "n", "tstat", "fstat", "pvalue", "df", "fdr")
    }
  } else if (modeling.method %in% c("pearson", "spearman")) {
    nc <- c("estimate", "n", "df", "significant", "pvalue", "lower", "upper")
  }

  for (ll in seq_len(length(ltype))) {
    iix <- !is.na(type) & is.element(type, ltype[[ll]])
    # ccix <- complete.cases(data[iix, , drop=FALSE], drugpheno[iix,,drop=FALSE], type[iix], batch[iix]) ### HACK???

    ccix <- rowSums(!is.na(data)) > 0 |
      rowSums(!is.na(drugpheno)) > 0 |
      is.na(type) |
      is.na(batch)
    ccix <- ccix[iix]
    # ccix <- !vapply(seq_len(NROW(data[iix,,drop=FALSE])), function(x) {
    #   return(all(is.na(data[iix,,drop=FALSE][x,])) || all(is.na(drugpheno[iix,,drop=FALSE][x,])) || all(is.na(type[iix][x])) || all(is.na(batch[iix][x])))
    # }, FUN.VALUE=logical(1))

    if (sum(ccix) < 3) {
      ## not enough experiments
      rest <- list(matrix(
        NA,
        nrow = ncol(data),
        ncol = length(nc),
        dimnames = list(colnames(data), nc)
      ))
      res <- c(res, rest)
    } else {
      feature_names <- colnames(data[iix, , drop = FALSE])
      apply_over_features <- if (nthread > 1) {
        function(X, FUN, ...) {
          parallel::mclapply(
            X,
            FUN,
            ...,
            mc.cores = nthread,
            mc.preschedule = TRUE
          )
        }
      } else {
        lapply
      }
      mcres <- apply_over_features(
        seq_len(ncol(data)),
        fit_feature,
        data = data[iix, , drop = FALSE],
        type = type[iix],
        batch = batch[iix],
        drugpheno = drugpheno[iix, , drop = FALSE]
      )
      rest <- do.call(
        rbind,
        lapply(seq_along(mcres), function(i) {
          matrix(
            mcres[[i]],
            nrow = 1,
            dimnames = list(feature_names[[i]], names(mcres[[i]]))
          )
        })
      )
      rest <- cbind(rest, "fdr" = p.adjust(rest[, "pvalue"], method = "fdr"))
      # rest <- rest[ , nc, drop=FALSE]
      res <- c(res, list(rest))
    }
  }
  names(res) <- names(ltype)
  return(res)
}

## End
