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

rankGeneDrugSensitivity <- function (data, drugpheno, type, batch,
                                     single.type=FALSE, standardize = "SD",
                                     nthread=1, verbose=FALSE,
                                     modeling.method = c("anova", "pearson"),
                                     inference.method = c("analytic", "resampling"), req_alpha = 0.05) {

  if (nthread != 1) {
    availcore <- parallel::detectCores()
    if ((missing(nthread) || nthread < 1 || nthread > availcore) && verbose) {
      warning("nthread undefined, negative or larger than available cores. Resetting to maximum number of cores.")
      nthread <- availcore
    }
  }

  modeling.method <- match.arg(modeling.method)
  inference.method <- match.arg(inference.method)

  if(modeling.method == "anova" && inference.method == "resampling") {
    stop("Resampling based inference for anova model is not yet implemented.")
  }

  if(is.null(dim(drugpheno))){

    drugpheno <- data.frame(drugpheno)

  } else if(!is(drugpheno, "data.frame")) {
    drugpheno <- as.data.frame(drugpheno)

  }

  if (missing(type) || all(is.na(type))) {
    type <- array("other", dim=nrow(data), dimnames=list(rownames(data)))
  }
  if (missing(batch) || all(is.na(batch))) {
    batch <- array(1, dim=nrow(data), dimnames=list(rownames(data)))
  }
  if (any(c(nrow(drugpheno), length(type), length(batch)) != nrow(data))) {
    stop("length of drugpheno, type, duration, and batch should be equal to the number of rows of data!")
  }
  rownames(drugpheno) <- names(type) <- names(batch) <- rownames(data)

  res <- NULL
  utype <- sort(unique(as.character(type)))
  ltype <- list("all"=utype)
  if (single.type) {
    ltype <- c(ltype, as.list(utype))
    names(ltype)[-1] <- utype
  }
  res <- NULL
  ccix <- complete.cases(data, type, batch, drugpheno)
  nn <- sum(ccix)


  if(modeling.method == "anova"){
    if(!any(unlist(lapply(drugpheno,is.factor)))){
      if(ncol(drugpheno)>1){
        ##### FIX NAMES!!! This is important
        nc <- lapply(seq_len(ncol(drugpheno)), function(i){

          est <- paste("estimate", i, sep=".")
          se <-  paste("se", i, sep=".")
          tstat <- paste("tstat", i, sep=".")

          nc <- c(est, se, tstat)
          return(nc)

        })
        nc  <- c(nc, n=nn, "fstat"=NA, "pvalue"=NA, "fdr")
      } else {
        nc  <- c("estimate", "se", "n", "tstat", "fstat", "pvalue", "df", "fdr")
      }
    } else {
      nc <- c("estimate", "se", "n", "tstat", "fstat", "pvalue", "df", "fdr")
    }
  } else if (modeling.method == "pearson") {
    nc <- c("estimate", "n", "df", "significant", "pvalue", "lower", "upper")
  }

  for (ll in seq_len(length(ltype))) {
    iix <- !is.na(type) & is.element(type, ltype[[ll]])
    # ccix <- complete.cases(data[iix, , drop=FALSE], drugpheno[iix,,drop=FALSE], type[iix], batch[iix]) ### HACK???

    ccix <- rowSums(!is.na(data)) > 0 | rowSums(!is.na(drugpheno)) > 0 | is.na(type) | is.na(batch)
    ccix <- ccix[iix]
    # ccix <- !vapply(seq_len(NROW(data[iix,,drop=FALSE])), function(x) {
    #   return(all(is.na(data[iix,,drop=FALSE][x,])) || all(is.na(drugpheno[iix,,drop=FALSE][x,])) || all(is.na(type[iix][x])) || all(is.na(batch[iix][x])))
    # }, FUN.VALUE=logical(1))

    if (sum(ccix) < 3) {
      ## not enough experiments
      rest <- list(matrix(NA, nrow=ncol(data), ncol=length(nc), dimnames=list(colnames(data), nc)))
      res <- c(res, rest)
    } else {
      # splitix <- parallel::splitIndices(nx=ncol(data), ncl=nthread)
      # splitix <- splitix[vapply(splitix, length, FUN.VALUE=numeric(1)) > 0]
      mcres <- parallel::mclapply(seq_len(ncol(data)), function(x, data, type, batch, drugpheno, standardize, modeling.method, inference.method, req_alpha) {
        if(modeling.method == "anova"){
          res <- t(apply(data[ , x, drop=FALSE], 2, geneDrugSensitivity, type=type, batch=batch, drugpheno=drugpheno, verbose=verbose, standardize=standardize))
        } else if(modeling.method == "pearson") {
          if(!is.character(data)){
            res <- t(apply(data[ , x, drop=FALSE], 2, geneDrugSensitivityPCorr,
                                                      type=type,
                                                      batch=batch,
                                                      drugpheno=drugpheno,
                                                      verbose=verbose,
                                                      test=inference.method,
                                                      req_alpha = req_alpha))
          } else {
            res <- t(apply(data[ , x, drop=FALSE], 2, function(dataIn) {
              geneDrugSensitivityPBCorr(as.factor(dataIn),
                                                      type=type,
                                                      batch=batch,
                                                      drugpheno=drugpheno,
                                                      verbose=verbose,
                                                      test=inference.method,
                                                      req_alpha = req_alpha)}))
          }

        }


        return(res)
      }, data=data[iix, , drop=FALSE],
         type=type[iix], batch=batch[iix],
         drugpheno=drugpheno[iix,,drop=FALSE],
         standardize=standardize,
         modeling.method = modeling.method,
         inference.method = inference.method,
         req_alpha = req_alpha, mc.cores = nthread, mc.preschedule = TRUE)
      rest <- do.call(rbind, mcres)
      rest <- cbind(rest, "fdr"=p.adjust(rest[ , "pvalue"], method="fdr"))
      # rest <- rest[ , nc, drop=FALSE]
      res <- c(res, list(rest))
    }
  }
  names(res) <- names(ltype)
  return(res)
}

## End
