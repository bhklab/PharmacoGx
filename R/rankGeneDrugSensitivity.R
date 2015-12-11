
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

rankGeneDrugSensitivity <- function (data, drugpheno, type, batch, single.type=FALSE, nthread=1, verbose=FALSE) {
  if (nthread != 1) {
    availcore <- parallel::detectCores()
    if (missing(nthread) || nthread < 1 || nthread > availcore) {
      nthread <- availcore
    }
  }
  if (missing(type) || all(is.na(type))) {
    type <- array("other", dim=nrow(data), dimnames=list(rownames(data)))
  }
  if (missing(batch) || all(is.na(batch))) {
    batch <- array(1, dim=nrow(data), dimnames=list(rownames(data)))
  }
    if (any(c(length(drugpheno), length(type), length(batch)) != nrow(data))) {
    stop("length of drugpheno, type, duration, and batch should be equal to the number of rows of data!")
  }
    names(drugpheno) <- names(type) <- names(batch) <- rownames(data)
  
  res <- NULL
  utype <- sort(unique(as.character(type)))
  ltype <- list("all"=utype)
  if (single.type) {
    ltype <- c(ltype, as.list(utype))
    names(ltype)[-1] <- utype
  }
  res <- NULL
#  nc <- c("estimate", "se", "n", "tstat", "fstat", "pvalue", "fdr")
  nc <- c("estimate", "se", "n", "pvalue", "fdr")
  
  for (ll in 1:length(ltype)) {
    iix <- !is.na(type) & is.element(type, ltype[[ll]])
    ccix <- complete.cases(data[iix, , drop=FALSE], drugpheno[iix], type[iix], batch[iix])
    if (sum(ccix) < 3) {
      ## not enough experiments
      rest <- list(matrix(NA, nrow=ncol(data), ncol=length(nc), dimnames=list(colnames(data), nc)))
      res <- c(res, rest)
    } else {
      splitix <- parallel::splitIndices(nx=ncol(data), ncl=nthread)
      splitix <- splitix[sapply(splitix, length) > 0]
      mcres <- parallel::mclapply(splitix, function(x, data, type, batch, drugpheno) {
        res <- t(apply(data[ , x, drop=FALSE], 2, geneDrugSensitivity, type=type, batch=batch, drugpheno=drugpheno, verbose=verbose))
        return(res)
      }, data=data[iix, , drop=FALSE], type=type[iix], batch=batch[iix], drugpheno=drugpheno[iix])
      rest <- do.call(rbind, mcres)
      rest <- cbind(rest, "fdr"=p.adjust(rest[ , "pvalue"], method="fdr"))
      rest <- rest[ , nc, drop=FALSE]
      res <- c(res, list(rest))
    }
  }
  names(res) <- names(ltype)
    return(res)
}

## End
