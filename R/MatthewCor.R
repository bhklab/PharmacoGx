## Matthews correlatipon coefficient
#' Compute a Mathews Correlation Coefficient 
#' 
#' The function computes a Matthews correlation coefficient for two factors provided to the function. It assumes each factor is a factor of class labels, 
#' and the enteries are paired in order of the vectors.
#' 
#' @examples
#' x <- factor(c(1,2,1,2,3,1))
#' y <- factor(c(2,1,1,1,2,2))
#' mcc(x,y)
#' 
#' @param x,y factor of the same length with the same number of levels
#' @param nperm number of permutations for significance estimation. If 0, no permutation testing is done
#' @param setseed seed for permutation testing
#' @param nthread can parallelize permutation texting using parallel's mclapply 
#' @return A list with the MCC as the $estimate, and p value as $p.value
#' @export
mcc <- 
  function(x, y, nperm=1000, setseed=12345, nthread=1) {
    set.seed(setseed)
    if ((length(x) != length(y)) || (!is.factor(x) || length(levels(x)) < 2) || (!is.factor(y) || length(levels(y)) < 2)) { stop("x and y must be factors of the same length with at least two levels") }
    if(length(levels(x))!= length(levels(y))){

      warning("The number of levels x and y was different. Taking the union of all class labels.")
      levels(x) <- union(levels(x), levels(y))
      levels(y) <- union(levels(x), levels(y))

    }
    res <- list("estimate"=NA, "p.value"=NA)
    ## compute MCC
    res[["estimate"]] <- .mcc(ct=table(x, y))
    ## compute significance of MCC using a permutation test
    if (nperm > 0) {
      splitix <- parallel::splitIndices(nx=nperm, ncl=nthread)
      splitix <- splitix[sapply(splitix, length) > 0]
      mcres <- parallel::mclapply(splitix, function(x, xx, yy) {
        res <- sapply(x, function(x, xx, yy) {
          xx <- sample(xx)
          yy <- sample(yy)
          return (.mcc(ct=table(xx, yy)))
        }, xx=xx, yy=yy)
        return (res)
      }, xx=x, yy=y)
      mcres <- unlist(mcres)
      res[["p.value"]] <- sum(mcres > res["estimate"]) / sum(!is.na(mcres))
      if (res["p.value"] == 0) { res["p.value"] <- 1 / (nperm + 1) }
    }
    return (res)
  }
.mcc <- 
  function(ct, nbcat=nrow(ct)) {
    if(nrow(ct) != ncol(ct)) { stop("the confusion table should be square!") }
    if(!(sum(ct)==sum(diag(ct))) &&  (length(which(apply(ct, 1, sum) == 0)) == (nbcat-1) & ((length(which(apply(ct, 2, sum) == 0)) != (nbcat-1)) | (length(which(apply(ct, 2, sum) == 0)) == (nbcat-1)))) || (length(which(apply(ct, 2, sum) == 0)) == (nbcat-1) & ((length(which(apply(ct, 1, sum) == 0)) != (nbcat-1)) | (length(which(apply(ct, 1, sum) == 0)) == (nbcat-1)) & sum(diag(ct)) == 0))) { ct <- ct + matrix(1, ncol=nbcat, nrow=nbcat) } ### add element to categories if nbcat-1 predictive categories do not contain elements. Not in case where all are correct!
    
    if(sum(ct, na.rm=TRUE) <= 0) { return(NA) }
    
    myx <- matrix(TRUE, nrow=nrow(ct), ncol=ncol(ct))
    diag(myx) <- FALSE
    if(sum(ct[myx]) == 0) { return(1) }
    myperf <- 0
    for(k in 1:nbcat) {
      for(m in 1:nbcat) {
        for(l in 1:nbcat) {
          myperf <- myperf + ((ct[k, k] * ct[m, l]) - (ct[l, k] * ct[k, m]))
        }
      }
    }
    aa <- 0
    for(k in 1:nbcat) {
      cc <- 0
      for(l in 1:nbcat) { cc <- cc + ct[l, k] }
      dd <- 0
      for(f in 1:nbcat) {
        for(g in 1:nbcat) { if(f != k) { dd <- dd + ct[g, f] } }
      }
      aa <- aa + (cc * dd)
    }
    bb <- 0
    for(k in 1:nbcat) {
      cc <- 0
      for(l in 1:nbcat) { cc <- cc + ct[k, l] }
      dd <- 0
      for(f in 1:nbcat) {
        for(g in 1:nbcat) { if(f != k) { dd <- dd + ct[f, g] } }
      }
      bb <- bb + (cc * dd)
    }
    
    myperf <- myperf / (sqrt(aa) * sqrt(bb))
    return(myperf)
  }


