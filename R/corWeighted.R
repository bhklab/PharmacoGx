#################################################
## Compute a weighted correlation coefficient
## inspired from package boot
##
## inputs:    
##      - 
##
## outputs, a list of items:
##            - 
##
#################################################

#' @importFrom stats cov.wt complete.cases

corWeighted <- 
function (x, y, w, method=c("pearson", "spearman"), alternative=c("two.sided", "greater", "less"), nperm=0, nthread=1, setseed, na.rm=FALSE) {
 
  ######################
  wcor <- function (d, w, na.rm=TRUE) {
### NOTE::: THIS FORMULA CAN SUFFER CATASTROPHIC CANCELATION AND SHOULD BE FIXED!!!
#     s <- sum(w, na.rm=na.rm)
#     m1 <- sum(d[ , 1L] * w, na.rm=na.rm) / s
#     m2 <- sum(d[ , 2L] * w, na.rm=na.rm) / s
#     res <- (sum(d[ , 1L] * d[ , 2L] * w, na.rm=na.rm) / s - m1 * m2) / sqrt((sum(d[ , 1L]^2 * w, na.rm=na.rm) / s - m1^2) * (sum(d[ , 2L]^2 * w, na.rm=na.rm) / s - m2^2))
    CovM <- cov.wt(d, wt=w)[["cov"]]
    res <- CovM[1,2]/sqrt(CovM[1,1]*CovM[2,2])
    return (res)
  }
  
  ######################
  
  if (missing(w)) { w <- rep(1, length(x)) / length(x) }
  if (length(x) != length(y) || length(x) != length(w)) { stop("x, y, and w must have the same length") }
  method <- match.arg(method)
  if (method == "spearman") {
    x <- rank(x)
    y <- rank(y)
  }
  alternative <- match.arg(alternative)
  
  res <- c("rho"=NA, "p"=NA)
  
  ## remove missing values
  ccix <- complete.cases(x, y, w)
  if(!all(ccix) && !na.rm) { warning("Missing values are present") }
  if(sum(ccix) < 3) {
    return(res)
  }
  x <- x[ccix]
  y <- y[ccix]
  w <- w[ccix]
  
  wc <- wcor(d=cbind(x, y), w=w)
  res["rho"] <- wc
  if (nperm > 1) {
    if (!missing(setseed)) { set.seed(setseed) }
    splitix <- parallel::splitIndices(nx=nperm, ncl=nthread)
    if (!is.list(splitix)) { splitix <- list(splitix) }
    splitix <- splitix[sapply(splitix, length) > 0]
    mcres <- parallel::mclapply(splitix, function(x, xx, yy, ww) {
      pres <- sapply(x, function(x, xx, yy, ww) {
        ## permute the data and the weights
        d2 <- cbind(xx[sample(1:length(xx))], yy[sample(1:length(yy))])
        w2 <- ww[sample(1:length(ww))]
        return(wcor(d=d2, w=w2))
      }, xx=xx, yy=yy, ww=ww)
      return(pres)
    }, xx=x, yy=y, ww=w)
    perms <- do.call(c, mcres)
    
    switch (alternative,
      "two.sided" = { 
        if (res["rho"] < 0) { p <- sum(perms <= res, na.rm=TRUE) } else { p <- sum(perms >= res, na.rm=TRUE) }
        if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
        p <- p * 2
      },
      "greater" = {
        p <- sum(perms >= res, na.rm=TRUE) 
        if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
      },
      "less" = {
        p <- sum(perms <= res, na.rm=TRUE) 
        if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
      })
    res["p"] <- p
  }
  return(res)
}

## End
