## adaptive Matthews correlation coefficient for binary classification
#' Calculate an Adaptive Matthews Correlation Coefficient
#' 
#' This function calculates an Adaptive Matthews Correlation Coefficient (AMCC) for two vectors of values identical length. It assumes the entries in the two vectors are paired.
#' The Adaptive Matthews Correlation Coefficient for two vectors of values is defined as the Maximum Matthews Coefficient over all possible 
#' binary splits of the ranks of the two vectors. In this way, it calculates the best possible agreement of a binary 
#' classifier on the two vectors of data. #If the AMCC is low, then it is impossible to find any binary classification of the two vectors with a high degree of concordance.   
#' 
#' @param x,y Two paired vectors of values. Could be replicates of observations for the same experiments for example.  
#' @param step.prct Instead of testing all possible splits of the data, it is possible to test steps of a percentage size of the total number of ranks in x/y. If this variable is 0, function defaults to testing all possible splits.
#' @param min.cat The minimum number of complete cases between x and y for which the function will attempt to calculate an AMCC. If the number of complete cases is smaller, the function will return NA for the AMCC and NA for the mcc table. 
#'                Useful for cases when the number of complete cases cannot be determined beforehand, but a minimum threshold should be enforced.
#' @param nperm The number of perumatation to use for estimating significance. If 0, then no p-value is calculated. 
#' @param setseed Allows setting a consitent seed for reproducibility of permutation testing results. Defaults to 12345.
#' @param nthread Number of threads to parallize over. Both the AMCC calculation and the permutation testing is done in parallel. 
#' @return Returns a list with two elements. $amcc contains the highest "mcc" value over all the splits, the p value, as well as the rank at which the split was done. 
#' @import parallel
#' @export
amcc <- 
  function(x, y, step.prct=0, min.cat=3, nperm=1000, setseed=12345, nthread=1) {
    ccix <- complete.cases(x, y)
    if (sum(ccix) >= (2 * min.cat)) {
      x2 <- rank(-x[ccix], ties.method="first")
      y2 <- rank(-y[ccix], ties.method="first")
      ## compute mcc for each rank
      iix <- 1:(min(max(x2), max(y2)) - 1)
      if (step.prct > 0) {
        iix <- round(quantile(iix, probs=seq(0, 1, by=step.prct)))
      }
      splitix <- parallel::splitIndices(nx=length(iix), ncl=nthread)
      splitix <- splitix[sapply(splitix, length) > 0]
      mcres <- parallel::mclapply(splitix, function(x, iix, x2, y2) {
        res <- t(sapply(iix[x], function(x, x2, y2) {
          x3 <- factor(ifelse (x2 <= x, "1", "0"))
          y3 <- factor(ifelse (y2 <= x, "1", "0"))
          res <- mcc(x=x3, y=y3, nperm=0, nthread=1)
          return (res)
        }, x2=x2, y2=y2))
        return (res)
      }, iix=iix, x2=x2, y2=y2)
      mm <- do.call(rbind, mcres)
      mode(mm) <- "numeric"
      ## remove extreme indices
      rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 2):nrow(mm))
      mccix <- max(which(mm[-rmix, "estimate", drop=FALSE] == max(mm[-rmix, "estimate", drop=FALSE], na.rm=TRUE))) + (min.cat - 1)
      ## compute significance only for the AMCC
      x3 <- factor(ifelse (x2 <= mccix, "1", "0"))
      y3 <- factor(ifelse (y2 <= mccix, "1", "0"))
      if (nperm > 0) {
        mm[mccix, "p.value"] <- mcc(x=x3, y=y3, nperm=nperm, nthread=nthread, setseed=setseed)[["p.value"]]
        ## bonferronni correction
        # mm[mccix, "p"] <- mm[mccix, "p"] * length(x3)
      }
      if (!is.na(mm[mccix, "p.value"]) && mm[mccix, "p.value"] > 1) { mm[mccix, "p.value"] <- 1 }
      res <- c("mcc"=mm[mccix, ], "n1"=mccix, "n2"=nrow(mm) - mccix, "n"=nrow(mm))
    } else {
      res <- c("mcc"=NA, "p"=NA, "n1"=0, "n2"=0, "n"=sum(ccix))
      mm <- NA
    }
    names(res) <- c("mcc", "p", "n1", "n2", "n")
    return(list("amcc"=res, "mcc"=mm))
  }


