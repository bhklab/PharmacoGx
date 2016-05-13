## adaptive Matthews correlation coefficient for binary classification
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
        mm[mccix, "p.value"] <- mcc(x=x3, y=y3, nperm=nperm, nthread=nthread)[["p.value"]]
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


