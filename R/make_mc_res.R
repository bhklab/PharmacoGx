make_mc_res <- function() {
    res <- NULL
    for(i in drugn) {
      ## using a linear model (x ~ concentration + cell + batch)
      dd <- lapply(drugpheno, function(x) x[ ,i]) #Can you write this down? Drugpheno is a list of dataframes, each one of which must be split columnwise you can stop writing.
      dd <- do.call(cbind, dd)
      colnames(dd) <- seq_len(ncol(dd))
      if(!is.na(sensitivity.cutoff)) {
        dd <- factor(ifelse(dd > sensitivity.cutoff, 1, 0), levels = c(0, 1))
      }
      rr <- rankGeneDrugSensitivity(data = expr,
                                    drugpheno = dd,
                                    type = type,
                                    batch = batch,
                                    single.type = FALSE,
                                    standardize = standardize,
                                    nthread = nthread,
                                    verbose = verbose)
      res <- c(res, list(rr$all))
    }
    names(res) <- drugn
    return(res)
}