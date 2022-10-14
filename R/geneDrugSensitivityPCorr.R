
cor.boot <- function(data, w){
  ddd <- data[w,]
        ## A question here is what to do when our bootstrap sample has 0 variance in one of
      ## the two variables (usually the molecular feature)
      ## If we return NA, we effectively condition on a sampling procedure that samples at least
      ## one expressor. If we return near 0 (the default of coop::pcor), then we effectively say that conditioned
      ## on not sampling any expressors, there is no association. Neither is correct, but the latter is certainly more
      ## conservative. We probably want to use a completely different model to discover "rare" biomarkers
      ## We will go with the later conservative option, however we will set it to zero exactly instead of relying on
      ## the behaviour of coop.

  if(length(unique(ddd[,1]))<2 || length(unique(ddd[,2]))<2){
    return(0)
  }

  return(coop::pcor(ddd[,1], ddd[,2]))
}


#' Calculate The Gene Drug Sensitivity
#'
#' This version of the function uses a partial correlation instead of standardized linear models.
#'
#' @param x A \code{numeric} vector of gene expression values
#' @param type A \code{vector} of factors specifying the cell lines or type types
#' @param batch A \code{vector} of factors specifying the batch
#' @param drugpheno A \code{numeric} vector of drug sensitivity values (e.g.,
#'   IC50 or AUC)
#' @param test A \code{character} string indicating whether resampling or analytic based tests should be used
#' @param req_alpha \code{numeric}, number of permutations for p value calculation
#' @param nBoot \code{numeric}, number of bootstrap resamplings for confidence interval estimation
#' @param conf.level \code{numeric}, between 0 and 1. Size of the confidence interval required
#' @param max_perm \code{numeric} the maximum number of permutations that QUICKSTOP can do before giving up and returning NA.
#'   Can be set globally by setting the option "PharmacoGx_Max_Perm", or left at the default of \code{ceiling(1/req_alpha*100)}.
#' @param verbose \code{boolean} Should the function display messages?
#'
#' @return A \code{vector} reporting the effect size (estimateof the coefficient
#'   of drug concentration), standard error (se), sample size (n), t statistic,
#'   and F statistics and its corresponding p-value.
#'
#' @examples
#' print("TODO::")
#'
#' @importFrom stats sd complete.cases lm glm anova pf formula var pt qnorm cor residuals runif
#' @importFrom boot boot boot.ci
#' @importFrom coop pcor
geneDrugSensitivityPCorr <- function(x, type, batch, drugpheno,
  test = c("resampling", "analytic"),
  req_alpha = 0.05,
  nBoot = 1e3,
  conf.level = 0.95,
  max_perm = getOption("PharmacoGx_Max_Perm", ceiling(1/req_alpha*100)),
  verbose=FALSE) {

  test <- match.arg(test)

  colnames(drugpheno) <- paste("drugpheno", seq_len(ncol(drugpheno)), sep=".")

  drugpheno <- data.frame(vapply(drugpheno, function(x) {
    if (!is.factor(x)) {
      x[is.infinite(x)] <- NA
    }
    return(list(x))
  }, USE.NAMES=TRUE,
  FUN.VALUE=list(1)), check.names=FALSE)


  ccix <- complete.cases(x, type, batch, drugpheno)
  nn <- sum(ccix)

  rest <- c("estimate"=NA_real_, "n"=as.numeric(nn), "df"=NA_real_, significant = NA_real_,"pvalue"=NA_real_, "lower" = NA_real_, "upper" = NA_real_)

  if(nn <= 3 || isTRUE(all.equal(var(x[ccix], na.rm=TRUE), 0))) {
    ## not enough samples with complete information or no variation in gene expression
    return(rest)
  }

  drugpheno <- drugpheno[ccix,,drop=FALSE]


  xx <- x[ccix]

  if(ncol(drugpheno)>1){
    stop("Partial Correlations not implemented for multiple output")
  } else {
    ffd <- "drugpheno.1 ~ . - x"
    ffx <- "x ~ . - drugpheno.1"
  }

  # ff1 <- sprintf("%s + x", ff0)

  dd <- data.frame(drugpheno, "x"=xx)

  ## control for tissue type
  if(length(sort(unique(type[ccix]))) > 1) {
    dd <- cbind(dd, type=type[ccix])
  }
  ## control for batch
  if(length(sort(unique(batch[ccix]))) > 1) {
    dd <- cbind(dd, batch=batch[ccix])
  }
  ## control for duration
  # if(length(sort(unique(duration))) > 1){
  #   ff0 <- sprintf("%s + duration", ff0)
  #   ff <- sprintf("%s + duration", ff)
  # }

  # if(is.factor(drugpheno[,1])){

  #   drugpheno <- drugpheno[,1]

  # } else {

  #   drugpheno <- as.matrix(drugpheno)

  # }
  if(any(unlist(lapply(drugpheno,is.factor)))){

    stop("Currently only continous output allowed for partial correlations")


  } else {

    if(ncol(dd) > 2){
      lm1 <- lm(formula(ffd), dd)
      var1 <- residuals(lm1)
      var2 <- residuals(lm(formula(ffx), dd))
      df <- lm1$df - 2L # taking the residual degrees of freedom minus 2 parameters estimated for pearson cor.
    } else { ## doing this if statement in the case there are some numerical differences between mean centred values and raw values
    var1 <- dd[,"drugpheno.1"]
    var2 <- dd[,"x"]
    df <- nn - 2L
  }

  obs.cor <- coop::pcor(var1, var2, use="complete.obs")


    ## NB: just permuting the residuals would leads to Type I error inflation,
    ## from an underestimation due to ignoring variance in the effects of the covariates.
    ## See: https://www.tandfonline.com/doi/abs/10.1080/00949650008812035
    ## Note that the above paper does not provide a single method recommended in all cases
    ## We apply the permutation of raw data method, as it is most robust to small sample sizes
  if(test == "resampling"){
      ## While the logic is equivalent regardless of if there are covariates for calculating the point estimate,
      ## (correlation is a subcase of partial correlation), for computational efficency in permuation testing we
      ## split here and don't do extranous calls to lm if it is unnecessay.

    if(ncol(dd) > 2){

      if(!getOption("PharmacoGx_useC")|| ncol(dd)!=3){ ## currently implementing c code only for 1 single grouping variable

        ## implementing a much more efficient method for the particular case where we have 3 columns with assumption that
        ## column 3 is the tissue.
        if(ncol(dd)==3){
          sample_function <- function(){

            partial.dp <- sample(dd[,1], nrow(dd))
            partial.x <- sample(dd[,2], nrow(dd))

            for(gp in unique(dd[,3])){
              partial.x[dd[,3]==gp] <- partial.x[dd[,3]==gp]-mean(partial.x[dd[,3]==gp])
              partial.dp[dd[,3]==gp] <- partial.dp[dd[,3]==gp]-mean(partial.dp[dd[,3]==gp])
            }

            perm.cor <- coop::pcor(partial.dp, partial.x, use="complete.obs")
            return(abs(obs.cor) < abs(perm.cor))
          }
        } else {
          sample_function <- function(){
            # browser()
            dd2 <- dd
            dd2[,1] <- sample(dd[,1], nrow(dd))
            dd2[,2] <- sample(dd[,2], nrow(dd))

            partial.dp <- residuals(lm(formula(ffd), dd2))
            partial.x <- residuals(lm(formula(ffx), dd2))

            perm.cor <- coop::pcor(partial.dp, partial.x, use="complete.obs")
            return(abs(obs.cor) < abs(perm.cor))
          }
        }

        p.value <- corPermute(sample_function, req_alpha = req_alpha, max_iter=max_perm)
        significant <- p.value$significant
        p.value <- p.value$p.value


      } else {

        x <- dd[,1]
        y <- dd[,2]
        GR <- as.integer(factor(dd[,3]))-1L
        GS <- as.integer(table(factor(dd[,3])))
        NG <- length(table(factor(dd[,3])))
        N <- as.numeric(length(x))

        p.value <-PharmacoGx:::partialCorQUICKSTOP(x, y, obs.cor, GR, GS, NG, max_perm, N, req_alpha, req_alpha/100, 10L, runif(2))
        significant <- p.value[[1]]
        p.value <- p.value[[2]]
      }


    pcor.boot <- function(ddd, w){
      ddd <- ddd[w,]
          ## Taking care of an edge case where only one factor level is left after resampling
          ## However, we need to keep the first two numeric columns to properly return a value, otherwise
          ## if we remove gene expression because there were only non-detected samples, for example,
          ## we will try to take the correlation against a character vector.
      ddd <- ddd[,c(TRUE, TRUE, apply(ddd[,-c(1,2),drop=FALSE], 2, function(x) return(length(unique(x))))>=2)]


      ## A question here is what to do when our bootstrap sample has 0 variance in one of
      ## the two variables (usually the molecular feature)
      ## If we return NA, we effectively condition on a sampling procedure that samples at least
      ## one expressor. If we return near 0 (the default of coop::pcor), then we effectively say that conditioned
      ## on not sampling any expressors, there is no association. Neither is correct, but the latter is certainly more
      ## conservative. We probably want to use a completely different model to discover "rare" biomarkers
      ## We will go with the later conservative option, however we will set it to zero exactly instead of relying on
      ## the behaviour of coop.

      if(length(unique(ddd[,1]))<2 || length(unique(ddd[,2]))<2){
        return(0)
      }

      if(ncol(ddd)==3){
        partial.dp <- ddd[,1]
        partial.x <- ddd[,2]
        for(gp in unique(ddd[,3])){
          partial.x[ddd[,3]==gp] <- partial.x[ddd[,3]==gp]-mean(partial.x[ddd[,3]==gp])
          partial.dp[ddd[,3]==gp] <- partial.dp[ddd[,3]==gp]-mean(partial.dp[ddd[,3]==gp])
        }

      } else if(ncol(ddd)==2){
        partial.dp <- ddd[,1]
        partial.x <- ddd[,2]
      } else {

        partial.dp <- residuals(lm(formula(ffd), ddd))
        partial.x <- residuals(lm(formula(ffx), ddd))

      }

      return(coop::pcor(partial.dp, partial.x, use="complete.obs"))
    }

    boot.out <- boot(dd, pcor.boot, R=nBoot)
    cint <- tryCatch(boot.ci(boot.out, conf = conf.level, type="bca")$bca[,4:5],
      error = function(e) {
        if(e$message == "estimated adjustment 'w' is infinite"){
          warning("estimated adjustment 'w' is infinite for some features")
          return(c(NA_real_,NA_real_))
        } else {
          stop(e)
        }
      })
  } else {
    if(!getOption("PharmacoGx_useC")){
      sample_function <- function(){
        v1 <- sample(var1)
        return(abs(obs.cor) < abs(coop::pcor(v1, var2, use="complete.obs")))
      }

      p.value <- corPermute(sample_function, req_alpha = req_alpha, max_iter=max_perm)
      significant <- p.value$significant
      p.value <- p.value$p.value
    } else {

      x <- as.numeric(var1)
      y <- as.numeric(var2)
      GR <- rep(0L, length(x))
      GS <- as.integer(length(x))
      NG <- 1L
      N <- as.numeric(length(x))

      p.value <-PharmacoGx:::partialCorQUICKSTOP(x, y, obs.cor, GR, GS, NG, max_perm, N, req_alpha, req_alpha/100, 10L, runif(2))
      significant <- p.value[[1]]
      p.value <- p.value[[2]]
    }




    boot.out <- boot(data.frame(var1, var2), cor.boot, R=nBoot)
    cint <- tryCatch(boot.ci(boot.out, conf = conf.level, type="bca")$bca[,4:5],
      error = function(e) {
        if(e$message == "estimated adjustment 'w' is infinite" || e$message == "Error in if (const(t, min(1e-08, mean(t, na.rm = TRUE)/1e+06))) { : \n  missing value where TRUE/FALSE needed\n"){
          warning("estimated adjustment 'w' is infinite for some features")
          return(c(NA_real_,NA_real_))
        } else {
          stop(e)
        }
      })
  }


} else if(test == "analytic"){
      # if(ncol(dd) > 2){
      #   df <- nn - 2L - controlled.var

      # } else {
      #   df <- nn - 2L
      #   # cor.test.res <- cor.test(dd[,"drugpheno.1"], dd[,"x"], method="pearson", use="complete.obs")
      # }
  stat <- sqrt(df) * obs.cor/sqrt(1-obs.cor^2) ## Note, this is implemented in same order of operations as cor.test
  p.value <- 2*min(pt(stat, df=df), pt(stat, df=df, lower.tail = FALSE))
      ## Implementing with fisher transform and normal dist for consistency with R's cor.test
  z <- atanh(obs.cor)
  sigma <- 1/sqrt(df - 1)
  cint <- tanh(z + c(-1, 1) * sigma * qnorm((1 + conf.level)/2))
  significant <- p.value < req_alpha
}

}

rest <- c("estimate"=obs.cor, "n"=nn, df=df, significant = as.numeric(significant), "pvalue"=p.value, "lower" = cint[1], "upper" = cint[2])


return(rest)
}
