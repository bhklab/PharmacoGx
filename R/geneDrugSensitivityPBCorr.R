log_denom <- function(suc, total, p){
  tmp <- 0;
    # if(log10(p)>-32)
    # {
  tmp <- (suc)*log(p)
    # } else {
    #  warning("p reached precision threshold")
    # }
    # if(log10(1.0-p)>-32)
    # {
  tmp <- tmp + (total - suc)*log(1 - p)
    # } else {
    #    warning("1-p reached precision threshold")
    # }
  return(tmp)

}
log_denom <- function(suc, total, p){
  return((suc)*log(p) + (total - suc)*log(1 - p))

}

### Implementing algorithm from quick stop paper

## decision boundary is inverse of error prob
## do everything in log scale because of numerical precision.
corPermute <- function(sample_function, req_alpha=0.05, tolerance_par = req_alpha*0.001, log_decision_boundary = 10, max_iter = 1/req_alpha*100){

  num.larger <- 0

  cur_success <- 0
  cur_iter <- 1

  log_cur_PiN <- log(1) # Initialization so the logic can stay the same through all loops

  p1 <- req_alpha
  p2 <- p1 + tolerance_par

  pr_min_1 <- 1/2

  while(cur_iter < max_iter){
    # vec1 <- sample(vec1)
    # perm.cor <- cor(vec1, vec2, use="complete.obs")

    # success <- abs(perm.cor) > abs(obs.cor)

    success <- sample_function()

    if(success){
      cur_success <- cur_success + 1
      log_cur_PiN <- log_cur_PiN + log_denom(1,1,pr_min_1)
    } else {
      log_cur_PiN <- log_cur_PiN + log_denom(0,1,pr_min_1)
    }
    # if(pr_min_1 >= p2){
    #   log_cur_suph1 <- log_denom(cur_success, cur_iter, p1)
    #   log_cur_suph2 <- log_denom(cur_success, cur_iter, pr_min_1)
    # } else if(pr_min_1 <= p1){
    #   log_cur_suph1 <- log_denom(cur_success, cur_iter, pr_min_1)
    #   log_cur_suph2 <- log_denom(cur_success, cur_iter, p2)
    # } else {
    #   log_cur_suph1 <- log_denom(cur_success, cur_iter, p1)
    #   log_cur_suph2 <- log_denom(cur_success, cur_iter, p2)
    # }
    if(pr_min_1<p1) {
      log_cur_suph1 <- log_denom(cur_success, cur_iter, pr_min_1)
    } else {
      log_cur_suph1 <- log_denom(cur_success, cur_iter, p1)
    }
    if(pr_min_1>p2) {
      log_cur_suph2 <- log_denom(cur_success, cur_iter, pr_min_1)
    } else {
      log_cur_suph2 <- log_denom(cur_success, cur_iter, p2)
    }

    cur_iter <- cur_iter + 1
    pr_min_1 <- (cur_success + 1/2)/cur_iter

    # if(cur_success == 0){
    #   next
    # }

    if(log_cur_PiN - log_cur_suph2 > log_decision_boundary){
      return(list(significant = TRUE, "p.value" = pr_min_1, num_iter=cur_iter, num_larger=cur_success))
    }
    if(log_cur_PiN - log_cur_suph1 > log_decision_boundary){
      return(list(significant = FALSE, "p.value" = pr_min_1, num_iter=cur_iter, num_larger=cur_success))
    }
  }
  return(list(significant = NA, "p.value" = pr_min_1, num_iter=cur_iter, num_larger=cur_success))
}



## Helper Functions
##TODO:: Add  function documentation
#' @importFrom stats quantile
.rescale <- function(x, na.rm=FALSE, q=0)
{
  if(q == 0) {
    ma <- max(x, na.rm=na.rm)
    mi <- min(x, na.rm=na.rm)
  } else {
    ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
    mi <- quantile(x, probs=q/2, na.rm=na.rm)
  }
  xx <- (x - mi) / (ma - mi)
  return(xx)
}

## TODO: decide better what to do with no variance cases.
cor.boot <- function(data, w){
  return(cor(data[w,1], data[w,2]))
}


#' Calculate The Gene Drug Sensitivity
#'
#' This version of the function uses a partial correlation instead of standardized linear models, for discrete predictive features
#' Requires at least 3 observations per group.
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
#' @importFrom stats sd complete.cases lm glm anova pf formula var
#' @importFrom boot boot boot.ci
geneDrugSensitivityPBCorr <- function(x, type, batch, drugpheno,
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

  if(nn <= 3 || all(duplicated(x)[-1L])) {
    ## not enough samples with complete information or no variation in gene expression
    return(rest)
  }

  ## taking at least 5 times the number of boot samples as number of observations, as with binary data many boot samples
  ## tend to be NA, and if the number of non-NA drops below number of obs, the emperical influence cannot be determined
  ## for BCA interval calculation

  if(test=="resampling"){
    nBoot <- max(nBoot, nn*5)
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


  if(!is.factor(dd[[2]])){
    stop("Molecular Feature is not discrete, but point biserial correlation was requested")
  } else if(length(unique(dd[[2]]))>2) {

    stop('More than two discrete settings for moleuclar feature not currently supported')

  } else if(length(unique(dd[[2]]))==1 || min(table(dd[[2]]))<3) {
    warning("Some features had less than 3 observations per category, returning NA.")
    return(rest)

  } else{
    dd[[2]] <- as.numeric(dd[[2]]) ## converting to numeric codings for downstream modeling
  }


  if(any(unlist(lapply(drugpheno,is.factor)))){

    stop("Currently only continous output allowed for point biserial correlations")


  } else{

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

  obs.cor <- cor(var1, var2, use="complete.obs")


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

      # if(!getOption("PharmacoGx_useC")|| ncol(dd)!=3){ ## not yet implemented

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

          perm.cor <- cor(partial.dp, partial.x, use="complete.obs")
          return(abs(obs.cor) < abs(perm.cor))
        }
      } else {
        sample_function <- function(){

          dd2 <- dd
          dd2[,1] <- sample(dd[,1], nrow(dd))
          dd2[,2] <- sample(dd[,2], nrow(dd))

          partial.dp <- residuals(lm(formula(ffd), dd2))
          partial.x <- residuals(lm(formula(ffx), dd2))

          perm.cor <- cor(partial.dp, partial.x, use="complete.obs")
          return(abs(obs.cor) < abs(perm.cor))
        }
      }

      p.value <- corPermute(sample_function, req_alpha = req_alpha, max_iter=max_perm)
      significant <- p.value$significant
      p.value <- p.value$p.value


      # } else {

      #   x <- dd[,1]
      #   y <- dd[,2]
      #   GR <- as.integer(factor(dd[,3]))-1L
      #   GS <- as.integer(table(factor(dd[,3])))
      #   NG <- length(table(factor(dd[,3])))
      #   N <- as.numeric(length(x))

      #   p.value <-PharmacoGx:::partialCorQUICKSTOP(x, y, obs.cor, GR, GS, NG, 1e7,N, req_alpha, req_alpha/100, 10L, runif(2))
      #   significant <- p.value[[1]]
      #   p.value <- p.value[[2]]
      # }


      pcor.boot <- function(ddd, w){
        ddd <- ddd[w,]
          ## Taking care of an edge case where only one covariate factor level is left after resampling
        ddd[,-c(1,2)] <- ddd[,-c(1,2),drop=FALSE][,apply(ddd[,-c(1,2),drop=FALSE], 2, function(x) return(length(unique(x))))>=2]

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

        return(cor(partial.dp, partial.x, use="complete.obs"))
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
    # if(!getOption("PharmacoGx_useC")){

      ## At this point we have verified that we are doing the normal (nor partial) PBCC,
      ## and we also verified that only 2 unique values of var2 exist. Therefore, diff
      ## should return a single result.
      ## Note that the PBCC permutation only depends on the mean differences, the
      ## denominator is proprtional to the total variance
      ## in var1 and inverse of the sqrt of the proportions between groups,
      ## both of which stay constant through the permutation. Therefore, we skip
      ## the needless normalization step in this permutation procedure.
      ## Note that this does not apply to bootstrapping.

      obs.mean.diff <- diff(tapply(var1, var2, mean))
      sample_function <- function(){
        v1 <- sample(var1)
        return(abs(obs.mean.diff) < abs(diff(tapply(v1, var2, mean))))
      }

      p.value <- corPermute(sample_function, req_alpha = req_alpha, max_iter=max_perm)
      significant <- p.value$significant
      p.value <- p.value$p.value
    # } else {

    #   x <- as.numeric(var1)
    #   y <- as.numeric(var2)
    #   GR <- rep(0L, length(x))
    #   GS <- as.integer(length(x))
    #   NG <- 1
    #   N <- as.numeric(length(x))

    #   p.value <-PharmacoGx:::partialCorQUICKSTOP(x, y, obs.cor, GR, GS, NG, 1e7,N, req_alpha, req_alpha/100, 10L, runif(2))
    #   significant <- p.value[[1]]
    #   p.value <- p.value[[2]]
    # }




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

      ## Think about if the partial cor should also be refit for each (Probably, different lines would be fit if points are missing...)

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
