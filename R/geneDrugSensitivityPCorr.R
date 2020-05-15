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
    # perm.cor <- coop::pcor(vec1, vec2, use="complete.obs")

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

partialPermute <- function(obs.cor, data, formula1, formula2, req_alpha = 0.05, tolerance_par = req_alpha*0.1, error_prob = 1e-10){
  
  num.larger <- 0


  for(i in seq_len(req_alpha)){
    data[,1] <- sample(data[,1], nrow(data))
    data[,2] <- sample(data[,2], nrow(data))

    partial.dp <- residuals(lm(formula(formula1), data))
    partial.x <- residuals(lm(formula(formula2), data))

    perm.cor <- coop::pcor(partial.dp, partial.x, use="complete.obs")

    if(abs(perm.cor) >= abs(obs.cor)){
      num.larger <- num.larger + 1
    }
  }

  return(num.larger/req_alpha)
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


cor.boot <- function(data, w){
  return(coop::pcor(data[w,1], data[w,2]))
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
#' @param conf.leve \code{numeric}, between 0 and 1. Size of the confidence interval required
#' @param verbose \code{boolean} Should the function display messages?
#'  
#' @return A \code{vector} reporting the effect size (estimateof the coefficient 
#'   of drug concentration), standard error (se), sample size (n), t statistic, 
#'   and F statistics and its corresponding p-value.
#'
#' @importFrom stats sd complete.cases lm glm anova pf formula var
#' @importFrom boot boot boot.ci
geneDrugSensitivityPCorr <- function(x, type, batch, drugpheno, 
  test = c("resampling", "analytic"),
  req_alpha = 0.05,
  nBoot = 1e3,
  conf.level = 0.95,
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

  if(nn < 3 || isTRUE(all.equal(var(x[ccix], na.rm=TRUE), 0))) {
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
  if(length(sort(unique(type))) > 1) { 
    dd <- cbind(dd, type=type[ccix])
  }
  ## control for batch
  if(length(sort(unique(batch))) > 1) {
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
      df <- lm1$df[2] - 2L # taking the residual degrees of freedom minus 2 parameters estimated for pearson cor. 
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



        p.value <- corPermute(sample_function, req_alpha = req_alpha)
        significant <- p.value$significant
        p.value <- p.value$p.value
        pcor.boot <- function(ddd, w){
          ddd <- ddd[w,] 

          partial.dp <- residuals(lm(formula(ffd), ddd))
          partial.x <- residuals(lm(formula(ffx), ddd))

          return(coop::pcor(partial.dp, partial.x, use="complete.obs"))
        }
        boot.out <- boot(dd, pcor.boot, R=nBoot)
        cint <- boot.ci(boot.out, conf = conf.level, type="bca")$bca[,4:5]
      } else {

        sample_function <- function(){
          v1 <- sample(var1)
          return(abs(obs.cor) < abs(coop::pcor(v1, var2, use="complete.obs")))
        }



        p.value <- corPermute(sample_function, req_alpha = req_alpha)
        significant <- p.value$significant
        p.value <- p.value$p.value
        boot.out <- boot(data.frame(var1, var2), cor.boot, R=nBoot)
        cint <- boot.ci(boot.out, conf = conf.level, type="bca")$bca[,4:5]
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
