corPermute <- function(obs.cor, vec1, vec2, nPerm=1e3){

  num.larger <- 0


  for(i in seq_len(nPerm)){
    vec1 <- sample(vec1, length(vec1))
    
    perm.cor <- coop::pcor(vec1, vec2, use="complete.obs")

    if(abs(perm.cor) >= abs(obs.cor)){
      num.larger <- num.larger + 1
    }
  }

  return(num.larger/nPerm)

}

partialPermute <- function(obs.cor, data, formula1, formula2, nPerm = 100){
  
  num.larger <- 0


  for(i in seq_len(nPerm)){
    data[,1] <- sample(data[,1], nrow(data))
    data[,2] <- sample(data[,2], nrow(data))

    partial.dp <- residuals(lm(formula(formula1), data))
    partial.x <- residuals(lm(formula(formula2), data))

    perm.cor <- coop::pcor(partial.dp, partial.x, use="complete.obs")

    if(abs(perm.cor) >= abs(obs.cor)){
      num.larger <- num.larger + 1
    }
  }

  return(num.larger/nPerm)
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


#' Calcualte The Gene Drug Sensitivity
#' 
#' This version of the function uses a partial correlation instead of standardized linear models. 
#' 
#' @param x A \code{numeric} vector of gene expression values
#' @param type A \code{vector} of factors specifying the cell lines or type types
#' @param batch A \code{vector} of factors specifying the batch
#' @param drugpheno A \code{numeric} vector of drug sensitivity values (e.g., 
#'   IC50 or AUC)
#' @param test A \code{character} string indicating whether resampling or analytic based tests should be used
#' @param nPerm \code{numeric}, number of permutations for p value calculation
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
  nPerm = 1e3,
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

  rest <- c("estimate"=NA_real_, "n"=as.numeric(nn), "df"=NA_real_, "pvalue"=NA_real_, "lower" = NA_real_, "upper" = NA_real_)

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
      var1 <- residuals(lm(formula(ffd), dd))
      var2 <- residuals(lm(formula(ffx), dd))
      ## TODO:: Fix degrees of freedom calculation for partial cor - binary variables use more degrees of freedom!
      controlled.var <- ncol(dd[,!colnames(dd) %in% c("drugpheno.1", "x"), drop=F])
      df <- nn - 2L - controlled.var
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
        p.value <- partialPermute(obs.cor, data=dd, formula1 = ffd, formula2 = ffx, nPerm = nPerm)

      } else {
        p.value <- corPermute(obs.cor, var1, var2, nPerm = nPerm)

      }

      boot.out <- boot(data.frame(var1, var2), cor.boot, R=nBoot)
      cint <- boot.ci(boot.out, conf = conf.level, type="bca")$bca[,4:5]
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

    }

  }

  rest <- c("estimate"=obs.cor, "n"=nn, df=df, "pvalue"=p.value, "lower" = cint[1], "upper" = cint[2])


  return(rest)
}
