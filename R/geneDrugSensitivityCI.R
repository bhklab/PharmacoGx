########################
## Benjamin Haibe-Kains & Petr Smirnov
## October 23, 2013
########################

#' @importFrom stats sd
#' @importFrom stats complete.cases
#' @importFrom stats lm
#' @importFrom stats glm
#' @importFrom stats anova
#' @importFrom stats pf
#' @importFrom stats formula
#' @importFrom stats var

geneDrugSensitivityCI <- function(x,  drugpheno, verbose=FALSE) {
## input:
##  x: numeric vector of gene expression values
##  type: vector of factors specifying the cell lines or type types
##  batch: vector of factors specifying the batch
##  drugpheno: numeric vector of drug sensitivity values (e.g., IC50 or AUC)
##  duration: numeric vector of experiment duration in hours
##  interaction.typexgene: Should interaction between gene expression and cell/type type be computed? Default set to FALSE
##  model: Should the full linear model be returned? Default set to FALSE
##
## output:
##  vector reporting the effect size (estimateof the coefficient of drug concentration), standard error (se), sample size (n), t statistic, and F statistics and its corresponding p-value

standardize <- match.arg(standardize)

colnames(drugpheno) <- paste("drugpheno", 1:ncol(drugpheno), sep=".")  

drugpheno <- data.frame(sapply(drugpheno, function(x) {
  if (!is.factor(x)) {
    x[is.infinite(x)] <- NA
  }
  return(list(x))
  }, USE.NAMES=FALSE), check.names=FALSE)


ccix <- complete.cases(x, type, batch, drugpheno)
nn <- sum(ccix)

if(length(table(drugpheno)) > 2){
 if(ncol(drugpheno)>1){
      ##### FIX NAMES!!!
      rest <- lapply(1:ncol(drugpheno), function(i){

        est <- paste("estimate", i, sep=".")
        se <-  paste("se", i, sep=".")
        tstat <- paste("tstat", i, sep=".")

        rest <- rep(NA, 3)
        names(rest) <- c(est, se, tstat)
        return(rest)

        })
      rest <- do.call(c, rest)
      rest <- c(rest, n=nn, "fstat"=NA, "pvalue"=NA)
      } else {
        rest <- c("estimate"=NA, "se"=NA, "n"=nn, "tstat"=NA, "fstat"=NA, "pvalue"=NA, "df"=NA)
      }
      } else {
    # rest <- c("estimate"=NA, "se"=NA, "n"=nn, "pvalue"=NA)
    if (is.factor(drugpheno[,1])){
      rest <- c("estimate"=NA_real_, "se"=NA_real_, "n"=nn, "pvalue"=NA_real_, df=NA_real_)
      } else {
        rest <- c("estimate" = NA, "se" = NA , "n" = nn, "tstat" = NA , "fstat" = NA , "pvalue" = NA , "df" = NA , "fdr" = NA)
      }
    }  

    if(nn < 3 || isTRUE(all.equal(var(x[ccix], na.rm=TRUE), 0))) {
    ## not enough samples with complete information or no variation in gene expression
    return(rest)
  }

  ## standardized coefficient in linear model 
  if(length(table(drugpheno)) > 2 & standardize!= "none") {
    switch(standardize, 
      "SD" = drugpheno <- apply(drugpheno, 2, function(x){
        return(x[ccix]/sd(as.numeric(x[ccix])))}) ,
      "rescale" = drugpheno <- apply(drugpheno, 2, function(x){
        return(rescale(as.numeric(x[ccix]), q=0.05, na.rm=TRUE))    })
      )

    }else{
      drugpheno <- drugpheno[ccix,,drop=FALSE]
    }
    if(length(table(x)) > 2  & standardize!= "none"){
      switch(standardize, 
        "SD" = xx <- x[ccix]/sd(as.numeric(x[ccix])) ,
        "rescale" = xx <- rescale(as.numeric(x[ccix]), q=0.05, na.rm=TRUE)
        )
      }else{
        xx <- x[ccix]
      }
      if(ncol(drugpheno)>1){
        ff0 <- paste("cbind(", paste(paste("drugpheno", 1:ncol(drugpheno), sep="."), collapse=","), ")", sep="")
        } else {
          ff0 <- sprintf("drugpheno.1")
        }

  # ff1 <- sprintf("%s + x", ff0)

  dd <- data.frame(drugpheno, "x"=xx)
  # , "x"=xx, "type"=type[ccix], "batch"=batch[ccix])
  
  ## control for tissue type
  if(length(sort(unique(type))) > 1) { 
    dd <- cbind(dd, type=type[ccix])
  }
  ## control for batch
  if(length(sort(unique(batch))) > 1) {
        dd <- cbind(dd, batch=batch[ccix])
  }

  CI.value <- fastCI(observations = dd[,"x"], predictions = dd[,"drugpheno.1"], outx = FALSE)$cindex
  permutation_done <- FALSE
  b.perm.par <- choose_b(alpha/n.tests, p.confidence)
  r.perm.par <- choose_r(alpha/n.tests, p.confidence)
  cur.i <- 1
  cur.stat.better <- 0
  while(!permutation_done){
    dd[,"x"] <- sample(dd[,"x"], NROW(dd))
    CI.perm <- fastCI(observations = dd[,"drugpheno.1"], predictions = dd[,"x"], outx = FALSE)$cindex
    if(abs(CI.perm - 0.5) > abs(CI.value - 0.5)){
      cur.stat.better <- cur.stat.better + 1
    }
    if(cur.stat.better >= r.perm.par){
      CI.p.val <- cur.stat.better / cur.i
      permutation_done <- TRUE
    }
    if(cur.i >= b.perm.par){
      CI.p.val <- cur.stat.better / cur.i
      permutation_done <- TRUE
    }
    cur.i <- cur.i + 1
  }
  if(CI.p.val == 0){
    CI.p.val <- 1/(b.perm.par + 1)
  }
  rest.CI <- c("CI" = CI.value, "CI.perm.p" = CI.p.val)
  rest <- rest.CI
  if(length(rest) != 2){
    rest <- c("CI" = NA_real_, "CI.perm.p" = NA_real_)
  } 
  return(rest)
}





## End
