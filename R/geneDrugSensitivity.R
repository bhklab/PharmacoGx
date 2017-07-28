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

geneDrugSensitivity <- function(x, type, batch, drugpheno, interaction.typexgene=FALSE, model=FALSE,  standardize=c("SD", "rescale", "none"), verbose=FALSE) {
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
    rest <- c("estimate"=NA, "se"=NA, "n"=nn, "pvalue"=NA)
  }  
  
  if(nn < 3 || var(x[ccix], na.rm=TRUE) == 0) {
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

rr0 <- tryCatch(try(glm(formula(drugpheno.1 ~ . - x), data=dd, model=FALSE, x=FALSE, y=FALSE, family="binomial")), 
    warning=function(w) {
      if(verbose) {
        ww <- "Null model did not convrge"
        print(ww)
        if("type" %in% colnames(dd)) {
          tt <- table(dd[,"type"])
          print(tt)
        }
      }
    })
  rr1 <- tryCatch(try(glm(formula(drugpheno.1 ~ .), data=dd, model=FALSE, x=FALSE, y=FALSE, family="binomial")), 
    warning=function(w) {
      if(verbose) {
        ww <- "Model did not converge"
        tt <- table(dd[,"drugpheno.1"])
        print(ww)
        print(tt)
      }
      return(ww)
    })


} else{

rr0 <- tryCatch(try(lm(formula(paste(ff0, "~ . -x", sep=" ")), data=dd)), 
    warning=function(w) {
      if(verbose) {
        ww <- "Null model did not converge"
        print(ww)
        if("type" %in% colnames(dd)) {
          tt <- table(dd[,"type"])
          print(tt)
        }      
      }
    })
  rr1 <- tryCatch(try(lm(formula(paste(ff0, "~ . ", sep=" ")), data=dd)), 
    warning=function(w) {
      if(verbose) {
        ww <- "Model did not converge"
        tt <- table(dd[,"drugpheno.1"])
        print(ww)
        print(tt)
      }
      return(ww)
    })


}
  
  
  if (class(rr0) != "try-error" && class(rr1) != "try-error" & class(rr0) != "character" && class(rr1) != "character") {
    rr <- summary(rr1)

    if(any(unlist(lapply(drugpheno,is.factor)))){
      rrc <- stats::anova(rr0, rr1, test="Chisq")
      rest <- c("estimate"=rr$coefficients[grep("^x", rownames(rr$coefficients)), "Estimate"], "se"=rr$coefficients[grep("^x", rownames(rr$coefficients)), "Std. Error"], "n"=nn, "pvalue"=rrc$'Pr(>Chi)'[2])
      names(rest) <- c("estimate", "se", "n", "pvalue")

    } else {
      if(ncol(drugpheno)>1){
        rrc <- summary(stats::manova(rr1))
        rest <- lapply(1:ncol(drugpheno), function(i) {
          est <- paste("estimate", i, sep=".")
          se <-  paste("se", i, sep=".")
          tstat <- paste("tstat", i, sep=".")
          rest <- c(rr[[i]]$coefficients[grep("^x", rownames(rr[[i]]$coefficients)), "Estimate"], rr[[i]]$coefficients[grep("^x", rownames(rr[[i]]$coefficients)), "Std. Error"], rr[[i]]$coefficients[grep("^x", rownames(rr[[i]]$coefficients)), "t value"])
          names(rest) <- c(est, se, tstat)
          return(rest)
        })
        rest <- do.call(c, rest)
        rest <- c(rest,"n"=nn, "fstat"=rrc$stats[grep("^x", rownames(rrc$stats)), "approx F"], "pvalue"=rrc$stats[grep("^x", rownames(rrc$stats)), "Pr(>F)"])
      } else {
        rrc <- stats::anova(rr0, rr1, test = "F") 
        rest <- c("estimate"=rr$coefficients[grep("^x", rownames(rr$coefficients)), "Estimate"], "se"=rr$coefficients[grep("^x", rownames(rr$coefficients)), "Std. Error"],"n"=nn, "tstat"=rr$coefficients[grep("^x", rownames(rr$coefficients)), "t value"], "fstat"=rrc$F[2], "pvalue"=rrc$'Pr(>F)'[2], "df"=rr1$df.residual)
        names(rest) <- c("estimate", "se", "n", "tstat", "fstat", "pvalue", "df")
      }
    }
    
    
#    rest <- c("estimate"=rr$coefficients["x", "Estimate"], "se"=rr$coefficients["x", "Std. Error"], "n"=nn, "tsat"=rr$coefficients["x", "t value"], "fstat"=rrc$F[2], "pvalue"=rrc$'Pr(>F)'[2])
    
#   names(rest) <- c("estimate", "se", "n", "tstat", "fstat", "pvalue")

## add tissue type/cell line statistics
#     if(length(sort(unique(type))) > 1) {
#       rr <- summary(rr0)
#       ttype <- c("type.fstat"=rr$fstatistic["value"], "type.pvalue"=pf(q=rr$fstatistic["value"], df1=rr$fstatistic["numdf"], df2=rr$fstatistic["dendf"], lower.tail=FALSE))
#       names(ttype) <- c("type.fstat", "type.pvalue")
#     } else { ttype <- c("type.fstat"=NA, "type.pvalue"=NA) }
#     rest <- c(rest, ttype)
    ## add model
    if(model) { rest <- list("stats"=rest, "model"=rr1) }
  }
  return(rest)
}





## End
