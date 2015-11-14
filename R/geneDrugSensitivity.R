########################
## Benjamin Haibe-Kains
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

geneDrugSensitivity <- function(x, type, batch, drugpheno, interaction.typexgene=FALSE, model=FALSE, verbose=FALSE) {
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

  ccix <- complete.cases(x, type, batch, drugpheno)
  nn <- sum(ccix)
#  rest <- c("estimate"=NA, "se"=NA, "n"=nn, "tstat"=NA, "fstat"=NA, "pvalue"=NA, "type.fstat"=NA, "type.pvalue"=NA)
  rest <- c("estimate"=NA, "se"=NA, "n"=nn, "pvalue"=NA)
  
  if(nn < 3 || var(x, na.rm=TRUE) == 0) {
    ## not enough samples with complete information or no variation in gene expression
    return(rest)
  }
  ## standardized coefficient in linear model
  if(length(table(drugpheno)) > 2) {
    drugpheno <- drugpheno[ccix] / sd(drugpheno[ccix], na.rm=TRUE)
  }else{
    drugpheno <- drugpheno[ccix]
  }
  if(length(table(x)) > 2){
    xx <- x[ccix] / sd(x[ccix], na.rm=TRUE)
  }else{
    xx <- x[ccix]
  }
  
  dd <- data.frame("drugpheno"=drugpheno, "x"=xx, "type"=type[ccix], "batch"=batch[ccix])
  ff0 <- sprintf("drugpheno ~ 1")
  ff1 <- sprintf("%s + x", ff0)
  ## control for tissue type
  if(length(sort(unique(type))) > 1) { 
    ff0 <- sprintf("%s + type", ff0)
    ff1 <- sprintf("%s + type", ff1)
  }
  ## control for batch
  if(length(sort(unique(batch))) > 1) {
    ff0 <- sprintf("%s + batch", ff0)
    ff1 <- sprintf("%s + batch", ff1)
  }
  ## control for duration
  # if(length(sort(unique(duration))) > 1){
  #   ff0 <- sprintf("%s + duration", ff0)
  #   ff <- sprintf("%s + duration", ff)
  # }
  rr0 <- tryCatch(try(glm(formula(ff0), data=dd, model=FALSE, x=FALSE, y=FALSE, family=ifelse(is.factor(drugpheno), "binomial", "gaussian"))), 
                  warning=function(w) {
                    if(verbose) {
                      ww <- "Null model did not convrge"
                      tt <- table(dd[,"type"])
                      print(ww)
                      print(tt)
                    }
                  })
  rr1 <- tryCatch(try(glm(formula(ff1), data=dd, model=FALSE, x=FALSE, y=FALSE, family=ifelse(is.factor(drugpheno), "binomial", "gaussian"))), 
                  warning=function(w) {
                    if(verbose) {
                      ww <- "Model did not convrge"
                      tt <- table(dd[,"drugpheno"])
                      print(ww)
                      print(tt)
                      return(ww)
                    }
                  })

if (class(rr0) != "try-error" && class(rr1) != "try-error" & class(rr0) != "character" && class(rr1) != "character") {
    rrc <- stats::anova(rr0, rr1, test="Chisq")
    
    rr <- summary(rr1)
#    rest <- c("estimate"=rr$coefficients["x", "Estimate"], "se"=rr$coefficients["x", "Std. Error"], "n"=nn, "tsat"=rr$coefficients["x", "t value"], "fstat"=rrc$F[2], "pvalue"=rrc$'Pr(>F)'[2])
    rest <- c("estimate"=rr$coefficients[grep("^x", rownames(rr$coefficients)), "Estimate"], "se"=rr$coefficients[grep("^x", rownames(rr$coefficients)), "Std. Error"], "n"=nn, "pvalue"=rrc$'Pr(>Chi)'[2])
    
#   names(rest) <- c("estimate", "se", "n", "tstat", "fstat", "pvalue")
    names(rest) <- c("estimate", "se", "n", "pvalue")

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
