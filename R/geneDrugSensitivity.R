########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

geneDrugSensitivity <- function(x, type, batch, drugpheno, interaction.typexgene=FALSE, model=FALSE) {
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

  nc <- c("estimate", "se", "n", "tstat", "fstat", "pvalue")
  ccix <- complete.cases(x, type, batch, drugpheno)
  nn <- sum(ccix)
  if(nn < 3 || var(x, na.rm=TRUE) == 0) {
    ## not enough samples with complete information or no variation in gene expression
    return(c("estimate"=NA, "se"=NA, "n"=nn, "tstat"=NA, "fstat"=NA, "pvalue"=NA, "type.fstat"=NA, "type.pvalue"=NA))
  }
  ## standardized coefficient in linear model
  drugpheno <- drugpheno[ccix] / sd(drugpheno[ccix], na.rm=TRUE)
  xx <- x[ccix] / sd(x[ccix], na.rm=TRUE)
  dd <- data.frame("drugpheno"=drugpheno, "x"=xx, "type"=type[ccix], "batch"=batch[ccix])
  ff0 <- sprintf("drugpheno ~ 1")
  ff1 <- sprintf("%s + x", ff0)
  if(length(sort(unique(type))) > 1) { 
    ff0 <- sprintf("%s + type", ff0)
    ff1 <- sprintf("%s + type", ff1)
  }
  if(length(sort(unique(batch))) > 1) {
    ff0 <- sprintf("%s + batch", ff0)
    ff1 <- sprintf("%s + batch", ff1)
  }
  # if(length(sort(unique(duration))) > 1){
#   	ff0 <- sprintf("%s + duration", ff0)
#   	ff <- sprintf("%s + duration", ff)
#   }
  rr0 <- lm(formula(ff0), data=dd, model=FALSE, x=FALSE, y=FALSE, qr=TRUE)
  rr1 <- lm(formula(ff1), data=dd, model=FALSE, x=FALSE, y=FALSE, qr=TRUE)
  rrc <- anova(rr0, rr1)
  rr <- summary(rr1)
  tt <- c("estimate"=rr$coefficients["x", "Estimate"], "se"=rr$coefficients["x", "Std. Error"], "n"=nn, "tsat"=rr$coefficients["x", "t value"], "fstat"=rrc$F[2], "pvalue"=rrc$'Pr(>F)'[2])
  names(tt) <- nc
  ## add tissue type/cell line statistics
  if(length(sort(unique(type))) > 1) {
    rr <- summary(rr0)
    ttype <- c("type.fstat"=rr$fstatistic["value"], "type.pvalue"=pf(q=rr$fstatistic["value"], df1=rr$fstatistic["numdf"], df2=rr$fstatistic["dendf"], lower.tail=FALSE))
    names(ttype) <- c("type.fstat", "type.pvalue")
  } else { ttype <- c("type.fstat"=NA, "type.pvalue"=NA) }
  tt <- c(tt, ttype)
  ## add model
  if(model) { tt <- list("stats"=tt, "model"=rr1) }
  return(tt)
}





## End
