#' @importFrom stats median
#' @importFrom stats complete.cases
#' @importFrom stats lm
#' @importFrom stats anova
#' @importFrom stats pf
#' 

## function computing gene-drug associations from perturbation data (CMAP)
geneDrugPerturbation <- function(x, concentration, type, batch, duration, model=FALSE) {
## input:
##  x: numeric vector of gene expression values
##  concentration: numeric vector with drug concentrations/doses
##  type: vector of factors specifying the cell lines or type types
##  batch: vector of factors specifying the batch
##  duration: numeric vector of measurement times (in hours)
##  model: Should the full linear model be returned? Default set to FALSE
##
## output:
##  vector reporting the effect size (estimateof the coefficient of drug concentration), standard error (se), sample size (n), t statistic, and F statistics and its corresponding p-value
    
  
  
    nc <- c("estimate", "se", "n", "tstat", "fstat", "pvalue")
    if (length(sort(unique(concentration))) < 2) {
        warning("No drug concentrations tested")
        tt <- rep(NA, length(nc))
        names(tt) <- nc
        return(tt)
    }
    ff0 <- sprintf("x ~ 1")
    ff <- sprintf("%s + concentration", ff0)
    

    if (length(sort(unique(type))) > 1) { 
        ff0 <- sprintf("%s + type", ff0)
        ff <- sprintf("%s + type", ff)
    }
    if (length(sort(unique(batch))) > 1) {
        ff0 <- sprintf("%s + batch", ff0)
        ff <- sprintf("%s + batch", ff)
    }
  
### add experiment duration if the vector consists of more than one different value

  if(length(sort(unique(duration))) > 2){
      ff0 <- sprintf("%s + duration", ff0)
      ff <- sprintf("%s + duration", ff)
  }

    dd <- data.frame("x"=x, "concentration"=concentration, "duration"=duration, "type"=type, "batch"=batch)
    nn <- sum(complete.cases(dd))
    if(nn < 3) {
        tt <- c("estimate"=NA, "se"=NA, "n"=nn, "tsat"=NA, "fstat"=NA, "pvalue"=NA)
    } else {
        names(dd)[1]<-"x"
        mm0 <- lm(formula=ff0, data=dd, model=FALSE, x=FALSE, y=FALSE, qr=TRUE)
        mm <- lm(formula=ff, data=dd, model=model, x=FALSE, y=FALSE, qr=TRUE)

        mmc <- stats::anova(mm0, mm)
        mm <- summary(mm)
## extract statistics
        tt <- c("estimate"=mm$coefficients["concentration", "Estimate"], "se"=mm$coefficients["concentration", "Std. Error"], "n"=nn, "tsat"=mm$coefficients["concentration", "t value"], "fstat"=mmc$F[2], "pvalue"=mmc$'Pr(>F)'[2])
    }
    names(tt) <- nc
## add tissue type/cell line statistics
    if(length(sort(unique(type))) > 1) {
        rr <- summary(mm0)
        ttype <- c("type.fstat"=rr$fstatistic["value"], "type.pvalue"=pf(q=rr$fstatistic["value"], df1=rr$fstatistic["numdf"], df2=rr$fstatistic["dendf"], lower.tail=FALSE))
        names(ttype) <- c("type.fstat", "type.pvalue")
    } else { ttype <- c("type.fstat"=NA, "type.pvalue"=NA) }
    tt <- c(tt, ttype)
## add model
    if (model) { tt <- list("stats"=tt, "model"=mm)}
    return(tt)
}


## End
