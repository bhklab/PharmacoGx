########################
## Benjamin Haibe-Kains
## October 23, 2013
########################

#################################################
#'Calculate the gwc score between two vectors, using either a weighted spearman or pearson correlation
#'
#'@examples
#'data(CCLEsmall)
#'x <- molecularProfiles(CCLEsmall,"rna")[,1]
#'y <- molecularProfiles(CCLEsmall,"rna")[,2]
#'x_p <- rep(0.05, times=length(x))
#'y_p <- rep(0.05, times=length(y))
#'names(x_p) <- names(x)
#'names(y_p) <- names(y)
#'gwc(x,x_p,y,y_p, nperm=100)
#'
#'@param x1 \code{numeric} vector of effect sizes (e.g., fold change or t statitsics) for the first experiment
#'@param p1 \code{numeric} vector of p-values for each corresponding effect size for the first experiment
#'@param x2 \code{numeric} effect size (e.g., fold change or t statitsics) for the second experiment
#'@param p2 \code{numeric} vector of p-values for each corresponding effect size for the second experiment
#'@param method.cor \code{character} string identifying if a \code{pearson} or
#'\code{spearman} correlation should be used
#'@param nperm \code{numeric} how many permutations should be done to determine
#'@param truncate.p \code{numeric} Truncation value for extremely low p-values
#'@param ... Other passed down to internal functions
#'
#'@return \code{numeric} a vector of two values, the correlation and associated p-value.
#'@export
##            -
##
## http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient#Calculating_a_weighted_correlation
## http://www.mathworks.com/matlabcentral/fileexchange/20846-weighted-correlation-matrix
## F. Pozzi, T. Di Matteo, T. Aste, "Exponential smoothing weighted correlations", The European Physical Journal B, Vol. 85, No 6, 2012. DOI: 10.1140/epjb/e2012-20697-x
#################################################

gwc <-
function (x1, p1, x2, p2, method.cor=c("pearson", "spearman"), nperm=1e4, truncate.p=1e-16, ...) {
    
    method.cor <- match.arg(method.cor)
    ## intersection between x and y
    ii <- intersectList(names(x1), names(p1), names(x2), names(p2))
    if(length(ii) < 10) {
        stop ("Less than 10 probes/genes in common between x and y")
    }
    x1 <- x1[ii]
    p1 <- p1[ii]
    x2 <- x2[ii]
    p2 <- p2[ii]
    ## truncate extremely low p-values
    p1[!is.na(p1) & p1 < truncate.p] <- truncate.p
    p2[!is.na(p2) & p2 < truncate.p] <- truncate.p
    ## scaled weights
    p1 <- -log10(p1)
    p1 <- p1 / sum(p1, na.rm=TRUE)
    p2 <- -log10(p2)
    p2 <- p2 / sum(p2, na.rm=TRUE)
    w <- p1 + p2
    ## compute genome-wide connectivity score
    res <- corWeighted(x=x1, y=x2, w=w, method=method.cor, nperm=nperm, ...)
    return(res)
}
