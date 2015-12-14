#' @importFrom stats pchisq
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @importFrom stats pt

combineTest <-
function(p, weight, method=c("fisher", "z.transform", "logit"), hetero=FALSE, na.rm=FALSE) {
    if(hetero) { stop("function to deal with heterogeneity is not implemented yet!") }
    method <- match.arg(method)
    na.ix <- is.na(p)
    if(any(na.ix) && !na.rm) { stop("missing values are present!") }
    if(all(na.ix)) { return(NA) } ## all p-values are missing
    p <- p[!na.ix]
    k <- length(p)
    if(k == 1) { return(p) }
    if(missing(weight)) { weight <- rep(1, k); }
    switch(method,  
    "fisher"={
        cp <- pchisq(-2 * sum(log(p)), df=2*k, lower.tail=FALSE)
    }, 
    "z.transform"={
        z <- qnorm(p, lower.tail=FALSE)
        cp <- pnorm(sum(weight * z) / sqrt(sum(weight^2)), lower.tail=FALSE)
    }, 
    "logit"={
        tt <- (- sum(log(p / (1 - p)))) / sqrt(k * pi^2 * (5 * k + 2) / (3 * (5 * k + 4)))
        cp <- pt(tt,df=5*k+4, lower.tail=FALSE)
    })
    return(cp)
}