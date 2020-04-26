#' GWC Score
#' 
#' @inherit CoreGx::gwc
#' 
#' @examples
#' data(CCLEsmall)
#' x <- molecularProfiles(CCLEsmall,"rna")[,1]
#' y <- molecularProfiles(CCLEsmall,"rna")[,2]
#' x_p <- rep(0.05, times=length(x))
#' y_p <- rep(0.05, times=length(y))
#' names(x_p) <- names(x)
#' names(y_p) <- names(y)
#' gwc(x,x_p,y,y_p, nperm=100)
#' 
#' @export
gwc <-
function (x1, p1, x2, p2, method.cor=c("pearson", "spearman"), nperm=1e4, 
          truncate.p=1e-16, ...) {
  CoreGx::gwc(x1, p1, x2, p2, method.cor, nperm, truncate.p, ...)
}