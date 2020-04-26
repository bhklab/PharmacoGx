#' Adaptive Matthews Correlation Coefficient
#' 
#' @inherit CoreGx::amcc
#' 
#' @export
amcc <- 
  function(x, y, step.prct=0, min.cat=3, nperm=1000, nthread=1) {
    CoreGx::amcc(x, y, step.prct, min.cat, nperm, nthread)
  }
