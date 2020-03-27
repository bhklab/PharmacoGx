## Matthews correlatipon coefficient
#' Compute a Mathews Correlation Coefficient 
#' 
#' @inherit CoreGx::mcc
#' 
#' @export
mcc <- 
  function(x, y, nperm=1000, nthread=1) {
    CoreGx::mcc(x, y, nperm, nthread)
  }