## Matthews correlation coefficient
#' Compute a Matthews Correlation Coefficient
#'
#' @inherit CoreGx::mcc
#'
#' @export
mcc <- function(x, y, nperm = 1000, nthread = 1) {
  CoreGx::mcc(x, y, nperm, nthread)
}
