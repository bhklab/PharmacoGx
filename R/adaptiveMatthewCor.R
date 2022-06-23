#' Adaptive Matthews Correlation Coefficient
#'
#' @inherit CoreGx::amcc
#'
#' @examples
#' amcc(0.6^(1:5), 0.5^(1:5))
#'
#' @export
amcc <- function(x, y, step.prct=0, min.cat=3, nperm=1000, nthread=1) {
    CoreGx::amcc(x, y, step.prct, min.cat, nperm, nthread)
}