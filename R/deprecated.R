#### Aliases for deprecated methods in PharmacoGx

#' Deprecated method wanring for subsetTo method
#'
#' Replace by subset method
#'
#' @importFrom crayon magenta cyan
#' @export
#' @noRd
subsetTo <- function(...) {
    warning(cyan$bold("The subsetTo method is now deprecated, please use subset instead"))
    eval(substitute(subset(...)))
}

##' Deprecated method warning for intersectPSet method
##'
##' Replaced by intersect method
##'
##' @importFrom crayon cyan magenta
##' @export
##' @noRd
#intersectPSet <- function(...) {
#    warning(cyan$bold("The intersectPSet method is not depreacted, please use intersect instead"))
#    eval(substitute(intersect(...)))
#}