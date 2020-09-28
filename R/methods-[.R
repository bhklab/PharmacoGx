# ==== PharmacoSet Class

#'`[`
#'
#' @examples
#' data(CCLEsmall)
#' CCLEsmall["WM1799", "Sorafenib"]
#'
#' @param x object
#' @param i Cell lines to keep in object
#' @param j Drugs to keep in object
#' @param ... further arguments
#' @param drop A boolean flag of whether to drop single dimensions or not
#'
#' @return Returns the subsetted object
#'
#' @export
setMethod(`[`, 'PharmacoSet', function(x, i, j, ..., drop = FALSE){
    eval(substitute(subset(x, i, j, i)))  # eval(substitute()) idiom allows
                                          # direct pass through of unevaluated
                                          # function arguments. Stops attempts
                                          # to evaluate arguments in current
                                          # scope.
})

# ==== LongTable Class

#' `[`
#'
#' @param x object
#' @param i row query
#' @param j col query
#' @param ... ad args
#' @param drop does nothing
#'
#' @export
setMethod('[', signature('LongTable'), function(x, i, j, ..., drop=FALSE) {
    eval(substitute(subset(x, i, j)))
})