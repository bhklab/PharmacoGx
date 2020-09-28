#' Get the dimensions of a `LongTable` object.
#'
#' @param x A [`LongTable`] object to retrieve dimensions for.
#'
#' @return [`numeric`] Vector of object dimensions.
#'
#' @importMethodsFrom SummarizedExperiment dim
#' @export
setMethod('dim', signature(x='LongTable'), function(x) {
    return(c(nrow(rowData(x)), nrow(colData(x))))
})