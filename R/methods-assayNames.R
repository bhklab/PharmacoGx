#' Retrieve the assay names from a `LongTable` object.
#'
#' @param x A [`LongTable`] object to retrieve the assay names from.
#'
#' @return [`character`] Names of the assays contained in the `LongTable`.
#'
#' @importMethodsFrom SummarizedExperiment assayNames
#' @export
setMethod('assayNames', signature(x='LongTable'), function(x) {
    return(names(assays(x)))
})