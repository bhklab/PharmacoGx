#' Return a list of `data.table` objects with the assay measurements from a
#'  `LongTable` object.
#'
#' @param x [`LongTable`] What to extract the assay data from.
#' @param withDimnames [`logical`] Should the returned assays be joined to
#'   the row and column identifiers (i.e., the pseudo dimnames of the object).
#'
#' @return
#'
#' @importMethodsFrom SummarizedExperiment assays
#' @export
setMethod('assays', signature(x='LongTable'), function(x, withDimnames=FALSE) {
    if (withDimnames)
        return(structure(
            lapply(assayNames(x), assay, x=x, withDimnames=withDimnames),
            .Names=assayNames(x))
            )
    else
        return(x@assays)
})