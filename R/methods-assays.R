#' Return a list of `data.table` objects with the assay measurements from a
#'  `LongTable` object.
#'
#' @param x [`LongTable`] What to extract the assay data from.
#' @param withDimnames [`logical`] Should the returned assays be joined to
#'   the row and column identifiers (i.e., the pseudo dimnames of the object).
#' @param metadata [`logical`] Should row and column metadata also be joined
#'   to the returned assays. This is useful for modifying assays before
#'   reconstructing a new LongTable.
#'
#' @return
#'
#' @importMethodsFrom SummarizedExperiment assays
#' @export
setMethod('assays', signature(x='LongTable'), function(x, withDimnames=FALSE, metadata=FALSE) {

    if (withDimnames)
        return(structure(
            lapply(assayNames(x),
                   FUN=assay,
                   x=x, withDimnames=withDimnames, metadata=metadata),
            .Names=assayNames(x)))
    else
        return(x@assays)
})