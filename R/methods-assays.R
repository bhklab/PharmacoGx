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


#'
#'
#'
#' @importMethodsFrom SummarizedExperiment assays<-
#' @export
setReplaceMethod('assays', signature(x='LongTable', value='list'), function(x, value) {

    # check input is correct
    isDF <- is.items(value, 'data.frame')
    isDT <- is.items(value, 'data.table')
    if (!all(isDF))
        stop(.errorMsg('Items ', which(!isDT)), ' in value are not data.tables or
            data.frames. These are the only types allowed in the value argument!',
            collapse=', ')

    if (!all(isDT))
        for (i in which(!isDT)) setDT(values[[i]])

    # extract the row and column values


})