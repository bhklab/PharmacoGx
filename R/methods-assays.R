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
##TODO:: Add key argument with default to FALSE to remove rowKey and colKey
setMethod('assays', signature(x='LongTable'),
    function(x, withDimnames=FALSE, metadata=FALSE, key=TRUE) {

    return(structure(
        lapply(assayNames(x),
               FUN=assay,
               x=x, withDimnames=withDimnames, metadata=metadata, key=key),
               .Names=assayNames(x)))
})


#'
#'
#' @param x A [`LongTable`] to modify the assays in.
#' @param value A [`list`] of `data.frame` or `data.table` objects, all of which
#'   contain the row and column identifiers and metadata.
#'
#' @return A copy of the [`LongTable`] with the assays modified
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

    # check new assay names
    if (is.null(names(value))) {
        warning(.warnMsg('\nThe list being assigned to assays has no names.',
            'Defaulting to numbered assays. You can correct his with',
            'assayNames(x) <- value.'))
        names(value) <- paste0('assay', seq_along(value))
    }

    # extract the row and column values
    rowIDCols <- colnames(.rowIDData(x)[, -'rowKey'])
    colIDCols <- colnames(.colIDData(x)[, -'colKey'])
    rowMetaCols <- setdiff(colnames(rowData(x)), rowIDCols)
    colMetaCols <- setdiff(colnames(colData(x)), colIDCols)

    # get the rowData and colData column mappings
    rowDataCols <- if (length(rowMetaCols) > 0) list(rowIDCols, rowMetaCols) else list(rowIDCols)
    colDataCols <- if (length(colMetaCols) > 0) list(colIDCols, colMetaCols) else list(colIDCols)

    # get assay column names
    allCols <- c(unlist(rowDataCols), unlist(colDataCols))
    assayCols <- lapply(value, colnames)
    assayCols <- lapply(assayCols, setdiff, y=allCols)
    names(assayCols) <- names(value)

    # reconstruct a new LongTable
    buildLongTable(from=value, rowDataCols, colDataCols, assayCols)

})