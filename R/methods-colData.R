#' Retrieve the column metadata table from a LongTable object
#'
#' @param x A [`LongTable`] to retrieve column metadata from.
#' @param key [`logical`] She the colKey column also be returned? Defaults to
#'     FALSE.
#'
#' @return A [`data.table`] containing row identifiers and metadata.
#'
#' @import data.table
#' @export
setMethod('colData', signature(x='LongTable'), function(x, key=FALSE) {
    return(if (key) x@colData[, -'.colnames'] else
        x@colData[, -c('.colnames', 'colKey')])
})

#' Updates the `colData` slot as long as the ID columns are not changed.
#'
#' @param x A [`LongTable`] object to modify.
#' @param value A [`data.table`] or [`data.frame`] to update with.
#'
#' @return A copy of the [`LongTable`] object with the `colData`
#'   slot updated.
#'
#' @importFrom crayon cyan magenta
#' @importFrom SummarizedExperiment colData<-
#' @export
setReplaceMethod('colData', signature(x='LongTable'), function(x, value) {

    # type check input
    if (is(value, 'data.frame'))
        setDT(value)
    if (!is(value, 'data.table'))
        stop(magenta$bold('Please pass a data.frame or data.table to update
            the colData slot. We recommend modifying the object returned by
            colData(x) then reassigning it with colData(x) <- newColData'))

    # remove key column
    ## TODO: Do I need this warning?
    if ('colKey' %in% colnames(value)) {
        value[, colKey := NULL]
        warning(cyan$bold('Dropping colKey from replacement value, this
            function will deal with mapping the colKey automatically.'))
    }

    # assemble information to select proper update method
    colIDCols <- colIDs(x)
    sharedColIDCols <- intersect(colIDCols, colnames(value))

    metadataCols <- setdiff(sharedColIDCols, colnames(colData(x)))
    sharedMetadataCols <- intersect(metadataCols, colnames(value))

    ## TODO:: Add message indicating if metadata columns change
    ## TODO:: Add parameter that allows modification of colID cols which errors if they change when FALSE

    # error if all the colID columns are not present in the new colData
    equalColIDs <- sharedColIDCols %in% colIDCols
    if (!all(equalColIDs)) stop(.errorMsg('The ID columns ',
        sharedIDCols[!equalColIDs]), 'are not present in value. Currently
            this function only supports updates with the same
            colID columns as the current colData!')

    colIDs <- colIDs(x, data=TRUE, key=TRUE)

    colData <- colIDs[value, on=sharedColIDCols]

    ## TODO:: Refactor .pasteColons into utilities.R or make a function to update .colnames
    .pasteColons <- function(...) paste(..., collapse=':')
    colData[, `:=`(.colnames=mapply(.pasteColons, transpose(.SD))), .SDcols=sharedColIDCols]
    setkey(colData, 'colKey')

    x@colData <- colData
    x
})