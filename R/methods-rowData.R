# ==== LongTable Class

#' Retrieve the row metadata table from a LongTable object
#'
#' @param x A [`LongTable`] object to retrieve the row metadata from.
#' @param key [`logical`] Should the rowKey column also be returned? Defaults
#'     to FALSE.
#' @param use.names [`logical`] This parameter is just here to stop matching
#'     they key argument to use.names from the rowData generic. It doesn't do
#'     anything.
#'
#' @return A [`data.table`] containing rowID, row identifiers, and row metadata.
#'
#' @import data.table
#' @export
setMethod('rowData', signature(x='LongTable'), function(x, key=FALSE, use.names=FALSE) {
    return(if (key) x@rowData[, -'.rownames'] else
        x@rowData[, -c('.rownames', 'rowKey')])
})

#' Updates the `rowData` slot as long as the ID columns are not changed.
#'
#' @param x A [`LongTable`] object to modify.
#' @param value A [`data.table`] or [`data.frame`] to update
#'
#' @return A copy of the [`LongTable`] object with the `rowData`
#'   slot updated.
#'
#' @importFrom crayon cyan magenta
#' @importFrom SummarizedExperiment `rowData<-`
#' @export
setReplaceMethod('rowData', signature(x='LongTable'), function(x, value) {

    # type check input
    if (is(value, 'data.frame'))
        value <- data.table(value, keep.rownames=FALSE)
    if (!is(value, 'data.table'))
        stop(magenta$bold('Please pass a data.frame or data.table to update
            the rowData slot. We recommend modifying the object returned by
            rowData(x) then reassigning it with rowData(x) <- newRowData'))

    # remove key column
    if ('rowKey' %in% colnames(value)) {
        value[, rowKey := NULL]
        warning(cyan$bold('Dropping rowKey from replacement value, this
            function will deal with mapping the rowKey automatically.'))
    }

    # assemble information to select proper update method
    rowIDCols <- colnames(.rowIDData(x))
    sharedRowIDCols <- intersect(rowIDCols, colnames(value))

    metadataCols <- colnames(rowData(x)[, -c(rowIDCols, 'rowKey')])
    sharedMetadataCols <- intersect(metadataCols, colnames(value))

    # case where no row ids are in update
    if (length(sharedRowIDCols) > 0) {
        if (length(sharedMetadataCols) > 0) {

        } else {

        }
    }

    # case where row ids are in the updated table
    if (all(rowIDCols %in% sharedRowIDCols)) {

    } else {

    }

})



#'
#'
#' @param rowIDs [`LongTable`]
#' @param rowData [`data.table`]
#'
#' @export
#' @keywords internal
.joinOnRowIDs <- function(rowIDs, rowData) {
    if (!('rowKey' %in% colnames(rowIDs)))
        stop(magenta$bold("No rowKey column in rowIDs?"))
    row_data <- rowIDs[rowData, on=colnames(rowIDs)]
    setkeyv(row_data, 'rowKey')
    return(row_data)
}

#' Helper to determine if the rowData rowIDs are identical between a new
#'   rowData data.table object and an existing `LongTable`
#'
#' @export
#' @keywords internal
.row_IDs_are_identical <- function(rowIDs, newRowData) {
    if ('rowKey' %in% colnames(rowIDs)) rowIDs[, rowKey := NULL]
    return(all(colnames(rowIDs) %in% colnames(newRowData)) &&
        all.equal(rowIDs, newRowData[, colnames(rowIDs), with=FALSE]))
}