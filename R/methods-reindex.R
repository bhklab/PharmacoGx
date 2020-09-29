# ==== LongTable Class

#' Redo indexing for a LongTable object to remove any gaps in integer indexes
#'
#' After subsetting a LongTable, it is possible that values of rowKey or colKey
#'   could no longer be present in the object. As a result there the indexes
#'   will no longer be contiguous integers. This method will calcualte a new
#'   set of rowKey and colKey values such that integer indexes are the smallest
#'   set of contiguous integers possible for the data.
#'
#' @param object The [`LongTable`] object to recalcualte indexes (rowKey and
#'     colKey values) for.
#'
#' @return A copy of the [`LongTable`] with all keys as the smallest set of
#'     contiguous integers possible given the current data.
#'
#' @export
setMethod('reindex', signature(object='LongTable'), function(object) {

    # extract assays joined to row/colData
    assayDataList <- assays(object, withDimnames=TRUE, metadata=TRUE)

    # find names of ID columns
    rowIDCols <- colnames(rowData(object))
    colIDCols <- colnames(colData(object))

    # extract the ID columns from the assay data
    newRowData <- unique(rbindlist(lapply(assayDataList,
                                          FUN=`[`,
                                          i=TRUE, j=get('rowIDCols'))
                                          ))[, rowKey := seq_len(.N)]
    setkeyv(newRowData, 'rowKey')
    newColData <- unique(rbindlist(lapply(assayDataList,
                                          FUN=`[`,
                                          i=TRUE, j=get('colIDCols'))
                                          ))[, colKey := seq_len(.N)]
    setkeyv(newColData, 'colKey')

    # remap the rowKey and colKey columns to the assays
    newAssayData <- lapply(assayDataList,
                           FUN=.joinDropOn,
                           DT2=newRowData, on=rowIDCols)
    newAssayData <- lapply(newAssayData,
                           FUN=.joinDropOn,
                           DT2=newColData, on=colIDCols)
    newAssayData <- lapply(newAssayData, setkeyv, cols=c('rowKey', 'colKey'))

    return(LongTable(rowData=newRowData, rowIDs=.rowIDs(object),
                     colData=newColData, colIDs=.colIDs(object),
                     assays=newAssayData, metadata=metadata(object)))

})

#' @keywords interal
#' @noRd
.joinDropOn <- function(DT1, DT2, on) {
    DT1[DT2, on=on][, -get('on')]
}