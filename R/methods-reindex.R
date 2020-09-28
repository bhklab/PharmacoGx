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

    assayDataList <- assays(object, withDimnames=TRUE)

})